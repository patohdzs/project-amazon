import math
import time
from functools import partial
from itertools import product

import numpy as np
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    NonNegativeReals,
    Objective,
    Param,
    RangeSet,
    Var,
    maximize,
)
from pyomo.opt import SolverFactory

from pysrc.optimization import np_to_dict


def solve_planner_problem(
    T,
    theta,
    gamma,
    x0,
    z0,
    zbar,
    pa_paths,
    pa_path_probs,
    dt=1,
    pe=20.76,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,
):
    model = ConcreteModel()

    # Indexing sets for time and sites
    model.T = RangeSet(T + 1)
    model.S = RangeSet(gamma.size)
    model.P = RangeSet(pa_paths.shape[0])

    # Parameters
    model.x0 = Param(model.S, model.P, initialize=lambda model, s, p: x0[s - 1, p - 1])
    model.z0 = Param(model.S, model.P, initialize=lambda model, s, p: z0[s - 1, p - 1])
    model.zbar = Param(model.S, initialize=np_to_dict(zbar))
    model.gamma = Param(model.S, initialize=np_to_dict(gamma))
    model.theta = Param(model.S, initialize=np_to_dict(theta))
    model.delta = Param(initialize=delta)
    model.pe = Param(initialize=pe)
    model.alpha = Param(initialize=alpha)
    model.kappa = Param(initialize=kappa)
    model.zeta = Param(initialize=zeta)
    model.dt = Param(initialize=dt)

    # Cattle price parameters
    model.pa = Param(
        model.P, model.T, initialize=lambda model, p, t: pa_paths[p - 1, t - 1]
    )
    model.pa_prob = Param(model.P, initialize=lambda model, p: pa_path_probs[p - 1])

    # Variables
    model.w = Var(model.T, model.P)
    model.x = Var(model.T, model.S, model.P)
    model.z = Var(model.T, model.S, model.P, within=NonNegativeReals)
    model.u = Var(model.T, model.S, model.P, within=NonNegativeReals)
    model.v = Var(model.T, model.S, model.P, within=NonNegativeReals)

    # Law of motion constraints
    model.zdot_def = Constraint(model.T, model.S, model.P, rule=_zdot_const)
    model.xdot_def = Constraint(model.T, model.S, model.P, rule=_xdot_const)

    # Adjustment costs constraint
    model.w_def = Constraint(model.T, model.P, rule=_w_const)

    # Price path constraints
    model.zp_def = Constraint(
        model.T, model.S, model.P, rule=partial(_zp_const, tau=max(model.P))
    )
    model.xp_def = Constraint(
        model.T, model.S, model.P, rule=partial(_xp_const, tau=max(model.P))
    )
    model.up_def = Constraint(
        model.T, model.S, model.P, rule=partial(_up_const, tau=max(model.P))
    )
    model.vp_def = Constraint(
        model.T, model.S, model.P, rule=partial(_vp_const, tau=max(model.P))
    )
    model.wp_def = Constraint(
        model.T, model.P, rule=partial(_wp_const, tau=max(model.P))
    )

    # Define the objective
    model.obj = Objective(rule=_planner_obj, sense=maximize)

    # Initial and terminal conditions
    for s in model.S:
        for p in model.P:
            model.z[:, s, p].setub(model.zbar[s])
            model.x[min(model.T), s, p].fix(model.x0[s, p])
            model.z[min(model.T), s, p].fix(model.z0[s, p])
            model.u[max(model.T), s, p].fix(0)
            model.v[max(model.T), s, p].fix(0)
            model.w[max(model.T), p].fix(0)

    # Solve the model
    solver = SolverFactory("gurobi")

    print("Solving the optimization problem...")
    start_time = time.time()
    solver.solve(model, tee=True)
    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")

    Z = np.array(
        [[[model.z[t, s, p].value for p in model.P] for s in model.S] for t in model.T]
    )
    X = np.array(
        [[[model.x[t, s, p].value for p in model.P] for s in model.S] for t in model.T]
    )
    U = np.array(
        [[[model.u[t, s, p].value for p in model.P] for s in model.S] for t in model.T]
    )
    V = np.array(
        [[[model.v[t, s, p].value for p in model.P] for s in model.S] for t in model.T]
    )
    w = np.array([[model.w[t, p].value for p in model.P] for t in model.T])

    return {
        "Z": Z,
        "X": X,
        "U": U,
        "V": V,
        "w": w,
    }


def _planner_obj(model):
    return sum(
        model.pa_prob[p]
        * sum(
            math.exp(-model.delta * (t * model.dt - model.dt))
            * (
                -model.pe
                * pyo.quicksum(
                    model.kappa * model.z[t, s, p]
                    - (model.x[t + 1, s, p] - model.x[t, s, p]) / model.dt
                    for s in model.S
                )
                + model.pa[p, t]
                * sum(model.theta[s] * model.z[t, s, p] for s in model.S)
                - model.zeta / 2 * (model.w[t, p] ** 2)
            )
            * model.dt
            for t in model.T
            if t < max(model.T)
        )
        for p in model.P
    )


def _zdot_const(model, t, s, p):
    if t < max(model.T):
        return (model.z[t + 1, s, p] - model.z[t, s, p]) / model.dt == (
            model.u[t, s, p] - model.v[t, s, p]
        )
    else:
        return Constraint.Skip


def _xdot_const(model, t, s, p):
    if t < max(model.T):
        return (model.x[t + 1, s, p] - model.x[t, s, p]) / model.dt == (
            -model.gamma[s] * model.u[t, s, p]
            - model.alpha * model.x[t, s, p]
            + model.alpha * model.gamma[s] * (model.zbar[s] - model.z[t, s, p])
        )
    else:
        return Constraint.Skip


def _w_const(model, t, p):
    if t < max(model.T):
        return model.w[t, p] == pyo.quicksum(
            model.u[t, s, p] + model.v[t, s, p] for s in model.S
        )
    else:
        return Constraint.Skip


def _zp_const(model, t, s, p, tau):
    if t <= tau:
        pivot = _pivot_point(t, p, tau)
        return model.z[t, s, p] == model.z[t, s, pivot]
    else:
        return Constraint.Skip


def _xp_const(model, t, s, p, tau):
    if t <= tau:
        pivot = _pivot_point(t, p, tau)
        return model.x[t, s, p] == model.x[t, s, pivot]
    else:
        return Constraint.Skip


def _up_const(model, t, s, p, tau):
    if t <= tau:
        pivot = _pivot_point(t, p, tau)
        return model.u[t, s, p] == model.u[t, s, pivot]
    else:
        return Constraint.Skip


def _vp_const(model, t, s, p, tau):
    if t <= tau:
        pivot = _pivot_point(t, p, tau)
        return model.v[t, s, p] == model.v[t, s, pivot]
    else:
        return Constraint.Skip


def _wp_const(model, t, p, tau):
    if t <= tau:
        pivot = _pivot_point(t, p, tau)
        return model.w[t, p] == model.w[t, pivot]
    else:
        return Constraint.Skip


def _pivot_point(t, p, tau):
    # Go back to zero-based indexing
    t = t - 1

    # Compute segment size
    segment_size = 2 ** (tau - t)

    # Find split points
    points = [segment_size * i + 1 for i in range(2**t)]

    # Find greatest upper bound among split points
    for point in reversed(points):
        if point <= p:
            return point
    return None


def price_paths(T, tau, states, start_high=True):
    """
    Get matrix of all possible price paths
    """
    # Compute initial price state
    p0 = np.ones(2**tau, np.int8)
    p0 = states[1] * p0 if start_high else states[0] * p0

    # Get stochastic and deterministic horizon paths
    rand_hzn = np.array(list(product(states, repeat=tau)))
    det_hzn = np.tile(rand_hzn[:, -1][:, np.newaxis], T - tau)

    # Merge into full price path
    price_paths = np.column_stack((p0, rand_hzn, det_hzn))
    return price_paths


def price_path_probs(M, tau, paths):
    """
    Compute likelihood of a price path realization
    """
    # Turn into binary matrix (1 for high, 0 for low)
    binary_paths = np.where(paths == paths.max(), 1, 0)

    # Find the probability of each path
    return [_price_path_prob(M, tau, path) for path in binary_paths]


def _price_path_prob(M, tau, price_path):
    prob = 1.0

    for t in range(tau):
        current_state = price_path[t]
        next_state = price_path[t + 1]
        prob *= M[current_state, next_state]

    return prob
