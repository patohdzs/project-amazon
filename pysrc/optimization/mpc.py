import math
import time
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
    model.x0 = Param(model.S, initialize=np_to_dict(x0))
    model.z0 = Param(model.S, initialize=np_to_dict(z0))
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
    model.w = Var(model.T)
    model.x = Var(model.T, model.S)
    model.z = Var(model.T, model.S, within=NonNegativeReals)
    model.u = Var(model.T, model.S, within=NonNegativeReals)
    model.v = Var(model.T, model.S, within=NonNegativeReals)

    # Constraints
    model.zdot_def = Constraint(model.T, model.S, rule=_zdot_const)
    model.xdot_def = Constraint(model.T, model.S, rule=_xdot_const)
    model.w_def = Constraint(model.T, rule=_w_const)

    # Define the objective
    model.obj = Objective(rule=_planner_obj, sense=maximize)

    # Initial and terminal conditions
    for s in model.S:
        model.z[:, s].setub(model.zbar[s])
        model.x[min(model.T), s].fix(model.x0[s])
        model.z[min(model.T), s].fix(model.z0[s])
        model.u[max(model.T), s].fix(0)
        model.v[max(model.T), s].fix(0)
        model.w[max(model.T)].fix(0)

    # Solve the model
    solver = SolverFactory("gurobi")

    print("Solving the optimization problem...")
    start_time = time.time()
    solver.solve(model, tee=True)
    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")

    Z = np.array([[model.z[t, r].value for r in model.S] for t in model.T])
    X = np.array([[model.x[t, r].value for r in model.S] for t in model.T])
    U = np.array([[model.u[t, r].value for r in model.S] for t in model.T])
    V = np.array([[model.v[t, r].value for r in model.S] for t in model.T])
    w = np.array([model.w[t].value for t in model.T])

    X_agg = X.sum(axis=1)
    X_agg = X_agg.reshape(X_agg.size, 1)

    sol_val_Ua = (w[:-1] ** 2).T.flatten()
    sol_val_X = np.concatenate((Z.T, X_agg.T, np.ones((1, Z.T.shape[1]))))
    sol_val_Up = U[:-1, :].T
    sol_val_Um = V[:-1, :].T
    sol_val_Z = sol_val_Up - sol_val_Um
    return (sol_val_X, sol_val_Up, sol_val_Um, sol_val_Z, sol_val_Ua)


def _planner_obj(model):
    return sum(
        model.pa_prob[p]
        * sum(
            math.exp(-model.delta * (t * model.dt - model.dt))
            * (
                -model.pe
                * pyo.quicksum(
                    model.kappa * model.z[t, s]
                    - (model.x[t + 1, s] - model.x[t, s]) / model.dt
                    for s in model.S
                )
                + model.pa[p, t] * sum(model.theta[s] * model.z[t, s] for s in model.S)
                - model.zeta / 2 * (model.w[t] ** 2)
            )
            * model.dt
            for t in model.T
            if t < max(model.T)
        )
        for p in model.P
    )


def _zdot_const(model, t, s):
    if t < max(model.T):
        return (model.z[t + 1, s] - model.z[t, s]) / model.dt == (
            model.u[t, s] - model.v[t, s]
        )
    else:
        return Constraint.Skip


def _xdot_const(model, t, s):
    if t < max(model.T):
        return (model.x[t + 1, s] - model.x[t, s]) / model.dt == (
            -model.gamma[s] * model.u[t, s]
            - model.alpha * model.x[t, s]
            + model.alpha * model.gamma[s] * (model.zbar[s] - model.z[t, s])
        )
    else:
        return Constraint.Skip


def _w_const(model, t):
    if t < max(model.T):
        return model.w[t] == pyo.quicksum(
            model.u[t, s] + model.v[t, s] for s in model.S
        )
    else:
        return Constraint.Skip


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
    print(binary_paths)
    # Find the probability of each path
    return [_price_path_prob(M, tau, path) for path in binary_paths]


def _price_path_prob(M, tau, price_path):
    prob = 1.0

    for t in range(tau):
        current_state = price_path[t]
        next_state = price_path[t + 1]
        prob *= M[current_state, next_state]

    return prob
