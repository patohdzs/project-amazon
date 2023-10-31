import math
import time

import numpy as np
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    NonNegativeReals,
    Objective,
    Param,
    Set,
    Var,
    maximize,
)
from pyomo.opt import SolverFactory


def solve_planner_problem(
    T,
    theta,
    gamma,
    x0_vals,
    zbar_2017,
    z_2017,
    pe=20.76,
    pa=44.75,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    zeta=1.66e-4 * 1e9,
    dt=1,
):
    model = ConcreteModel()

    # Indexing sets for time and sites
    model.T = Set(initialize=list(range(T + 1)), ordered=True)
    model.S = Set(initialize=list(range(gamma.size)), ordered=True)

    # Parameters
    model.x0 = Param(model.S, initialize=x0_vals)
    model.z0 = Param(model.S, initialize=z_2017)
    model.zbar = Param(model.S, initialize=zbar_2017)
    model.gamma = Param(model.S, initialize=gamma)
    model.theta = Param(model.S, initialize=theta)
    model.delta = Param(initialize=delta)
    model.p_e = Param(initialize=pe)
    model.p_a = Param(initialize=pa)
    model.alpha = Param(initialize=alpha)
    model.kappa = Param(initialize=kappa)
    model.zeta = Param(initialize=zeta)
    model.dt = Param(initialize=dt)

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
    return (Z, X, U, V, w)


def _planner_obj(model):
    return sum(
        math.exp(-model.delta * t)
        * (
            -model.p_e
            * pyo.quicksum(
                model.kappa * model.z[t, s] - (model.x[t + 1, s] - model.x[t, s])
                for s in model.S
            )
            + model.p_a * sum(model.theta[s] * model.z[t, s] for s in model.S)
            - model.zeta / 2 * (model.w[t] ** 2)
        )
        for t in model.T
        if t < max(model.T)
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
