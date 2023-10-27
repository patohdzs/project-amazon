import math
import time

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


def solve_outer_optimization_problem(
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

    model.T = Set(initialize=list(range(T + 1)), ordered=True)
    model.R = Set(initialize=list(range(gamma.size)), ordered=True)

    # Parameters
    model.x0 = Param(model.R, initialize=x0_vals)
    model.z0 = Param(model.R, initialize=z_2017)
    model.zbar = Param(model.R, initialize=zbar_2017)
    model.gamma = Param(model.R, initialize=gamma)
    model.theta = Param(model.R, initialize=theta)
    model.delta = Param(initialize=delta)
    model.p_e = Param(initialize=pe)
    model.p_a = Param(initialize=pa)
    model.alpha = Param(initialize=alpha)
    model.kappa = Param(initialize=kappa)
    model.zeta = Param(initialize=zeta)
    model.dt = Param(initialize=dt)

    # Variables
    model.w = Var(model.T)
    model.x = Var(model.T, model.R)
    model.z = Var(model.T, model.R, within=NonNegativeReals)
    model.u = Var(model.T, model.R, within=NonNegativeReals)
    model.v = Var(model.T, model.R, within=NonNegativeReals)

    # Constraints
    model.zdot_def = Constraint(model.T, model.R, rule=_zdot_rule)
    model.xdot_def = Constraint(model.T, model.R, rule=_xdot_rule)
    model.w_def = Constraint(model.T, rule=_w_rule)

    # Define the objective
    model.obj = Objective(expr=_objective(model), sense=maximize)

    # Initial and terminal conditions
    for r in model.R:
        model.z[:, r].setub(model.zbar[r])
        model.x[0, r].fix(model.x0[r])
        model.z[0, r].fix(model.z0[r])
        model.u[T, r].fix(0)
        model.v[T, r].fix(0)
        model.w[T].fix(0)

    # Solve the model
    solver = SolverFactory("gurobi")

    print("Solving the optimization problem...")
    start_time = time.time()
    solver.solve(model)
    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")


def _objective(model):
    return sum(
        math.exp(-model.delta * t)
        * (
            -model.p_e
            * sum(
                model.kappa * model.z[t, r] - (model.x[t + 1, r] - model.x[t, r])
                for r in model.R
            )
            + model.p_a * sum(model.theta[r] * model.z[t, r] for r in model.R)
            - model.zeta / 2 * (model.w[t] ** 2)
        )
        for t in model.T
        if t < max(model.T)
    )


def _zdot_rule(model, t, r):
    if t < max(model.T):
        return (model.z[t + 1, r] - model.z[t, r]) / model.dt == (
            model.u[t, r] - model.v[t, r]
        )
    else:
        return Constraint.Skip


def _xdot_rule(model, t, r):
    if t < max(model.T):
        return (model.x[t + 1, r] - model.x[t, r]) / model.dt == (
            -model.gamma[r] * model.u[t, r]
            - model.alpha * model.x[t, r]
            + model.alpha * model.gamma[r] * (model.zbar[r] - model.z[t, r])
        )
    else:
        return Constraint.Skip


def _w_rule(model, t):
    if t < max(model.T):
        return model.w[t] == sum(model.u[t, r] + model.v[t, r] for r in model.R)
    else:
        return Constraint.Skip
