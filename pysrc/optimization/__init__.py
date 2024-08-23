import math
import time
from dataclasses import dataclass

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


@dataclass
class PlannerSolution:
    Z: np.ndarray
    X: np.ndarray
    U: np.ndarray
    V: np.ndarray
    w: np.ndarray


def solve_planner_problem(
    x0,
    z0,
    zbar,
    gamma,
    theta,
    low_pq,
    dt=1,
    time_horizon=200,
    price_emissions=20.76,
    price_cattle=44.75,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    zeta_u=1.66e-4 * 1e9,
    zeta_v_1=1.66e-4 * 1e9,
    zeta_v_2=1.66e-4 * 1e9,
    solver="gurobi",
):
    model = ConcreteModel()

    # Indexing sets for time and sites
    model.T = RangeSet(time_horizon + 1)
    model.S = RangeSet(gamma.size)

    # Parameters
    model.x0 = Param(model.S, initialize=_np_to_dict(x0))
    model.z0 = Param(model.S, initialize=_np_to_dict(z0))
    model.zbar = Param(model.S, initialize=_np_to_dict(zbar))
    model.gamma = Param(model.S, initialize=_np_to_dict(gamma))
    model.theta = Param(model.S, initialize=_np_to_dict(theta))
    model.delta = Param(initialize=delta)
    model.pe = Param(initialize=price_emissions)

    # Set cattle price as series
    if isinstance(price_cattle, float):
        price_cattle = {t + 1: price_cattle for t in range(time_horizon)}
    else:
        price_cattle = {t + 1: price_cattle[t] for t in range(time_horizon)}

    model.pa = Param(model.T, initialize=price_cattle)

    model.low_pq = Param(model.S, initialize=_np_to_dict(low_pq))
    model.high_pq = Param(model.S, initialize=_np_to_dict(1 - low_pq))

    # Asymmetric adj. costs
    model.zeta_u = Param(initialize=zeta_u)
    model.zeta_v_low = Param(initialize=zeta_v_1)
    model.zeta_v_high = Param(initialize=zeta_v_2)

    model.alpha = Param(initialize=alpha)
    model.kappa = Param(initialize=kappa)
    model.dt = Param(initialize=dt)

    # Variables
    model.x = Var(model.T, model.S)
    model.z = Var(model.T, model.S, within=NonNegativeReals)
    model.u = Var(model.T, model.S, within=NonNegativeReals)
    model.v = Var(model.T, model.S, within=NonNegativeReals)

    # Constraints
    model.zdot_def = Constraint(model.T, model.S, rule=_zdot_const)
    model.xdot_def = Constraint(model.T, model.S, rule=_xdot_const)

    # Define the objective
    model.obj = Objective(rule=_planner_obj, sense=maximize)

    # Initial and terminal conditions
    for s in model.S:
        model.z[:, s].setub(model.zbar[s])
        model.x[min(model.T), s].fix(model.x0[s])
        model.z[min(model.T), s].fix(model.z0[s])
        model.u[max(model.T), s].fix(0)
        model.v[max(model.T), s].fix(0)

    # Solve the model
    opt = SolverFactory(solver)
    print("Solving the optimization problem...")
    start_time = time.time()
    if solver == "gams":
        opt.solve(model, tee=True, solver="cplex", mtype="qcp")
    else:
        opt.solve(model, tee=True)

    print(f"Done! Time elapsed: {time.time()-start_time} seconds.")

    Z = np.array([[model.z[t, r].value for r in model.S] for t in model.T])
    X = np.array([[model.x[t, r].value for r in model.S] for t in model.T])
    U = np.array([[model.u[t, r].value for r in model.S] for t in model.T])
    V = np.array([[model.v[t, r].value for r in model.S] for t in model.T])
    w = np.array([model.w[t].value for t in model.T])

    return PlannerSolution(Z, X, U, V, w)


def vectorize_trajectories(traj: PlannerSolution):
    X_agg = traj.X.sum(axis=1)
    X_agg = X_agg.reshape(X_agg.size, 1)

    sol_val_Ua = (traj.w[:-1] ** 2).T.flatten()
    sol_val_X = np.concatenate((traj.Z.T, X_agg.T, np.ones((1, traj.Z.T.shape[1]))))
    sol_val_Up = traj.U[:-1, :].T
    sol_val_Um = traj.V[:-1, :].T
    sol_val_Z = sol_val_Up - sol_val_Um
    return {
        "sol_val_X": sol_val_X,
        "sol_val_Up": sol_val_Up,
        "sol_val_Um": sol_val_Um,
        "sol_val_Z": sol_val_Z,
        "sol_val_Ua": sol_val_Ua,
    }


def _planner_obj(model):
    return sum(
        math.exp(-model.delta * (t * model.dt - model.dt))
        * (
            -model.pe
            * pyo.quicksum(
                model.kappa * model.z[t + 1, s]
                - (model.x[t + 1, s] - model.x[t, s]) / model.dt
                for s in model.S
            )
            + model.pa[t]
            * pyo.quicksum(model.theta[s] * model.z[t + 1, s] for s in model.S)
            - (model.zeta_u / 2) * (pyo.quicksum(model.u[t, s] for s in model.S) ** 2)
            - (model.zeta_v_low / 2)
            * (pyo.quicksum(model.low_pq[s] * model.v[t, s] for s in model.S) ** 2)
            - (model.zeta_v_high / 2)
            * (pyo.quicksum(model.high_pq[s] * model.v[t, s] for s in model.S) ** 2)
        )
        * model.dt
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


def _np_to_dict(x):
    return dict(enumerate(x.flatten(), 1))
