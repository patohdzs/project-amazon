import numpy as np

from pysrc.optimization import gams, gurobi
from pysrc.sampling import adjusted
from pysrc.services.data_service import load_site_data_1995


def shadow_price_opt(
    zbar_1995,
    z_1995,
    forest_area_1995,
    z_2008,
    theta,
    gamma,
    sitenum=78,
    opt="gams",
    timehzn=200,
    pa=41.11,
    pe=7.1,
    model="det",
):
    # Computing carbon absorbed in start period
    x0_vals_1995 = gamma * forest_area_1995

    # Choose optimizer
    if opt == "gurobi":
        solve_planner_problem = gurobi.solve_planner_problem

    elif opt == "gams":
        solve_planner_problem = gams.solve_planner_problem

    else:
        raise ValueError("Optimizer must be one of ['gurobi', 'gams']")

    if model == "mpc":
        solve_planner_problem = gams.mpc_shadow_price

    results = solve_planner_problem(
        T=timehzn,
        theta=theta,
        gamma=gamma,
        x0=x0_vals_1995,
        zbar=zbar_1995,
        z0=z_1995,
        pe=pe,
        pa=pa,
        model="shadow_price",
    )
    Z = results["Z"]
    z_2008_agg = np.sum(z_2008) / 1e9
    ratio = np.abs((np.sum(Z[13]) - z_2008_agg) / z_2008_agg)

    return ratio


def shadow_price_cal(sitenum=78, pa=41.11, opt="gams", model="det", xi=1):
    if model == "det":
        (
            zbar_1995,
            z_1995,
            forest_area_1995,
            _,
            _,
            _,
            _,
            z_2008,
            theta,
            gamma,
        ) = load_site_data_1995(sitenum)

        pe_values = np.arange(5, 8, 0.1)

        results = np.array(
            [
                shadow_price_opt(
                    zbar_1995,
                    z_1995,
                    forest_area_1995,
                    z_2008,
                    theta,
                    gamma,
                    sitenum=sitenum,
                    opt=opt,
                    timehzn=200,
                    pa=pa,
                    pe=pe,
                )
                for pe in pe_values
            ]
        )

    elif model == "hmc":
        (
            zbar_1995,
            z_1995,
            forest_area_1995,
            _,
            _,
            _,
            _,
            z_2008,
            theta,
            gamma,
        ) = load_site_data_1995(sitenum)

        pe_values = np.arange(5, 7, 0.1)
        results = []
        for pe in pe_values:
            samples = adjusted.sample(
                xi=xi,
                pe=pe,
                pa=pa,
                weight=0.25,
                num_sites=sitenum,
                T=200,
                solver=opt,
                max_iter=100,
                final_sample_size=5_000,
                iter_sampling=1000,
                iter_warmup=500,
                show_progress=True,
                seed=1,
                inits=0.2,
            )

            theta = np.mean(samples["final_sample"][:, :78], axis=0)
            gamma = np.mean(samples["final_sample"][:, 78:], axis=0)
            result = shadow_price_opt(
                zbar_1995,
                z_1995,
                forest_area_1995,
                z_2008,
                theta,
                gamma,
                sitenum=sitenum,
                opt=opt,
                timehzn=200,
                pa=pa,
                pe=pe,
            )
            results.append(result)
        results = np.array(results)

    if model == "mpc":
        (
            zbar_1995,
            z_1995,
            forest_area_1995,
            _,
            _,
            _,
            _,
            z_2008,
            theta,
            gamma,
        ) = load_site_data_1995(sitenum)

        pe_values = np.arange(6.5, 7.1, 0.05)

        results = []
        for pe in pe_values:
            result = shadow_price_opt(
                zbar_1995,
                z_1995,
                forest_area_1995,
                z_2008,
                theta,
                gamma,
                sitenum=sitenum,
                opt=opt,
                timehzn=200,
                pa=pa,
                pe=pe,
                model="mpc",
            )
            results.append(result)
            print(f"pe: {pe}, result: {result}")

        results = np.array(results)

    min_index = np.argmin(results)
    min_result = results[min_index]
    min_pe = pe_values[min_index]

    return min_result, min_pe


## det 1043 sites
min_result, det_1043_pe = shadow_price_cal(sitenum=1043, model="det")
print("min_result", min_result, "min_pe", det_1043_pe)
min_result, det_78_pe = shadow_price_cal(sitenum=78, model="det")
print("min_result", min_result, "min_pe", det_78_pe)


## hmc 78 sites
min_result, hmc_78_pe = shadow_price_cal(sitenum=78, model="hmc", xi=1)
print("min_result", min_result, "min_pe", hmc_78_pe)

## mpc 78 sites
min_result, mpc_78_pe = shadow_price_cal(sitenum=78, model="mpc")
print("min_result", min_result, "min_pe", mpc_78_pe)
