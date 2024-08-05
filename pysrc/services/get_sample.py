import os
import pickle

from pysrc.sampling import adjusted, baseline
from pysrc.services.file_service import get_path


def get_sampling(opt="gams", num_sites=78, pa=41.11, pee=5, xi=1):
    b = [0, 10, 15, 20, 25]
    pe_values = [pee + bi for bi in b]
    for pe in pe_values:
        results = adjusted.sample(
            xi=xi,
            pe=pe,
            pa=pa,
            weight=0.25,
            num_sites=num_sites,
            T=200,
            optimizer=opt,
            max_iter=100,
            final_sample_size=5_000,
            iter_sampling=1000,
            iter_warmup=500,
            show_progress=True,
            seed=1,
            inits=0.2,
        )
        output_base_path = os.path.join(
            str(get_path("output")),
            "sampling",
            opt,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"xi_{xi}",
            f"pe_{pe}",
        )
        if not os.path.exists(output_base_path):
            os.makedirs(output_base_path)
        outfile_path = output_base_path + "/results.pcl"
        with open(outfile_path, "wb") as outfile:
            pickle.dump(results, outfile)
            print(f"Results saved to {outfile_path}")

    return


def get_prior(num_sites=78):
    baseline_fit = baseline.sample(
        num_sites=num_sites, iter_sampling=10**4, chains=5, seed=1
    )

    theta = baseline_fit.stan_variable("theta")
    gamma = baseline_fit.stan_variable("gamma")

    results_dir = os.path.join(
        str(get_path("output")), "sampling", "prior", f"{num_sites}sites"
    )
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    filename = "prior.pcl"
    file_path = os.path.join(results_dir, filename)

    data = {"theta": theta, "gamma": gamma}

    with open(file_path, "wb") as f:
        pickle.dump(data, f)

    print(f"Results (theta and gamma) saved to {file_path}")

    return
