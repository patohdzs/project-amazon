import argparse
import pickle

from pysrc.sampling import adjusted
from pysrc.services.file_service import logs_dir_path, output_dir_path, plots_dir_path

# Read arguments from stdin
parser = argparse.ArgumentParser(description="parameter settings")

parser.add_argument("--model", type=str, default="sigma_correction_model")
parser.add_argument("--opt", type=str, default="gurobi")
parser.add_argument("--xi", type=float, default=2.0)
parser.add_argument("--pf", type=float, default=21.5)
parser.add_argument("--pa", type=float, default=42.03)
parser.add_argument("--weight", type=float, default=0.25)
parser.add_argument("--sitenum", type=int, default=78)
parser.add_argument("--timehzn", type=int, default=200)

# Parse arguments
args = parser.parse_args()

# Create output and plots directories
output_dir = output_dir_path(**vars(args))
plots_dir = plots_dir_path(**vars(args))
logs_dir = logs_dir_path(**vars(args))

# Sample from adjusted distribution
results = adjusted.sample(
    model_name=args.model,
    xi=args.xi,
    pe=args.pf,
    pa=args.pa,
    weight=args.weight,
    num_sites=args.sitenum,
    T=args.timehzn,
    optimizer=args.opt,
    max_iter=100,
    final_sample_size=5_000,
    iter_sampling=1000,
    iter_warmup=500,
    show_progress=True,
    seed=1,
    inits=0.2,
)

# Save results
outfile_path = output_dir / "results.pcl"
with open(outfile_path, "wb") as outfile:
    pickle.dump(results, outfile)
    print(f"Results saved to {outfile_path}")
