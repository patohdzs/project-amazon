from pysrc.services.get_opt import get_optimization
from pysrc.services.get_sample import get_sampling
from pysrc.analysis.figures import land_allocation,density,trajectory_diff
from pysrc.analysis.tables import value_decom,transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation



# density(num_sites=1043,pee=5.3,xi=10.0,solver="gams")

density(num_sites=1043,pee=4.5,xi=5.0,solver="gams")


# density(num_sites=1043,pee=2.2,xi=1.0,solver="gams")
