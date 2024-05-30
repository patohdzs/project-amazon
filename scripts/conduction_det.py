from pysrc.services.get_opt import get_optimization
from pysrc.services.get_sample import get_sampling,get_prior
from pysrc.analysis.figures import land_allocation,density,trajectory_diff
from pysrc.analysis.tables import value_decom,transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation

### Section 7.2 Results for case without stochasticity or ambiguity aversion


get_optimization(num_sites=1043,pee=7.6,model="det")
land_allocation(num_sites=1043) # Figure 5
value_decom(num_sites=1043,pee=7.6) # Table 2
transfer_cost(num_sites=1043,pee=7.6,y=30) # Table 3
transfer_cost(num_sites=1043,pee=7.6,y=15) # Table 3


get_optimization(num_sites=78,pee=7.1,model="det")
get_optimization(num_sites=78,pee=5.3,model="det")
value_decom(num_sites=78,pee=7.1) # Table 4

print("det All done!")