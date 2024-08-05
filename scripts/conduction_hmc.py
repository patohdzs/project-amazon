from pysrc.services.get_opt import get_optimization
from pysrc.services.get_sample import get_sampling,get_prior
from pysrc.analysis.figures import land_allocation,density,trajectory_diff
from pysrc.analysis.tables import value_decom,transfer_cost,ambiguity_decom
from pysrc.analysis.map import spatial_allocation



### Section 7.3 Results with robustness to parameter uncertainty

get_optimization(num_sites=78,pee=7.1,model="det")
get_sampling(num_sites=78,pee=5.3,xi=1)
get_optimization(num_sites=78,pee=5.3,model="hmc")
get_prior(num_sites=78)
density(num_sites=78,pee=5.3,xi=1)
get_sampling(num_sites=78,pee=7.1,xi=1)
get_optimization(num_sites=78,pee=7.1,model="hmc")

spatial_allocation(pe_hmc=7.1,pe_det=7.1,b=0) # Figure 10
spatial_allocation(pe_hmc=5.3,pe_det=7.1,b=0) # Figure 12
spatial_allocation(pe_hmc=5.3,pe_det=7.1,b=15) # Figure 13
trajectory_diff(pe_hmc=7.1,pe_det=7.1,b=0,opt="gams",pa=41.11) # Figure 11
trajectory_diff(pe_hmc=5.3,pe_det=7.1,b=0,opt="gams",pa=41.11) # Figure 11
trajectory_diff(pe_hmc=5.3,pe_det=7.1,b=15,opt="gams",pa=41.11) # Figure 11

value_decom(num_sites=78,pee=5.3,model='hmc')
ambiguity_decom(num_sites=78,pe_det=7.1,pe_hmc=5.3) # Table 5


print("hmc All done!")