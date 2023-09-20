dataname="tests"
sitenum=24
xi=0.1
pf=25
pa=44.75
theta=1.0
gamma=1.0
T=200
weight=0.1
mix_in=2
mass_matrix_theta_scale=1.0
mass_matrix_gamma_scale=1.0
mass_matrix_weight=0.0
symplectic_integrator_num_steps=10
stepsize=0.05
scale=0.0
mode=1.0

python scripts/sampler.py \
    --dataname ${dataname} \
    --weight ${weight} \
    --xi ${xi} \
    --pf ${pf} \
    --pa ${pa} \
    --sitenum ${sitenum} \
    --time ${T} \
    --mix_in ${mix_in} \
    --mass_matrix_theta_scale ${mass_matrix_theta_scale} \
    --mass_matrix_gamma_scale ${mass_matrix_gamma_scale} \
    --mass_matrix_weight ${mass_matrix_weight} \
    --symplectic_integrator_num_steps ${symplectic_integrator_num_steps} \
    --stepsize ${stepsize} \
    --scale ${stepsize} \
    --mode ${mode}
