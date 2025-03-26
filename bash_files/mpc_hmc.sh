
pe=30.0
# idarray=($(seq 51 100))
idarray=($(seq 1 50))
# idarray=($(seq 997 998))
xi=0.2
trig=1
type="constrained"

for id in "${idarray[@]}"; do
   
                            count=0
                                        
                            action_name="mpc"

                            dataname="${action_name}"

                            mkdir -p ./job-outs/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}

                            if [ -f ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh ]; then
                                rm ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh
                            fi

                            mkdir -p ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}

                            touch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh

                            tee -a ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=id_${id}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake	
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load python/anaconda-2022.05  
module load gurobi/11.0

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 -u /project/lhansen/HMC_rp/project-amazon/pysrc/mpc/mpc_hmc.py --id ${id} --pe ${pe} --xi ${xi} --trig ${trig} --type ${type}
echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
    count=$(($count + 1))
    sbatch ./bash/${action_name}/xi_${xi}/pe_${pe}/id_${id}/trig_${trig}/type_${type}/run.sh


done
