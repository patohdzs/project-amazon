sitenumarray=(78)
xiarray=(2.0)
pfarray=(16.5)
paarray=(42.03)
timehznarray=(200)
weightarray=(0.25)


hmc_python_name="adjusted_sampler.py"


for sitenum in "${sitenumarray[@]}"; do
    for xi in "${xiarray[@]}"; do
        for pf in "${pfarray[@]}"; do
            for pa in "${paarray[@]}"; do
                for timehzn in "${timehznarray[@]}"; do
                    for weight in "${weightarray[@]}"; do

                        count=0
                                    
                        action_name="test4"
                        opt="gams"

                        dataname="${action_name}"

                        mkdir -p ./job-outs/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/
                        
                        if [ -f ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.sh ]; then
                            rm ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.sh
                        fi

                        mkdir -p ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/
                        touch ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.sh

                        tee -a ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=pi-lhansen
#SBATCH --job-name=sn_${sitenum}_xi_${xi}_w_${weight}_${action_name}
#SBATCH --output=./job-outs/$job_name/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.out
#SBATCH --error=./job-outs/$job_name/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.err
#SBATCH --time=1-11:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=21G


module load python/anaconda-2022.05

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_timehzn=\$(date +%s)


python3 -u /project/lhansen/HMC/project-amazon/scripts/$hmc_python_name  --pf ${pf} --pa ${pa} --timehzn ${timehzn} --sitenum ${sitenum} --xi ${xi} --weight ${weight} --opt ${opt} 
echo "Program ends \$(date)"
end_timehzn=\$(date +%s)
elapsed=\$((end_timehzn - start_timehzn))

eval "echo Elapsed timehzn: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"
EOF
                        count=$(($count + 1))
                        sbatch ./bash/${action_name}/pf_${pf}_pa_${pa}_timehzn_${timehzn}/sitenum_${sitenum}_xi_${xi}/weight_${weight}/run.sh
                    done
                done
            done
        done
    done
done