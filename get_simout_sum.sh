#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="get_simout"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1-36
#SBATCH --partition=all


module load R/3.6.0-foss-2019a


## file with three parameters per line, one line per array id:
## simparam    nrep    fragparam
## (fragparam can be NA or missing for runs without genome fragmentation)
params=$HOME/Rearrangements/Simulations/simparams.txt

mysim=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
nrep=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
myfrag=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')

simoutdir=$HOME/Rearrangements/Simulations/simout
scriptdir=$HOME/Rearrangements/Simulations/R

cd ${simoutdir}/${mysim}



## run R script to summarize replicates (saved as .RData)

if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
    ## without genome fragmentation
    Rscript --vanilla ${scriptdir}/get_simout_sum.R \
            ${mysim} ${nrep}
else
    ## with genome fragmentation
    Rscript --vanilla ${scriptdir}/get_simout_sum.R \
            ${mysim} ${nrep} ${myfrag}
fi


echo "-----------------------------------"
echo "Finished summarizing ${mysim} ${myfrag}"
echo "-----------------------------------"


exit
