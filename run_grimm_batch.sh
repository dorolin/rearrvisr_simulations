#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="run_grimm"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:50:00
#SBATCH --mem-per-cpu=4G
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
grimmoutdir=$HOME/Rearrangements/Simulations/GRIMM
scriptdir=$HOME/Rearrangements/Simulations/R

mkdir -p ${grimmoutdir}/${mysim}
cd ${grimmoutdir}/${mysim}

## go through replicates
for (( i=1; i<=${nrep}; i++ ))
do
    ## get correct simulation file
    if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
	myrun=${mysim}_${i}
    else
	myrun=${mysim}_${i}_${myfrag}
    fi

    ## make GRIMM input file
    Rscript --vanilla ${scriptdir}/genomes2grimm.R \
	    ${simoutdir}/${mysim}/${myrun}.RData ${myrun}

    ## run GRIMM
    $HOME/bin/GRIMM-2.01/grimm \
	-f ${myrun}_grimm.txt \
	-o ${myrun}_run1.txt -g 1,2 -v \
	-d -c -z -s


    echo "-----------------------------------"
    echo "Finished ${myrun}"
    echo "-----------------------------------"

    ## clean up
    unset myrun

done

echo "-----------------------------------"
echo "Finished ${nrep} runs for ${mysim} ${myfrag}"
echo "-----------------------------------"

exit
