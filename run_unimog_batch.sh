#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="run_unimog"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:50:00
#SBATCH --mem-per-cpu=6G
#SBATCH --array=1-36
#SBATCH --partition=all

module load Java/11.0.2

## file with three parameters per line, one line per array id:
## simparam    nrep    fragparam
## (fragparam can be NA or missing for runs without genome fragmentation)
params=$HOME/Rearrangements/Simulations/simparams.txt

mysim=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
nrep=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
myfrag=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')


grimmoutdir=$HOME/Rearrangements/Simulations/GRIMM
unimogoutdir=$HOME/Rearrangements/Simulations/UniMoG

mkdir -p ${unimogoutdir}/${mysim}
cd ${unimogoutdir}/${mysim}

## gop through replicates
for (( i=1; i<=${nrep}; i++ ))
do
    ## get correct genome file
    if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
	myrun=${mysim}_${i}
    else
	myrun=${mysim}_${i}_${myfrag}
    fi

    ## make UniMoG input file
    cat ${grimmoutdir}/${mysim}/${myrun}_grimm.txt \
	| sed 's/\$/|/g' >${myrun}_unimog.txt


    # ## run UniMoG (uniform sampling DCJ)
    # java -jar $HOME/bin/UniMoG/update/UniMoG.jar \
    # 	 -m=1 -s ${myrun}_unimog.txt -p >${myrun}_p.out

    ## run UniMoG (standard DCJ)
    java -jar $HOME/bin/UniMoG/update/UniMoG.jar \
	 -m=1 ${myrun}_unimog.txt -p >${myrun}_p_st.out


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
