#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="get_unimog"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1-36
#SBATCH --partition=all



## file with three parameters per line, one line per array id:
## simparam    nrep    fragparam
## (fragparam can be NA or missing for runs without genome fragmentation)
params=$HOME/Rearrangements/Simulations/simparams.txt

mysim=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
nrep=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
myfrag=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')

unimogoutdir=$HOME/Rearrangements/Simulations/UniMoG

cd ${unimogoutdir}/${mysim}


## make outfile name
if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
    # myout=${mysim}_summary.txt
    myout2=${mysim}_st_summary.txt
else
    # myout=${mysim}_${myfrag}_summary.txt
    myout2=${mysim}_${myfrag}_st_summary.txt
fi

## prepare output
line="Run Distance"
# echo ${line} > ${myout}
echo ${line} > ${myout2}


## go through each replicate

for (( i=1; i<=${nrep}; i++ ))
do

    if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
	myrun=${mysim}_${i}
    else
	myrun=${mysim}_${i}_${myfrag}
    fi

    # ## uniform sampling DCJ distance
    # ## -----------------------------
    # dcj=$(grep '"compgenome" & "focalgenome" (DCJ)' ${myrun}_p.out | cut -d ":" -f2 | sed -E 's/\s+//g')

    # ## as a few runs failed, assign NA instead
    # if [[ -z "${dcj// }" ]]; then
    # 	dcj=NA
    # fi

    # ## print out
    # line="${myrun} ${dcj}"
    # echo ${line} >> ${myout}

    ## standard DCJ distance
    ## -----------------------------
    dcj=$(grep '"compgenome" & "focalgenome" (DCJ)' ${myrun}_p_st.out | cut -d ":" -f2 | sed -E 's/\s+//g')

    ## if runs failed, assign NA instead (should not happen)
    if [[ -z "${dcj// }" ]]; then
	dcj=NA
    fi

    ## print out
    line="${myrun} ${dcj}"
    echo ${line} >> ${myout2}



    ## clean up
    unset myrun

done


echo "-----------------------------------"
echo "Finished writing ${myout2}"
echo "-----------------------------------"


exit
