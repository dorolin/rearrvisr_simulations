#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="get_grimm"
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

grimmoutdir=$HOME/Rearrangements/Simulations/GRIMM

cd ${grimmoutdir}/${mysim}


events=events_${SLURM_ARRAY_TASK_ID}.txt
tmpevents=tmpevents_${SLURM_ARRAY_TASK_ID}.txt

## first, check out which rearrangement events exist
##  (and make outfile name on the way)

for (( i=1; i<=${nrep}; i++ ))
do

    if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
	myrun=${mysim}_${i}
	myout=${mysim}_summary.txt
    else
	myrun=${mysim}_${i}_${myfrag}
	myout=${mysim}_${myfrag}_summary.txt
    fi

    grep "^Step" ${myrun}_run1.txt | grep -v "(Source)" | cut -d ":" -f3 | sed -E 's/\s+//g' | sed 's/(Destination)//' | sort | uniq >>${tmpevents}

done


cat ${tmpevents} | sort | uniq > ${events}
rm -f ${tmpevents}

## prepare output

neve=$(wc -l ${events} | cut -d ' ' -f1)

myeves=$(cat ${events} | tr '\n' ' ')

line="Run Int_bpts Ext_bpts Distance ${myeves}"
echo ${line} > ${myout}


## go through each replicate

for (( i=1; i<=${nrep}; i++ ))
do

    if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
	myrun=${mysim}_${i}
    else
	myrun=${mysim}_${i}_${myfrag}
    fi

    ## number of breakpoints
    ib=$(grep "Number of internal breakpoints" ${myrun}_run1.txt | cut -d ":" -f2 | sed -E 's/\s+//g')
    eb=$(grep "Number of external breakpoints" ${myrun}_run1.txt | cut -d ":" -f2 | sed -E 's/\s+//g')
    ## distance (i.e., total events)
    md=$(grep "Multichromosomal Distance" ${myrun}_run1.txt | cut -d ":" -f2 | sed -E 's/\s+//g')

    ## number of events per event type
    myeves=""
    for (( j=1; j<=${neve}; j++  ))
    do
	myeve=$(cat $events | awk -v var=$j 'NR==var {print $1}')
	me=$(grep "^Step" ${myrun}_run1.txt | grep -v "(Source)" | cut -d ":" -f3 | sed -E 's/\s+//g' | sed 's/(Destination)//' | grep -c ${myeve})
	myeves="${myeves} ${me}"
    done

    ## print out
    line="${myrun} ${ib} ${eb} ${md} ${myeves}"
    echo ${line} >> ${myout}

    ## clean up
    unset myrun

done

rm -f ${events}

echo "-----------------------------------"
echo "Finished writing ${myout}"
echo "-----------------------------------"


exit
