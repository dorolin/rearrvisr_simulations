#!/bin/bash

#SBATCH --mail-user=email@host.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="runsims"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:50:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-4 ## 5-36 !! when doing simulations the first time, runs 1-4 need to finish before 13-24 start
#SBATCH --partition=all

module load R/3.6.0-foss-2019a

## file with three parameters per line, one line per array id:
## simparam    nrep    fragparam
## (fragparam can be NA or missing for runs without genome fragmentation)
params=$HOME/Rearrangements/Simulations/simparams.txt

mysim=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
nrep=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
myfrag=$(cat $params | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')


paramdir=$HOME/Rearrangements/Simulations/simparams
outdir=$HOME/Rearrangements/Simulations/simout
scriptdir=$HOME/Rearrangements/Simulations/R

mkdir -p ${outdir}/${mysim}


if [[ ! -f "${outdir}/${mysim}/${mysim}_${nrep}.RData" ]]; then
    ## run simulations
    echo "-----------------------------------"
    echo "Running ${nrep} simulations ${mysim}"
    echo "-----------------------------------"

    Rscript --vanilla ${scriptdir}/runSims.R \
	    ${paramdir}/${mysim}_params.RData \
	    ${outdir}/${mysim}/${mysim}.RData \
	    ${nrep}
    ## requires correct source location in runSims.R
fi


if [[ -z "${myfrag// }" || "${myfrag}" == "NA" ]]; then
    ## run rearrvisr without genome fragmentation
    echo "-----------------------------------"
    echo "Running rearrvisr for ${mysim}"
    echo "-----------------------------------"

    Rscript --vanilla ${scriptdir}/runRearrvisr.R \
	    ${outdir}/${mysim}/${mysim}.RData \
	    ${nrep}
    ## requires correct source and library location in runRearrvisr.R
else
    ## run rearrvisr with genome fragmentation
    echo "-----------------------------------"
    echo "Running rearrvisr for ${mysim} with ${myfrag}"
    echo "-----------------------------------"

    Rscript --vanilla ${scriptdir}/runRearrvisr.R \
	    ${outdir}/${mysim}/${mysim}.RData \
	    ${nrep} \
	    ${paramdir}/${myfrag}_params.RData
    ## requires correct source and library location in runRearrvisr.R
fi


echo "-----------------------------------"
echo "Finished ${nrep} simulations ${mysim} ${myfrag}"
echo "-----------------------------------"

exit



## to install and load rearrvisr:
## ------------------------------
## library(devtools)
## withr::with_libpaths(new = "~/R/site-library/", install(pkg = "~/Source/R/rearrvisr"))
## library("rearrvisr", lib.loc="~/R/site-library/")
## ------------------------------
