# rearrvisr_simulations

Steps and scripts to simulate rearrangements for assessing the performance of *rearrvisr*


* Requires R package [`rearrvisr`](https://github.com/dorolin/rearrvisr) to be installed.
* Requires directory `~/Rearrangements/Simulations/` and R scripts in `~/Rearrangements/Simulations/R/`
* Requires file `simparams.txt` in `~/Rearrangements/Simulations/`
* Bash scripts were prepared for running on a Linux cluster with a [slurm](https://slurm.schedmd.com/documentation.html) batch-queueing system
* The calculations below were performed on [UBELIX](http://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern

## Make simulation settings

    ```bash
    cd ~/Rearrangements/Simulations/
    mkdir simparams
    Rscript --vanilla R/makeSimParams.R
    Rscript --vanilla R/makeFragParams.R 
    ```

## Run simulations and *rearrvisr*

* Run `runit_batch.sh` for first set of simulations (`#SBATCH --array=1-4`)
* After finishing, continue with remaining simulations (set `#SBATCH --array=5-36`)
* After all runs finished, summarize results by running `get_simout_sum.sh`
* Optionally, plot results with `Rscript --vanilla R/summarizePrecRec.R`


## Optionally, run the software [GRIMM](http://grimm.ucsd.edu/GRIMM/index.html) on the same set of simulations

* Requires *rearrvisr* runs above to be finished
* Run `run_grimm_batch.sh`
* After all runs finished, summarize results by running `get_grimm_sum.sh`


## Optionally, run the software [UniMoG](https://bibiserv.cebitec.uni-bielefeld.de/dcj) on the same set of simulations

* Requires GRIMM runs above to be finished
* Run `run_unimog_batch.sh`
* After all runs finished, summarize results by running `get_unimog_sum.sh`


