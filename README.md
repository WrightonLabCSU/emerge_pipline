# emerge_pipline
The pipline for the folks at emerge, they will want it eventually. I know this README is short, and it would be nice to fix that, but this Snake make pipeline is mostly self explanatory. 

## How to Run

 * It will help if you know Snakemake, so look up the "Read The Docks" for Snakemake, download and install it.
 * A conda environment file is included for requirements, do something like `conda env create --name emerge -f environment.yaml` and then `conda activate emerge` to hopefully get all the python dependency's.
 * Open the config in `config/config.yaml` and add in the paths To your DRAM annotations you want to distill.
 * Run something like `snakemake -j 30` and you are good to go.

Sorry this pipeline is not setup for slurm but it should be fast anyway.

## TODO
 * Finish the read me
 * add other things to the todo list
