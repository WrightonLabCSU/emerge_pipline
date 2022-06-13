# DRAM product refined pipeline

## Description

* Input: DRAM annotations file.
* Output: distillate (gene presence/absence), product (reaction % presence) and product refined (pathway calls).

Each pathway (e.g. Calvin-Benson carbon fixation) is defined as a combination of reactions (e.g. RubisCO).
For a genome to have a positive pathway call, it must have >60% of the total pathway (accounting for
alternate reaction paths), it must have some presence in >70% of the reactions (at least one gene) and,
it must have the signature genes, where applicable (used to distinguish pathways with high overlap, e.g.
the Wood-Ljungdahl methanogenic pathway and hydrogenotrophic methanogenesis).

## How to Run

* Create conda environment: `conda env create --name DRAM_product_refined -f environment.yaml`.
* Activate environment: `conda activate DRAM_product_refined`.
* Add path to DRAM annotations to the config in `config/config.yaml`.
* Run pipeline: `snakemake -j 30`.
