# PETRI Project Repository
This repository contains the raw data and code for running simulations and processing data for the PETRI (Plasmid Evolution of Transfer Rate Investigation) Project.
The corresponding manuscript is currently available on bioRxiv: [Molecular pet or parasite? Exploring selection for vertical and horizontal plasmid transfer] (https://doi.org/10.64898/2026.06.02.729688)

## Project description
This project aims to explore the environmental conditions that select for vertical and horizontal plasmid transfer. We use a dynamical model to run simulations of mutant plasmid invasion, and pair it with experimental data to validate simulations.

## Repo structure
Main folders:
*src: contains scripts for running simulations, processing data, and plotting data
*input_data: contains raw experimental data and simulation input parameters
*results: contains processed data and preliminary plots
*figures: contains final figures for the manuscript
*snakemake files:
  - snakefile: snakemake pipeline script
  - PETRI_env.yaml: conda environment for snakemake and R runtime
Others:
*.snakemake: snakemake usage folder, logs are not synced due to large file quantity
*.hyak: slurm logs from UW Hyak computing cluster use

## Snakemake pipeline usage
The entire analysis pipeline (raw data -> figures) is executed using snakemake. The pipeline is as follows:

To execute the pipeline:
```
conda env create -f PETRI_env.yaml
conda activate PETRI
snakemake --cores 7 --use-conda
```
