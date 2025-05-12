# Snakemake file

TREATMENTS = ["T1", "T2", "T3", "T4"]

rule all:
  input:
    expand("results/case_study_sims/{treatment}/{treatment}_data.csv", treatment = TREATMENTS),
    expand("results/case_study_sims/{treatment}/{treatment}_data_long.csv", treatment = TREATMENTS),
    expand("results/case_study_sims/{treatment}/{treatment}_plot.pdf", treatment = TREATMENTS)

# Define rule for running case study invasion simulations
rule case_study_sims:
  input:
    "src/simrun_case_study.R",
    "input_data/Treatments.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv",
    "results/case_study_sims/{treatment}/{treatment}_plot.pdf"
  shell:
    "Rscript src/simrun_case_study.R --treatment {wildcards.treatment}"
