# Snakemake file

TREATMENTS = ["T1", "T2", "T3", "T4", "T5", "T6"]

rule all:
  input:
    expand("results/case_study_sims/{treatment}/{treatment}_density_plot.pdf", treatment = TREATMENTS),
    expand("results/case_study_sims/{treatment}/{treatment}_ratio_plot.pdf", treatment = TREATMENTS)

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
    
# Define rule for processing and plotting case study invasion simulations
rule plot_case_study_sims:
  input:
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_density_plot.pdf",
    "results/case_study_sims/{treatment}/{treatment}_ratio_plot.pdf",
    "results/case_study_sims/{treatment}/{treatment}_density_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_ratio_plot_df.csv"
  shell:
    "Rscript src/process_simrun.R --treatment {wildcards.treatment}"
