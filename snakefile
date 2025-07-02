# Snakemake file

# Snakemake Configuration
import json

## Specify treatments
TREATMENTS = ["A","B","C","D"]

## Specify global variables
### Colors
plot_colors = {
  "p_Anc": "#8394F6",
  "p_Mut": "#8A407A",
  "p_F": "gray40",
  "p_growth": "#FFFBEF",
  "p_conj": "#D0D9FF",
  "p_tselect": "#FBE9FF"
}

## Linetypes
plot_lines = {
  "rifR_l": "solid",
  "nalR_l": "dashed"
}


rule all:
  input:
    expand("results/case_study_sims/{treatment}/{treatment}_density_plot.pdf", treatment = TREATMENTS),
    expand("results/case_study_sims/{treatment}/{treatment}_frequency_plot.pdf", treatment = TREATMENTS)

# Define rule for running case study invasion simulations
rule case_study_sims:
  input:
    "src/simrun_case_study.R",
    "input_data/Treatments_case_study.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv",
    "results/case_study_sims/{treatment}/{treatment}_sim_plot.pdf"
  shell:
    "Rscript src/simrun_case_study.R --treatment {wildcards.treatment}"
    
# Define rule for processing case study invasion simulations
rule process_case_study_sims:
  input:
    "src/process_simrun.R",
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_density_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_phases_plot_df.csv"
  shell:
    "Rscript src/process_simrun.R --treatment {wildcards.treatment}"


# Define rule for plotting case study invasion simulations
rule plot_case_study_sims:
  input:
    "src/plot_simrun_snakemake.R",
    "results/case_study_sims/{treatment}/{treatment}_density_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_phases_plot_df.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_density_plot.pdf",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot.pdf"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_lines = json.dumps(plot_lines)
  shell:
    """
    Rscript src/plot_simrun_snakemake.R \
      --treatment {wildcards.treatment} \
      --colors '{params.plot_colors}' \
      --lines '{params.plot_lines}'
    """
