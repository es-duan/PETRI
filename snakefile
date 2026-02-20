# Snakemake file

# Snakemake Configuration
import json

## Specify treatments
#TREATMENTS = ["A.2","C.2","O","E","F"]
TREATMENTS = ["DG_pB10", "HFC_pB10", "LFC_pB10-A"]
PSWEEPS = ["pHFC_S.pB10","pLFC_S.pB10","pDim90_E.R1","pDim90_E.RP4","pHFCi_S.pB10","pLFCi_S.pB10","pHFC_S.pB10-A","pLFC_S.pB10-A"]


## Specify global variables
### Colors
plot_colors = {
  "p_Anc": "#8394F6",
  "p_Mut": "#8A407A",
  "p_F": "gray40",
  "p_growth": "#FFFBEF",
  "p_conj": "#D0D9FF",
  "p_tselect": "#FBE9FF",
  "p_imm": "#d0ffeb",
  "p_select": "#ffe9fc"
}

plot_colors_psweep = {
  "p_Exc": "gray95",
  "p_Dis": "#140433",
  "p_lowI": "#FFC2C2",
  "p_highI": "#C20000",
  "p_lowD": "#BABAFF",
  "p_highD": "#0000AB"
}

### Linetypes
plot_lines = {
  "rifR_l": "solid",
  "nalR_l": "dashed"
}


rule all:
  input:
    expand("results/case_study_sims/{treatment}/{treatment}_density_plot.pdf", treatment = TREATMENTS),
    expand("results/case_study_sims/{treatment}/{treatment}_frequency_plot.pdf", treatment = TREATMENTS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_change_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_change_strain_plot.pdf", psweep = PSWEEPS)

# Define rule for running case study invasion simulations
rule case_study_sims:
  input:
    "src/simrun_case_study.R",
    "input_data/case_study_sims/{treatment}_inv_settings.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv",
    "results/case_study_sims/{treatment}/{treatment}_sim_plot.pdf"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/simrun_case_study.R \
      --treatment {wildcards.treatment} \
      --colors '{params.plot_colors}'
    """
    
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
    "src/plot_simrun.R",
    "results/case_study_sims/{treatment}/{treatment}_density_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot_df.csv",
    "results/case_study_sims/{treatment}/{treatment}_phases_plot_df.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_density_plot.pdf",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot.pdf",
    "results/case_study_sims/{treatment}/{treatment}_density_plot.rds",
    "results/case_study_sims/{treatment}/{treatment}_frequency_plot.rds"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_lines = json.dumps(plot_lines)
  shell:
    """
    Rscript src/plot_simrun.R \
      --treatment {wildcards.treatment} \
      --colors '{params.plot_colors}' \
      --lines '{params.plot_lines}'
    """

# Define rule for generating parameter sets for sweeps
rule psweep_parameter:
  input:
    "src/generate_parameter_sets.R",
    "input_data/parameter_sweeps/{psweep}_psweep_settings.csv"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_params.csv",
    "results/parameter_sweeps/{psweep}/{psweep}_settings.rds"
  shell:
    """
    Rscript src/generate_parameter_sets.R \
      --psweepsetting {wildcards.psweep} 
    """

# Define rule for running parameter sweep simulations
rule psweep_sims:
  input:
    "src/simrun_psweep.R",
    "results/parameter_sweeps/{psweep}/{psweep}_params.csv",
    "results/parameter_sweeps/{psweep}/{psweep}_settings.rds"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_out.csv"
  threads:
    6
  shell:
    """
    Rscript src/simrun_psweep.R \
      --psweepsetting {wildcards.psweep} \
      --threads {threads}
    """

# Define rule for processing parameter sweep simulations
rule process_psweep:
  input:
    "src/process_psweep.R",
    "results/parameter_sweeps/{psweep}/{psweep}_out.csv",
    "results/parameter_sweeps/{psweep}/{psweep}_settings.rds"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_plot.csv"
  shell:
    """
    Rscript src/process_psweep.R \
      --psweepsetting {wildcards.psweep}
    """
    
# Define rule for plotting parameter sweeps
rule plot_psweep:
  input:
    "src/plot_psweep.R",
    "results/parameter_sweeps/{psweep}/{psweep}_plot.csv"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_change_plot.pdf",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.pdf",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.rds"
  params:
    plot_colors_psweep = json.dumps(plot_colors_psweep)
  shell:
    """
    Rscript src/plot_psweep.R \
      --psweepsetting {wildcards.psweep} \
      --colors '{params.plot_colors_psweep}'
    """
    
# Define rule for plotting parameter sweeps with strains
rule plot_psweep_strain:
  input:
    "input_data/strain_phenotypes.csv",
    "src/plot_psweep_strain.R",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.rds"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_strain_plot.pdf"
  shell:
    """
    Rscript src/plot_psweep_strain.R \
      --psweepsetting {wildcards.psweep}
    """
