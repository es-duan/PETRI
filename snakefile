# Snakemake file

# Snakemake Configuration
import json

## Specify treatments
#TREATMENTS = ["A.2","C.2","O","E","F"]
TREATMENTS = ["E.R1"]
SETTINGS = ["test", "testmut", "DimR1", "DimR1finO"]
#PSWEEPS = ["pHFC.l_DimR1","pHFC.l_DimR1finO","pHFC.l_test"]
#PSWEEPS = ["pLFC.l_DimR1finO","pLFC.l_DimR1","pLFC.l_test"]
#PSWEEPS = ["pF_test","pF_DimR1finO","pF_DimR1","pE_test","pE_DimR1finO","pE_DimR1"]
PSWEEPS = ["pE_DimR1"]
#PSWEEPS = ["pHFC.l_DimR1","pHFC.l_DimR1finO","pHFC.l_test","pLFC.l_DimR1finO","pLFC.l_DimR1","pLFC.l_test","pF_test","pF_DimR1finO","pF_DimR1"]


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
  "p_lowI": "#d2c5eb",
  "p_highI": "#2b1457",
  "p_parI": "#b7bde8",
  "p_parD": "gray80",
  "p_lowD": "#a38cd1"
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
    expand("results/parameter_sets/{parameters}_sweep_param.csv", parameters = SETTINGS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_change_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_time_plot.pdf", psweep = PSWEEPS)

# Define rule for running case study invasion simulations
rule case_study_sims:
  input:
    "src/simrun_case_study_snakemake.R",
    "input_data/Treatments_case_study.csv"
  output:
    "results/case_study_sims/{treatment}/{treatment}_data.csv",
    "results/case_study_sims/{treatment}/{treatment}_data_long.csv",
    "results/case_study_sims/{treatment}/{treatment}_sim_plot.pdf"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/simrun_case_study_snakemake.R \
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

# Define rule for generating parameter sets for sweeps
rule psweep_parameter:
  input:
    "src/generate_parameter_sets.R",
    "input_data/Parameter_sweep_settings.csv"
  output:
    "results/parameter_sets/{parameters}_sweep_param.csv"
  shell:
    """
    Rscript src/generate_parameter_sets.R \
      --setting {wildcards.parameters} 
    """
# GPT thingy that is supposed to help me with the wildcard issue but idk
# def get_setting_from_psweep(wc):
#     return wc.psweep.split("_")[-1]

# Define rule for running parameter sweep simulations
rule psweep_sims:
  input:
    "src/simrun_psweep.R",
    "input_data/Treatments_parameter_sweep.csv"
  output:
    "results/parameter_sweeps/{psweep}/{psweep}_out.csv"
  threads:
    5
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
    "results/parameter_sweeps/{psweep}/{psweep}_out.csv"
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
    "results/parameter_sweeps/{psweep}/{psweep}_inv_plot.pdf",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_time_plot.pdf"
  params:
    plot_colors_psweep = json.dumps(plot_colors_psweep)
  shell:
    """
    Rscript src/plot_psweep.R \
      --psweepsetting {wildcards.psweep} \
      --colors '{params.plot_colors_psweep}'
    """
