# Snakemake file

# Snakemake Configuration
conda: "PETRI_config.yaml"

import json

## Specify treatments
TREATMENTS = ["DG_S.pB10", "HFC_S.pB10", "LFC_S.pB10-A", "HFC_S.pB10_full", "LFC_S.pB10-A_full"]
PSWEEPS = ["pHFC_S.pB10","pLFC_S.pB10","pDim90_E.R1","pDim90_E.RP4","pHFC_S.pB10-A","pLFC_S.pB10-A","pHFCi_S.pB10","pLFCi_S.pB10","pHFCi2_S.pB10","pLFCi2_S.pB10"]

## Specify global variables
### Colors
plot_colors = {
  "p_Anc": "#8A407A",
  "p_Mut": "#8394F6",
  "p_F": "gray40",
  "p_growth": "#e7e9fd",
  "p_conj": "#d8aace",
  "p_tselect": "#c4bbeb",
  "p_imm": "#d0ffeb",
  "p_hc": "#e7e9fd",
  "p_syn": "#A69AE0",
  "p_par": "#d8aace"
}

plot_colors_psweep = {
  "p_Exc": "gray95",
  "p_Dis": "#140433",
  "p_lowI": "#FFC2C2",
  "p_highI": "#c93251",
  "p_lowD": "#BABAFF",
  "p_highD": "#6c78a8",
  "p_axes": "white",
  "p_mid": "white",
  "p_invline": "gray90"
}

### Linetypes
plot_lines = {
  "rifR_l": "solid",
  "nalR_l": "22",
  "exp_l": "solid",
  "sim_l": "12",
  "plot_lw": 2
}

### Points
plot_points = {
  "psweep_point_size": 2,
  "exp_point_size": 3,
  "ph_point_size": 4,
  "sh_Anc": 16,
  "sh_Mut": 17,
  "sh_R1": 16,
  "sh_copA": 18,
  "sh_finO": 15
}

rule all:
  input:
    expand("results/case_study_sims/{treatment}/{treatment}_frequency_plot.pdf", treatment = TREATMENTS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_change_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.pdf", psweep = PSWEEPS),
    expand("results/parameter_sweeps/{psweep}/{psweep}_inv_change_strain_plot.pdf", psweep = PSWEEPS),
#    "figures/panels/fig1a_axes.pdf",
    "figures/panels/fig3b_phenotyping.pdf",
    "figures/fig4_DG_invasion.pdf",
    "figures/fig5_validation2.pdf",
    "figures/fig5_validation3.pdf",
    "figures/fig6_psweep.pdf",
    "figures/fig7_Dim_psweep.pdf"

# Download cran packages
rule setup_r_environment:
    output:
        touch("results/.r_setup_complete.flag")
    conda:
        "PETRI_config.yaml"
    script:
        "src/install_cran_pkgs.R"

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
    19
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
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot.rds",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot2.pdf",
    "results/parameter_sweeps/{psweep}/{psweep}_inv_change_plot2.rds"
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
  params:
    plot_points = json.dumps(plot_points)
  shell:
    """
    Rscript src/plot_psweep_strain.R \
      --psweepsetting {wildcards.psweep}\
      --points '{params.plot_points}'
    """

# Define rule for processing extinction and colony size data
rule process_extinction:
  input:
    "src/process_extinction.R",
    "input_data/experimental_data/2026-01-27_T_extinction_plating.csv",
    "input_data/experimental_data/2026-01-27_DR_extinction_plating.csv"
  output:
    "results/experimental_validation/extinction_cell_counts/extinction_size_av.csv",
    "results/experimental_validation/extinction_cell_counts/colony_size_T_out.pdf",
    "results/experimental_validation/extinction_cell_counts/extinction_all_out.pdf"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/process_extinction.R \
      --colors '{params.plot_colors}'
    """
    
# Define rule for processing invasion experiments
rule process_experiments:
  input:
    "src/process_experiments.R",
    "results/experimental_validation/extinction_cell_counts/extinction_size_av.csv",
    "input_data/experimental_data/2026-01-27_HFC_plating.csv",
    "input_data/experimental_data/2026-01-27_LFC_plating.csv"
  output:
    "results/experimental_validation/HFC/HFC_plating_processed_av.csv",
    "results/experimental_validation/HFC/HFC_frequency_processed_av.csv",
    "results/experimental_validation/HFC/HFC_ttest.csv",
    "results/experimental_validation/LFC/LFC_plating_processed_av.csv",
    "results/experimental_validation/LFC/LFC_frequency_processed_av.csv",
    "results/experimental_validation/LFC/LFC_ttest.csv"
  shell:
    """
    Rscript src/process_experiments.R
    """

# Define rule for plotting invasion experiments
rule plot_experiments:
  input:
    "src/plot_experiments.R",
    "results/experimental_validation/HFC/HFC_plating_processed_av.csv",
    "results/experimental_validation/HFC/HFC_frequency_processed_av.csv",
    "results/experimental_validation/LFC/LFC_plating_processed_av.csv",
    "results/experimental_validation/LFC/LFC_frequency_processed_av.csv"
  output:
    "results/experimental_validation/HFC/HFC_density_plot.pdf",
    "results/experimental_validation/HFC/HFC_density_plot.rds",
    "results/experimental_validation/HFC/HFC_frequency_plot.pdf",
    "results/experimental_validation/HFC/HFC_frequency_plot.rds",
    "results/experimental_validation/LFC/LFC_density_plot.pdf",
    "results/experimental_validation/LFC/LFC_density_plot.rds",
    "results/experimental_validation/LFC/LFC_frequency_plot.pdf",
    "results/experimental_validation/LFC/LFC_frequency_plot.rds"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_lines = json.dumps(plot_lines)
  shell:
    """
    Rscript src/plot_experiments.R \
      --colors '{params.plot_colors}' \
      --lines '{params.plot_lines}'
    """

# Define rule for processing growth rate data
rule process_growth:
  input:
    "src/process_growth.R",
    "input_data/experimental_data/2026-04-16_OD_growth_curve.xlsx"
  output:
    "results/phenotyping/growth_rate/growthratemax_av.csv",
    "results/phenotyping/growth_rate/growthcurve_allrep.pdf",
    "results/phenotyping/growth_rate/growthcurve_av.pdf",
    "results/phenotyping/growth_rate/growthcurve_av.rds",
    "results/phenotyping/growth_rate/growthrate_av.pdf",
    "results/phenotyping/growth_rate/growthrate_av.rds",
    "results/phenotyping/growth_rate/growthratemax_av.pdf",
    "results/phenotyping/growth_rate/growthratemax_av.rds"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/process_growth.R \
      --colors '{params.plot_colors}'
    """
    
# Define rule for processing phenotyping data
rule process_phenotyping:
  input:
    "src/process_phenotyping.R",
    "results/phenotyping/growth_rate/growthratemax_av.csv"
  output:
    "results/phenotyping/phenotyping_av.csv"
  shell:
    """
    Rscript src/process_phenotyping.R
    """

# Define rule for plotting phenotyping data
rule plot_phenotyping:
  input:
    "src/plot_phenotyping.R",
    "results/phenotyping/phenotyping_av.csv"
  output:
    "results/phenotyping/pB10_phenotyping.pdf",
    "results/phenotyping/pB10_phenotyping.rds",
    "results/phenotyping/pB10_mut_inv.pdf",
    "results/phenotyping/pB10_mut_inv.rds",
    "results/phenotyping/pB10_anc_inv.pdf",
    "results/phenotyping/pB10_anc_inv.rds"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_points = json.dumps(plot_points)
  shell:
    """
    Rscript src/plot_phenotyping.R \
      --colors '{params.plot_colors}'\
      --points '{params.plot_points}'
    """

# FIGURE SCRIPTS 
# Define rule for generating fig 1a: phenotypic axes
rule fig1a:
  input:
    "src/fig_axes.R",
    setup = "results/.r_setup_complete.flag"
  output:
    "figures/panels/fig1a_axes.pdf",
    "figures/panels/fig1a_axes.rds"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/fig_axes.R \
      --colors '{params.plot_colors}'
    """

# Define rule for generating fig 3b: phenotyping annotated
rule fig3b:
  input:
    "src/fig_phenotyping_labeled.R",
    "results/phenotyping/pB10_phenotyping.rds",
    "results/phenotyping/phenotyping_av.csv",
    "results/phenotyping/growth_rate/growthratemax_tt.csv",
    "results/phenotyping/LDM_conjugation/LDM_conjugation_tt.csv"
  output:
    "figures/panels/fig3b_phenotyping.pdf",
    "figures/panels/fig3b_phenotyping.rds"
  params:
    plot_colors = json.dumps(plot_colors)
  shell:
    """
    Rscript src/fig_phenotyping_labeled.R \
      --colors '{params.plot_colors}'
    """

# Define rule for generating fig 5bd: exp validation plots
rule fig5bd:
  input:
    "src/fig_exp_sim.R",
    "src/ggplot_theme.R",
    "results/experimental_validation/HFC/HFC_frequency_processed_av.csv",
    "results/experimental_validation/LFC/LFC_frequency_processed_av.csv",
    "results/case_study_sims/HFC_S.pB10/HFC_S.pB10_frequency_plot_df.csv",
    "results/case_study_sims/LFC_S.pB10-A/LFC_S.pB10-A_frequency_plot_df.csv"
  output:
    "figures/panels/fig5d_HFC.pdf",
    "figures/panels/fig5d_HFC.rds",
    "figures/panels/fig5b_LFC.pdf",
    "figures/panels/fig5b_LFC.rds"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_lines = json.dumps(plot_lines),
    plot_points = json.dumps(plot_points)
  shell:
    """
    Rscript src/fig_exp_sim.R \
      --colors '{params.plot_colors}'\
      --points '{params.plot_points}'\
      --lines '{params.plot_lines}'
    """


# Define rule for generating fig 4: De Gelder invasion simulation
rule fig4:
  input:
    "src/fig_DG.R",
    "results/case_study_sims/DG_S.pB10/DG_S.pB10_density_plot.rds",
    "results/case_study_sims/DG_S.pB10/DG_S.pB10_frequency_plot.rds"
  output:
    "figures/fig4_DG_invasion.pdf"
  shell:
    """
    Rscript src/fig_DG.R 
    """

# Define rule for generating fig 5: selection protocol simulations & experiments
rule fig5:
  input:
    "src/fig_validation.R",
    "results/case_study_sims/HFC_S.pB10/HFC_S.pB10_frequency_plot.rds",
    "results/case_study_sims/LFC_S.pB10-A/LFC_S.pB10-A_frequency_plot.rds",
    "results/experimental_validation/HFC/HFC_frequency_plot.rds",
    "results/experimental_validation/LFC/LFC_frequency_plot.rds",
    "results/phenotyping/phenotyping_av.csv",
    "results/phenotyping/pB10_anc_inv.rds",
    "results/phenotyping/pB10_mut_inv.rds",
    "figures/panels/fig5d_HFC.rds",
    "figures/panels/fig5b_LFC.rds"
  output:
    "figures/fig5_validation.pdf",
    "figures/fig5_validation2.pdf",
    "figures/fig5_validation3.pdf"
  shell:
    """
    Rscript src/fig_validation.R 
    """

# Define rule for generating fig 6: parameter sweeps
rule fig6:
  input:
    "src/fig_psweep.R",
    "results/parameter_sweeps/pHFC_S.pB10/pHFC_S.pB10_inv_change_plot2.rds",
    "results/parameter_sweeps/pLFC_S.pB10/pLFC_S.pB10_inv_change_plot2.rds",
    "results/parameter_sweeps/pHFC_S.pB10-A/pHFC_S.pB10-A_inv_change_plot2.rds",
    "results/parameter_sweeps/pLFC_S.pB10-A/pLFC_S.pB10-A_inv_change_plot2.rds",
    "input_data/strain_phenotypes.csv"
  output:
    "figures/fig6_psweep.pdf"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_points = json.dumps(plot_points)
  shell:
    """
    Rscript src/fig_psweep.R \
      --colors '{params.plot_colors}'\
      --points '{params.plot_points}'
    """

# Define rule for generating fig 7: Dimitriu et al. parameter sweeps
rule fig7:
  input:
    "src/fig_Dim_psweep.R",
    "input_data/strain_phenotypes.csv",
    "results/parameter_sweeps/pDim90_E.R1/pDim90_E.R1_inv_change_plot.rds"
  output:
    "figures/fig7_Dim_psweep.pdf"
  params:
    plot_colors = json.dumps(plot_colors),
    plot_points = json.dumps(plot_points)
  shell:
    """
    Rscript src/fig_Dim_psweep.R \
      --colors '{params.plot_colors}'\
      --points '{params.plot_points}'
    """


