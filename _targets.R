
# set up compute environment ------------------------
library(targets)
tar_option_set(packages = c("tidyverse", "data.table", "janitor", "terra",
                            "raster", "cmdstanr", "wCorr", "furrr", "patchwork",
                            "ggridges"))
sapply(list.files("R", full.names = TRUE), source)
options(future.globals.maxSize = 2000*1024^2)


# define analysis pipeline -------------------------
list(

      # prep data
      tar_target(fia, assemble_fia()),
      tar_target(fia_env, add_environment(fia)),
      tar_target(species, focal_species(fia_env$annual)),
      tar_target(trends, scale_trends(fia_env$trend, species)),

      # load stan models
      tar_target(recr_model, "stan/recruitment.stan", format = "file"),
      tar_target(grow_model, "stan/growth.stan", format = "file"),
      tar_target(mort_model, "stan/mortality.stan", format = "file"),

      # format model inputs and related data
      tar_target(recr_ts, prep_recruitment_data(fia, fia_env$annual, species, summarize = FALSE)),
      tar_target(recr_data, prep_recruitment_data(fia, fia_env$annual, species)),
      tar_target(grow_data, prep_growth_data(fia_env$annual, species)),
      tar_target(mort_data, prep_mortality_data(fia_env$annual, species)),
      tar_target(out_scl, outcome_scale(recr_data, grow_data, mort_data)),

      # fit models
      tar_target(recr_draws, fit_models(recr_data, recr_model)),
      tar_target(grow_draws, fit_models(grow_data, grow_model)),
      tar_target(mort_draws, fit_models(mort_data, mort_model)),

      # calculate predicted values
      tar_target(recr_pred, predict_recruitment(recr_draws, recr_data, trends, recr_ts)),
      tar_target(grow_pred, predict_growth(grow_draws, grow_data, trends)),
      tar_target(mort_pred, predict_mortality(mort_draws, mort_data, trends)),

      # distill results
      tar_target(esr, compile_esr(recr_pred, grow_pred, mort_pred)),
      tar_target(eval, evaluate(recr_pred, grow_pred, mort_pred)),

      # create plots
      tar_target(scatter_plt, combined_scatter(eval)),
      tar_target(exp_plt, exposure_plots(esr, trends, species)),
      tar_target(sens_plt, sensitivity_plots(esr, species, recr_draws, grow_draws, mort_draws, out_scl)),
      tar_target(resp_plt, response_plots(esr, trends)),
      tar_target(imp_plt, importance_plots(esr, species, out_scl))
)

