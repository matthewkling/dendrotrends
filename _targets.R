# set up compute environment ------------------------
library(targets)
library(tarchetypes)
tar_option_set(packages = c("tidyverse", "data.table", "janitor", "terra",
                            "raster", "cmdstanr", "wCorr", "furrr", "patchwork",
                            "ggridges"))
sapply(list.files("R", full.names = TRUE), source)
options(future.globals.maxSize = 2000*1024^2)


# specify climate parameterizations ------------------------
#
# "baseline" = mean annual temperature + annual precipitation (the primary
#   analysis reported in the main text)
# "alt" = annual max VPD + annual min CMI (robustness analysis)
clim_values <- list(
      clim        = c("baseline", "alt"),
      env         = rlang::syms(c("fia_env_baseline", "fia_env_alt")),
      clim_labels = list(c("temperature", "precipitation"),
                         c("max VPD", "min CMI")),
      suffix      = c("", "_altclim")
)


# define analysis pipeline -------------------------
list(

      # raw FIA data: shared across parameterizations
      tar_target(fia, assemble_fia()),

      # stan models: shared across parameterizations
      tar_target(recr_model, "stan/recruitment.stan", format = "file"),
      tar_target(grow_model, "stan/growth.stan", format = "file"),
      tar_target(mort_model, "stan/mortality.stan", format = "file"),

      # environmental data, one build per climate parameterization
      tar_map(
            values = clim_values[c("clim")],
            names = clim,
            tar_target(fia_env, add_environment(fia, clim = clim))
      ),

      # focal species list: derived ONCE from the baseline environment and shared,
      # so both parameterizations analyze an identical species set
      tar_target(spp, focal_spp(fia_env_baseline$annual)),

      # everything downstream, branched by parameterization
      tar_map(
            values = clim_values,
            names = clim,

            # scaling constants must be recomputed per parameterization
            tar_target(species, species_scaling(env$annual, spp)),
            tar_target(trends,  scale_trends(env$trend, species)),

            # format model inputs and related data
            tar_target(recr_ts,   prep_recruitment_data(fia, env$annual, species, summarize = FALSE)),
            tar_target(recr_data, prep_recruitment_data(fia, env$annual, species)),
            tar_target(grow_data, prep_growth_data(env$annual, species)),
            tar_target(mort_data, prep_mortality_data(env$annual, species)),
            tar_target(out_scl,   outcome_scale(recr_data, grow_data, mort_data)),

            # fit models
            tar_target(recr_draws, fit_models(recr_data, recr_model)),
            tar_target(grow_draws, fit_models(grow_data, grow_model)),
            tar_target(mort_draws, fit_models(mort_data, mort_model)),

            # calculate predicted values
            tar_target(recr_pred, predict_recruitment(recr_draws, recr_data, trends, recr_ts)),
            tar_target(grow_pred, predict_growth(grow_draws, grow_data, trends)),
            tar_target(mort_pred, predict_mortality(mort_draws, mort_data, trends)),

            # distill results
            tar_target(esr,  compile_esr(recr_pred, grow_pred, mort_pred)),
            tar_target(eval, evaluate(recr_pred, grow_pred, mort_pred)),

            # variable importance (Fig. 5) for both parameterizations
            tar_target(imp_plt, importance_plots(esr, species, out_scl,
                                                 clim_labels = clim_labels,
                                                 suffix = suffix))
      ),

      # east-west stratified variance partitioning (baseline parameterization),
      # addressing whether sulfur's importance holds within regions or reflects a
      # between-region contrast
      tar_target(plot_lon, {
            sub <- fia_env_baseline$annual %>%
                  dplyr::select(plot_id, lon) %>%
                  distinct()
            plt <- sub %>%
                  mutate(plot_id = str_sub(plot_id, 1, -3)) %>%
                  group_by(plot_id) %>%
                  summarize(lon = mean(lon), .groups = "drop")
            bind_rows(sub, plt) %>% distinct(plot_id, .keep_all = TRUE)
      }),
      tar_target(strat_imp, stratified_importance(esr_baseline, species_baseline,
                                                  out_scl_baseline, plot_lon)),

      # compile numeric values cited in results section
      tar_target(reported, reported_values(esr_baseline, species_baseline,
                                           out_scl_baseline, trends_baseline,
                                           eval_baseline,
                                           strat_imp = strat_imp,
                                           esr_alt = esr_alt,
                                           species_alt = species_alt,
                                           oscl_alt = out_scl_alt)),

      # remaining figures: baseline only (these have climate-specific
      # back-transforms and axis breaks hardcoded; see note in 08_figures.R)
      tar_target(scatter_plt, combined_scatter(eval_baseline)),
      tar_target(exp_plt,     exposure_plots(esr_baseline, trends_baseline, species_baseline)),
      tar_target(sens_plt,    sensitivity_plots(esr_baseline, species_baseline,
                                                recr_draws_baseline, grow_draws_baseline,
                                                mort_draws_baseline, out_scl_baseline)),
      tar_target(resp_plt,    response_plots(esr_baseline, trends_baseline))
)
