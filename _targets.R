
library(targets)
tar_option_set(packages = c("tidyverse", "data.table", "janitor", "terra",
                            "raster", "cmdstanr", "wCorr", "furrr"))
sapply(list.files("R", full.names = TRUE), source)
tt <- tar_target
options(future.globals.maxSize = 2000*1024^2)


list(
      # data prep
      tt(fia, assemble_fia()),
      tt(fia_env, extract_fia(fia)),
      tt(species, focal_species(fia_env$annual)),
      tt(trends, scale_trends(fia_env$trend, species)),

      # stan models
      tt(recr_model, "stan/recruit_hurdle.stan", format = "file"),
      tt(grow_model, "stan/growth_lm.stan", format = "file"),
      tt(mort_model, "stan/mortality.stan", format = "file"),

      # model input data
      tt(recr_data, prep_recruitment_data(fia, fia_env$annual, species)),
      tt(grow_data, prep_growth_data(fia_env$annual, species)),
      tt(mort_data, prep_mortality_data(fia_env$annual, species)),

      # model fitting
      tt(recr_draws, fit_models(recr_data, recr_model, dbh = F)),
      tt(grow_draws, fit_models(grow_data, grow_model)),
      tt(mort_draws, fit_models(mort_data, mort_model)),

      # # model predictions
      tt(recr_pred, predict_recruitment(recr_draws, recr_data, trends, ndraws = 5)),
      tt(grow_pred, predict_growth(grow_draws, grow_data, trends, ndraws = 5)),
      tt(mort_pred, predict_mortality(mort_draws, mort_data, trends, ndraws = 5)),

      # evaluation
      tt(eval, evaluate(recr_pred, grow_pred, mort_pred))

)

