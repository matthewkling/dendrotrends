
library(targets)
tar_option_set(packages = c("tidyverse", "data.table", "janitor", "terra",
                            "raster", "cmdstanr"))
sapply(list.files("R", full.names = TRUE), source)
tt <- tar_target


list(
      # data prep
      tt(fia, assemble_fia()),
      tt(fia_env, extract_fia(fia)),
      tt(species, focal_species(fia_env)),
      tt(trend_env, env_trends(fia, species)),

      # stan models
      tt(recruitment_model, "stan/recruit_hurdle.stan", format = "file"),
      tt(growth_model, "stan/growth_lm.stan", format = "file"),
      tt(mortality_model, "stan/mortality.stan", format = "file"),

      # model input data
      tt(recruitment_data, prep_recruitment_data(fia, fia_env, species)),
      tt(growth_data, prep_growth_data(fia_env, species)),
      tt(mortality_data, prep_mortality_data(fia_env, species)),

      # model fitting
      tt(recruitment_draws, fit_recruitment_model(recruitment_data, species, recruitment_model)),
      tt(growth_draws, fit_growth_model(growth_data, species, growth_model)),
      tt(mortality_draws, fit_mortality_model(mortality_data, species, mortality_model)),

      # model predictions
      tt(recruitment_pred, predict_recruitment(recruitment_draws, recruitment_data, trend_env, ndraws = 5)),
      tt(growth_pred, predict_growth(growth_draws, growth_data, trend_env, ndraws = 5)),
      tt(mortality_pred, predict_mortality(mortality_draws, mortality_data, trend_env, ndraws = 5)),

      # evaluation
      tt(eval, evaluate(recruitment_pred, growth_pred, mortality_pred))

)

