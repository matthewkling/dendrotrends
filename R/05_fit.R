
fit_models <- function(d, model_file, n_cores = 9){

      # whether to include DBH as predictor
      dbh <- !str_detect(model_file, "recr")

      # function to fit model for a single species
      fit_model <- function(ds, dbh, init, model_file){

            # independent variable
            y <- ds$outcome

            # predictors
            x <- cbind(ds$bacon,
                       ds$bahet,
                       ds$sulfur,
                       ds$nitrogen,
                       ds$bio1,
                       ds$bio12)
            if(dbh) x <- cbind(x, ds$dia)

            # years between inventories
            t <- ds$t

            # structured model data
            data <- list(N = nrow(x),
                         P = ncol(x),
                         x = x,
                         y = y,
                         t = t)

            # MCMC initial values
            npreds <- ncol(x)
            init <- switch(model_file,
                           "stan/mortality.stan" = function() list(beta = c(runif(1, -4, -3), # intercept
                                                                            rnorm(npreds, 0, .1))), # betas
                           "stan/growth.stan" = function() list(beta = c(runif(1, .09, .14),
                                                                            rnorm(npreds, 0, .01)),
                                                                   sigma = runif(1, .03, .09)),
                           "stan/recruitment.stan" = function() list(zeta = c(runif(1, -6, -4),
                                                                                 rnorm(npreds, 0, .1)),
                                                                        mu = c(runif(1, 1, 3),
                                                                               rnorm(npreds, 0, .1)),
                                                                        sigma = runif(1, .03, .09)))

            # fit model
            fit <- model$sample(data = data, init = init,
                                iter_warmup = 500, iter_sampling = 500, chains = 5,
                                refresh = 0, show_messages = FALSE, show_exceptions = FALSE)

            # return fitted model info
            list(species = ds$species[1],
                 draws = fit$draws(format = "draws_df"),
                 diagnostics = fit$sampler_diagnostics(format = "draws_df"),
                 summary = fit$summary())
      }

      model <- cmdstan_model(model_file)

      future::plan(multisession, workers = n_cores)

      d <- split(d, d$species)
      d[sample(length(d))] %>%
            future_map(possibly(fit_model), dbh = dbh, init = init, model_file = model_file)
}
