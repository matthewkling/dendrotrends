
fit_models <- function(d, model_file, n_cores = 9){

      # whether to include DBH as predictor
      dbh <- !str_detect(model_file, "recr")

      # MCMCM initial values
      npreds <- ifelse(dbh, 7, 6)
      init <- switch(model_file,
                     "stan/mortality.stan" = function() list(beta = c(runif(1, -4, -3), # intercept
                                                                      rnorm(npreds, 0, .1))), # betas
                     "stan/growth_lm.stan" = function() list(beta = c(runif(1, .09, .14),
                                                                      rnorm(npreds, 0, .01)),
                                                             sigma = runif(1, .03, .09)),
                     "stan/recruit_hurdle.stan" = function() list(zeta = c(runif(1, -6, -4),
                                                                           rnorm(npreds, 0, .1)),
                                                                  mu = c(runif(1, 1, 3),
                                                                         rnorm(npreds, 0, .1)),
                                                                  sigma = runif(1, .03, .09)))


      d <- split(d, d$species)

      model <- cmdstan_model(model_file)

      fit_model <- function(ds, dbh, init){

            y <- ds$outcome

            x <- cbind(ds$bacon,
                       ds$bahet,
                       ds$sulfur,
                       ds$nitrogen,
                       ds$bio1,
                       ds$bio12)
            if(dbh) x <- cbind(x, ds$dia)

            t <- ds$t

            data <- list(N = nrow(x),
                         P = ncol(x),
                         x = x,
                         y = y,
                         t = t)

            fit <- model$sample(data = data, init = init,
                                iter_warmup = 400, iter_sampling = 300,
                                refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
                                chains = 5)

            list(species = ds$species[1],
                 draws = fit$draws(format = "draws_df"),
                 diagnostics = fit$sampler_diagnostics(format = "draws_df"),
                 summary = fit$summary())
      }

      future::plan(multisession, workers = n_cores)
      d[sample(length(d))] %>% future_map(possibly(fit_model), dbh = dbh, init = init)
}
