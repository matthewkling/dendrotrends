
fit_models <- function(d, model_file, n_cores = 9, dbh = TRUE){

      d <- split(d, d$species)

      model <- cmdstan_model(model_file)

      fit_model <- function(ds, dbh){

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

            fit <- model$sample(data = data,
                                iter_warmup = 500, iter_sampling = 500,
                                chains = 3)

            fit$draws(format = "draws_df") %>%
                  as.data.frame() %>% as_tibble() %>%
                  mutate(sp = ds$species[1])
      }

      future::plan(multisession, workers = n_cores)
      d %>% future_map_dfr(possibly(fit_model), dbh = dbh)
}
