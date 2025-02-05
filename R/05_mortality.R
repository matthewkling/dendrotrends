
prep_mortality_data <- function(d, species){
      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

      # d <- read_csv("data/formatted_fia_data.csv")
      d <- d %>%
            mutate(bio12 = log(bio12)) %>%
            na.omit() %>%
            filter(species %in% spp) %>%
            filter(class != "d")

      dd <- d %>%
            mutate(mortality = as.integer(class_next == "d"),
                   dia = log((dia + dia_next) / 2),
                   t = year_next - year) %>%
            left_join(scl) %>%
            mutate(dia = (dia - dia_mean) / dia_sd,
                   bacon = (sqrt(ba_j_con + ba_a_con) - bacon_mean) / bacon_sd,
                   bahet = (sqrt(ba_j_het + ba_a_het) - bahet_mean) / bahet_sd,
                   nitrogen = (nitrogen - nitrogen_mean) / nitrogen_sd,
                   sulfur = (sulfur - sulfur_mean) / sulfur_sd,
                   bio1 = (bio1 - bio1_mean) / bio1_sd,
                   bio12 = (bio12 - bio12_mean) / bio12_sd)
      dd
}

# note: this could be combined w growth and recr versions, which differ very little
fit_mortality_model <- function(dd, species, model_file){

      select <- dplyr::select
      spp <- species$spp

      model <- cmdstan_model(model_file)

      fit_model <- function(sp){

            ds <- filter(dd, species == sp)

            y <- ds$mortality

            x <- cbind(ds$dia,
                       ds$bacon,
                       ds$bahet,
                       ds$sulfur,
                       ds$nitrogen,
                       ds$bio1,
                       ds$bio12)

            t <- ds$t

            data <- list(N = nrow(x),
                         P = ncol(x),
                         x = x,
                         y = y,
                         t = t)

            fit <- model$sample(data = data,
                                iter_warmup = 500, iter_sampling = 500,
                                chains = 4, parallel_chains = 4)

            fit$draws(format = "draws_df") %>% as.data.frame() %>% as_tibble() %>%
                  mutate(sp = sp)
      }

      spp %>% map_dfr(possibly(fit_model))
}
