
prep_growth_data <- function(d, species){

      spp <- species$spp
      scl <- species$scl

      # d <- read_csv("data/formatted_fia_data.csv")
      d <- d %>%
            mutate(bio12 = log(bio12)) %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit()

      dd <- d %>%
            mutate(t = year_next - year,
                   agr = (dia_next / dia) ^ (1/(t)) - 1,
                   agr2 = (dia_next^2 / dia^2) ^ (1/(t)) - 1,
                   dia_next_2 = dia * (1 + agr) ^ t,
                   dens = ba_j_tot + ba_a_tot) %>%
            filter(species %in% spp,
                   class_next != "d") %>%
            group_by(species) %>%
            mutate(dia = log((dia + dia_next) / 2), # average diameter over survey interval
                   ba_con = ba_j_con + ba_a_con,
                   ba_het = ba_j_het + ba_a_het) %>%
            # mutate(bacon = log(ba_con),
            #        bacon = ifelse(!is.finite(bacon), NA, bacon),
            #        bahet = log(ba_het + 1)) %>%
            mutate(bacon = sqrt(ba_con),
                   bahet = sqrt(ba_het)) %>%
            left_join(scl) %>%
            mutate(dia = (dia - dia_mean) / dia_sd,
                   bacon = (bacon - bacon_mean) / bacon_sd,
                   bahet = (bahet - bahet_mean) / bahet_sd,
                   nitrogen = (nitrogen - nitrogen_mean) / nitrogen_sd,
                   sulfur = (sulfur - sulfur_mean) / sulfur_sd,
                   bio1 = (bio1 - bio1_mean) / bio1_sd,
                   bio12 = (bio12 - bio12_mean) / bio12_sd) %>%
            sample_n(min(1e5, length(species))) %>%
            ungroup() %>%
            mutate(spid = as.integer(factor(species)),
                   agr = ifelse(agr < 0, 0, agr)) # clamp dependent var to nonnegative
      # agr = log(agr)) %>%
      # filter(is.finite(agr)) # could instead replace -Inf w smallest finite value -1
      dd
}


fit_growth_model <- function(dd, species, model_file){

      select <- dplyr::select
      spp <- species$spp

      model <- cmdstan_model(model_file)

      fit_model <- function(sp){

            ds <- filter(dd, species == sp) %>%
                  mutate(agr = sqrt(agr)) # sqrt transform

            x <- cbind(ds$dia,
                       ds$bacon,
                       ds$bahet,
                       ds$sulfur,
                       ds$nitrogen,
                       ds$bio1,
                       ds$bio12)

            data <- list(N = nrow(x),
                         P = ncol(x),
                         x = x,
                         y = ds$agr)

            fit <- model$sample(data = data,
                                iter_warmup = 500, iter_sampling = 500,
                                chains = 4, parallel_chains = 4)

            fit$draws(format = "draws_df") %>% as.data.frame() %>% as_tibble() %>%
                  mutate(sp = sp)
      }

      spp %>% map_dfr(fit_model)
}
