# Identify the focal species set.
#
# This is deliberately computed ONCE (from the baseline environmental data) and
# shared across climate parameterizations, so that every parameterization
# analyzes the same 100 species. Deriving it separately per branch would risk
# slightly different species sets (via na.omit() on differing NA coverage in the
# climate layers), which would confound comparisons between parameterizations.
focal_spp <- function(d, n_spp = 100){

      select <- dplyr::select

      d %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit() %>%
            select(species, plot_id) %>%
            distinct() %>%
            count(species) %>%
            arrange(desc(n)) %>%
            slice(1:n_spp) %>%
            pull(species)
}


# Compute per-species and global predictor scaling constants.
#
# Unlike the species list, this MUST be recomputed for each climate
# parameterization, because it standardizes the climate predictors (bio1 /
# bio12) against their own spatial variation.
species_scaling <- function(d, spp){

      select <- dplyr::select

      d <- d %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit()

      scale_data <- function(d){
            d %>%
                  filter(species %in% spp) %>%
                  mutate(dia = log(dia/2 + dia_next/2)) %>%
                  rename(ba_con = ba_a_con,
                         ba_het = ba_a_het) %>%
                  mutate(bacon = sqrt(ba_con),
                         bahet = sqrt(ba_het)) %>%
                  summarize_at(.vars = vars(dia, # log native units
                                            bacon, # sqrt
                                            bahet, # sqrt
                                            sulfur, # sqrt native units
                                            nitrogen, # sqrt native units
                                            bio1,
                                            bio12), # units depend on parameterization
                               .funs = list(mean = mean, sd = sd), na.rm = T)
      }
      scl <- d %>% group_by(species) %>% scale_data()
      gscl <- d %>% scale_data()

      list(spp = spp,
           scl = scl,
           gscl = gscl)
}


scale_trends <- function(e, species){

      scl <- species$scl

      scl2 <- scl %>%
            gather(var, value, -species) %>%
            mutate(var = str_replace(var, "ba_", "ba")) %>%
            separate(var, c("var", "stat"), sep = "_") %>%
            spread(stat, value)

      # NOTE: "temperature" / "precipitation" are internal slot names for the two
      # climate predictors throughout the downstream pipeline, regardless of which
      # variables actually occupy those slots. Display labels are applied at the
      # figure stage via the `clim_labels` argument.
      e %>%
            left_join(scl2) %>%
            mutate(value = (value - mean) / sd) %>%
            dplyr::select(-mean, -sd) %>%
            mutate(var = case_when(var == "bio1" ~ "temperature",
                                   var == "bio12" ~ "precipitation",
                                   TRUE ~ var))
}
