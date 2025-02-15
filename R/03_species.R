
focal_species <- function(d){

      select <- dplyr::select

      d <- d %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit()

      spp <- d %>%
            select(species, plot_id) %>% distinct() %>%
            count(species) %>%
            arrange(desc(n)) %>% slice(1:100) %>% pull(species)

      scl <- d %>%
            filter(species %in% spp) %>%
            group_by(species) %>%
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
                                      bio12), # log native units
                         .funs = list(mean = mean, sd = sd), na.rm = T)

      gscl <- d %>%
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
                                      bio12), # log native units
                         .funs = list(mean = mean, sd = sd), na.rm = T)

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

      e %>%
            left_join(scl2) %>%
            mutate(value = (value - mean) / sd) %>%
            dplyr::select(-mean, -sd) %>%
            mutate(var = case_when(var == "bio1" ~ "temperature",
                                   var == "bio12" ~ "precipitation",
                                   TRUE ~ var))
}

