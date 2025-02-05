
focal_species <- function(d){

      select <- dplyr::select

      # d <- read_csv("data/formatted_fia_data.csv")
      d <- d %>%
            mutate(bio12 = log(bio12)) %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit()

      spp <- d %>%
            select(species, plot_id) %>% distinct() %>%
            count(species) %>%
            arrange(desc(n)) %>% slice(1:100) %>% pull(species)
      # saveRDS(spp, "data/focal_species.rds")

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
      # write_csv(scl, "data/var_scl.rds")

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
      # write_csv(gscl, "data/var_scl_global.rds")

      list(spp = spp,
           scl = scl,
           gscl = gscl)
}
