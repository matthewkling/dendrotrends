

prep_recruitment_data <- function(trees, fia_env, species, summarize = TRUE){

      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

      r <- trees %>%
            mutate(sp = paste(genus, species)) %>%
            mutate(status = recode(statuscd,
                                   "1" = "live", "2" = "dead", "3" = "harvested",
                                   .default = NA_character_),
                   class = ifelse(dia >= 5, "a", "j")) %>%
            group_by(plot_id) %>%
            filter(!is.na(status),
                   !is.na(sp)) %>%
            filter(!any(status == "harvested")) %>% # exclude plots with any harvesting
            ungroup() %>%
            filter(is.finite(yr),
                   !is.na(status)) %>%
            mutate(inv = paste(plot_id, yr))

      r <- r %>%
            group_by(plot_id) %>%
            mutate(plot_first_yr = min(yr)) %>%
            group_by(plot_id, tree_id) %>%
            mutate(tree_first_yr = min(yr)) %>%
            filter(length(unique(sp)) == 1) %>%
            ungroup() %>%
            mutate(recruit = yr == tree_first_yr & yr != plot_first_yr
                   & prevcond == 1 & condid == 1 # exclude appearance due to nonforest-forest conversion
                   & dia < 5
                   & reconcilecd == 1 & !is.na(reconcilecd)
                   & status == "live"
            )

      # FIA subplot land area, in ha:
      area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000
      area <- area * 4 # aggregate to plot level

      rr <- r %>%
            mutate(adult = dia >= 5,
                   ba = pi * (dia / 2 * 2.54 / 100) ^ 2, # basal area, in m2, e
                   macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
            group_by(plot_id, sp, yr) %>%
            summarize(recruits = sum(recruit, na.rm = T) / area["micro"], # recruits / ha
                      n = sum(adult & dia < macro_bpd) / area["sub"] + # adult popn / ha
                            sum(adult & dia >= macro_bpd) / area["macro"],
                      ba = sum(ba[adult & dia < macro_bpd]) / area["sub"] + # adult ba / ha
                            sum(ba[adult & dia >= macro_bpd]) / area["macro"]) %>%
            group_by(plot_id, sp) %>%
            arrange(yr) %>%
            reframe(nn = n() - 1,
                    ba = ba[1:nn],
                    n = n[1:nn],
                    recruits = recruits[(1:nn)+1],
                    yr0 = yr[1:nn],
                    yr1 = yr[(1:nn)+1]) %>%
            mutate(t = yr1 - yr0) %>%
            filter(sp %in% spp,
                   !is.na(t))

      d <- fia_env %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit() %>%
            mutate(plot_id = str_sub(plot_id, 1, -3))

      # average predictor vars across subplots and years
      d <- d %>%
            group_by(plot_id, species) %>%
            filter(!is.na(year_next),
                   year_next == max(year_next, na.rm = T)) %>%
            summarize(ba_a_tot = mean(ba_a_tot, na.rm = T),
                      n_a_tot = mean(n_a_tot, na.rm = T),
                      ba_a_het = mean(ba_a_het, na.rm = T),
                      n_a_het = mean(n_a_het, na.rm = T),
                      ba_a_con = mean(ba_a_con, na.rm = T),
                      n_a_con = mean(n_a_con, na.rm = T),
                      bio1 = mean(bio1, na.rm = T),
                      bio12 = mean(bio12, na.rm = T),
                      nitrogen = mean(nitrogen, na.rm = T),
                      sulfur = mean(sulfur, na.rm = T)) %>%
            ungroup()

      # standardize predictors
      d <- d %>%
            left_join(scl) %>%
            mutate(bacon = (sqrt(ba_a_con) - bacon_mean) / bacon_sd,
                   bahet = (sqrt(ba_a_het) - bahet_mean) / bahet_sd,
                   nitrogen = (nitrogen - nitrogen_mean) / nitrogen_sd,
                   sulfur = (sulfur - sulfur_mean) / sulfur_sd,
                   bio1 = (bio1 - bio1_mean) / bio1_sd,
                   bio12 = (bio12 - bio12_mean) / bio12_sd)

      d <- rr %>%
            rename(species = sp) %>%
            left_join(d, by = join_by(plot_id, species))


      # summarize over years (T for fitting, F for prediction)
      if(summarize){
            d <- d %>%
                  group_by(plot_id, species) %>%
                  summarize(
                        recruits = sum(recruits),
                        ba = mean(ba),
                        t = sum(t),
                        n = sum(n),
                        yr0 = min(yr0), yr1 = max(yr1),
                        bacon = mean(bacon),
                        bahet = mean(bahet),
                        bio1 = mean(bio1),
                        bio12 = mean(bio12),
                        nitrogen = mean(nitrogen),
                        sulfur = mean(sulfur)) %>%
                  ungroup()
      }


      # targeted smoothing to avoid undefined recruitment rates where adult BA is zero
      pool_ba <- function(sp){
            s <- d %>%
                  filter(species == sp,
                         t != 0,
                         !is.na(ba),
                         !is.na(n),
                         !is.na(bio1))

            # identify orphan and adult pools
            orphans <- s %>% filter(recruits > 0, n == 0)
            if(nrow(orphans) == 0) return(s)
            adults <- s %>% filter(n > 0)
            others <- s %>% filter(recruits == 0 & n == 0)

            # identify the adult plot most similar to each orphan plot
            vars <- c("bio1", "bio12", "bahet", "nitrogen", "sulfur") # exclude bacon
            nn <- FNN::get.knnx(adults %>% select(all_of(vars)) %>% as.matrix(),
                                orphans %>% select(all_of(vars)) %>% as.matrix(),
                                k = 1)$nn.index
            adults$group <- 1:nrow(adults)
            orphans$group <- adults$group[nn]
            others$group <- -(1:nrow(others))

            # average adult basal area within orphan-adult pairs
            bind_rows(adults, orphans, others) %>%
                  group_by(group) %>%
                  mutate(ba_pooled = mean(ba),
                         n_ba = n()) %>%
                  ungroup()
      }
      d <- map_dfr(spp, pool_ba)

      d <- d %>%
            mutate(outcome = recruits/ba_pooled/t) %>%
            filter(is.finite(outcome)) %>%
            na.omit() %>%
            select(plot_id, species, yr0, yr1, t, ba_pooled, n_ba, outcome,
                   bacon, bahet, sulfur, nitrogen, bio1, bio12)
      d
}

prep_growth_data <- function(d, species){

      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

      d <- d %>%
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

      dd %>% mutate(outcome = sqrt(agr)) %>% filter(species %in% spp) %>%
            select(plot_id, lon, lat, species, tree_id, year, year_next, t, outcome,
                   bacon, bahet, sulfur, nitrogen, bio1, bio12, dia)
}

prep_mortality_data <- function(d, species){

      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

      d <- d %>%
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

      dd %>% mutate(outcome = mortality) %>% filter(species %in% spp) %>%
            select(plot_id, lon, lat, species, tree_id, year, year_next, t, outcome,
                   bacon, bahet, sulfur, nitrogen, bio1, bio12, dia)
}
