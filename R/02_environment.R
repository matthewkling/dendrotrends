
# add environmental variables (forest density, climate, pollution) to FIA records

add_environment <- function(trees){
      select <- dplyr::select

      # operating at the subplot level
      trees <- trees %>% mutate(plot_id = subplot_id)

      # initial cleaning
      t <- trees %>%
            filter(!is.na(dia)) %>%
            mutate(sp = paste(genus, species)) %>%
            mutate(status = recode(statuscd,
                                   "1" = "live", "2" = "dead", "3" = "harvested",
                                   .default = NA_character_),
                   class = ifelse(dia >= 5, "a", "j")) %>%
            group_by(plot_id) %>%
            filter(!any(status == "harvested")) %>% # exclude plots with any harvesting
            ungroup() %>%
            filter(is.finite(yr),
                   !is.na(status)) %>%
            mutate(inv = paste(plot_id, yr))

      # inventory summaries
      invn <- t %>%
            select(inv, plot_id, yr) %>%
            distinct() %>%
            arrange(plot_id, yr) %>%
            group_by(plot_id) %>%
            arrange(yr) %>%
            mutate(n_surveys = n(),
                   yrs = lead(yr) - yr,
                   final_survey = is.na(yrs)) %>%
            ungroup()

      # encode state and transition variables
      t <- t %>%
            select(inv, yr, plot_id, sp, tree_id, status, class, dia, macroplot, macro_bpd) %>%
            filter(inv %in% invn$inv) %>%
            group_by(tree_id) %>%

            # fix resurrections
            arrange(desc(yr)) %>%
            mutate(status01 = as.integer(factor(status, levels = c("dead", "live"))) - 1,
                   lives_on = cummax(status01) == 1,
                   status = ifelse(status == "dead" & lives_on, "live", status)) %>%

            # fix reversions
            arrange(desc(yr)) %>%
            mutate(class01 = as.integer(factor(class, levels = c("j", "a"))) - 1,
                   future_j = cummin(class01) == 0,
                   class = ifelse(class == "a" & future_j, "j", class)) %>%

            arrange(yr) %>%
            mutate(yr_next = lead(yr),
                   dia_next = lead(dia),
                   class_next = lead(class),
                   status_next = lead(status)) %>%
            ungroup() %>%

            full_join(invn) %>%
            mutate(trans = case_when(status == "dead" & status_next %in% c("dead", NA) ~ "already_dead",
                                     status == status_next & class == class_next ~ "survived",
                                     status == "live" & status_next == "dead" ~ "died",
                                     status == "dead" & status_next == "live" ~ "resurrected", # shouldn't be necessary given reassignments above
                                     class == "j" & class_next == "a" ~ "graduated",
                                     class == "a" & class_next == "j" ~ "reverted", # shouldn't be necessary given reassignments above
                                     is.na(status_next) & final_survey ~ "not_resurveyed",
                                     is.na(status_next) & !final_survey ~ "died", # implicit death
                                     TRUE ~ "other"))

      # remove cases where trees and plots disagree about sampling interval
      # (i.e. some trees were skipped for an inventory, or sampled twice in one year)
      t <- t %>%
            mutate(years = yr_next - yr,
                   years = ifelse(is.na(years), yrs, years)) %>%
            group_by(inv) %>%
            mutate(year_mismatch = any(years != yrs, na.rm = T)) %>%
            filter(!year_mismatch) %>%
            ungroup() %>%
            select(-year_mismatch)

      # check that dataset does not contain excluded transition types
      if(any(t$trans %in% c("resurrected", "reverted"))) stop("reversions/resurrections detected")

      d <- t %>%
            filter(status != "dead",
                   ! final_survey) %>%
            mutate(yr_next = ifelse(is.na(yr_next), yr + years, yr_next)) %>% # populate this field for implicit deaths
            mutate(class_next = case_when(class == "j" & trans == "survived" ~ "j",
                                          class == "j" & trans == "graduated" ~ "a",
                                          class == "j" & trans == "died" ~ "d",
                                          class == "a" & trans == "survived" ~ "a",
                                          class == "a" & trans == "died" ~ "d")) %>%
            select(species = sp, plot_id, tree_id, dia, dia_next,
                   year = yr, year_next = yr_next,
                   macroplot, macro_bpd,
                   class, class_next)


      ### forest density ===================

      # FIA subplot land area, in ha:
      area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000

      ### baseline ###

      # basal area (of neighbors, NOT including self), in m2/ha
      d <- d %>%
            mutate(dia_next = ifelse(is.na(dia_next), dia, dia_next),
                   diam = (dia + dia_next) / 2, # mean dbh over survey interval
                   ba = pi * (diam / 2 * 2.54 / 100) ^ 2, # basal area, in m2, e
                   macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
            group_by(plot_id, year) %>%
            mutate(ba_j_tot = (sum(ba[class == "j"]) - ifelse(class == "j", ba, 0)) / area["micro"], # sum neighbors' basal area
                   ba_a_tot = (sum(ba[class == "a" & diam < macro_bpd]) -
                                     ifelse(class == "a" & diam < macro_bpd, ba, 0)) / area["sub"] +
                         (sum(ba[class == "a" & diam >= macro_bpd]) -
                                ifelse(class == "a" & diam >= macro_bpd, ba, 0)) / area["macro"]) %>%
            group_by(plot_id, year, species) %>%
            mutate(ba_j_con = (sum(ba[class == "j"]) - ifelse(class == "j", ba, 0)) / area["micro"],
                   ba_a_con = (sum(ba[class == "a" & diam < macro_bpd]) -
                                     ifelse(class == "a" & diam < macro_bpd, ba, 0)) / area["sub"] +
                         (sum(ba[class == "a" & diam >= macro_bpd]) -
                                ifelse(class == "a" & diam >= macro_bpd, ba, 0)) / area["macro"]) %>%
            ungroup() %>%
            mutate(ba_j_het = ba_j_tot - ba_j_con,
                   ba_a_het = ba_a_tot - ba_a_con)

      # population (including self)
      d <- d %>%
            group_by(plot_id, year) %>%
            mutate(n_j_tot = sum(class == "j") / area["micro"],
                   n_a_tot = sum(class == "a" & diam < macro_bpd) / area["sub"] +
                         sum(class == "a" & diam >= macro_bpd) / area["macro"]) %>%
            group_by(plot_id, year, species) %>%
            mutate(n_j_con = sum(class == "j") / area["micro"],
                   n_a_con = sum(class == "a" & diam < macro_bpd) / area["sub"] +
                         sum(class == "a" & diam >= macro_bpd) / area["macro"]) %>%
            ungroup() %>%
            mutate(n_j_het = n_j_tot - n_j_con,
                   n_a_het = n_a_tot - n_a_con)



      ### trends ###

      area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000

      # for each species, for each plot, we want trends for con_ba and tot_ba
      dd <- trees %>%
            mutate(species = paste(genus, species),
                   # year = yr,
                   adult = dia >= 5,
                   ba = pi * (dia / 2 * 2.54 / 100) ^ 2, # basal area, in m2
                   macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
            filter(is.finite(ba), ba > 0) %>%
            group_by(lon, lat, plot_id, yr) %>%
            mutate(ba_tot = sum(ba[adult & dia < macro_bpd], na.rm = T) / area["sub"] + # adult + juvenile ba / ha
                         sum(ba[adult & dia >= macro_bpd], na.rm = T) / area["macro"] +
                         sum(ba[!adult], na.rm = T) / area["sub"]) %>%
            group_by(lon, lat, plot_id, yr, species) %>%
            mutate(ba_con = sum(ba[adult & dia < macro_bpd], na.rm = T) / area["sub"] + # adult + juvenile ba / ha
                         sum(ba[adult & dia >= macro_bpd], na.rm = T) / area["macro"] +
                         sum(ba[!adult], na.rm = T) / area["sub"]) %>%
            group_by(species) %>%
            mutate(ba_max = quantile(ba_tot, .95)) %>% # max density across all plots a species occurs in
            group_by(plot_id) %>%
            filter(min(ba_tot) > 0) %>%
            ungroup() %>%
            select(lon, lat, plot_id, species, yr,
                   ba_tot, ba_con, ba_max) %>%
            distinct()

      bb <- dd %>%
            select(lon, lat, plot_id, species, yr) %>%
            na.omit() %>%
            distinct() %>%
            group_by(lon, lat, plot_id) %>%
            filter(length(unique(yr)) >= 2) %>%
            expand(species, yr) %>%
            ungroup()

      ba_trends <- left_join(bb, dd) %>%
            mutate(ba_con = ifelse(is.na(ba_con), 0, ba_con)) %>%
            group_by(plot_id, yr) %>%
            mutate(ba_tot = ifelse(is.na(ba_tot), mean(ba_tot, na.rm = T), ba_tot)) %>%
            mutate(ba_het = ba_tot - ba_con) %>%
            group_by(lon, lat, plot_id, species) %>%
            summarize(ba_con_2009 = pmax(0, project(yr, ba_con, 2009)),
                      ba_con_2010 = pmax(0, project(yr, ba_con, 2010)),
                      ba_het_2009 = pmax(0, project(yr, ba_het, 2009)),
                      ba_het_2010 = pmax(0, project(yr, ba_het, 2010)),
                      .groups = "drop")


      ### climate =======================

      ll <- trees %>%
            filter(plot_id %in% unique(d$plot_id)) %>%
            dplyr::select(plot_id, lon, lat) %>%
            distinct()

      # trends
      pr_ann <- list.files("/Volumes/T7/CHELSA/v2/derived/annual/", full.names = T, pattern = "_pr_") %>%
            rast() %>%
            setNames(paste0("y", 1:nlyr(.))) %>%
            extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
            "/"(100) %>% log() %>%
            apply(1, function(x) c(project(2000:2018, x, 2009),
                                   project(2000:2018, x, 2010))) %>%
            t() %>%
            as.data.frame() %>%
            setNames(c("bio12_2009", "bio12_2010")) %>%
            as_tibble() %>%
            mutate(plot_id = ll$plot_id) %>%
            distinct()
      tas_ann <- list.files("/Volumes/T7/CHELSA/v2/derived/annual/", full.names = T, pattern = "_tas_") %>%
            rast() %>%
            setNames(paste0("y", 1:nlyr(.))) %>%
            extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
            "*"(0.1) %>% "+"(-273.15) %>%
            apply(1, function(x) c(project(2000:2018, x, 2009),
                                   project(2000:2018, x, 2010))) %>%
            t() %>%
            as.data.frame() %>%
            setNames(c("bio1_2009", "bio1_2010")) %>%
            as_tibble() %>%
            mutate(plot_id = ll$plot_id) %>%
            distinct()
      clim_trends <- left_join(pr_ann, tas_ann)

      # means
      clim_means <- clim_trends %>%
            mutate(bio1 = (bio1_2009 + bio1_2010) / 2,
                   bio12 = (bio12_2009 + bio12_2010) / 2) %>%
            select(plot_id, bio1, bio12)
      d <- d %>%
            left_join(clim_means) %>%
            left_join(ll)


      ### pollution ============================

      # CMAQ deposition data, 2002-2019
      f <- list.files("~/data/CMAQ/annual_total_adjusted_wet_deposition/",
                      full.names = T, pattern = "\\.tif", recursive = T)

      # average
      poll <- c(f[grepl("ADJ_TOTDEP_N_CONUS", f)] %>% rast() %>% log() %>% mean(),
                f[grepl("ADJ_TOTDEP_S_CONUS", f)] %>% rast() %>% log() %>% mean()) %>%
            setNames(c("nitrogen", "sulfur"))

      # trend
      slope <- function(x, ...){
            d <- x - mean(x)
            w <- ((1:length(x) - 1) - (length(x)-1)/2)
            weighted.mean(d/w, abs(w))
      }
      trend <- c(f[grepl("ADJ_TOTDEP_N_CONUS", f)] %>% rast() %>% log() %>% app(slope),
                 f[grepl("ADJ_TOTDEP_S_CONUS", f)] %>% rast() %>% log() %>% app(slope)) %>%
            setNames(c("nitrogen", "sulfur"))
      py <- mean(2002:2019)
      p2009 <- poll + trend * (2009 - py)
      p2010 <- poll + trend * (2010 - py)

      poll <- poll %>%
            setNames(c("nitrogen", "sulfur")) %>%
            exp() %>% sqrt()
      poll_trends <- c(p2009, p2010) %>%
            setNames(c("nitrogen_2009", "sulfur_2009", "nitrogen_2010", "sulfur_2010")) %>%
            exp() %>% sqrt()

      pts <- dplyr::select(d, lon, lat)
      coordinates(pts) <- c("lon", "lat")
      crs(pts) <- "+proj=longlat"
      pts <- spTransform(pts, crs(poll))
      poll <- terra::extract(poll, coordinates(pts))
      poll_trends <- terra::extract(poll_trends, coordinates(pts))

      d <- cbind(d, poll)


      # join and format trend data ===========

      e <- d %>%
            select(plot_id) %>%
            bind_cols(poll_trends) %>%
            distinct() %>%
            left_join(ll) %>%
            left_join(clim_trends) %>%
            left_join(ba_trends) %>%
            na.omit()

      e <- e %>%
            mutate(ba_con_2009 = sqrt(ba_con_2009),
                   ba_con_2010 = sqrt(ba_con_2010),
                   ba_het_2009 = sqrt(ba_het_2009),
                   ba_het_2010 = sqrt(ba_het_2010)) %>%
            gather(var, value, -plot_id, -lon, -lat, -species) %>%
            mutate(var = str_replace(var, "ba_", "ba")) %>%
            separate(var, c("var", "year"), sep = "_") %>%
            mutate(year = as.integer(year)) %>%
            filter(year %in% c(2009, 2010),
                   var != "batot")


      d <- filter(d, is.finite(nitrogen))

      list(annual = d,
           trend = e)
}

