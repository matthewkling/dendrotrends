
prep_recruitment_data <- function(trees, fia_env, species){
      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl


      # ja <- trees %>%
      #       mutate(species = paste(genus, species),
      #              yr = measyear) %>%
      #       filter(species %in% spp,
      #              !is.na(species),
      #              is.finite(lon),
      #              is.finite(dia)) %>%
      #       group_by(species, plot_id, lon, lat) %>%
      #       filter(yr == max(yr)) %>%
      #       summarize(adults = sum(dia >= 5),
      #                 juveniles = sum(dia < 5))

      r <- trees %>%
            mutate(yr = measyear) %>%
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
            ##### filter(yr %in% rev(sort(unique(yr)))[1:2]) %>% # only use data from the two most recent surveys [recently updated to remove this]
            mutate(plot_first_yr = min(yr)) %>%
            group_by(plot_id, tree_id) %>%
            mutate(tree_first_yr = min(yr)) %>%
            # tree_min_dia = min(dia)
            filter(length(unique(sp)) == 1) %>%
            ungroup() %>%
            mutate(recruit = yr == tree_first_yr & yr != plot_first_yr
                   & prevcond == 1 & condid == 1 # exclude appearance due to nonforest-forest conversion
                   & dia < 5
                   & reconcilecd == 1 & !is.na(reconcilecd)
                   & status == "live"
            )


      # identify macroplots
      library(data.table)
      plot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv",
                    stringsAsFactors = F, colClasses = c(CN = "character")) %>%
            select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD, MACRO_BREAKPOINT_DIA,
                   STATECD, UNITCD, COUNTYCD, PLOT,
                   LAT, LON, ELEV) %>% as_tibble()
      subplot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv", stringsAsFactors = F,
                       colClasses = c(CN = "character", PLT_CN = "character")) %>%
            select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
                   SLOPE, ASPECT) %>% as_tibble()
      macro <- plot %>%
            left_join(subplot) %>%
            mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
                   SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
            select(plot_id = PLOT_ID,
                   subplot_id = SUBPLOT_ID,
                   macro_bpd = MACRO_BREAKPOINT_DIA) %>%
            mutate(macroplot = is.finite(macro_bpd)) %>%
            distinct()

      # FIA subplot land area, in ha:
      area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000

      area <- area * 4 # if aggregating to plot level

      rr <- r %>%
            left_join(ungroup(macro)) %>%
            mutate(yr1 = yr == plot_first_yr,
                   adult = dia >= 5,
                   ba = pi * (dia / 2 * 2.54 / 100) ^ 2, # basal area, in m2, e
                   macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
            ##### group_by(plot_id, sp, yr, yr1) %>%
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
            mutate(t = yr1 - yr0,
                   rec_ba = recruits/ba/t,
                   rec_n = recruits/n/t)



      d <- fia_env %>%
            mutate(bio12 = log(bio12)) %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit() %>%
            mutate(plot_id = str_sub(plot_id, 1, -3))

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
            ungroup() %>%
            rename(sp = species)

      d <- left_join(rr, d) %>%
            rename(species = sp)

      d <- d %>%
            left_join(scl) %>%
            mutate(bacon = (sqrt(ba) - bacon_mean) / bacon_sd,
                   bahet = (sqrt(ba_a_het) - bahet_mean) / bahet_sd,
                   nitrogen = (nitrogen - nitrogen_mean) / nitrogen_sd,
                   sulfur = (sulfur - sulfur_mean) / sulfur_sd,
                   bio1 = (bio1 - bio1_mean) / bio1_sd,
                   bio12 = (bio12 - bio12_mean) / bio12_sd)
      d
}


fit_recruitment_model <- function(d, species, model_file){

      select <- dplyr::select
      spp <- species$spp

      model <- cmdstan_model(model_file)

      fit_model <- function(sp){

            dd <- d %>% filter(species == sp, is.finite(rec_ba)) %>% na.omit()

            x <- bind_cols(bacon = dd$bacon,
                           bahet = dd$bahet,
                           sulfur = dd$sulfur,
                           nitrogen = dd$nitrogen,
                           bio1 = dd$bio1,
                           bio12 = dd$bio12)

            data <- list(N = nrow(dd),
                         y = dd$rec_ba,
                         P = ncol(x),
                         x = x,
                         t = dd$t)

            fit <- model$sample(data = data,
                                iter_warmup = 500, iter_sampling = 500,
                                chains = 4, parallel_chains = 4)

            fit$draws(format = "draws_df") %>% as.data.frame() %>% as_tibble() %>%
                  mutate(sp = sp)
      }

      spp %>% map_dfr(fit_model)
}
