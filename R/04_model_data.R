

prep_recruitment_data <- function(trees, fia_env, species){

      # note: recruitment is plot-level, rather than subplot- or tree-level

      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

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


      # # identify macroplots
      # library(data.table)
      # plot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv",
      #               stringsAsFactors = F, colClasses = c(CN = "character")) %>%
      #       select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD, MACRO_BREAKPOINT_DIA,
      #              STATECD, UNITCD, COUNTYCD, PLOT,
      #              LAT, LON, ELEV) %>% as_tibble()
      # subplot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv", stringsAsFactors = F,
      #                  colClasses = c(CN = "character", PLT_CN = "character")) %>%
      #       select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
      #              SLOPE, ASPECT) %>% as_tibble()
      # macro <- plot %>%
      #       left_join(subplot) %>%
      #       mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
      #              SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
      #       select(plot_id = PLOT_ID,
      #              subplot_id = SUBPLOT_ID,
      #              macro_bpd = MACRO_BREAKPOINT_DIA) %>%
      #       mutate(macroplot = is.finite(macro_bpd)) %>%
      #       distinct()

      # FIA subplot land area, in ha:
      area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000
      area <- area * 4 # if aggregating to plot level

      rr <- r %>%
            # left_join(ungroup(macro)) %>%
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
            # mutate(bio12 = log(bio12)) %>%
            filter(!str_detect(species, "spp")) %>%
            na.omit() %>%
            mutate(plot_id = str_sub(plot_id, 1, -3))

      # average environmental vars across subplots and trees
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
                   bio12 = (bio12 - bio12_mean) / bio12_sd) %>%
            filter(is.finite(rec_ba)) %>% na.omit()

      d %>% mutate(outcome = rec_ba) %>% filter(species %in% spp)
}

prep_growth_data <- function(d, species){

      spp <- species$spp
      scl <- species$scl

      # d <- read_csv("data/formatted_fia_data.csv")
      d <- d %>%
            # mutate(bio12 = log(bio12)) %>%
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

      dd %>% mutate(outcome = sqrt(agr)) %>% filter(species %in% spp)
}

prep_mortality_data <- function(d, species){
      select <- dplyr::select
      spp <- species$spp
      scl <- species$scl

      # d <- read_csv("data/formatted_fia_data.csv")
      d <- d %>%
            # mutate(bio12 = log(bio12)) %>%
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

      dd %>% mutate(outcome = mortality) %>% filter(species %in% spp)
}
