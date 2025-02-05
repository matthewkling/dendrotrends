
# Assemble relevant FIA data from original CSV files.
assemble_fia <- function(dir = "~/data/FIA/2024_03"){

      select <- dplyr::select

      ## load and join FIA tables ============================

      ## load FIA tables ############

      # PLOT table
      # primary key: CN
      # unique key: STATECD, INVYR, UNITCD, COUNTYCD, PLOT
      plot <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv"), stringsAsFactors = F,
                    colClasses = c(CN = "character")) %>%
            select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD,
                   # KINDCD, QA_STATUS, ##### these two are added ######
                   STATECD, INVYR, MEASYEAR, UNITCD, COUNTYCD, PLOT,
                   LAT, LON, ELEV) %>% as_tibble()

      # SUBPLOT table
      # primary key: CN
      # unique key: PLT_CN, SUBP
      # When PLOT.MANUAL <1.0, the field crew measured slope at a condition level
      # When PLOT.MANUAL >=1.0, slope is collected at a subplot level, and then the aspect
      # from the subplot representing the greatest proportion of the condition is assigned as a surrogate
      subplot <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv"), stringsAsFactors = F,
                       colClasses = c(CN = "character", PLT_CN = "character")) %>%
            select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN,
                   #MICRCOND, SUBPCOND, # these two were included in dendrodemog but not transitions
                   SUBP_STATUS_CD, POINT_NONSAMPLE_REASN_CD,
                   SLOPE, ASPECT) %>% as_tibble()

      # CONDITION table (used to identify and exclude planted plots)
      # primary key: CN
      # unique key: PLT_CN, CONDID
      cond <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_COND.csv"), stringsAsFactors = F,
                    colClasses = c(CN = "character", PLT_CN = "character")) %>%
            select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
                   SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>% as_tibble()
      planted <- cond %>%
            filter(STDORGCD == 1) %>%
            select(PLT_CN) %>%
            distinct() %>%
            mutate(planted = TRUE)

      # TREE table
      # primary key: CN
      # unique key: PLT_CN, SUBP, TREE
      # Tree numbers (TREE) can be used to track trees when PLOT.DESIGNCD is the same between inventories.
      tree <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_TREE.csv"),
                    stringsAsFactors = F,
                    select = c("PLT_CN", "PLOT", "SUBP", "INVYR", "TREE", "SPCD",
                               "PREVCOND", "CONDID", #"DIST", #### these three were added in dendrodemog ####
                               "DIA", "DIAHTCD", "STATUSCD", "DIACHECK",
                               "RECONCILECD"#, "PREVDIA" #### these two were added in dendrodemog ####
                    ),
                    colClasses = c(PLT_CN = "character", DIST = "numeric")) %>% as_tibble()

      # SPECIES table
      species <- fread(paste0(dir, "/FIADB_REFERENCE/REF_SPECIES.csv"), stringsAsFactors = F) %>%
            select(SPCD, COMMON_NAME, GENUS, SPECIES) %>% as_tibble()


      ## join tables ########

      d <- plot %>%
            left_join(planted) %>%
            left_join(subplot) %>%
            left_join(tree) %>%
            left_join(species) %>%
            mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
                   SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%

            #### this chunk differs slightly from transitions ######
      # filter(is.na(planted),
      #        !is.na(SUBP_STATUS_CD), SUBP_STATUS_CD %in% c(1, 2),
      #        DESIGNCD %in% c(1, 501:506), # 4-subplot design
      #        KINDCD %in% 1:2) %>%
      #       group_by(PLOT_ID) %>%
      #       filter(max(SUBP) < 5) %>% ungroup() %>%
      filter(is.na(planted),
             DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
             SUBP <= 4) %>%



            clean_names() %>%
            mutate(lon = ifelse(lon>0, lon-360, lon)) %>%
            filter(lat > 24) %>%


            mutate(tree_id = paste(subplot_id, designcd, tree)) %>%
            filter(manual >= 2.0,
                   is.finite(lon), is.finite(lat))

      # # measurment year (not invyr) is the correct sampling date
      # d <- d %>%
      #       mutate(yr = measyear) %>%
      #       select(-invyr, -measyear)
      #
      # # operating at the subplot level
      # d <- d %>% mutate(plot_id = subplot_id)

      d
}






# format_fia <- function(trees){
#       # trees <- readRDS("data/fia_trees.rds")
#       select <- dplyr::select
#
#       # initial cleaning
#       t <- trees %>%
#             filter(!is.na(dia)) %>%
#             mutate(sp = paste(genus, species)) %>%
#             mutate(status = recode(statuscd,
#                                    "1" = "live", "2" = "dead", "3" = "harvested",
#                                    .default = NA_character_),
#                    class = ifelse(dia >= 5, "a", "j")) %>%
#             group_by(plot_id) %>%
#             filter(!any(status == "harvested")) %>% # exclude plots with any harvesting
#             ungroup() %>%
#             filter(is.finite(yr),
#                    !is.na(status)) %>%
#             mutate(inv = paste(plot_id, yr))
#       # 8,356,100
#
#       # inventory summaries
#       invn <- t %>%
#             select(inv, plot_id, yr) %>%
#             distinct() %>%
#             arrange(plot_id, yr) %>%
#             group_by(plot_id) %>%
#             arrange(yr) %>%
#             mutate(n_surveys = n(),
#                    yrs = lead(yr) - yr,
#                    final_survey = is.na(yrs)) %>%
#             ungroup()
#
#       # encode state and transition variables
#       t <- t %>%
#             select(inv, yr, plot_id, sp, tree_id, status, class, dia) %>%
#             filter(inv %in% invn$inv) %>%
#             group_by(tree_id) %>%
#
#             # fix resurrections
#             arrange(desc(yr)) %>%
#             mutate(status01 = as.integer(factor(status, levels = c("dead", "live"))) - 1,
#                    lives_on = cummax(status01) == 1,
#                    status = ifelse(status == "dead" & lives_on, "live", status)) %>%
#
#             # fix reversions
#             arrange(desc(yr)) %>%
#             mutate(class01 = as.integer(factor(class, levels = c("j", "a"))) - 1,
#                    future_j = cummin(class01) == 0,
#                    class = ifelse(class == "a" & future_j, "j", class)) %>%
#
#             arrange(yr) %>%
#             mutate(yr_next = lead(yr),
#                    dia_next = lead(dia),
#                    class_next = lead(class),
#                    status_next = lead(status)) %>%
#             ungroup() %>%
#
#             full_join(invn) %>%
#             mutate(trans = case_when(status == "dead" & status_next %in% c("dead", NA) ~ "already_dead",
#                                      status == status_next & class == class_next ~ "survived",
#                                      status == "live" & status_next == "dead" ~ "died",
#                                      status == "dead" & status_next == "live" ~ "resurrected", # shouldn't be necessary given reassignments above
#                                      class == "j" & class_next == "a" ~ "graduated",
#                                      class == "a" & class_next == "j" ~ "reverted", # shouldn't be necessary given reassignments above
#                                      is.na(status_next) & final_survey ~ "not_resurveyed",
#                                      is.na(status_next) & !final_survey ~ "died", # implicit death
#                                      TRUE ~ "other"))
#       # 8,356,100
#
#       # remove cases where trees and plots disagree about sampling interval
#       # (i.e. some trees were skipped for an inventory, or sampled twice in one year)
#       t <- t %>%
#             mutate(years = yr_next - yr,
#                    years = ifelse(is.na(years), yrs, years)) %>%
#             group_by(inv) %>%
#             mutate(year_mismatch = any(years != yrs, na.rm = T)) %>%
#             filter(!year_mismatch) %>%
#             ungroup() %>%
#             select(-year_mismatch)
#       # 8,355,217
#
#       # check that dataset does not contain excluded transition types
#       if(any(t$trans %in% c("resurrected", "reverted"))) stop("reversions/resurrections detected")
#
#
#
#       d <- t %>%
#             filter(status != "dead",
#                    ! final_survey) %>%
#             mutate(yr_next = ifelse(is.na(yr_next), yr + years, yr_next)) %>% # populate this field for implicit deaths
#             mutate(class_next = case_when(class == "j" & trans == "survived" ~ "j",
#                                           class == "j" & trans == "graduated" ~ "a",
#                                           class == "j" & trans == "died" ~ "d",
#                                           class == "a" & trans == "survived" ~ "a",
#                                           class == "a" & trans == "died" ~ "d")) %>%
#             select(species = sp, plot_id, tree_id, dia, dia_next,
#                    year = yr, year_next = yr_next,
#                    class, class_next)
#       # 3,982,087
#
#       ### add densities
#
#       # identify macroplots
#       plot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv",
#                     stringsAsFactors = F, colClasses = c(CN = "character")) %>%
#             select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD, MACRO_BREAKPOINT_DIA,
#                    STATECD, UNITCD, COUNTYCD, PLOT,
#                    LAT, LON, ELEV) %>% as_tibble()
#       subplot <- fread("~/data/FIA/2024_03/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv", stringsAsFactors = F,
#                        colClasses = c(CN = "character", PLT_CN = "character")) %>%
#             select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
#                    SLOPE, ASPECT) %>% as_tibble()
#       macro <- plot %>%
#             left_join(subplot) %>%
#             mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
#                    SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
#             select(plot_id = SUBPLOT_ID, macro_bpd = MACRO_BREAKPOINT_DIA) %>%
#             mutate(macroplot = is.finite(macro_bpd)) %>%
#             distinct()
#
#       # FIA subplot land area, in ha:
#       area <- base::pi * (c(micro = 6.8, sub = 24, macro = 58.9) / 3.28084) ^ 2 / 10000
#
#       # basal area (of neighbors, NOT including self), in m2/ha
#       d <- d %>%
#             left_join(macro) %>%
#             mutate(dia_next = ifelse(is.na(dia_next), dia, dia_next),
#                    diam = (dia + dia_next) / 2, # mean dbh over survey interval
#                    ba = pi * (diam / 2 * 2.54 / 100) ^ 2, # basal area, in m2, e
#                    macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
#             group_by(plot_id, year) %>%
#             mutate(ba_j_tot = (sum(ba[class == "j"]) - ifelse(class == "j", ba, 0)) / area["micro"], # sum neighbors' basal area
#                    ba_a_tot = (sum(ba[class == "a" & diam < macro_bpd]) -
#                                      ifelse(class == "a" & diam < macro_bpd, ba, 0)) / area["sub"] +
#                          (sum(ba[class == "a" & diam >= macro_bpd]) -
#                                 ifelse(class == "a" & diam >= macro_bpd, ba, 0)) / area["macro"]) %>%
#             group_by(plot_id, year, species) %>%
#             mutate(ba_j_con = (sum(ba[class == "j"]) - ifelse(class == "j", ba, 0)) / area["micro"],
#                    ba_a_con = (sum(ba[class == "a" & diam < macro_bpd]) -
#                                      ifelse(class == "a" & diam < macro_bpd, ba, 0)) / area["sub"] +
#                          (sum(ba[class == "a" & diam >= macro_bpd]) -
#                                 ifelse(class == "a" & diam >= macro_bpd, ba, 0)) / area["macro"]) %>%
#             ungroup() %>%
#             mutate(ba_j_het = ba_j_tot - ba_j_con,
#                    ba_a_het = ba_a_tot - ba_a_con)
#
#       # population (including self)
#       d <- d %>%
#             group_by(plot_id, year) %>%
#             mutate(n_j_tot = sum(class == "j") / area["micro"],
#                    n_a_tot = sum(class == "a" & diam < macro_bpd) / area["sub"] +
#                          sum(class == "a" & diam >= macro_bpd) / area["macro"]) %>%
#             group_by(plot_id, year, species) %>%
#             mutate(n_j_con = sum(class == "j") / area["micro"],
#                    n_a_con = sum(class == "a" & diam < macro_bpd) / area["sub"] +
#                          sum(class == "a" & diam >= macro_bpd) / area["macro"]) %>%
#             ungroup() %>%
#             mutate(n_j_het = n_j_tot - n_j_con,
#                    n_a_het = n_a_tot - n_a_con)
#       # 3,982,087
#
#       d
# }

# add_env <- function(trees, d){
#       select <- dplyr::select
#
#       #### add climate
#       ll <- trees %>%
#             filter(plot_id %in% unique(d$plot_id)) %>%
#             dplyr::select(plot_id, lon, lat) %>%
#             distinct()
#       clim <- list.files("/Volumes/T7/CHELSA/v2/raw", full.names = T,
#                          pattern = "bio1_|bio5_|bio6_|bio12_") %>%
#             c("/Volumes/T7/CHELSA/v2/derived/CHELSA_hydro_annual_1981-2010_V.2.1.tif") %>%
#             rast() %>%
#             setNames(c("bio1", "bio12", "bio5", "bio6", "ppt", "pet", "aet", "cwd", "rnr")) %>%
#             extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
#             as.data.frame() %>%
#             as_tibble() %>%
#             mutate(plot_id = ll$plot_id) %>%
#             distinct()
#
#       d <- d %>%
#             left_join(clim) %>%
#             left_join(ll)
#       if(any(is.na(d$bio1))) stop("problem: NA climate values")
#
#
#
#       #### add pollution
#       f <- list.files("~/data/CMAQ/annual_total_adjusted_wet_deposition/",
#                       full.names = T, pattern = "\\.tif", recursive = T)
#       poll <- c(f[grepl("ADJ_TOTDEP_N_CONUS", f)] %>% rast() %>% mean(),
#                 f[grepl("ADJ_TOTDEP_S_CONUS", f)] %>% rast() %>% mean()) %>%
#             setNames(c("nitrogen", "sulfur")) %>%
#             sqrt() ### note this ###
#       pts <- dplyr::select(d, lon, lat)
#       coordinates(pts) <- c("lon", "lat")
#       crs(pts) <- "+proj=longlat"
#       pts <- spTransform(pts, crs(poll))
#       d <- cbind(d, terra::extract(poll, coordinates(pts)))
#       d <- filter(d, is.finite(nitrogen))
#
#
#       d <- d %>%
#             mutate(bio12 = log(bio12)) %>%
#             filter(!str_detect(species, "spp")) %>%
#             na.omit()
#
#       d
# }
