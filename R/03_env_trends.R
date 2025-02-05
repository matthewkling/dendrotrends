
# env_trends.R from dendrodemography
env_trends <- function(trees, species){

      select <- dplyr::select
      # trees <- readRDS("data/fia_trees.rds")



      # forest density ################

      # identify macroplots
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


      # for each species, for each plot, we want trends for con_ba and tot_ba

      d <- trees %>%
            left_join(macro) %>%
            mutate(species = paste(genus, species),
                   year = measyear,
                   adult = dia >= 5,
                   ba = pi * (dia / 2 * 2.54 / 100) ^ 2, # basal area, in m2
                   macro_bpd = ifelse(macroplot, macro_bpd, Inf)) %>%
            filter(is.finite(ba), ba > 0) %>%
            # filter(plot_id == "27 2 135 20936", species == "Acer spicatum") %>% arrange(year) %>% view()
            group_by(lon, lat, plot_id, year) %>%
            mutate(ba_tot = sum(ba[adult & dia < macro_bpd], na.rm = T) / area["sub"] + # adult + juvenile ba / ha
                         sum(ba[adult & dia >= macro_bpd], na.rm = T) / area["macro"] +
                         sum(ba[!adult], na.rm = T) / area["sub"]) %>%
            group_by(lon, lat, plot_id, year, species) %>%
            mutate(ba_con = sum(ba[adult & dia < macro_bpd], na.rm = T) / area["sub"] + # adult + juvenile ba / ha
                         sum(ba[adult & dia >= macro_bpd], na.rm = T) / area["macro"] +
                         sum(ba[!adult], na.rm = T) / area["sub"]) %>%
            group_by(species) %>%
            mutate(ba_max = quantile(ba_tot, .95)) %>% # max density across all plots a species occurs in
            group_by(plot_id) %>%
            filter(min(ba_tot) > 0) %>%
            ungroup() %>%
            select(lon, lat, plot_id, species, year,
                   ba_tot, ba_con, ba_max) %>%
            distinct()

      # p <- d %>%
      #       filter(species %in% sample(species, 100)) %>%
      #       select(species, plot_id, year, ba_tot) %>%
      #       distinct() %>%
      #       arrange(species, plot_id, year) %>%
      #       group_by(species, plot_id) %>%
      #       mutate(dd = lead(ba_tot) / ba_tot) %>%
      #       ggplot(aes(ba_tot, dd, group = species)) +
      #       geom_smooth(se = F, linewidth = .3) +
      #       scale_y_log10() +
      #       coord_cartesian(xlim = c(0, 100))
      # ggsave("figures/density/density_delta_curves.png", p, width = 6, height = 6, units = "in")


      project <- function(x, y, p){
            mx <- mean(x)
            x <- x - mx
            beta <- sum(1/sum(x^2) * x * y) # formula for regression coef if mean(x) == 0
            mean(y) + beta * (p - mx)
      }

      project_logistic <- function(x, y, xp){ # not using because there are too many data points outside valid range
            mx <- mean(x)
            x <- x - mx
            y <- log(y / (1-y))
            beta <- sum(1/sum(x^2) * x * y) # formula for regression coef if mean(x) == 0
            z <- mean(y) + beta * (xp - mx)
            1 / (1 + exp(-z))
      }


      # demonstrate that lm code is legit
      # x <- rnorm(3)
      # x <- x - mean(x)
      # y <- rnorm(3, 10)
      #
      # coef(lm(y ~ x))
      #
      # sum(1/sum(x^2) * x * y)
      #
      # x <- matrix(x, ncol = 1)
      # solve(t(x) %*% x) %*% t(x) %*% y



      bb <- d %>%
            select(lon, lat, plot_id, species, year) %>%
            na.omit() %>%
            distinct() %>%
            group_by(lon, lat, plot_id) %>%
            filter(length(unique(year)) >= 2) %>%
            expand(species, year) %>%
            ungroup()

      dens <- left_join(bb, d) %>% #filter(plot_id == "27 2 135 20936") %>%
            mutate(ba_con = ifelse(is.na(ba_con), 0, ba_con)) %>%
            group_by(plot_id, year) %>%
            mutate(ba_tot = ifelse(is.na(ba_tot), mean(ba_tot, na.rm = T), ba_tot)) %>%
            group_by(plot_id) %>%
            mutate(ba_max = weighted.mean(ba_max, ba_con, na.rm = T)) %>%
            group_by(lon, lat, plot_id, species) %>%
            summarize(ba_max = ba_max[1],
                      ba_con_2000 = pmax(0, project(year, ba_con, 2000)),
                      ba_con_2020 = pmax(0, project(year, ba_con, 2020)),
                      ba_con_2010 = pmax(0, project(year, ba_con, 2010)),
                      ba_con_2080 = pmax(0, project(year, ba_con, 2080))) %>%
            group_by(plot_id) %>%
            mutate(ba_tot_2000 = sum(ba_con_2000),
                   ba_con_2000 = ifelse(ba_tot_2000 > ba_max, ba_con_2000 * ba_max / ba_tot_2000, ba_con_2000),
                   ba_tot_2000 = pmin(ba_max, ba_tot_2000),

                   ba_tot_2020 = sum(ba_con_2020),
                   ba_con_2020 = ifelse(ba_tot_2020 > ba_max, ba_con_2020 * ba_max / ba_tot_2020, ba_con_2020),
                   ba_tot_2020 = pmin(ba_max, ba_tot_2020),

                   ba_tot_2010 = sum(ba_con_2010),
                   ba_con_2010 = ifelse(ba_tot_2010 > ba_max, ba_con_2010 * ba_max / ba_tot_2010, ba_con_2010),
                   ba_tot_2010 = pmin(ba_max, ba_tot_2010),

                   ba_tot_2080 = sum(ba_con_2080),
                   ba_con_2080 = ifelse(ba_tot_2080 > ba_max, ba_con_2080 * ba_max / ba_tot_2080, ba_con_2080),
                   ba_tot_2080 = pmin(ba_max, ba_tot_2080)) %>%
            ungroup()

      # p <- dens %>%
      #       ggplot(aes(ba_tot_2080 / ba_tot_2010)) +
      #       geom_histogram(boundary = 0) +
      #       geom_vline(aes(xintercept = 1)) +
      #       scale_x_log10(limits = c(.05, 20),
      #                     breaks = c(.1, .2, .5, 1, 2, 5, 10))
      # ggsave("figures/density/projected_density_change.png", p, width = 6, height = 6, units = "in")


      # # lat trends
      # dens %>%
      #       select(species, lon, lat, plot_id, ba_con_2010, ba_con_2080) %>%
      #       group_by(species) %>%
      #       filter(n() > 100) %>%
      #       distinct() %>%
      #       mutate(lat = scales::rescale(lat),
      #              ratio = ba_con_2080 / ba_con_2010,
      #              ratio = pmin(10, pmax(1/10, ratio))) %>%
      #       ggplot(aes(lat, ratio, group = species)) +
      #       geom_smooth(method = lm, se = F, size = .5) +
      #       geom_smooth(aes(group = NULL), method = lm, se = F, color = "red") +
      #       scale_y_log10() +
      #       geom_hline(yintercept = 1)


      # p <- dens %>%
      #       filter(lat < 51) %>%
      #       select(lon, lat, ba_tot_2010, ba_tot_2080) %>%
      #       distinct() %>%
      #       mutate(ratio = ba_tot_2080 / ba_tot_2010,
      #              ratio = pmin(4, pmax(1/4, ratio))) %>%
      #       ggplot(aes(lon, lat, color = ratio)) +
      #       geom_point(size = .5) +
      #       scale_color_gradientn(colors = c("darkred", "orange", "gray80", "dodgerblue", "darkblue"),
      #                             trans = "log10") +
      #       theme_void() +
      #       labs(color = "density ratio\n(2080 / 2010)")
      # ggsave("figures/density/projected_density_change_map.png", p, width = 10, height = 5, units = "in")
      #
      # c("Pinus", "Acer", "Quercus") %>%
      #       map(function(taxon){
      #             p <- dens %>%
      #                   filter(grepl(taxon, species)) %>%
      #                   select(species, lon, lat, plot_id, ba_con_2010, ba_con_2080) %>%
      #                   group_by(species) %>%
      #                   filter(n() > 100) %>%
      #                   distinct() %>%
      #                   mutate(ratio = ba_con_2080 / ba_con_2010,
      #                          ratio = pmin(4, pmax(1/4, ratio))) %>%
      #                   ggplot(aes(lon, lat, color = ratio)) +
      #                   facet_wrap(~species, scales = "free") +
      #                   geom_point(size = .5) +
      #                   scale_color_gradientn(colors = c("darkred", "orange", "gray80", "dodgerblue", "darkblue"),
      #                                         trans = "log10") +
      #                   theme_bw() +
      #                   theme(panel.grid = element_blank()) +
      #                   labs(color = "density ratio\n(2080 / 2010)")
      #             ggsave(paste0("figures/density/projected_density_change_maps_", taxon, ".png"),
      #                    p, width = 15, height = 10, units = "in")
      #       })








      # climate #####################

      # library(terra)

      ll <- d %>%
            ungroup() %>%
            dplyr::select(plot_id, lon, lat) %>%
            distinct()


      # trends based on recent annual climate
      pr_ann <- list.files("/Volumes/T7/CHELSA/v2/derived/annual/", full.names = T, pattern = "_pr_") %>%
            rast() %>%
            setNames(paste0("y", 1:nlyr(.))) %>%
            extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
            "/"(100) %>% log() %>%
            apply(1, function(x) c(project(2000:2018, x, 2000),
                                   project(2000:2018, x, 2020))) %>%
            t() %>%
            as.data.frame() %>%
            setNames(c("bio12_2000", "bio12_2020")) %>%
            as_tibble() %>%
            mutate(plot_id = ll$plot_id) %>%
            distinct()
      tas_ann <- list.files("/Volumes/T7/CHELSA/v2/derived/annual/", full.names = T, pattern = "_tas_") %>%
            rast() %>%
            setNames(paste0("y", 1:nlyr(.))) %>%
            extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
            "*"(0.1) %>% "+"(-273.15) %>%
            apply(1, function(x) c(project(2000:2018, x, 2000),
                                   project(2000:2018, x, 2020))) %>%
            t() %>%
            as.data.frame() %>%
            setNames(c("bio1_2000", "bio1_2020")) %>%
            as_tibble() %>%
            mutate(plot_id = ll$plot_id) %>%
            distinct()


      # # trends based on recent and future climatologies
      # clim <- c(list.files("/Volumes/T7/CHELSA/v2/raw", full.names = T, pattern = "bio1_|bio12_"),
      #           c("/Volumes/T7/CHELSA/v2/cmip/CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1.tif",
      #             "/Volumes/T7/CHELSA/v2/cmip/CHELSA_bio12_2071-2100_gfdl-esm4_ssp585_V.2.1.tif")) %>%
      #       rast() %>%
      #       setNames(c("bio1_2010", "bio12_2010", "bio1_2080", "bio12_2080")) %>%
      #       extract(ll %>% select(lon, lat) %>% as.matrix()) %>%
      #       as.data.frame() %>%
      #       as_tibble() %>%
      #       mutate(plot_id = ll$plot_id) %>%
      #       distinct() %>%
      #       mutate(bio12_2010 = log(bio12_2010),
      #              bio12_2080 = log(bio12_2080))
      #
      # clim <- clim %>%
      #       left_join(pr_ann) %>%
      #       left_join(tas_ann)
      clim <- left_join(pr_ann, tas_ann)



      # pollution ##############

      # library(raster)
      slope <- function(x, ...){
            d <- x - mean(x)
            w <- ((1:length(x) - 1) - (length(x)-1)/2)
            weighted.mean(d/w, abs(w))
      }

      # CMAQ deposition data, 2002-2019
      f <- list.files("~/data/CMAQ/annual_total_adjusted_wet_deposition/",
                      full.names = T, pattern = "\\.tif", recursive = T)

      poll <- c(f[grepl("ADJ_TOTDEP_N_CONUS", f)] %>% rast() %>% log() %>% mean(),
                f[grepl("ADJ_TOTDEP_S_CONUS", f)] %>% rast() %>% log() %>% mean()) %>%
            setNames(c("nitrogen", "sulfur"))

      trend <- c(f[grepl("ADJ_TOTDEP_N_CONUS", f)] %>% rast() %>% log() %>% app(slope),
                 f[grepl("ADJ_TOTDEP_S_CONUS", f)] %>% rast() %>% log() %>% app(slope)) %>%
            setNames(c("nitrogen", "sulfur"))

      py <- mean(2002:2019)
      p2000 <- poll + trend * (2000 - py)
      # p2010 <- poll + trend * (2010 - py)
      p2020 <- poll + trend * (2020 - py)
      # p2080 <- poll + trend * (2080 - py)

      pts <- dplyr::select(ll, lon, lat)
      coordinates(pts) <- c("lon", "lat")
      crs(pts) <- "+proj=longlat"
      pts <- spTransform(pts, crs(poll))

      poll <- c(p2000, #p2010,
                p2020#, p2080
                ) %>%
            setNames(c("nitrogen_2000", "sulfur_2000",
                       # "nitrogen_2010", "sulfur_2010",
                       "nitrogen_2020", "sulfur_2020"#,
                       # "nitrogen_2080", "sulfur_2080"
                       )
                     ) %>%
            exp() %>%
            sqrt() %>% # sqrt pollution is what went into the models
            terra::extract(coordinates(pts)) %>%
            as.data.frame() %>%
            as_tibble() %>%
            mutate(plot_id = ll$plot_id) %>%
            distinct()

      pd <- left_join(poll, ll) %>%
            # gather(var, value, nitrogen_2000:sulfur_2080) %>%
            gather(var, value, nitrogen_2000:sulfur_2020) %>%
            mutate(value = value^2) %>% # back to native units
            separate(var, c("var", "year")) %>%
            mutate(year = paste0("y", year)) %>%
            spread(year, value) %>%
            filter(lat < 51)

      # p <- pd %>%
      #       ggplot(aes(lon, lat, color = y2010)) +
      #       facet_wrap(~var) +
      #       geom_point(size = .25) +
      #       scale_color_viridis_c(trans = "sqrt") +
      #       guides(color = guide_colorbar(barwidth = 15)) +
      #       theme_void() +
      #       theme(legend.position = "top") +
      #       labs(color = "2010 deposition rate (kg/ha/yr) ")
      # ggsave("figures/pollution/pollution_maps_plots.png", p, width = 12, height = 5, units = "in")
      #
      # p <- pd %>%
      #       ggplot(aes(lon, lat, color = y2080 / y2010)) +
      #       facet_wrap(~var) +
      #       geom_point(size = .25) +
      #       scale_color_gradientn(colors = c("darkblue", "dodgerblue", "gray80", "orange", "darkred"),
      #                             limits = 10^(max(abs(log10(pd$y2080/pd$y2010)), na.rm = T) * c(-1, 1)),
      #                             breaks = c(.0001, .01, 1, 100, 10000),
      #                             trans = "log10") +
      #       guides(color = guide_colorbar(barwidth = 15)) +
      #       theme_void() +
      #       theme(legend.position = "top") +
      #       labs(color = "2080/2010 deposition trends (ratio) ")
      # ggsave("figures/pollution/pollution_maps_plots_trend.png", p, width = 12, height = 5, units = "in")



      # combine and export ###############

      e <- d %>%
            select(plot_id) %>%
            distinct() %>%
            left_join(ll) %>%
            left_join(clim) %>%
            left_join(poll) %>%
            left_join(dens) %>%
            na.omit()
      # write_csv(e, "data/env_trends.csv")

      # from rgm_evaluation.R:
      # projected environmental data for each plot
      # e <- read_csv("data/env_trends.csv") %>%
      e <- e %>%
            mutate(ba_het_2000 = ba_tot_2000 - ba_con_2000,
                   ba_het_2020 = ba_tot_2020 - ba_con_2020) %>%
            mutate(ba_con_2000 = sqrt(ba_con_2000),
                   ba_con_2020 = sqrt(ba_con_2020),
                   ba_het_2000 = sqrt(ba_het_2000),
                   ba_het_2020 = sqrt(ba_het_2020)) %>%
            select(-ba_max) %>%
            gather(var, value, -plot_id, -lon, -lat, -species) %>%
            mutate(var = str_replace(var, "ba_", "ba")) %>%
            separate(var, c("var", "year"), sep = "_") %>%
            mutate(year = as.integer(year)) %>%
            filter(year %in% c(2000, 2020),
                   var != "batot")

      # standardize
      # scl <- read_csv("data/var_scl.rds")
      scl <- species$scl
      scl2 <- scl %>%
            gather(var, value, -species) %>%
            mutate(var = str_replace(var, "ba_", "ba")) %>%
            separate(var, c("var", "stat"), sep = "_") %>%
            spread(stat, value)
      e <- e %>%
            left_join(scl2) %>%
            mutate(value = (value - mean) / sd) %>%
            select(-mean, -sd)

      e <- e %>%
            mutate(var = case_when(var == "bio1" ~ "temperature",
                                   var == "bio12" ~ "precipitation",
                                   TRUE ~ var))
      # write_csv(e, "data/e.csv")
      e
}
