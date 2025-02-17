
# this fx needs to be renamed, at a minimum:
compile_esr <- function(recr_pred, grow_pred, mort_pred){
      f <- bind_rows(recr_pred$esr %>% mutate(model = "recruitment",
                                              tree_id = NA_character_),
                     grow_pred$esr %>% mutate(model = "growth"),
                     mort_pred$esr %>% mutate(model = "mortality")) %>%
            mutate(stat = factor(stat, levels = c("exposure", "sensitivity", "response")),
                   model = paste(model, param),
                   model = factor(model,
                                  levels = c("growth beta", "mortality beta", "recruitment zeta", "recruitment mu"),
                                  labels = c("growth", "mortality", "recruitment\nprobability", "recruitment\nrate")),
                   value = ifelse(stat == "sensitivity", y2009, delta))

      f <- f %>%
            filter(stat == "response") %>%
            group_by(i, species, model, plot_id, stat, tree_id) %>%
            summarize(value = sum(value), .groups = "drop") %>%
            mutate(var = "combined") %>%
            bind_rows(f) %>%
            mutate(var = factor(var, levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet", "combined")))

      f
}


format_scl <- function(x){
      x %>%
            gather(var, value, contains("_mean"), contains("_sd")) %>%
            separate(var, c("var", "stat")) %>%
            spread(stat, value) %>%
            mutate(var = recode(var, "bio1" = "temperature", "bio12" = "precipitation", "dia" = "dbh"))
}

exposure_plots <- function(f, e, species){
      # f <- compile_esr(tar_read(recr_pred), tar_read(grow_pred), tar_read(mort_pred))
      # e <- tar_read(trends)

      select <- dplyr::select

      ll <- select(e, plot_id, lon, lat) %>% distinct()



      # e <- read_csv("data/e.csv")


      # f <- bind_rows(read_csv("data/sensitivity_recruitment.csv") %>% mutate(model = "recruitment"),
      #                read_csv("data/sensitivity_growth.csv") %>% mutate(model = "growth"),
      #                read_csv("data/sensitivity_mortality.csv") %>% mutate(model = "mortality")) %>%
      #       mutate(stat = factor(stat, levels = c("exposure", "sensitivity", "response")),
      #              model = paste(model, param),
      #              model = factor(model,
      #                             levels = c("growth beta", "mortality beta", "recruitment zeta", "recruitment mu"),
      #                             labels = c("growth", "mortality", "recruitment\nprobability", "recruitment\nrate")),
      #              value = ifelse(stat == "sensitivity", y2009, delta))

      # f <- f %>%
      #       filter(stat == "response") %>%
      #       group_by(i, species, model, plot_id, stat, tree_id) %>%
      #       summarize(value = sum(value)) %>%
      #       mutate(var = "combined") %>%
      #       bind_rows(f) %>%
      #       mutate(var = factor(var, levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet", "combined")))

      factor_var <- function(x){
            mutate(x,
                   var = factor(var,
                                levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet"),
                                labels = c("temperature", "precipitation", "nitrogen\ndeposition", "sulfur\ndeposition",
                                           "conspecific\nforest density", "heterospecific\nforest density")))
      }




      ### scaled by species -----------

      pd <- f %>%
            filter(stat == "exposure") %>%
            left_join(ll) %>%
            mutate(lon = plyr::round_any(lon, .5),
                   lat = plyr::round_any(lat, .5))
      pdll <- pd %>%
            factor_var() %>%
            group_by(var, lon, lat) %>%
            summarize(value = mean(value, na.rm = T)) %>%
            group_by(var) %>%
            mutate(value = value / sd(value),
                   value = pmax(quantile(value, .0025), pmin(quantile(value, .9975), value)))

      maps <- pdll %>%
            ggplot(aes(lon, lat, fill = value)) +
            facet_grid(var ~ .) +
            geom_raster() +
            scale_fill_gradientn(colors = c("black", "darkred", "orangered", "gray", "dodgerblue", "darkblue", "black"),
                                 values = c(0, .3, .45, .5, .55, .7, 1), limits = max(abs(c(pdll$value, pdll$value))) * c(-1, 1)) +
            theme_bw() +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black"))

      dens <- pdll %>%
            ggplot(aes(value, 0, fill = stat(x), height = after_stat(scaled))) +
            geom_density_ridges_gradient(adjust = .5, stat = "density", trim = F) +
            geom_vline(xintercept = 0, color = "black", alpha = .35) +
            facet_grid(var ~ ., scales = "free") +
            scale_fill_gradientn(colors = c("black", "darkred", "orangered", "gray", "dodgerblue", "darkblue", "black"),
                                 values = c(0, .3, .45, .5, .55, .7, 1), limits = max(abs(c(pdll$value, pdll$value))) * c(-1, 1)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 1.9)) +
            coord_flip() +
            theme_bw() +
            theme(legend.position = "none",
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black"),
                  panel.grid = element_blank(),
                  strip.text.y = element_blank()) +
            labs(x = "mean standardized rate of change in environmental variable (species range stdev / yr)",
                 y = "relative frequency across\ncommunities (grid cells)")

      p <- dens + maps + plot_layout(nrow = 1, widths = c(1, 2))
      ggsave("figures/exposure/exposure_maps_dens.pdf", p, width = 4, height = 8, units = "in")



      exposure_plot <- function(v){

            ev <- filter(pdll, var == v)

            peak <- .6 # desired location of peak
            dens <- density(ev$value, adjust = 2)
            mode <- dens$x[dens$y == max(dens$y)]
            xmin <- min(ev$value)
            xmax <- max(ev$value)
            pp <- (mode - xmin) / (xmax - xmin)
            if(pp > peak) xmax <- xmin + (xmax - xmin) * pp / peak
            if(pp < peak) xmin <- xmax - (xmax - xmin) * peak / pp

            pal <- c("black", "darkred", "orangered", "gray", "dodgerblue", "darkblue", "black")
            if(v != "precipitation") pal <- rev(pal)
            colors <- scale_fill_gradientn(colors = pal,
                                           values = c(0, .3, .4, .5, .6, .7, 1),
                                           limits = max(abs(ev$value)) * c(-1, 1))

            map <- ev %>%
                  ggplot(aes(lon, lat, fill = value)) +
                  geom_raster() +
                  colors +
                  scale_y_continuous(limits = c(min(ev$lat) - 8, NA)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            dens <- ev %>%
                  ggplot(aes(value, 0, fill = stat(x), height = after_stat(scaled))) +
                  geom_density_ridges_gradient(adjust = .5, stat = "density", trim = F) +
                  geom_vline(xintercept = 0, color = "black", alpha = .35) +
                  geom_hline(yintercept = 0, color = "black", alpha = .35) +
                  colors +
                  scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.9)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.background = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank())

            map + inset_element(dens, left = 0, right = 1, bottom = .07, top = .35)
      }

      p <- levels(pdll$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot) %>%
            Reduce("+", .)
      ggsave("figures/exposure/exposure_maps_dens_v2.png", p, width = 8, height = 4.5, units = "in")


      ### raw units not scaled by species -----------

      # remove scaling
      scl <- species$scl %>%
            gather(var, value, -species) %>%
            separate(var, c("var", "stat"), sep = "_") %>%
            spread(stat, value) %>%
            mutate(var = case_when(var == "bio1" ~ "temperature",
                                   var == "bio12" ~ "precipitation",
                                   TRUE ~ var))
      ee <- e %>%
            left_join(scl) %>%
            mutate(value = value * sd + mean)

      dt <- 1

      ell <- ee %>%
            mutate(lon = plyr::round_any(lon, .5),
                   lat = plyr::round_any(lat, .5)) %>%
            factor_var() %>%
            group_by(var, lon, lat, year) %>%
            summarize(value = mean(value, na.rm = T)) %>%
            mutate(year = paste0("y", year)) %>%
            mutate(value = ifelse(var %in% c("bacon", "bahet", "nitrogen", "sulfur"),
                                  value^2, value),
                   value = ifelse(var == "precipitation", exp(value), value)) %>% # back to native mm units
            spread(year, value) %>%
            mutate(mean = y2010/2 + y2009/2,
                   delta = y2010 - y2009,
                   rate = delta / dt,
                   rate_prop = ifelse(var == "temperature", rate, (y2010 / y2009) ^ (1/dt) - 1)) %>%
            group_by(var) %>%
            mutate(rate_std = rate / sd(mean, na.rm = T))

      exposure_plot_3 <- function(v = "temperature", variable = "rate", limits = "global"){

            ev <- ell %>% filter(var == v)
            ev$value <- ev[[variable]]
            ev <- ev %>%
                  filter(is.finite(value)) %>%
                  mutate(value = pmax(quantile(value, .005), pmin(quantile(value, .995), value)))


            lim <- max(abs(ev$value)) * c(-1, 1)
            if(limits == "global") lim <- max(abs(ell[[variable]]), na.rm = T) * c(-1, 1)


            pal <- c("black", "darkred", "orangered", "gray", "dodgerblue", "darkblue", "black")
            if(v != "precipitation") pal <- rev(pal)
            colors <- scale_fill_gradientn(colors = pal,
                                           values = c(0, .3, .4, .5, .6, .7, 1),
                                           limits = lim)

            peak <- .6 # desired location of peak
            dens <- density(ev$value, adjust = 3)
            mode <- dens$x[dens$y == max(dens$y)]
            xmin <- min(ev$value)
            xmax <- max(ev$value)
            pp <- (mode - xmin) / (xmax - xmin)
            if(pp > peak) xmax <- xmin + (xmax - xmin) * pp / peak
            if(pp < peak) xmin <- xmax - (xmax - xmin) * peak / pp

            brk <- switch(v,
                          "temperature" = c(0, .04, .08),
                          "precipitation" = c(-20, 0, 20),
                          "nitrogen\ndeposition" = c(-.06, -.03, 0, .03),
                          "sulfur\ndeposition" = c(-.2, -.1, 0),
                          "conspecific\nforest density" = c(-.04, 0, .04),
                          "heterospecific\nforest density" = c(-.05, 0, .05))
            xscl <- scale_x_continuous(expand = c(0, 0), breaks = brk, limits = c(xmin, xmax))
            if(limits == "global") xscl <- scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax))
            if(variable == "rate_prop" & v != "temperature") xscl <- scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax),
                                                                                        labels = scales::percent)
            if(variable == "rate_prop" & v == "temperature") xscl <- scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax),
                                                                                        breaks = brk)

            map <- ev %>%
                  ggplot(aes(lon, lat, fill = value)) +
                  geom_raster() +
                  colors +
                  scale_y_continuous(limits = c(min(ev$lat) - 8, NA)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            dens <- ev %>%
                  ggplot(aes(value, 0, fill = stat(x), height = after_stat(scaled))) +
                  geom_density_ridges_gradient(adjust = .5, linewidth = .25,
                                               stat = "density", trim = F) +
                  geom_vline(xintercept = 0, color = "black", alpha = .35) +
                  geom_hline(yintercept = 0, color = "black", alpha = .35) +
                  colors +
                  xscl +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.background = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank())

            # map + inset_element(dens, left = 0.01, right = .5, bottom = -0.02, top = .45)
            map + inset_element(dens, left = 0, right = 1, bottom = .07, top = .35)
      }


      p <- levels(ell$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot_3, limits = "specific") %>%
            Reduce("+", .)
      ggsave("figures/exposure/exposure_maps_dens_v3_rate.png", p, width = 8, height = 4.5, units = "in")

      p <- levels(ell$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot_3, variable = "rate_std") %>%
            Reduce("+", .)
      ggsave("figures/exposure/exposure_maps_dens_v3_rate_std.png", p, width = 8, height = 4.5, units = "in")

      p <- levels(ell$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot_3, variable = "rate_prop", limits = "specific") %>%
            Reduce("+", .)
      ggsave("figures/exposure/exposure_maps_dens_v3_rate_prop.png", p, width = 8, height = 4.5, units = "in")




      # color by rate_std, label by rate
      exposure_plot_4 <- function(v = "temperature"){

            variable <- "rate_std"

            s <- 10

            limits <- "global"

            ev <- ell %>% filter(var == v)

            ev$rate <- ev$rate * s
            ev$rate_std <- ev$rate_std * s

            ev$value <- ev[[variable]]
            ev <- ev %>%
                  filter(is.finite(value)) %>%
                  mutate(value = pmax(quantile(value, .005), pmin(quantile(value, .995), value)))

            lim <- max(abs(ell[[variable]] * s), na.rm = T) * c(-1, 1)


            pal <- c("black", "darkred", "orangered", "gray", "dodgerblue", "darkblue", "black")
            if(v != "precipitation") pal <- rev(pal)
            colors <- scale_fill_gradientn(colors = pal,
                                           values = c(0, .44, .48, .5, .52, .56, 1),
                                           limits = lim)
            # colors2 <- scale_color_gradientn(colors = pal,
            #                                values = c(0, .44, .48, .5, .52, .56, 1),
            #                                limits = lim)

            peak <- .6 # desired location of peak
            dens <- density(ev$value, adjust = 3)
            mode <- dens$x[dens$y == max(dens$y)]
            xmin <- min(ev$value)
            xmax <- max(ev$value)
            pp <- (mode - xmin) / (xmax - xmin)
            if(pp > peak) xmax <- xmin + (xmax - xmin) * pp / peak
            if(pp < peak) xmin <- xmax - (xmax - xmin) * peak / pp

            brk <- switch(v,
                          "temperature" = c(0, .04, .08),
                          "precipitation" = c(-20, 0, 20),
                          # "nitrogen\ndeposition" = c(-.06, -.03, 0, .03),
                          "nitrogen\ndeposition" = c(-.04, -.02, 0, .02),
                          "sulfur\ndeposition" = c(-.2, -.1, 0, .1),
                          "conspecific\nforest density" = c(-.04, 0, .04),
                          "heterospecific\nforest density" = c(-.05, 0, .05)) * 10

            stdev <- mean(ev$rate / ev$rate_std, na.rm = T)
            xscl <- scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax),
                                       breaks = brk / stdev, labels = brk)
            # brk1 <- brk[which(brk == 0) + 1]
            brk1 <- tail(brk, 1)

            map <- ev %>%
                  ggplot(aes(lon, lat, fill = value)) +
                  geom_raster() +
                  colors +
                  scale_y_continuous(limits = c(min(ev$lat) - 8, NA)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            # ticks <- tibble(x = seq(-10, 10, .1)) %>%
            #       filter(between(x, lim[1], lim[2]))

            dens <- ev %>%
                  ggplot(aes(value, 0, fill = after_stat(x), height = after_stat(scaled))) +
                  geom_density_ridges_gradient(adjust = .5, linewidth = .15,
                                               stat = "density", trim = F) +
                  # geom_vline(data = ticks, aes(color = x, xintercept = x), linewidth = .15) +
                  annotate(geom = "text", x = brk1/stdev, y = 0,
                           vjust = -.8, size = 2.25, lineheight = .8, fontface = "bold", color = "gray30",
                           label = paste0("+", signif(brk1/stdev, 2), "\nsd / dec")) +
                  annotate(geom = "text", x = brk1/stdev, y = 0, vjust = 0.25, hjust = .4, size = 4,
                           label = paste0("|")) +
                  # annotate(geom = "text", x = brk1/stdev, y = 0, label = "0.1 sd/dec\n|", vjust = 0) +
                  geom_vline(xintercept = 0, color = "black", linewidth = .45) +
                  geom_hline(yintercept = 0, color = "black", linewidth = .15) +
                  colors + #colors2 +
                  xscl +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.background = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank())

            # map + inset_element(dens, left = 0.01, right = .5, bottom = -0.02, top = .45)
            map + inset_element(dens, left = 0, right = 1, bottom = .07, top = .38, on_top = F, clip = F)
      }


      p <- levels(ell$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot_4) %>%
            Reduce("+", .)
      ggsave("figures/exposure/exposure_maps_v4.png", p, width = 8, height = 4.5, units = "in")

}


sensitivity_plots <- function(f, species, recr_draws, grow_draws, mort_draws){

      select <- dplyr::select

      # get_fits <- function(x, recruitment = F){
      #
      #       if(recruitment){
      #             x <- x %>%
      #                   read_csv() %>%
      #                   select(-lp__) %>%
      #                   gather(param, value, `zeta[1]`:`mu[7]`) %>%
      #                   mutate(param = str_remove(param, "\\]"),
      #                          param = str_replace(param, "\\[", "_")) %>%
      #                   separate(param, c("param", "var")) %>%
      #                   mutate(var = recode(var, "1" = "intercept",
      #                                       "2" = "bacon", "3" = "bahet",
      #                                       "4" = "sulfur", "5" = "nitrogen",
      #                                       "6" = "temperature", "7" = "precipitation")) %>%
      #                   select(i = .draw, species = sp, var, param, param_value = value)
      #       }else{
      #             x <- x %>%
      #                   read_csv() %>%
      #                   select(-lp__) %>%
      #                   gather(param, value, `beta[1]`:`beta[8]`) %>%
      #                   mutate(param = str_remove(param, "\\]"),
      #                          param = str_replace(param, "\\[", "_")) %>%
      #                   separate(param, c("param", "var")) %>%
      #                   select(i = .draw, species = sp, var, param, param_value = value) %>%
      #                   mutate(var = recode(var, "1" = "intercept", "2" = "dbh",
      #                                       "3" = "bacon", "4" = "bahet",
      #                                       "5" = "sulfur", "6" = "nitrogen",
      #                                       "7" = "temperature", "8" = "precipitation"))
      #       }
      #       x
      # }
      #
      # s <- bind_rows(get_fits("data/draws_growth.csv") %>% mutate(model = "growth"),
      #                get_fits("data/draws_mortality.csv") %>% mutate(model = "mortality"),
      #                get_fits("data/draws_recruitment.csv", recruitment = TRUE) %>% mutate(model = "recruitment"))



      s <- bind_rows(load_fits(get_draws(grow_draws)) %>% mutate(model = "growth"),
                     load_fits(get_draws(mort_draws)) %>% mutate(model = "mortality"),
                     load_fits(get_draws(recr_draws)) %>% mutate(model = "recruitment"))


      nn <- f %>% ungroup() %>% count(species)


      scl <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)

      # scl <- scl %>% #read_csv("data/var_scl.rds") %>%
      #       gather(var, value, -species) %>%
      #       separate(var, c("var", "stat")) %>%
      #       spread(stat, value) %>%
      #       mutate(var = recode(var, "bio1" = "temperature", "bio12" = "precipitation", "dia" = "dbh"))
      # gscl <- gscl %>% #read_csv("data/var_scl_global.rds") %>%
      #       gather(var, value) %>%
      #       separate(var, c("var", "stat")) %>%
      #       spread(stat, value) %>%
      #       mutate(var = recode(var, "bio1" = "temperature", "bio12" = "precipitation", "dia" = "dbh"))
      s <- s %>%
            filter(var != "intercept") %>%
            left_join(scl) %>%
            rename(spec = param_value) %>%
            mutate(glob = spec * sd) %>% # remove species scale
            select(-mean, -sd) %>%
            left_join(gscl) %>%
            mutate(glob = glob / sd) %>% # add global scale
            select(-mean, -sd) %>%
            gather(scaling, param_value, spec, glob)


      ds <- s %>%
            group_by(model, species, param, var, scaling) %>%
            summarize(credible = sign(quantile(param_value, .05)) == sign(quantile(param_value, .95)),
                      param_value = median(param_value, na.rm = T)) %>%
            mutate(param = paste(model, param)) %>%
            ungroup() %>%
            filter(var != "intercept") %>%
            mutate(param = str_replace(param, "mu", "rate"),
                   param = str_replace(param, "zeta", "probability"),
                   param = str_replace(param, " beta", ""),
                   param = ifelse(param == "mortality", "survival", param),
                   param_value = ifelse(param == "survival", -param_value, param_value)) %>%
            mutate(var_group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                         var %in% c("bahet", "bacon") ~ "forest density",
                                         var %in% c("sulfur", "nitrogen") ~ "pollution",
                                         TRUE ~ "other"),
                   var_group = factor(var_group, levels = c("climate", "pollution", "forest density", "other")),
                   param = factor(param, levels = c("recruitment probability", "recruitment rate", "growth", "survival")),
                   var = factor(var, levels = c("dbh", "bacon", "bahet", "nitrogen", "sulfur", "temperature", "precipitation"),
                                labels = c("DBH", "conspecific BA", "heterospecific BA", "nitrogen", "sulfur", "temperature", "precipitation"))) %>%
            group_by(var, param, scaling) %>%
            mutate(ppos = mean(param_value > 0)) %>%
            left_join(nn)

      p <- ds %>% filter(scaling == "spec") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free", space = "free_y") +
            geom_vline(xintercept = 0, color = "gray50", linewidth = .5) +
            geom_boxplot(aes(fill = ppos, weight = 1), alpha = .75, outliers = F) +
            geom_jitter(aes(alpha = credible, color = factor(sign(param_value)), size = n),
                        height = .2, width = 0) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid", "orchid", "gray90", "limegreen", "forestgreen"),
                                 limits = 0:1) +
            scale_color_manual(values = c("darkorchid4", "darkgreen")) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  # strip.text.y = element_text(angle = 0, vjust = 0),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "effect on parameter (predictor-standardized; 1 point per species)",
                 y = NULL)
      ggsave("figures/sensitivity/param_boxplots_species.pdf", p, width = 9, height = 5, units = "in")


      p <- ds %>% filter(scaling == "glob") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free", space = "free_y") +
            geom_vline(xintercept = 0, color = "gray50", size = .5) +
            geom_boxplot(aes(fill = ppos, weight = 1), alpha = .75, outliers = F) +
            geom_jitter(aes(alpha = credible, color = factor(sign(param_value)), size = n),
                        height = .2, width = 0) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid", "orchid", "gray90", "limegreen", "forestgreen"),
                                 limits = 0:1) +
            scale_color_manual(values = c("darkorchid4", "darkgreen")) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  # strip.text.y = element_text(angle = 0, vjust = 0),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "effect on parameter (predictor-standardized; 1 point per species)",
                 y = NULL)
      ggsave("figures/sensitivity/param_boxplots_global.pdf", p, width = 9, height = 5, units = "in")

      p <- ds %>% filter(scaling == "glob",
                         var != "DBH") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free", space = "free_y") +
            geom_vline(xintercept = 0, color = "gray50", size = .5) +
            geom_boxplot(aes(fill = ppos, weight = 1), alpha = .75, outliers = F) +
            geom_jitter(aes(alpha = credible, color = factor(sign(param_value)), size = n),
                        height = .2, width = 0) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid", "orchid", "gray90", "limegreen", "forestgreen"),
                                 limits = 0:1) +
            scale_color_manual(values = c("darkorchid4", "darkgreen")) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  # strip.text.y = element_text(angle = 0, vjust = 0),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "effect of environmental variable on demographic rate (predictor-standardized)",
                 y = NULL)
      ggsave("figures/sensitivity/param_boxplots_global_noDBH.pdf", p, width = 9, height = 5, units = "in")



      # fully standardized including outcome
      oscl <- f %>%
            filter(stat == "response") %>%
            mutate(y2095 = y2009/2 + y2010/2) %>%
            group_by(model, param) %>%
            summarize(mean = mean(y2095, na.rm = T),
                      sd = sd(y2095, na.rm = T))
      oscl2 <- oscl %>%
            mutate(model = str_remove(as.character(model), "\nprobability|\nrate"))
      ds <- s %>%
            group_by(model, species, param, var, scaling) %>%
            summarize(credible = sign(quantile(param_value, .05)) == sign(quantile(param_value, .95)),
                      param_value = median(param_value, na.rm = T)) %>%
            left_join(oscl2) %>%
            mutate(param_value = param_value / sd) %>% # outcome-standardize
            mutate(param = paste(model, param)) %>%
            ungroup() %>%
            filter(var != "intercept") %>%
            mutate(param = str_replace(param, "mu", "rate"),
                   param = str_replace(param, "zeta", "probability"),
                   param = str_replace(param, " beta", ""),
                   param = ifelse(param == "mortality", "survival", param),
                   param_value = ifelse(param == "survival", -param_value, param_value)) %>%
            mutate(var_group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                         var %in% c("bahet", "bacon") ~ "forest density",
                                         var %in% c("sulfur", "nitrogen") ~ "pollution",
                                         TRUE ~ "other"),
                   var_group = factor(var_group, levels = c("climate", "pollution", "forest density", "other")),
                   param = factor(param, levels = c("recruitment probability", "recruitment rate", "growth", "survival")),
                   var = factor(var, levels = c("dbh", "bacon", "bahet", "nitrogen", "sulfur", "temperature", "precipitation"),
                                labels = c("DBH", "conspecific BA", "heterospecific BA", "nitrogen", "sulfur", "temperature", "precipitation"))) %>%
            group_by(var, param, scaling) %>%
            mutate(ppos = mean(param_value > 0)) %>%
            left_join(nn)


      p <- ds %>% filter(scaling == "glob") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free", space = "free") +
            geom_vline(xintercept = 0, color = "gray50", size = .5) +
            geom_boxplot(aes(fill = ppos, weight = 1), alpha = .75, outliers = F) +
            geom_jitter(aes(alpha = credible, color = factor(sign(param_value)), size = n),
                        height = .2, width = 0) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid", "orchid", "gray90", "limegreen", "forestgreen"),
                                 limits = 0:1) +
            scale_color_manual(values = c("darkorchid4", "darkgreen")) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "effect of environmental variable on demographic rate (predictor-standardized)",
                 y = NULL)
      ggsave("figures/sensitivity/param_boxplots_global_std.pdf", p, width = 9, height = 5, units = "in")

      p <- ds %>% filter(scaling == "glob",
                         var != "DBH") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free", space = "free") +
            geom_vline(xintercept = 0, color = "gray50", size = .5) +
            geom_boxplot(aes(fill = ppos, weight = 1), alpha = .75, outliers = F) +
            geom_jitter(aes(alpha = credible, color = factor(sign(param_value)), size = n),
                        height = .2, width = 0) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid", "orchid", "gray90", "limegreen", "forestgreen"),
                                 limits = 0:1) +
            scale_color_manual(values = c("darkorchid4", "darkgreen")) +
            scale_x_continuous(breaks = seq(-6, 6, 2)) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  # strip.text.y = element_text(angle = 0, vjust = 0),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "standardized effect of environmental variable on demographic rate",
                 y = NULL)
      ggsave("figures/sensitivity/param_boxplots_global_std_noDBH.pdf", p, width = 9, height = 5, units = "in")

}


response_plots <- function(f, e){

      select <- dplyr::select

      ll <- e %>% dplyr::select(plot_id, lon, lat) %>% distinct()
      ll <- ll %>% mutate(plot_id = str_sub(plot_id, 1, -3)) %>% bind_rows(ll) %>% distinct()

      pd <- f %>%
            filter(var == "combined") %>%
            left_join(ll, by = join_by(plot_id)) %>%
            mutate(lon = plyr::round_any(lon, .5),
                   lat = plyr::round_any(lat, .5)) %>%
            mutate(model = str_replace(model, "\n", " "),
                   value = ifelse(model == "mortality", -value, value),
                   model = ifelse(model == "mortality", "survival", model))

      pdll <- pd %>%
            group_by(model, lon, lat) %>%
            summarize(pp = mean(value > 0, na.rm = T),
                      value = median(value, na.rm = T), .groups = "drop")
      pdsp <- pd %>%
            group_by(model, species) %>%
            summarize(pp = mean(value > 0, na.rm = T),
                      value = median(value, na.rm = T),
                      lon = mean(lon, na.rm = T),
                      lat = mean(lat, na.rm = T), .groups = "drop")

      maps <- pdll %>%
            ggplot(aes(lon, lat, fill = pp)) +
            facet_grid(model ~ .) +
            geom_raster() +
            scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                 values = c(0, .1, .5, .9, 1), limits = 0:1) +
            theme_bw() +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black")) +
            labs(fill = "modeled rate\nof change\n(units vary\nby parameter)")

      dens <- ggplot() +
            geom_density_ridges_gradient(data = pdsp, aes(pp, 0, fill = stat(x), height = after_stat(density)),
                                         adjust = .5, stat = "density") +
            geom_density_ridges_gradient(data = pdll, aes(pp, 0, fill = stat(x), height = after_stat(density)),
                                         adjust = .5, stat = "density") +
            geom_density(data = pdsp, aes(pp, y = after_stat(density)*1.8), adjust = .5, color = "white") +
            geom_density(data = pdsp, aes(pp, y = after_stat(density)*1.8), adjust = .5, linetype = "dashed") +
            facet_grid(model ~ .) +
            scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                 values = c(0, .1, .5, .9, 1), limits = 0:1) +
            scale_x_continuous(expand = c(0, 0), limits = 0:1, breaks = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1")) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 6), breaks = c(0, 2, 4)) +
            theme_bw() +
            theme(legend.position = "none",
                  # axis.text.y = element_blank(),
                  # axis.ticks.y = element_blank(),
                  panel.grid = element_blank(),
                  strip.text = element_blank()) +
            labs(x = "proportion of trees with modeled\nincrease in demographic rate",
                 y = "relative frequency across species (dashed) and grid cells (solid)")

      p <- dens + maps +
            plot_layout(nrow = 1)
      ggsave("figures/response/ppos_maps_dens.pdf", p, width = 6, height = 7, units = "in")



      # v2

      maps <- pdll %>%
            ggplot(aes(lon, lat, fill = pp)) +
            facet_grid(model ~ .) +
            geom_raster() +
            geom_point(data = pdsp, shape = 21, color = "black") +
            scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                 values = c(0, .1, .5, .9, 1), limits = 0:1) +
            theme_bw() +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black")) +
            labs(fill = "modeled rate\nof change\n(units vary\nby parameter)")

      dens <- bind_rows(pdll %>% mutate(group = "communities"),
                        pdsp %>% mutate(group = "species")) %>%
            mutate(group = factor(group, levels = c("species", "communities"))) %>%
            ggplot(aes(pp, 0, fill = stat(x), height = after_stat(density))) +
            geom_density_ridges_gradient(adjust = .5, stat = "density") +
            geom_vline(xintercept = 0.5, color = "black", alpha = .35) +
            facet_grid(model ~ group) +
            scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                 values = c(0, .1, .5, .9, 1), limits = 0:1) +
            scale_x_continuous(expand = c(0, 0), limits = 0:1, breaks = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1")) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 6), breaks = c(0, 2, 4)) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black"),
                  panel.grid = element_blank(),
                  strip.text.y = element_blank()) +
            labs(x = "proportion of trees with modeled\nincrease in demographic rate",
                 y = "relative frequency across species or communities (grid cells)")

      p <- dens + maps + plot_layout(nrow = 1)
      ggsave("figures/response/ppos_maps_dens_v2.pdf", p, width = 6, height = 7, units = "in")
      ggsave("figures/response/ppos_maps_dens_v2.png", p, width = 6, height = 7, units = "in")




      ## alt version ##

      pd <- f %>% ungroup() %>%
            filter(stat == "response", var != "combined") %>%
            mutate(model = ifelse(str_detect(model, "recruit"), "recruitment", as.character(model))) %>%
            group_by(i, species, model, param, plot_id, tree_id) %>%
            summarize(y2010 = sum(y2010),
                      y2009 = sum(y2009), .groups = "drop")

      # inv_logit <- function(x) 1 / (1 + exp(-x))
      dt <- 1

      pd <- pd %>%
            mutate(y2009 = case_when(param == "mu" ~ exp(y2009),
                                     param == "zeta" ~ inv_logit(y2009),
                                     param == "beta" ~ y2009),
                   y2010 = case_when(param == "mu" ~ exp(y2010),
                                     param == "zeta" ~ inv_logit(y2010),
                                     param == "beta" ~ y2010)) %>%
            group_by(i, species, model, plot_id, tree_id) %>%
            summarize(y2009 = prod(y2009), # multiply mu and zeta to compute overall recruitment
                      y2010 = prod(y2010),
                      value = (y2010 - y2009) / dt,
                      n = n(),
                      .groups = "drop")

      pd <- pd %>%
            ungroup() %>%
            left_join(ll) %>%
            mutate(lon = plyr::round_any(lon, .5),
                   lat = plyr::round_any(lat, .5)) %>%
            mutate(model = str_replace(model, "\n", " "),
                   value = ifelse(model == "mortality", -value, value),
                   model = ifelse(model == "mortality", "survival", model))

      pdll <- pd %>%
            group_by(model, lon, lat) %>%
            summarize(pp = mean(value > 0, na.rm = T),
                      value = median(value, na.rm = T))
      pdsp <- pd %>%
            group_by(model, species) %>%
            summarize(pp = mean(value > 0, na.rm = T),
                      value = median(value, na.rm = T),
                      lon = mean(lon, na.rm = T),
                      lat = mean(lat, na.rm = T))

      response_plot <- function(m = "growth"){

            fill <- scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                         values = c(0, .1, .5, .9, 1), limits = 0:1)

            cmap <- pdll %>% filter(str_detect(model, m)) %>%
                  ggplot(aes(lon, lat, fill = pp)) +
                  geom_tile() +
                  scale_y_continuous(limits = c(min(pdll$lat) - 6, NA)) +
                  fill +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            smap <- pdsp %>% filter(model == m) %>%
                  ggplot(aes(lon, lat, fill = pp)) +
                  geom_tile(data = pdll %>% filter(model == m), fill = "gray80") +
                  geom_point(shape = 21, color = "black", size = 2) +
                  scale_y_continuous(limits = c(min(pdll$lat) - 6, NA)) +
                  fill +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            cdens <- pdll %>% filter(model == m) %>%
                  ggplot(aes(pp, 0, fill = after_stat(x), height = after_stat(density))) +
                  geom_density_ridges_gradient(adjust = .5, stat = "density", linewidth = .25) +
                  # geom_vline(xintercept = 0.5, color = "black", alpha = .35) +
                  fill +
                  scale_x_continuous(expand = c(0, 0), limits = 0:1, breaks = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1")) +
                  # scale_y_continuous(expand = c(0, 0)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        strip.text.y = element_blank()) +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.background = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank())

            sdens <- pdsp %>% filter(model == m) %>%
                  ggplot(aes(pp, 0, fill = after_stat(x), height = after_stat(density))) +
                  geom_density_ridges_gradient(adjust = .5, stat = "density", linewidth = .25) +
                  # geom_vline(xintercept = 0.5, color = "black", alpha = .35) +
                  fill +
                  scale_x_continuous(expand = c(0, 0), limits = 0:1, breaks = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1")) +
                  # scale_y_continuous(expand = c(0, 0)) +
                  theme_bw() +
                  theme(legend.position = "none",
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        strip.text.y = element_blank()) +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.background = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title = element_blank())


            (cmap + inset_element(cdens, left = .02, right = .98, bottom = 0, top = .28)) /
                  (smap + inset_element(sdens, left = .02, right = .98, bottom = 0, top = .28))
      }

      p <- c("recruitment", "growth", "survival") %>%
            map(response_plot) %>%
            Reduce("|", .)
      ggsave("figures/response/response_maps_dens.pdf", p, width = 8, height = 4.5, units = "in")







      response_plot_native <- function(m = "growth"){

            clamp <- function(x, mn, mx) pmax(mn, pmin(mx, x))
            # y <- na.omit(pd$value[pd$model == m])
            # q <- .05
            # y <- clamp(y, quantile(y, q), quantile(y, 1-q))
            # lims <- max(abs(y))*c(-1,1)

            lims <- switch(m,
                           "recruitment" = .1,
                           "growth" = .0035,
                           "survival" = .15) * c(-1, 1)
            brks <- switch(m,
                           "recruitment" = .08,
                           "growth" = .0025,
                           "survival" = .1) * c(-1, 0, 1)

            fill <- scale_fill_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                         values = c(0, .35, .5, .65, 1),
                                         limits = lims)
            color <- scale_color_gradientn(colors = c("darkorchid4", "darkmagenta", "gray", "forestgreen", "darkgreen"),
                                           values = c(0, .35, .5, .65, 1),
                                           limits = lims)

            cmap <- pdll %>% filter(str_detect(model, m)) %>%
                  mutate(value = clamp(value, lims[1], lims[2])) %>%
                  ggplot(aes(lon, lat, fill = value)) +
                  geom_tile() +
                  scale_y_continuous(limits = c(min(pdll$lat) - 6, NA)) +
                  fill +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"))

            # sdens <- pd %>%
            #       filter(model == m) %>%
            #       mutate(value = clamp(value, lims[1], lims[2])) %>%
            #       group_by(species) %>%
            #       mutate(median = median(value, na.rm = T)) %>%
            #       ungroup() %>%
            #       arrange(median) %>%
            #       mutate(species = factor(species, levels = unique(species))) %>%
            #
            #       ggplot(aes(value, species, fill = after_stat(x), height = after_stat(scaled))) +
            #       geom_density_ridges_gradient(adjust = .5, stat = "density",
            #                                    color = "black", linewidth = .1) +
            #       geom_vline(xintercept = 0, color = "black") +
            #       fill +
            #       scale_x_continuous(expand = c(0, 0), limits = lims, breaks = brks) +
            #       scale_y_discrete(expand = expansion(mult = c(.025, .05))) +
            #       theme_bw() +
            #       theme(legend.position = "none",
            #             panel.grid = element_blank(),
            #             panel.background = element_blank(),
            #             plot.background = element_blank(),
            #             axis.title.y = element_blank(),
            #             axis.text.y = element_blank(),
            #             axis.ticks.y = element_blank(),
            #             strip.text = element_text(color = "white"),
            #             strip.background = element_rect(fill = "black")) +
            #       labs(x = switch(m,
            #                       "recruitment" = "N / m2 / ha / yr",
            #                       "growth" = "sqrt(%AGR) / yr",
            #                       "survival" = "logit(probability) / yr"))

            library(ggforce)

            sdens <- pd %>%
                  filter(model == m) %>%
                  mutate(value = clamp(value, lims[1], lims[2])) %>%
                  group_by(species) %>%
                  summarize(median = median(value, na.rm = T),
                            q95 = quantile(value, .95, na.rm = T),
                            q25 = quantile(value, .25, na.rm = T),
                            q75 = quantile(value, .75, na.rm = T),
                            q05 = quantile(value, .05, na.rm = T)) %>%
                  ungroup() %>%
                  arrange(median) %>%
                  mutate(species = factor(species, levels = unique(species))) %>%
                  ggplot(aes(color = after_stat(x))) +
                  # geom_link(aes(x = q05, xend = q95, y = species, yend = species),
                  #           lineend = 'round', n = 500, linewidth = .25) +
                  geom_link(aes(x = q05, xend = q95, y = species, yend = species),
                            n = 20, linewidth = .7, alpha = .3) +
                  geom_link(aes(x = q25, xend = q75, y = species, yend = species),
                            n = 10, linewidth = .7) +
                  geom_point(aes(median, species), size = .25, color = "black") +
                  geom_vline(xintercept = 0, color = "black") +

                  color +
                  scale_x_continuous(expand = c(0, 0), limits = lims, breaks = brks) +
                  # scale_y_discrete(expand = expansion(mult = c(.025, .05))) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black")) +
                  labs(x = switch(m,
                                  "recruitment" = "N / m2 / ha / yr",
                                  "growth" = "sqrt(%AGR) / yr",
                                  "survival" = "logit(probability) / yr"))

            cdens <- pd %>%
                  filter(model == m) %>%
                  mutate(value = clamp(value, lims[1], lims[2])) %>%
                  group_by(lon, lat) %>%
                  summarize(median = median(value, na.rm = T),
                            q95 = quantile(value, .95, na.rm = T),
                            q25 = quantile(value, .25, na.rm = T),
                            q75 = quantile(value, .75, na.rm = T),
                            q05 = quantile(value, .05, na.rm = T)) %>%
                  ungroup() %>%
                  arrange(median) %>%
                  mutate(cell = paste(lon, lat),
                         cell = factor(cell, levels = unique(cell))) %>%
                  ggplot(aes(color = after_stat(x))) +
                  # geom_link(aes(x = q05, xend = q95, y = cell, yend = cell),
                  #           n = 20, linewidth = .7, alpha = .3) +
                  geom_link(aes(x = q05, xend = q95, y = cell, yend = cell),
                            n = 100, linewidth = .15) +
                  geom_point(aes(median, cell), size = .1, color = "black") +
                  geom_vline(xintercept = 0, color = "black") +
                  color +
                  scale_x_continuous(expand = c(0, 0), limits = lims, breaks = brks) +
                  # scale_y_discrete(expand = expansion(mult = c(.025, .05))) +
                  theme_bw() +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black")) +
                  labs(x = switch(m,
                                  "recruitment" = "N / m2 / ha / yr",
                                  "growth" = "sqrt(%AGR) / yr",
                                  "survival" = "logit(probability) / yr"))

            cmap <- cmap + scale_y_continuous(limits = c(NA, NA))
            cmap + cdens + sdens + plot_layout(heights = c(1, 1.25, 1.25), ncol = 1)

      }

      p <- c("recruitment", "growth", "survival") %>%
            map(response_plot_native) %>%
            Reduce("|", .)
      ggsave("figures/response/response_maps_dens_native.pdf", p, width = 8, height = 7, units = "in")



      # compare changes in the three vital rates

      pdg <- bind_rows(pdll %>% ungroup() %>%
                             mutate(group = paste(lon, lat), level = "community") %>%
                             select(model, group, value, level),
                       pdsp %>% ungroup() %>%
                             mutate(group = species, level = "species") %>%
                             select(model, group, value, level))

      pdg1 <- pdg %>%
            group_by(model) %>%
            mutate(value = (rank(value)-1)/(length(value)-1)) %>%
            spread(model, value) %>%
            ecoclim::pairsData(xy_vars = c("growth", "recruitment", "survival"),
                               z_vars = c("group", "level"),
                               mirror = TRUE)
      pdg2 <- pdg %>%
            group_by(model, level) %>%
            summarize(value = ecdf(value)(0), .groups = "drop") %>%
            spread(model, value) %>%
            ecoclim::pairsData(xy_vars = c("growth", "recruitment", "survival"),
                               z_vars = c("level"),
                               mirror = TRUE)


      p <- ggplot(mapping = aes(x_value, y_value, color = level, fill = level,
                                size = level)) +
            facet_grid(y_var ~ x_var) +
            geom_hline(data = pdg2, aes(yintercept = y_value, color = level)) +
            geom_vline(data = pdg2, aes(xintercept = x_value, color = level)) +
            geom_point(data = pdg1) +
            geom_smooth(data = pdg1, linewidth = 1, method = lm) +
            scale_size_manual(values = c(.1, 1)) +
            scale_color_manual(values = c("orangered", "darkblue")) +
            scale_fill_manual(values = c("orangered", "darkblue")) +
            coord_fixed() +
            theme_bw() +
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black"),
                  legend.position = "bottom") +
            labs(x = "modeled rate of change in demographic rate listed at top (rank)",
                 y = "modeled rate of chagne in demographic rate listed at right (rank)",
                 color = NULL, fill = NULL, size = NULL)
      ggsave("figures/response/response_pairs.pdf", p, width = 6, height = 6.5, units = "in")

      # correlations at various levels (quantifying the pattern in the above plot)
      pdg <- pd %>%
            mutate(plot_id = ifelse(model == "recruitment", plot_id,
                                    str_sub(plot_id, 1, -3)),
                   cell = paste(plyr::round_any(lon, 3),
                                plyr::round_any(lat, 3))) %>%
            group_by(species, model, plot_id, cell) %>%
            summarize(value = weighted.mean(value, n, na.rm = T),
                      n = sum(n, na.rm = T),
                      .groups = "drop") %>%
            ungroup()
      pdg <- bind_rows(pdg %>% mutate(level = "spp-plt") %>%
                             group_by(species, plot_id) %>% mutate(n = sum(n, na.rm = T)),
                       pdg %>% group_by(species, model) %>%
                             summarize(value = weighted.mean(value, n, na.rm = T),
                                       n = sum(n, na.rm = T),
                                       level = "spp") %>%
                             group_by(species) %>% mutate(n = sum(n)),
                       pdg %>% group_by(cell, model) %>%
                             summarize(value = weighted.mean(value, n, na.rm = T),
                                       n = sum(n, na.rm = T),
                                       level = "cell") %>%
                             group_by(cell) %>% mutate(n = sum(n)),
                       pdg %>% group_by(plot_id, model) %>%
                             summarize(value = weighted.mean(value, n, na.rm = T),
                                       n = sum(n, na.rm = T),
                                       level = "plt") %>%
                             group_by(plot_id) %>% mutate(n = sum(n)))
      r <- pdg %>%
            spread(model, value) %>%
            filter(is.finite(recruitment),
                   is.finite(growth),
                   is.finite(survival)) %>%
            group_by(level) %>%
            summarize(pearson_gr = weightedCorr(growth, recruitment, method = "Pearson", weights = n),
                      pearson_gs = weightedCorr(growth, survival, method = "Pearson", weights = n),
                      pearson_sr = weightedCorr(survival, recruitment, method = "Pearson", weights = n),
                      spearman_gr = weightedCorr(growth, recruitment, method = "Spearman", weights = n),
                      spearman_gs = weightedCorr(growth, survival, method = "Spearman", weights = n),
                      spearman_sr = weightedCorr(survival, recruitment, method = "Spearman", weights = n))



      ### demographic compensation ------------------------

      ff <- f %>%
            ungroup() %>%
            filter(stat == "response") %>%
            mutate(value = ifelse(model == "mortality", -value, value)) %>%
            group_by(model) %>%
            mutate(value = value / sd(value, na.rm = T)) %>%
            group_by(i, species, plot_id, stat, tree_id, var) %>%
            summarize(sum = abs(sum(value)),
                      magnitude = sum(abs(value)),
                      stdev = sd(value),
                      n = n()) %>%
            mutate(additivity = sum / magnitude,
                   cv = stdev / (sum / n))

      p1 <- ff %>%
            ggplot(aes(additivity, color = var, fill = var)) +
            geom_density(alpha = .1, adjust = 2) +
            # scale_color_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            # scale_fill_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            scale_x_continuous(expand = c(0, 0), limits = 0:1) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
            labs(x = "response additivity: abs(sum(x)) / sum(abs(x))",
                 y = "relative frequency across trees") +
            theme_bw()
      p2 <- ff %>%
            filter(is.finite(cv)) %>%
            filter(between(cv, quantile(cv, .001), quantile(cv, .99))) %>%
            ggplot(aes(cv, color = var, fill = var)) +
            geom_density(alpha = .1, adjust = 2) +
            scale_x_log10(expand = c(0, 0)) +
            # scale_color_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            # scale_fill_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            labs(x = "response coefficient of variation: sd(x) / abs(mean(x))",
                 y = "relative frequency across trees") +
            theme_bw()

      p <- p1 / p2 + plot_layout(guides = "collect")
      ggsave("figures/response/rate_compensation_curves.pdf", p, width = 7, height = 6, units = "in")


      ### variable additivity --------------------

      ff <- f %>%
            filter(stat == "response",
                   var != "combined") %>%
            group_by(i, species, model, plot_id, stat, tree_id) %>%
            reframe(value = c(abs(sum(value)), sum(abs(value)), sd(value)),
                    var = c("sum", "magnitude", "stdev"),
                    n = n()) %>%
            select(i:n) %>%
            spread(var, value) %>%
            mutate(additivity = sum / magnitude,
                   cv = stdev / (sum / n),
                   model = str_replace(model, "\n", " "))

      p1 <- ff %>%
            ggplot(aes(additivity, color = model, fill = model)) +
            geom_density(alpha = .1) +
            scale_color_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            scale_fill_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            scale_x_continuous(expand = c(0, 0), limits = 0:1) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
            labs(x = "response additivity: abs(sum(x)) / sum(abs(x))",
                 y = "relative frequency across trees") +
            theme_bw()
      p2 <- ff %>%
            filter(is.finite(cv)) %>%
            filter(between(cv, quantile(cv, .001), quantile(cv, .99))) %>%
            ggplot(aes(cv, color = model, fill = model)) +
            geom_density(alpha = .1) +
            scale_x_log10(expand = c(0, 0)) +
            scale_color_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            scale_fill_manual(values = c("darkgreen", "orangered", "blue", "dodgerblue")) +
            labs(x = "response coefficient of variation: sd(x) / abs(mean(x))",
                 y = "relative frequency across trees") +
            theme_bw()

      library(patchwork)
      p <- p1 / p2 + plot_layout(guides = "collect")
      ggsave("figures/response/response_additivity_curves.pdf", p, width = 7, height = 6, units = "in")

}


importance_plots <- function(f, species){

      select <- dplyr::select
      scl <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)

      ### combined importance -----------------

      am <- f %>%
            filter(var != "combined") %>%
            group_by(species, model) %>%
            mutate(sp_wt = 1 / n()) %>%
            group_by(stat, model, var) %>%
            summarize(mav = mean(abs(value), na.rm = T),
                      mav_wt = weighted.mean(abs(value), w = sp_wt, na.rm = T),
                      amv = abs(mean(value, na.rm = T)))

      p <- am %>%
            mutate(category = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                        var %in% c("nitrogen", "sulfur") ~ "pollutant deposition",
                                        TRUE ~ "forest density"),
                   var = factor(var, levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet"),
                                labels = c("temperature", "precipitation", "nitrogen", "sulfur",
                                           "conspecific\nbasal area", "heterospecific\nbasal area"))) %>%
            gather(metric, value, mav, mav_wt, amv) %>%
            group_by(metric, model, stat) %>%
            mutate(value = value / max(value)) %>%
            spread(metric, value) %>%
            ggplot(aes(var, mav, fill = stat)) +
            facet_grid(model ~ category, scales = "free") +
            geom_bar(stat = "identity", position = "dodge",
                     width = .5, color = "white", alpha = .75) +
            geom_segment(aes(y = 0, yend = mav_wt), position = position_dodge(width = .5),
                         color = "black") +
            geom_point(aes(y = mav_wt), position = position_dodge(width = .5),
                       shape = 21, color = "black", size = 3) +
            scale_fill_manual(values = c("orangered", "dodgerblue", "darkorchid4")) +
            scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0),
                               breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
            theme_bw() +
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = "black"),
                  axis.ticks.x = element_blank(),
                  panel.grid = element_blank()) +
            labs(x = NULL,
                 y = "relative strength",
                 fill = NULL)
      ggsave("figures/importance/importance_bars_species.pdf", p, width = 8, height = 5, units = "in")



      p <- f %>%
            group_by(stat, model) %>%
            filter(is.finite(value)) %>%
            filter(between(value, quantile(value, .001), quantile(value, .999))) %>%
            group_by(stat, model, var) %>%
            mutate(mean = mean(value)) %>%
            ggplot() +
            facet_grid(var ~ model + stat, scales = "free") +
            geom_vline(xintercept = 0, color = "gray40", size = .5) +
            # geom_density(color = NA, alpha = .5, adjust = 1) +
            geom_histogram(aes(value,# after_stat(scaled),
                               after_stat(ncount),
                               fill = stat),
                           bins = 50, boundary = 0, alpha = .5) +
            geom_vline(aes(xintercept = mean, color = stat), size = .5) +
            scale_color_manual(values = c("orangered", "dodgerblue", "purple")) +
            scale_fill_manual(values = c("orangered", "dodgerblue", "purple")) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = NA)) +
            labs(y = "relative frequency",
                 x = "standardized value")
      ggsave("figures/importance/param_sens_densities.pdf", p, width = 13, height = 8, units = "in")



      # alternative version, with exposure and sensitivity scaled globally

      am <- f %>%
            select(-plot_id, -tree_id, -t, -y2009, -y2010, -delta, -param) %>%
            filter(var != "combined") %>%
            mutate(var = as.character(var)) %>%
            left_join(scl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"), value * sd, value)) %>%
            select(-sd, -mean) %>%
            left_join(gscl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"), value / sd, value)) %>%
            select(-sd, -mean) %>%
            group_by(species, model) %>%
            mutate(sp_wt = 1 / n()) %>%
            group_by(stat, model, var) %>%
            summarize(mav = mean(abs(value), na.rm = T),
                      mav_wt = weighted.mean(abs(value), w = sp_wt, na.rm = T),
                      amv = abs(mean(value, na.rm = T)))

      p <- am %>%
            mutate(category = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                        var %in% c("nitrogen", "sulfur") ~ "pollutant deposition",
                                        TRUE ~ "forest density"),
                   var = factor(var, levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet"),
                                labels = c("temperature", "precipitation", "nitrogen", "sulfur",
                                           "conspecific\nbasal area", "heterospecific\nbasal area"))) %>%
            gather(metric, value, mav, mav_wt, amv) %>%
            group_by(metric, model, stat) %>%
            mutate(value = value / max(value)) %>%
            spread(metric, value) %>%
            ggplot(aes(var, mav, fill = stat)) +
            facet_grid(model ~ category, scales = "free") +
            geom_bar(stat = "identity", position = "dodge",
                     width = .5, color = "white", alpha = .75) +
            geom_segment(aes(y = 0, yend = mav_wt), position = position_dodge(width = .5),
                         color = "black") +
            geom_point(aes(y = mav_wt), position = position_dodge(width = .5),
                       shape = 21, color = "black", size = 3) +
            scale_fill_manual(values = c("orangered", "dodgerblue", "darkorchid4")) +
            scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0),
                               breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
            theme_bw() +
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = "black"),
                  axis.ticks.x = element_blank(),
                  panel.grid = element_blank()) +
            labs(x = NULL,
                 y = "relative strength",
                 fill = NULL)
      ggsave("figures/importance/importance_bars_global.pdf", p, width = 8, height = 5, units = "in")



      # variance partitioning  -----------------------------------------------

      am <- f %>%
            select(#-plot_id, -tree_id,
                  -t, -y2009, -y2010, -delta, -param) %>%
            filter(var != "combined") %>%
            filter(stat %in% c("exposure", "sensitivity")) %>%
            mutate(var = as.character(var)) %>%
            left_join(scl) %>%
            mutate(value = value * sd) %>%
            select(-sd, -mean) %>%
            left_join(gscl) %>%
            mutate(value = value / sd) %>%
            select(-sd, -mean) %>%
            unite(var, stat, var)
      vars <- unique(am$var)
      z <- filter(am, var == vars[1]) %>%
            rename(!!vars[1] := value) %>%
            select(-var)
      for(v in vars[2:length(vars)]) z <- z %>%
            left_join(filter(am, var == v) %>%
                            rename(!!v := value) %>%
                            select(-var),
                      by = join_by(i, species, model, plot_id, tree_id))
      z <- f %>%
            select(-t, -y2009, -y2010, -delta, -param) %>%
            filter(var == "combined") %>%
            left_join(z) %>%
            rename(response = value) %>%
            select(-var)

      vp <- unique(z$model) %>%
            map(function(m){
                  zz <- z %>% filter(model == m) %>% sample_n(50000)
                  vp <- rdacca.hp::rdacca.hp(zz$response,
                                             select(zz, exposure_bacon:sensitivity_precipitation))
                  vp$Hier.part %>%
                        as.data.frame() %>%
                        mutate(total_explained = vp$Total_explained_variation,
                               model = m) %>%
                        rownames_to_column("var")
            })

      vps <- vp %>%
            bind_rows() %>%
            separate(var, c("stat", "var"), sep = "_") %>%
            janitor::clean_names() %>%
            mutate(percent = i_perc_percent * total_explained) %>%
            mutate(group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                     var %in% c("bacon", "bahet") ~ "forest density",
                                     var %in% c("nitrogen", "sulfur") ~ "pollution"))
      vps <- vps %>%
            filter(var == "bacon", stat == "exposure") %>%
            mutate(var = "unexplained", stat = NA, group = "",
                   percent = 100 * (1 - total_explained)) %>%
            bind_rows(vps)

      p <- vps %>%
            ggplot(aes(var, percent/100, fill = stat)) +
            facet_grid(group~model, scales = "free_y", space = "free_y") +
            geom_bar(stat = "identity", position = "stack") +
            coord_flip() +
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = c("orangered", "dodgerblue")) +
            theme_bw() +
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = "black"),
                  axis.title.y = element_blank()) +
            labs(y = "portion of response variation explained",
                 fill = NULL)
      ggsave("figures/importance/variance_paritioned.pdf",
             p, width = 8, height = 5, units = "in")

      p <- vps %>%
            arrange(model) %>%
            group_by(model) %>%
            mutate(var = factor(var, levels = c("temperature", "precipitation",
                                                "bacon", "bahet", "sulfur", "nitrogen",
                                                "unexplained")),
                   group = factor(group, levels = (c("climate", "forest density", "pollution", "")))) %>%
            arrange(model, var, stat) %>%
            mutate(pcum = cumsum(percent)) %>%
            ggplot(aes(x = var, y = (pcum - percent/2)/100,
                       width = 1, height = percent/100,
                       fill = stat)) +
            facet_grid(model ~ group,
                       scales = "free_x", space = "free_x") +
            geom_tile() +
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = c("orangered", "dodgerblue")) +
            theme_bw() +
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = "black"),
                  axis.title.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid = element_line(linewidth = .25),
                  legend.position = "right") +
            labs(y = "cumulative portion of modeled response variation explained",
                 fill = NULL)
      ggsave("figures/importance/variance_paritioned_cumulative.pdf",
             p, width = 8, height = 8, units = "in")

}
