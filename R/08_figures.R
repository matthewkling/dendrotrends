
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

# scales for outcome variables
outcome_scale <- function(recr_data, grow_data, mort_data){
      g <- bind_rows(
            grow_data %>% summarize(model = "growth", mean = mean(outcome), sd = sd(outcome)),
            mort_data %>% summarize(model = "mortality", mean = mean(outcome), sd = sd(outcome)),
            recr_data %>% reframe(model = c("recruitment\nprobability", "recruitment\nrate"),
                                  mean = c(mean(outcome > 0), mean(log(outcome[outcome > 0]))),
                                  sd = c(sd(outcome > 0), sd(log(outcome[outcome > 0]))))
      )

      s <- bind_rows(
            grow_data %>% group_by(species) %>% summarize(model = "growth", mean = mean(outcome), sd = sd(outcome)),
            mort_data %>% group_by(species) %>% summarize(model = "mortality", mean = mean(outcome), sd = sd(outcome)),
            recr_data %>% group_by(species) %>% reframe(model = c("recruitment\nprobability", "recruitment\nrate"),
                                                        mean = c(mean(outcome > 0), mean(log(outcome[outcome > 0]))),
                                                        sd = c(sd(outcome > 0), sd(log(outcome[outcome > 0]))))
      )

      return(list(global = g, species = s))
}


exposure_plots <- function(f, e, species){

      select <- dplyr::select
      ll <- select(e, plot_id, lon, lat) %>% distinct()

      factor_var <- function(x){
            mutate(x,
                   var = factor(var,
                                levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet"),
                                labels = c("temperature", "precipitation", "nitrogen\ndeposition", "sulfur\ndeposition",
                                           "conspecific\nforest density", "heterospecific\nforest density")))
      }


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


      exposure_plot <- function(v = "temperature"){

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
            colors2 <- scale_color_gradientn(colors = pal,
                                             values = c(0, .44, .48, .5, .52, .56, 1),
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
                          "nitrogen\ndeposition" = c(-.04, -.02, 0, .02),
                          "sulfur\ndeposition" = c(-.2, -.1, 0, .1),
                          "conspecific\nforest density" = c(-.1, 0, .1),
                          "heterospecific\nforest density" = c(-.1, 0, .1)) * 10

            stdev <- mean(ev$rate / ev$rate_std, na.rm = T)
            xscl <- scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax),
                                       breaks = brk / stdev, labels = brk)
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

            res <- .1
            x <- round(seq(plyr::round_any(lim[1], res), plyr::round_any(lim[2], res), res), 1)
            x <- x[x >= xmin & x <= xmax]
            sg_minor <- data.frame(x = x)
            sg_major <- data.frame(x = x[x == round(x)])

            dens <- ev %>%
                  ggplot(aes(value, 0, fill = after_stat(x), height = after_stat(scaled))) +
                  geom_vline(data = sg_minor, aes(xintercept = x, color = x), linewidth = .15) +
                  geom_vline(data = sg_major, aes(xintercept = x, color = x), linewidth = .5) +
                  geom_density_ridges_gradient(adjust = .5, linewidth = .15,
                                               stat = "density", trim = F) +
                  geom_vline(xintercept = 0, color = "black", linewidth = .75) +
                  geom_hline(yintercept = 0, color = "black", linewidth = .15) +
                  colors + colors2 +
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

            map + inset_element(dens, left = 0, right = 1, bottom = .07, top = .38, on_top = F, clip = F)
      }

      p <- levels(ell$var)[c(1, 3, 5, 2, 4, 6)] %>%
            map(exposure_plot) %>%
            Reduce("+", .)
      ggsave("figures/exposure_maps.png", p, width = 8, height = 4.5, units = "in")
}


sensitivity_plots <- function(f, species, recr_draws, grow_draws, mort_draws, oscl){

      require(ggdist)

      select <- dplyr::select

            s <- bind_rows(load_fits(get_draws(grow_draws)) %>% mutate(model = "growth"),
                     load_fits(get_draws(mort_draws)) %>% mutate(model = "mortality"),
                     load_fits(get_draws(recr_draws)) %>% mutate(model = "recruitment"))


      nn <- f %>% ungroup() %>% count(species)

      scl <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)

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

      oscl2 <- oscl$global %>%
            mutate(param = case_when(model %in% c("growth", "mortality") ~ "beta",
                                     model == "recruitment\nprobability" ~ "zeta",
                                     model == "recruitment\nrate" ~ "mu"),
                   model = ifelse(str_detect(model, "recruitment"), "recruitment", model))

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

      pd <- ds %>%
            ungroup() %>%
            filter(scaling == "glob") %>%
            mutate(ppos = (ppos * 2 - 1) * max(abs(param_value)))

      pdm <- pd %>% filter(var != "DBH") %>% group_by(var_group, var, param) %>%
            summarize(param_value = weighted.mean(param_value, n))

      p <- pd %>% filter(var != "DBH") %>%
            ggplot(aes(param_value, var)) +
            facet_grid(var_group ~ param, scales = "free") +
            geom_vline(xintercept = 0, color = "black", size = .25) +
            stat_dots(aes(fill = after_stat(x)), linewidth = NA, layout = "swarm", side = "both") +
            geom_point(data = pdm, color = "black", size = 2) +
            scale_alpha_manual(values = c(.25, 1)) +
            scale_size_continuous(range = c(.25, 2)) +
            scale_fill_gradientn(colors = c("darkorchid4", "darkorchid4", "orchid", "gray70", "limegreen", "darkgreen", "darkgreen"),
                                 values = c(0, .3, .45, .5, .55, .7, 1),
                                 limits = max(abs(pd$param_value)) * c(-1, 1)) +
            scale_y_discrete(expand = c(0, 0)) +
            theme_bw() +
            theme(legend.position = "none",
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  panel.grid = element_line(color = "gray95")) +
            labs(x = "standarized effect of environmental variable on demographic rate",
                 y = NULL)
      ggsave("figures/sensitivity_swarm_global_std_noDBH.pdf",
             p, width = 9, height = 5, units = "in")
}


response_plots <- function(f, e){

      select <- dplyr::select

      ll <- e %>% dplyr::select(plot_id, lon, lat) %>% distinct()
      ll <- ll %>% mutate(plot_id = str_sub(plot_id, 1, -3)) %>% bind_rows(ll) %>% distinct()

      pd <- f %>% ungroup() %>%
            filter(stat == "response", var != "combined") %>%
            mutate(model = ifelse(str_detect(model, "recruit"), "recruitment", as.character(model))) %>%
            group_by(i, species, model, param, plot_id, tree_id) %>%
            summarize(y2010 = sum(y2010),
                      y2009 = sum(y2009), .groups = "drop")

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

      response_plot <- function(m = "growth"){

            clamp <- function(x, mn, mx) pmax(mn, pmin(mx, x))

            lims <- switch(m,
                           "recruitment" = .1,
                           "growth" = .0035,
                           "survival" = .18) * c(-1, 1)
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

            cmap <- pdll %>%
                  filter(str_detect(model, m)) %>%
                  mutate(value = clamp(value, lims[1], lims[2])) %>%
                  ggplot(aes(lon, lat, fill = value)) +
                  geom_tile() +
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

            sd <- pd %>%
                  filter(model == m) %>%
                  mutate(value = clamp(value, lims[1], lims[2])) %>%
                  group_by(species) %>%
                  summarize(value = median(value, na.rm = T),
                            lon = mean(lon),
                            lat = mean(lat))
            smap <- ggplot() +
                  geom_raster(data = distinct(select(pdll, lon, lat)),
                              aes(lon, lat), fill = "gray90") +
                  geom_point(data = sd,
                             aes(lon, lat, color = value)) +
                  color +
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


            curves <- function(x){

                  x <- x %>%
                        mutate(value = ifelse(value %in% lims, NA, value)) %>%
                        group_by(group) %>%
                        summarize(q50 = median(value, na.rm = T),
                                  q25 = quantile(value, .25, na.rm = T),
                                  q75 = quantile(value, .75, na.rm = T)) %>%
                        gather(stat, value, -group)

                  ggplot() +
                        geom_density_ridges_gradient(
                              data = filter(x, stat == "q50"),
                              aes(value, 0, fill = after_stat(x), height = after_stat(scaled)),
                              stat = "density", linewidth = 0) +
                        geom_vline(xintercept = 0, color = "gray40", linewidth = .25) +
                        geom_density_ridges_gradient(
                              data = filter(x, stat != "q50"),
                              aes(value, 0, height = after_stat(scaled), group = stat,
                                  linetype = stat),
                              stat = "density", fill = NA, linewidth = .5) +
                        scale_linetype_manual(values = c(1, 6)) +
                        fill +
                        scale_x_continuous(expand = c(0, 0), limits = lims, breaks = brks) +
                        scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                        theme_bw() +
                        theme(legend.position = "none",
                              panel.grid = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              panel.background = element_blank(),
                              panel.border = element_blank(),
                              plot.background = element_blank(),
                              strip.text = element_text(color = "white"),
                              strip.background = element_rect(fill = "black")) +
                        labs(x = switch(m,
                                        "recruitment" = "N / m2 / ha / yr",
                                        "growth" = "sqrt(%AGR) / yr",
                                        "survival" = "logit(probability) / yr"))
            }

            cdens <- pd %>%
                  filter(model == m) %>%
                  mutate(value = clamp(value, lims[1], lims[2]),
                         group = paste(lon, lat)) %>%
                  curves()

            sdens <- pd %>%
                  filter(model == m) %>%
                  mutate(value = clamp(value, lims[1], lims[2]),
                         group = species) %>%
                  curves()

            (cmap + ylim(5, NA) + inset_element(cdens, 0, .03, 1, .5, align_to = 'full')) +
                  (smap + ylim(5, NA) + inset_element(sdens, 0, .03, 1, .5, align_to = 'full')) +
                  plot_layout(ncol = 1)
      }

      p <- response_plot("recruitment") | response_plot("growth") | response_plot("survival")
      ggsave("figures/response_maps_dens.pdf", p, width = 8, height = 6, units = "in")
}


importance_plots <- function(f, species, oscl){

      select <- dplyr::select
      scl <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)


      # global scaling
      goscl <- oscl$global
      amg <- f %>%
            select(-plot_id, -tree_id, -t, -y2009, -y2010, -delta, -param) %>%
            filter(var != "combined") %>%
            mutate(var = as.character(var)) %>%
            left_join(scl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"), value * sd, value)) %>%
            select(-sd, -mean) %>%
            left_join(gscl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"), value / sd, value)) %>%
            select(-sd, -mean) %>%
            left_join(goscl) %>%
            mutate(value = ifelse(stat %in% c("response", "sensitivity"), value / sd, value)) %>%
            group_by(species, model) %>%
            mutate(sp_wt = 1 / n()) %>%
            group_by(stat, model, var) %>%
            summarize(mav = mean(abs(value), na.rm = T),
                      mav_wt = weighted.mean(abs(value), w = sp_wt, na.rm = T),
                      amv = abs(mean(value, na.rm = T)))

      # species scaling
      goscl <- oscl$species
      ams <- f %>%
            select(-plot_id, -tree_id, -t, -y2009, -y2010, -delta, -param) %>%
            filter(var != "combined") %>%
            mutate(var = as.character(var)) %>%
            left_join(goscl) %>%
            mutate(value = ifelse(stat %in% c("response", "sensitivity"), value / sd, value)) %>%
            group_by(species, model) %>%
            mutate(sp_wt = 1 / n()) %>%
            group_by(stat, model, var) %>%
            summarize(mav = mean(abs(value), na.rm = T),
                      mav_wt = weighted.mean(abs(value), w = sp_wt, na.rm = T),
                      amv = abs(mean(value, na.rm = T)))

      am <- bind_rows(mutate(amg, scaling = "global"),
                      mutate(ams, scaling = "species")) %>%
            mutate(category = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                        var %in% c("nitrogen", "sulfur") ~ "pollutant deposition",
                                        TRUE ~ "forest density"),
                   model = factor(model, levels = c("recruitment\nprobability", "recruitment\nrate", "growth", "mortality"),
                                  labels = c("recruitment prob.", "recruitment rate", "growth", "survival")),
                   var = factor(var, levels = c("temperature", "precipitation", "nitrogen", "sulfur", "bacon", "bahet"),
                                labels = c("temperature", "precipitation", "nitrogen", "sulfur",
                                           "conspec. BA", "heterospec. BA"))) %>%
            gather(metric, value, mav, mav_wt, amv) %>%
            group_by(metric, stat, scaling) %>%
            mutate(value = value / max(value, na.rm = T),
                   value = pmax(value, .01))

      imp <- am %>%
            filter(metric == "mav", scaling == "global")


      # variance partitioning  -----------------------------------------------

      am <- f %>%
            select(-t, -y2009, -y2010, -delta, -param) %>%
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
            arrange(model) %>%
            group_by(model) %>%
            mutate(var = factor(var, levels = c("temperature", "precipitation",
                                                "bacon", "bahet", "sulfur", "nitrogen",
                                                "unexplained"),
                                labels = c("temperature", "precipitation",
                                           "conspec. BA", "heterospec. BA", "sulfur", "nitrogen",
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
                  legend.position = "none") +
            labs(y = "cumulative portion of modeled response variation explained",
                 fill = NULL)
      ggsave("figures/importance_vp_cumulative.pdf",
             p, width = 8, height = 8, units = "in")



      # combined plot ----------------------------

      v <- vps %>%
            arrange(model) %>%
            group_by(model) %>%
            mutate(var = factor(var,
                               levels = c("temperature", "precipitation", "bacon", "bahet", "nitrogen", "sulfur", "unexplained"),
                               labels = c("temperature", "precipitation", "conspec. BA", "heterospec. BA", "nitrogen", "sulfur", "unexplained")),
                  group = factor(group, levels = (c("climate", "pollution", "forest density", "")))) %>%
            select(model, stat, group, var, value = percent) %>%
            group_by(model) %>%
            mutate(value = value / sum(value),
                   value = pmax(value, .001)) %>%
            filter(var != "unexplained") %>%
            mutate(model = factor(model,
                                  levels = c("recruitment\nprobability", "recruitment\nrate", "growth", "mortality"),
                                  labels = c("recruitment\nprobability", "recruitment\nrate", "growth", "survival"))) %>%
            mutate(label = str_sub(round(value*100), 1, 2),
                   label_color = ifelse(value < .005, "white", "black"))

      i <- imp %>%
            ungroup() %>%
            mutate(var_group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                         str_detect(var, "BA") ~ "forest density",
                                         var %in% c("sulfur", "nitrogen") ~ "pollution",
                                         TRUE ~ "other"),
                   group = factor(var_group, levels = c("climate", "pollution", "forest density", "other"))) %>%
            select(model, stat, group, var, value) %>%
            mutate(stat = factor(stat, levels = c("exposure", "sensitivity", "response"),
                                 labels = c("\nexposure", "\nsensitivity", "\nresponse")),
                   model = str_replace(model, "recruitment", "recr.")) %>%
            mutate(model = factor(model,
                                  levels = c("recr. prob.", "recr. rate", "growth", "survival"))) %>%
            mutate(label = ifelse(value == 1, "1.0", str_sub(round(value, 2), 2, 4)),
                   label_color = ifelse(value < .05, "white", "black"))


      pi <- ggplot(i,
             aes(model, var, fill = value)) +
            facet_grid(group~stat, scales = "free") +
            geom_tile(color = "white") +
            geom_text(aes(label = label, color = label_color)) +
            scale_color_identity() +
            scale_fill_viridis_c(trans = "log",
                                 breaks = c(.01, .03, .1, .3, 1)) +
            scale_x_discrete(expand = c(0,0)) +
            scale_y_discrete(expand = c(0,0)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                  axis.ticks = element_blank(),
                  strip.background = element_rect(fill = "black", color = "black"),
                  panel.border = element_blank(),
                  strip.text = element_text(color = "white"),
                  strip.text.y = element_blank(),
                  legend.position = "top") +
            labs(x = NULL,
                 y = NULL,
                 fill = "relative importance") +
            guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                          barwidth = 10, barheight = .5))


      pv <- ggplot(v,
                   aes(stat, var, fill = value)) +
            facet_grid(group ~ model, scales = "free") +
            geom_tile(color = "white") +
            geom_text(aes(label = label, color = label_color)) +
            scale_color_identity() +
            scale_fill_viridis_c(trans = "log", option = "B",
                                 breaks = c(.003, .03,  .3),
                                 labels = c("0.3%", "3%", "30%")) +
            scale_x_discrete(expand = c(0,0)) +
            scale_y_discrete(expand = c(0,0)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  strip.background = element_rect(fill = "black", color = "black"),
                  panel.border = element_blank(),
                  strip.text = element_text(color = "white"),
                  strip.text.y = element_text(angle = 0, vjust = 0),
                  legend.position = "top") +
            labs(x = NULL,
                 y = NULL,
                 fill = "% of response variance explained") +
            guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                          barwidth = 10, barheight = .5))

      p <- pi + pv + plot_layout(nrow = 1, widths = c(12, 8))
      ggsave("figures/importance_combo.pdf",
             p, width = 10, height = 5, units = "in")
}


