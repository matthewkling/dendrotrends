# East-west stratified variance partitioning

# Assemble the wide predictor matrix used by rdacca.hp.
# This mirrors the variance partitioning section of importance_plots(), so that
# stratified and pooled results are directly comparable.
build_vp_data <- function(f, species){

      select <- dplyr::select
      scl <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)

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

      list(data = z, predictors = vars)
}


# Mean absolute standardized effects (MASE) within longitudinal strata.
#
# Companion to stratified_vp(). The variance partitioning expresses each variable
# as a share of the variance explained *within its own region*, so a lower share
# in one region can reflect either weaker absolute importance or a different
# denominator. MASE has no such denominator: it is in common standardized units,
# so values are directly comparable across regions.
#
# The response component is the most informative here, since it combines
# sensitivity with each region's actual exposure.
#
# NOTE: unlike Fig. 5a, values are NOT rescaled within region ??? they are rescaled
# by a single maximum per component across all regions, so that regions remain
# comparable. Raw (unscaled) values are also returned.
compute_stratified_mase <- function(f, species, oscl, plot_lon,
                                    brk = -100){

      select <- dplyr::select
      scl   <- format_scl(species$scl)
      gscl  <- format_scl(species$gscl)
      goscl <- oscl$global

      # global standardization, mirroring the `amg` block of importance_plots(),
      # but retaining plot_id so observations can be assigned to a region
      d <- f %>%
            select(-tree_id, -t, -y2009, -y2010, -delta, -param) %>%
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
            select(-sd, -mean) %>%
            left_join(plot_lon, by = "plot_id")

      miss <- d %>% group_by(model) %>% summarize(frac_na = mean(is.na(lon)), .groups = "drop")
      if(any(miss$frac_na > 0.01))
            stop("plot_lon does not cover all models; missing fractions: ",
                 paste(miss$model, round(miss$frac_na, 3), sep = "=", collapse = ", "))

      d <- filter(d, is.finite(lon)) %>%
            mutate(region = ifelse(lon < brk, "west", "east"))

      # abundance-weighted mean across trees, matching Fig. 5a
      mase <- d %>%
            group_by(region, stat, model, var) %>%
            summarize(mase = mean(abs(value), na.rm = T),
                      n = n(), .groups = "drop")

      res <- mase %>%
            group_by(stat) %>%                       # one denominator per component,
            mutate(rel = mase / max(mase, na.rm = T)) %>%   # shared across regions
            ungroup() %>%
            mutate(group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                     var %in% c("bacon", "bahet") ~ "forest density",
                                     var %in% c("nitrogen", "sulfur") ~ "pollution"))

      res
}


compute_stratified_vp <- function(f, species, plot_lon,
                                  brk = -100,
                                  n_sample = 50000,
                                  seed = 123){

      select <- dplyr::select
      set.seed(seed)

      vp_data <- build_vp_data(f, species)
      z <- vp_data$data %>% left_join(plot_lon, by = "plot_id")

      # Fail loudly on join problems: growth/mortality are keyed at subplot level
      # while recruitment is keyed at plot level, so an incomplete lookup silently
      # drops entire demographic rates rather than a few rows.
      miss <- z %>%
            group_by(model) %>%
            summarize(frac_na = mean(is.na(lon)), .groups = "drop")
      if(any(miss$frac_na > 0.01))
            stop("plot_lon does not cover all models; missing fractions: ",
                 paste(miss$model, round(miss$frac_na, 3), sep = "=", collapse = ", "))

      z <- filter(z, is.finite(lon)) %>%
            mutate(region = ifelse(lon < brk, "west", "east"))

      # each demographic rate is partitioned separately within each region
      grid <- expand_grid(.region = c("west", "east"),
                          .model = unique(z$model))

      res <- pmap(grid, function(.region, .model){

            zz <- z %>% filter(model == .model, region == .region)

            n_avail <- nrow(zz)
            message(.region, " / ", .model, ": n = ", n_avail)

            zz <- slice_sample(zz, n = min(n_sample, n_avail))

            vp <- rdacca.hp::rdacca.hp(zz$response,
                                       select(zz, all_of(vp_data$predictors)))

            vp$Hier.part %>%
                  as.data.frame() %>%
                  rownames_to_column("var") %>%
                  janitor::clean_names() %>%
                  mutate(total_explained = vp$Total_explained_variation,
                         percent = i_perc_percent * total_explained,
                         region = .region, model = .model, n = n_avail)
      }) %>%
            bind_rows() %>%
            separate(var, c("stat", "var"), sep = "_", fill = "right") %>%
            mutate(group = case_when(var %in% c("temperature", "precipitation") ~ "climate",
                                     var %in% c("bacon", "bahet") ~ "forest density",
                                     var %in% c("nitrogen", "sulfur") ~ "pollution"))

      res
}


# Combined stratified importance figure (SI companion to Fig. 5).

stratified_importance <- function(f, species, oscl, plot_lon,
                                  brk = -100, suffix = ""){

      tidy_axes <- function(x, brk){
            x %>%
                  mutate(var = factor(var,
                                      levels = rev(c("temperature", "precipitation",
                                                     "sulfur", "nitrogen",
                                                     "bacon", "bahet")),
                                      labels = rev(c("temperature", "precipitation",
                                                     "sulfur", "nitrogen",
                                                     "conspec. BA", "heterospec. BA"))),
                         model = str_replace(model, "recruitment", "recr.") %>%
                               str_replace("\n", " ") %>%
                               str_replace("mortality", "survival"),
                         model = factor(model, levels = c("recr. probability", "recr. rate",
                                                          "growth", "survival"),
                                        labels = c("recr. prob.", "recr. rate",
                                                   "growth", "survival")),
                         region = factor(region, levels = c("west", "east"),
                                         labels = c(paste0("West\n(lon < ", brk, ")"),
                                                    paste0("East\n(lon > ", brk, ")"))))
      }

      mase <- compute_stratified_mase(f, species, oscl, plot_lon, brk = brk)
      vp   <- compute_stratified_vp(f, species, plot_lon, brk = brk)

      i <- mase %>%
            mutate(rel = pmax(rel, .01),
                   label = ifelse(rel >= .995, "1.0", str_sub(round(rel, 2), 2, 4)),
                   label_color = ifelse(rel < .05, "white", "black"),
                   stat = factor(stat, levels = c("exposure", "sensitivity", "response"))) %>%
            tidy_axes(brk)

      # Match Fig. 5b
      v <- vp %>%
            group_by(region, model) %>%
            mutate(te = first(total_explained),
                   te = ifelse(te > 1, te / 100, te),
                   value = percent / sum(percent) * te,
                   value = pmax(value, .001)) %>%
            ungroup() %>%
            mutate(label = str_sub(round(value * 100), 1, 2),
                   label_color = ifelse(value < .005, "white", "black"),
                   stat = factor(stat, levels = c("exposure", "sensitivity"))) %>%
            tidy_axes(brk)

      base <- list(
            geom_tile(color = "white"),
            geom_text(aes(label = label, color = label_color), size = 2.7),
            scale_color_identity(),
            scale_x_discrete(expand = c(0, 0)),
            scale_y_discrete(expand = c(0, 0)),
            theme_bw(),
            theme(axis.ticks = element_blank(),
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  strip.background = element_rect(fill = "black", color = "black"),
                  strip.text = element_text(color = "white"),
                  legend.position = "top"),
            labs(x = NULL, y = NULL)
      )

      pi <- ggplot(i, aes(model, var, fill = rel)) +
            facet_grid(region ~ stat) +
            base +
            scale_fill_viridis_c(trans = "log", breaks = c(.01, .03, .1, .3, 1)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                  strip.text.y = element_blank()) +
            labs(fill = "(a)  relative importance") +
            guides(fill = guide_colourbar(title.position = "top", title.hjust = .5,
                                          barwidth = 10, barheight = .5))

      pv <- ggplot(v, aes(stat, var, fill = value)) +
            facet_grid(region ~ model) +
            base +
            scale_fill_viridis_c(trans = "log", option = "B",
                                 breaks = c(.003, .03, .3),
                                 labels = c("0.3%", "3%", "30%")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                  axis.text.y = element_blank()) +
            labs(fill = "(b)  % of response variance explained") +
            guides(fill = guide_colourbar(title.position = "top", title.hjust = .5,
                                          barwidth = 10, barheight = .5))

      p <- pi + pv + plot_layout(nrow = 1, widths = c(12, 8))
      ggsave(paste0("figures/importance_stratified", suffix, ".pdf"),
             p, width = 8, height = 6, units = "in")

      list(mase = mase, vp = vp)
}
