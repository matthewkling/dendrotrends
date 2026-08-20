# Compute numeric values cited in the results text, in one place

# Pooled (unstratified) variance partitioning, matching Fig. 5b.
importance_values <- function(esr, species, oscl, n_sample = 50000, seed = 123){
      set.seed(seed)
      vd <- build_vp_data(esr, species)
      unique(vd$data$model) %>%
            map_dfr(function(m){
                  zz <- vd$data %>% filter(model == m)
                  zz <- slice_sample(zz, n = min(n_sample, nrow(zz)))
                  vp <- rdacca.hp::rdacca.hp(zz$response,
                                             dplyr::select(zz, all_of(vd$predictors)))
                  vp$Hier.part %>%
                        as.data.frame() %>%
                        tibble::rownames_to_column("var") %>%
                        janitor::clean_names() %>%
                        mutate(total_explained = vp$Total_explained_variation,
                               model = m)
            }) %>%
            group_by(model) %>%
            mutate(percent = individual / sum(individual) *
                         ifelse(first(total_explained) > 1,
                                first(total_explained) / 100,
                                first(total_explained)) * 100) %>%
            ungroup() %>%
            separate(var, c("stat", "var"), sep = "_", fill = "right")
}

# Pooled mean absolute standardized effects, matching Fig. 5a but WITHOUT the
# rescaling by the per-component maximum, so that values are comparable between
# parameterizations rather than only within a panel.
mase_values <- function(esr, species, oscl){
      scl   <- format_scl(species$scl)
      gscl  <- format_scl(species$gscl)
      goscl <- oscl$global

      esr %>%
            dplyr::select(-tree_id, -t, -y2009, -y2010, -delta, -param) %>%
            filter(var != "combined") %>%
            mutate(var = as.character(var)) %>%
            left_join(scl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"),
                                  value * sd, value)) %>%
            dplyr::select(-sd, -mean) %>%
            left_join(gscl) %>%
            mutate(value = ifelse(stat %in% c("exposure", "sensitivity"),
                                  value / sd, value)) %>%
            dplyr::select(-sd, -mean) %>%
            left_join(goscl) %>%
            mutate(value = ifelse(stat %in% c("response", "sensitivity"),
                                  value / sd, value)) %>%
            dplyr::select(-sd, -mean) %>%
            group_by(stat, model, var) %>%
            summarize(mase = mean(abs(value), na.rm = TRUE), .groups = "drop")
}

# Weighted STS summary statistics, matching the annotations in Fig. S3.
eval_summary <- function(trend){
      norm <- function(n) n / sum(n)
      bin <- function(x){
            x %>% summarize(r = wCorr::weightedCorr(drdt_pred, drdt_obs, "Pearson",
                                                    weights = tree_yrs),
                            drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                            drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                            n = sum(tree_yrs), .groups = "drop") %>%
                  group_by(model) %>% mutate(n = norm(n))
      }
      pd <- bind_rows(
            trend %>% group_by(species, model) %>% bin() %>% mutate(group = "species"),
            trend %>% group_by(bin, model) %>% bin() %>% mutate(group = "landscapes")) %>%
            filter(n > .001)

      s <- pd %>% group_by(model, group) %>%
            summarize(ppos = weighted.mean(r > 0, n),
                      r = wCorr::weightedCorr(drdt_pred, drdt_obs,
                                              weights = n, method = "Pearson"),
                      .groups = "drop")

      bind_rows(
            s %>% group_by(group) %>%
                  summarize(label = paste0("sts_pct_positive_within_", first(group)),
                            value = paste0(round(min(ppos)*100), "-", round(max(ppos)*100)),
                            units = "%", .groups = "drop"),
            s %>% group_by(group) %>%
                  summarize(label = paste0("sts_r_among_", first(group)),
                            value = paste0(round(min(r), 2), "-", round(max(r), 2)),
                            units = "r", .groups = "drop"))
}


reported_values <- function(esr, species, oscl, trends, eval,
                            strat_imp = NULL, esr_alt = NULL, species_alt = NULL,
                            oscl_alt = NULL, brk = -100, suffix = ""){

      select <- dplyr::select
      out <- list()
      # value is stored as character because some entries are ranges ("23-37") or
      # paired comparisons rather than single numbers; this is a lookup table for
      # dropping into the text, not something to compute on.
      add <- function(label, value, units, note)
            out[[length(out) + 1]] <<- tibble(label = label,
                                              value = as.character(value),
                                              units = units, note = note)


      ## --- Exposure to environmental change ------------------------------------
      # Standardized rates of environmental change, in units of the global spatial
      # SD per decade, matching the Fig. 2 caption ("dividing rates of change by the
      # standard deviation in each variable's long-term mean across all trees").

      scl  <- format_scl(species$scl)
      gscl <- format_scl(species$gscl)

      exp_std <- esr %>% ungroup() %>%
            filter(stat == "exposure", var != "combined") %>%
            mutate(var = as.character(var)) %>%
            left_join(scl, by = join_by(species, var)) %>%
            mutate(delta = delta * sd) %>%          # per-species SD -> native units
            select(-sd, -mean) %>%
            left_join(gscl, by = join_by(var)) %>%
            mutate(delta = delta / sd) %>%          # native units -> global SD
            select(-sd, -mean) %>%
            group_by(var) %>%
            summarize(rate_std = mean(abs(delta), na.rm = TRUE) * 10, # per decade
                      n = n(), .groups = "drop")

      for(v in exp_std$var)
            add(paste0("exposure_sd_per_decade_", v),
                round(exp_std$rate_std[exp_std$var == v], 3),
                "sd/decade", "Results: Exposure to environmental change")

      add("exposure_ratio_sulfur_to_next",
          round(max(exp_std$rate_std) / sort(exp_std$rate_std, decreasing = TRUE)[2], 2),
          "ratio", "Results: sulfur vs. next fastest-changing variable")


      ## --- Modeled demographic responses: directional percentages --------------
      # Replicates response_plots(): per-tree response summed across variables,
      # link-transformed, mortality flipped to survival, then the share of 0.5-deg
      # grid cells and of species whose median trend is positive.

      ll <- trends %>%
            dplyr::select(plot_id, lon, lat) %>%
            distinct()
      ll <- ll %>%
            mutate(plot_id = str_sub(plot_id, 1, -3)) %>%
            group_by(plot_id) %>%
            summarize(lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
            bind_rows(ll) %>%
            distinct(plot_id, .keep_all = TRUE)

      pd <- esr %>% ungroup() %>%
            filter(stat == "response", var != "combined") %>%
            mutate(model = ifelse(str_detect(model, "recruit"), "recruitment",
                                  as.character(model))) %>%
            group_by(i, species, model, param, plot_id, tree_id) %>%
            summarize(y2010 = sum(y2010), y2009 = sum(y2009), .groups = "drop") %>%
            mutate(y2009 = case_when(param == "mu" ~ exp(y2009),
                                     param == "zeta" ~ inv_logit(y2009),
                                     param == "beta" ~ y2009),
                   y2010 = case_when(param == "mu" ~ exp(y2010),
                                     param == "zeta" ~ inv_logit(y2010),
                                     param == "beta" ~ y2010)) %>%
            group_by(i, species, model, plot_id, tree_id) %>%
            summarize(y2009 = prod(y2009), y2010 = prod(y2010),
                      value = y2010 - y2009, .groups = "drop") %>%
            left_join(ll, by = "plot_id") %>%
            mutate(value = ifelse(model == "mortality", -value, value),
                   model = ifelse(model == "mortality", "survival", model))

      cells <- pd %>%
            mutate(lon = plyr::round_any(lon, .5), lat = plyr::round_any(lat, .5)) %>%
            group_by(model, lon, lat) %>%
            summarize(value = median(value, na.rm = TRUE), .groups = "drop") %>%
            group_by(model) %>%
            summarize(pct = mean(value > 0, na.rm = TRUE) * 100, .groups = "drop")

      spp <- pd %>%
            group_by(model, species) %>%
            summarize(value = median(value, na.rm = TRUE), .groups = "drop") %>%
            group_by(model) %>%
            summarize(pct = mean(value > 0, na.rm = TRUE) * 100, .groups = "drop")

      for(m in cells$model){
            add(paste0("pct_gridcells_increasing_", m),
                round(cells$pct[cells$model == m]), "%",
                "Results: Modeled demographic responses (0.5-deg cells)")
            add(paste0("pct_species_increasing_", m),
                round(spp$pct[spp$model == m]), "%",
                "Results: Modeled demographic responses (species)")
      }


      ## --- Variance partitioning values ----------------------------------------
      # Pulled from the same object that produces Fig. 5b.

      vp <- importance_values(esr, species, oscl)

      vp_get <- function(d, s, v){
            x <- vp %>% filter(model == d, stat == s, var == v)
            round(x$percent[1])
      }

      vp_cite <- tribble(
            ~d,                        ~s,            ~v,
            "growth",                  "exposure",    "bacon",
            "recruitment\nrate",       "exposure",    "bacon",
            "growth",                  "exposure",    "bahet",
            "recruitment\nprobability","sensitivity", "sulfur",
            "recruitment\nrate",       "sensitivity", "sulfur",
            "growth",                  "sensitivity", "sulfur",
            "mortality",               "sensitivity", "sulfur",
            "mortality",               "sensitivity", "temperature")
      for(k in seq_len(nrow(vp_cite))){
            r <- vp_cite[k, ]
            add(paste0("vp_", r$s, "_", r$v, "_",
                       str_replace(r$d, "\n", "_")),
                vp_get(r$d, r$s, r$v), "% of variance",
                "Results: variance partitioning (Fig. 5b)")
      }

      add("vp_sensitivity_sulfur_range",
          paste0(min(filter(vp, stat == "sensitivity", var == "sulfur")$percent) %>% round(), "-",
                 max(filter(vp, stat == "sensitivity", var == "sulfur")$percent) %>% round()),
          "% of variance", "Results: sulfur sensitivity range across the four rates")


      ## --- Model evaluation (STS) ----------------------------------------------

      if(!is.null(eval)){
            ev <- eval_summary(eval)
            for(k in seq_len(nrow(ev)))
                  add(ev$label[k], ev$value[k], ev$units[k],
                      "Results: Model evaluation (Fig. S3)")
      }


      ## --- Alternative climate parameterization (Fig. S5) ----------------------

      if(!is.null(esr_alt)){
            m_alt <- mase_values(esr_alt, species_alt, oscl_alt)
            m_base <- mase_values(esr, species, oscl)

            vpd <- m_alt %>% filter(stat == "exposure", var == "temperature") %>%
                  summarize(v = mean(mase)) %>% pull(v)
            tas <- m_base %>% filter(stat == "exposure", var == "temperature") %>%
                  summarize(v = mean(mase)) %>% pull(v)
            add("altclim_vpd_vs_temp_exposure_ratio", round(vpd / tas, 2), "ratio",
                "Results: alt climate, VPD vs mean temperature exposure")

            # the two climate slots are named temperature/precipitation internally
            # in both parameterizations; label by what actually occupies them
            for(pair in list(c("temperature", "maxVPD_vs_temp"),
                             c("precipitation", "minCMI_vs_precip"))){
                  a <- m_alt %>% filter(stat == "sensitivity", var == pair[1],
                                        model == "recruitment\nprobability") %>% pull(mase)
                  b <- m_base %>% filter(stat == "sensitivity", var == pair[1],
                                         model == "recruitment\nprobability") %>% pull(mase)
                  add(paste0("altclim_recrprob_sensitivity_", pair[2]),
                      paste0(signif(a, 2), " vs ", signif(b, 2)), "MASE",
                      "Results: alt climate, recruitment probability sensitivity")
            }
      }


      ## --- Stratified analyses (Fig. S6) ---------------------------------------

      if(!is.null(strat_imp)){
            sv <- strat_imp$vp
            sg <- function(rg, d, s, v){
                  x <- sv %>% filter(region == rg, model == d, stat == s, var == v)
                  round(x$percent[1])
            }
            add("strat_vp_sulfur_sens_survival_west", sg("west", "mortality", "sensitivity", "sulfur"),
                "% of variance", "Results: stratified, western survival")
            add("strat_vp_sulfur_sens_survival_east", sg("east", "mortality", "sensitivity", "sulfur"),
                "% of variance", "Results: stratified, eastern survival")
            add("strat_vp_bacon_exp_growth_west", sg("west", "growth", "exposure", "bacon"),
                "% of variance", "Results: stratified, western growth")
            add("strat_vp_bacon_exp_recrrate_west",
                sg("west", "recruitment\nrate", "exposure", "bacon"),
                "% of variance", "Results: stratified, western recruitment rate")

            sm <- strat_imp$mase
            ratio <- function(s){
                  e <- sm %>% filter(region == "east", stat == s, var == "sulfur") %>% pull(mase)
                  w <- sm %>% filter(region == "west", stat == s, var == "sulfur") %>% pull(mase)
                  round(mean(e) / mean(w), 1)
            }
            add("strat_sulfur_exposure_east_west_ratio", ratio("exposure"), "ratio",
                "Results: stratified, sulfur exposure east:west")
            add("strat_sulfur_sensitivity_east_west_ratio", ratio("sensitivity"), "ratio",
                "Results: stratified, sulfur sensitivity east:west")

            # sample sizes and model fit, for figure captions
            for(rg in c("west", "east")){
                  n <- sv %>% filter(region == rg) %>% distinct(model, n)
                  for(k in seq_len(nrow(n)))
                        add(paste0("strat_n_", rg, "_", str_replace(n$model[k], "\n", "_")),
                            n$n[k], "rows", "Fig. S6 caption: sample size")
                  te <- sv %>% filter(region == rg) %>%
                        distinct(model, total_explained)
                  for(k in seq_len(nrow(te)))
                        add(paste0("strat_totexp_", rg, "_", str_replace(te$model[k], "\n", "_")),
                            round(te$total_explained[k] * 100), "%",
                            "Fig. S6 caption: total variance explained")
            }
      }


      res <- bind_rows(out)
      readr::write_csv(res, paste0("figures/reported_values", suffix, ".csv"))
      res
}
