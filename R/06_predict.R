get_draws <- function(x){
      x %>%
            map_dfr(function(y){
                  if(is.null(y)) return(NULL)
                  y$draws %>%
                        as.data.frame() %>% as_tibble() %>%
                        left_join(as.data.frame(y$diagnostics),
                                  by = join_by(.chain, .iteration, .draw)) %>%
                        mutate(sp = y$species)})
}

get_summary <- function(x){
      x %>%
            map_dfr(function(y){
                  if(is.null(y)) return(NULL)
                  y$summary %>%
                        mutate(sp = y$species)})
}

diagnose <- function(x){
      d <- get_draws(x) %>%
            group_by(sp) %>%
            summarize(divergence = max(divergent__) > 1,
                      n_chains = length(unique(.chain)),
                      .groups = "drop") %>%
            mutate(all_chains = n_chains == max(n_chains))

      s <- get_summary(x) %>%
            group_by(sp) %>%
            summarize(poor_mixing = max(rhat[variable != "sigma"]) > 1.05, .groups = "drop")
      v <- left_join(d, s, by = join_by(sp)) %>%
            mutate(valid = !divergence & !poor_mixing & all_chains)
      message(round(mean(v$valid) * 100), "% of species have acceptable MCMC diagnostics.")
      v
}

valid_draws <- function(x){
      v <- diagnose(x)
      if(any(!v$valid)) message("Dropping ", sum(!v$valid), " species with poor diagnostics: ",
                                paste(v$sp[!v$valid], collapse = ", "))
      get_draws(x) %>%
            filter(sp %in% v$sp[v$valid])
}

esr_data <- function(x){
      x %>%
            filter(i == i[1],
                   ! str_detect(var, "dbh|intercept")) %>%
            dplyr::select(-lon, -lat) %>%
            rename(exposure = value,
                   sensitivity = param_value) %>%
            mutate(response = exposure * sensitivity) %>%
            gather(stat, value, exposure, sensitivity, response) %>%
            mutate(year = paste0("y", year)) %>%
            spread(year, value) %>%
            mutate(delta = y2010 - y2009)
}


# fitted parameter values
load_fits <- function(x, draws = NULL){

      x <- x %>%
            dplyr::select(-lp__) %>%
            gather(param, value, contains("[")) %>%
            mutate(param = str_remove(param, "\\]"),
                   param = str_replace(param, "\\[", "_")) %>%
            separate(param, c("param", "var")) %>%
            mutate(var = recode(var, "1" = "intercept",
                                "2" = "bacon", "3" = "bahet",
                                "4" = "sulfur", "5" = "nitrogen",
                                "6" = "temperature", "7" = "precipitation",
                                "8" = "dbh"))

      if(!is.null(draws)) x <- x %>% filter(.draw %in% draws)

      x %>% dplyr::select(i = .draw, species = sp, var, param, param_value = value)
}

# Reproducible subset of posterior draw indices. Previously this selection was
# made with an unseeded sample(), so results shifted between runs.
choose_draws <- function(x, ndraws = NULL, seed = 123){
      d <- sort(unique(x$.draw))
      if(is.null(ndraws) || ndraws >= length(d)) return(d)
      set.seed(seed)
      sort(sample(d, ndraws))
}

# Posterior mean parameter values across ALL draws.
#
# Quantities that are linear in the parameters ??? the exposure-sensitivity-response
# decomposition, and growth predictions (identity link) ??? are computed exactly by
# plugging in posterior means, since the mean of a linear function equals the
# function of the mean. This is both exact and free, and removes Monte Carlo noise
# from these results entirely.
mean_fits <- function(x){
      load_fits(x) %>%
            group_by(species, var, param) %>%
            summarize(param_value = mean(param_value), .groups = "drop") %>%
            mutate(i = 1L)
}

# Compute the linear predictor for a batch of posterior draws, using data.table.
#
# This is the performance-critical step: for each draw, each tree-year's
# environmental values are multiplied by that species' coefficients and summed.
# The dplyr equivalent (many-to-many join, then grouped summarize over millions of
# groups) is both slow and memory-hungry, because grouped aggregation at this
# scale is dominated by the grouping machinery rather than the arithmetic.
#
#   ode  : long-format environmental data (one row per group x var)
#   fits : long-format coefficients (one row per draw x species x var [x param])
#   keys : grouping columns identifying a prediction unit
#   d    : draw ids in this batch
#
# Returns one row per draw x key (x param), with the summed linear predictor.
linpred_batch <- function(ode, fits, keys, d, by_param = FALSE){

      # `i` is data.table's own name for the row-subset argument, so the draw-index
      # column is renamed internally to avoid collisions in both the subset and the
      # `by=` grouping, then restored on the way out.
      fits <- fits[fits$i %in% d, , drop = FALSE]
      names(fits)[names(fits) == "i"] <- "draw_id"

      O <- data.table::as.data.table(ode)
      F <- data.table::as.data.table(fits)

      data.table::setkeyv(O, c("species", "var"))
      data.table::setkeyv(F, c("species", "var"))

      grp <- c("draw_id", keys, if(by_param) "param")

      # Join first, then aggregate. Doing both in one X[Y, j, by=] expression is
      # faster still, but the `i.` prefixing of the join argument's columns in the
      # j/by scope is easy to get wrong and fails silently; keeping them separate
      # means every column has a plain name at aggregation time. data.table's
      # grouped aggregation is the part that matters for speed here.
      J <- merge(O, F, by = c("species", "var"), allow.cartesian = TRUE)
      out <- J[, .(pred = sum(param_value * value)), by = grp]
      rm(J)

      data.table::setnames(out, "draw_id", "i")
      tibble::as_tibble(out)
}


# Average a per-draw computation over posterior draws in batches.
#
# Needed for mortality and recruitment, where a nonlinear link is applied before
# differencing, so posterior means are not exact. Draws are collapsed to one row
# per group within each batch, so peak memory scales with `batch` rather than
# `ndraws` and the accumulated result stays small.
batch_average <- function(draw_ids, batch, keys, fun){
      acc <- NULL
      for(d in split(draw_ids, ceiling(seq_along(draw_ids) / batch))){
            b <- data.table::as.data.table(fun(d))
            b <- b[, .(psum = sum(pred), nd = .N), by = keys]
            acc <- if(is.null(acc)) b else
                  rbind(acc, b)[, .(psum = sum(psum), nd = sum(nd)), by = keys]
            rm(b); gc(verbose = FALSE)
      }
      acc[, pred := psum / nd]
      acc[, psum := NULL][, nd := NULL]
      tibble::as_tibble(acc)
}

# observed data
load_obs <- function(x){
      x %>%
            mutate(intercept = 1,
                   year = (year + year_next) / 2) %>%
            rename(obs = outcome) %>%
            dplyr::select(species, plot_id, tree_id, year, obs, t,
                          intercept, dbh = dia, bacon, bahet, sulfur, nitrogen, temperature = bio1, precipitation = bio12) %>%
            gather(var, value, intercept:precipitation)
}

baseline_pred <- function(fi, od){
      left_join(fi, od, by = join_by(species, var)) %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(obs = mean(obs),
                      pred = sum(param_value * value),
                      .groups = "drop")
}

baseline_pred_rec <- function(fi, od){
      left_join(fi, od, by = join_by(species, var)) %>%
            group_by(i, species, plot_id, year, t, param) %>%
            summarize(obs = mean(obs),
                      pred = sum(param_value * value),
                      .groups = "drop")
}

expand_env_obs <- function(od, e){

      # special treatment for recruitment
      rec <- ! "tree_id" %in% names(od)
      if(rec) od$tree_id <- 1
      nterms <- ifelse(rec, 7, 8)

      # join env trends to obs trees
      ode <- od %>%
            dplyr::select(species, plot_id, tree_id) %>%
            distinct() %>%
            left_join(e, by = join_by(species, plot_id),
                      relationship = "many-to-many") %>%
            na.omit()

      # add trends in diameter and intercept
      ode <- od %>%
            filter(var %in% c("dbh", "intercept")) %>%
            group_by(species, plot_id, tree_id, var) %>%
            reframe(value1 = mean(value),
                    t = sum(t),
                    value = project(year, value, c(2009, 2010)),
                    year = c(2009, 2010)) %>%
            mutate(value = ifelse(is.finite(value), value, value1)) %>%
            dplyr::select(-value1) %>%
            bind_rows(ode)

      ode %>%
            group_by(species, plot_id, tree_id, year) %>%
            mutate(nrows = n()) %>% ungroup() %>%
            filter(nrows == nterms) %>%
            dplyr::select(-nrows)
}

predict_growth <- function(grow_draws, grow_data, trends){

      select <- dplyr::select
      grow_draws <- valid_draws(grow_draws)

      # Growth uses an identity link and all downstream quantities are linear in
      # the parameters, so posterior means give exactly the draw-averaged result.
      fi <- mean_fits(grow_draws)
      od <- load_obs(grow_data)


      ### baseline ================================================

      # combine with fitted parameters and calculate predicted growth
      fdp <- baseline_pred(fi, od)


      ### trends ================================

      # join FIA data, env trends, and fitted params
      ode <- expand_env_obs(od, trends)
      fe <- left_join(fi, ode, relationship = "many-to-many",  by = join_by(species, var))

      # exposure-sensitivity-response data
      esr <- esr_data(fe)

      # combine with fitted parameters and calculate predicted growth
      fe <- fe %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop") %>%
            select(-i)

      # observed trends for 2009 to 2010
      mt <- od %>%
            select(species, plot_id, tree_id, t, year, obs) %>%
            distinct() %>%
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    obs = project(year, obs, c(2009, 2010)),
                    year = c(2009, 2010))

      dt <- 1

      fm <- left_join(fe, mt, by = join_by(species, plot_id, year)) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2010 - pred_2009) / dt,
                   drdt_obs = (obs_2010 - obs_2009) / dt)

      fm0 <- fm

      fm <- filter(fm, nreps > 1, drdt_pred != 0)

      # geographic bins
      fm <- fm %>%
            left_join(grow_data %>% select(plot_id, lon, lat) %>% distinct(),
                      by = join_by(plot_id)) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))

      # figure showing n inventories and landscape grid cells
      cells <- fm %>% select(bin) %>% distinct() %>%
            separate(bin, c("x", "y"), sep = " ") %>%
            mutate(x = as.integer(x), y = as.integer(y))
      ninv <- fm0 %>% left_join(grow_data %>% select(plot_id, lon, lat) %>% distinct()) %>%
            select(lon, lat, nreps) %>% distinct() %>% arrange(nreps) %>% mutate(n = nreps + 1)
      p <- ggplot() +
            geom_point(data = ninv, aes(lon, lat, color = n), size = .5) +
            geom_tile(data = cells, aes(x, y), color = "black", linewidth = .5, fill = NA) +
            scale_color_gradientn(colors = c("palegreen", "forestgreen", "darkgreen", "black")) +
            guides(color = guide_legend()) +
            theme_void() +
            labs(color = "number of\ninventories")
      ggsave("figures/plot_map.png", p, width = 8, height = 4, units = "in")

      # return data
      list(baseline = fdp,
           esr = esr,
           trend = fm)
}

predict_mortality <- function(mort_draws, mort_data, trends, ndraws = 50, batch = 5, seed = 123){

      select <- dplyr::select
      mort_draws <- valid_draws(mort_draws)
      fi_all <- load_fits(mort_draws)
      dr <- choose_draws(mort_draws, ndraws, seed)
      od <- load_obs(mort_data)

      keys <- c("species", "plot_id", "tree_id", "year", "t")


      ### baseline ================================================

      # combine with fitted parameters and calculate predicted mortality.
      # inv_logit is nonlinear, so this is averaged over draws rather than
      # evaluated at the posterior mean.
      obs_lookup <- baseline_pred(filter(fi_all, i == dr[1]), od) %>%
            select(all_of(keys), obs)
      fdp <- batch_average(dr, batch, keys, function(d){
            baseline_pred(filter(fi_all, i %in% d), od) %>%
                  mutate(pred = ann2multi(inv_logit(pred), t)) # annual logit to multiyear probability
      }) %>%
            left_join(obs_lookup, by = keys)


      ### trends ================================

      # join FIA data, env trends, and fitted params
      ode <- expand_env_obs(od, trends)

      # exposure-sensitivity-response data. These are linear in the parameters,
      # so posterior means give the exact draw-averaged values.
      esr <- left_join(mean_fits(mort_draws), ode,
                       relationship = "many-to-many", by = join_by(species, var)) %>%
            esr_data()

      # Carry only the columns the prediction needs into the batched join. lon/lat
      # are used by esr_data() above but not here, and `param` is constant for this
      # model; at tens of millions of rows per draw each dropped column matters.
      ode_lean <- ode %>% select(species, plot_id, tree_id, var, value, year, t)
      fi_lean <- fi_all %>% select(i, species, var, param_value)

      # combine with fitted parameters and calculate predicted mortality
      fe <- batch_average(dr, batch, keys, function(d){
            linpred_batch(ode_lean, fi_lean, keys, d) %>%
                  mutate(pred = inv_logit(pred))
      })

      # observed trends for 2009 to 2010
      mt <- od %>%
            select(species, plot_id, tree_id, t, year, obs) %>%
            distinct() %>%
            mutate(obs = obs / t * 2) %>% # convert from multiyear to annual
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    obs = pmax(0, pmin(1, project(year, obs, c(2009, 2010)))),
                    year = c(2009, 2010))

      dt <- 1

      fm <- left_join(fe, mt, by = join_by(species, plot_id, year)) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2010 - pred_2009) / dt,
                   drdt_obs = (obs_2010 - obs_2009) / dt) %>%
            filter(nreps > 1)

      fm <- filter(fm, drdt_pred != 0)

      # geographic bins
      fm <- fm %>%
            left_join(mort_data %>% select(plot_id, lon, lat) %>% distinct(),
                      by = join_by(plot_id)) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))

      list(baseline = fdp,
           esr = esr,
           trend = fm)

}

predict_recruitment <- function(recr_draws, recr_data, trends, recr_ts,
                                ndraws = 50, batch = 5, seed = 123){

      select <- dplyr::select
      recr_draws <- valid_draws(recr_draws)
      fi_all <- load_fits(recr_draws)
      dr <- choose_draws(recr_draws, ndraws, seed)

      # observed data
      load_obs_rec <- function(x){
            x %>%
                  mutate(intercept = 1,
                         year = (yr0 + yr1) / 2) %>%
                  rename(obs = outcome) %>%
                  dplyr::select(species, plot_id, #tree_id,
                                year, obs, t,
                                intercept,
                                bacon, bahet, sulfur, nitrogen, temperature = bio1, precipitation = bio12) %>%
                  gather(var, value, intercept:precipitation)
      }
      od <- load_obs_rec(recr_data) %>%
            filter(is.finite(obs))
      odt <- load_obs_rec(recr_ts) %>%
            filter(is.finite(obs))


      ### baseline ================================================

      # Observed rates are draw-independent, so compute once.
      base_keys <- c("species", "plot_id", "year", "t")
      obs_rec <- baseline_pred_rec(filter(fi_all, i == dr[1]), od) %>%
            gather(stat, value, obs, pred) %>%
            unite(stat, param, stat) %>%
            spread(stat, value) %>%
            mutate(prob_obs = as.integer(zeta_obs > 0),
                   obs = prob_obs * zeta_obs) %>%
            select(all_of(base_keys), obs)

      # Predicted rates combine zeta and mu nonlinearly, so average over draws.
      fdp <- batch_average(dr, batch, base_keys, function(d){
            baseline_pred_rec(filter(fi_all, i %in% d), od) %>%
                  gather(stat, value, obs, pred) %>%
                  unite(stat, param, stat) %>%
                  spread(stat, value) %>%
                  mutate(pred = ann2multi(inv_logit(zeta_pred), t) * exp(mu_pred)) %>%
                  select(i, all_of(base_keys), pred) %>%
                  na.omit()
      }) %>%
            left_join(obs_rec, by = base_keys) %>%
            na.omit()


      ### trends ================================

      ptrends <- trends %>%
            mutate(plot_id = str_sub(plot_id, 1, -3)) %>%
            group_by(plot_id, lon, lat, species, year, var) %>%
            summarize(value = mean(value), .groups = "drop")

      # join FIA data, env trends, and fitted params
      ode <- expand_env_obs(odt, ptrends) %>%
            group_by(species, plot_id) %>%
            mutate(t = unique(na.omit(t))) %>% # make sure t is present for all params
            ungroup()

      # exposure-sensitivity-response data. These are linear in the parameters,
      # so posterior means give the exact draw-averaged values.
      esr <- left_join(mean_fits(recr_draws), ode,
                       relationship = "many-to-many", by = join_by(species, var)) %>%
            esr_data()

      # Carry only the columns the prediction needs into the batched join.
      # tree_id is a constant placeholder for recruitment; lon/lat are used by
      # esr_data() above but not here.
      ode_lean <- ode %>% select(species, plot_id, var, value, year, t)
      fi_lean <- fi_all %>% select(i, species, var, param, param_value)

      # calculate predicted recruitment rate. The zeta/mu combination below is
      # nonlinear, so this is averaged over draws in batches rather than
      # evaluated at the posterior mean.
      rec_keys <- c("species", "plot_id", "year", "t")
      fe <- batch_average(dr, batch, rec_keys, function(d){
            linpred_batch(ode_lean, fi_lean, rec_keys, d, by_param = TRUE) %>%
                  spread(param, pred) %>%
                  mutate(pred = ann2multi(inv_logit(zeta), t) * exp(mu)) %>%
                  select(-mu, -zeta)
      })

      # observed trends for 2009 to 2010
      mt <- odt %>%
            select(species, plot_id, t, year, obs) %>%
            distinct() %>%
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    obs = pmax(0, project(year, obs, c(2009, 2010))),
                    year = c(2009, 2010),
                    n = n())

      dt <- 1

      fm <- fe %>%
            left_join(mt, by = join_by(species, plot_id, year)) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2010 - pred_2009) / dt,
                   drdt_obs = (obs_2010 - obs_2009) / dt) %>%
            filter(nreps > 1)

      # geographic bins
      fm <- fm %>%
            left_join(ptrends %>% select(plot_id, lon, lat) %>% distinct(),
                      by = join_by(plot_id)) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))

      fm <- fm %>%
            filter(is.finite(drdt_pred), is.finite(drdt_obs))

      list(baseline = fdp,
           esr = esr,
           trend = fm)
}
