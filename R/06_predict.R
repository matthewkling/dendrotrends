
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
load_fits <- function(x, i = NULL){

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

      if(!is.null(i)) x <- x %>%  filter(.draw %in% sample(unique(.draw), i))

      x %>% dplyr::select(i = .draw, species = sp, var, param, param_value = value)
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

predict_growth <- function(grow_draws, grow_data, trends, ndraws = 5){

      select <- dplyr::select
      grow_draws <- valid_draws(grow_draws)
      fi <- load_fits(grow_draws, ndraws)
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
                      .groups = "drop")

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

predict_mortality <- function(mort_draws, mort_data, trends, ndraws = 5){

      select <- dplyr::select
      mort_draws <- valid_draws(mort_draws)
      fi <- load_fits(mort_draws, ndraws)
      od <- load_obs(mort_data)


      ### baseline ================================================

      # combine with fitted parameters and calculate predicted mortality
      fdp <- baseline_pred(fi, od) %>%
            mutate(pred = ann2multi(inv_logit(pred), t)) # convert from annual logit to multiyear probability


      ### trends ================================

      # join FIA data, env trends, and fitted params
      ode <- expand_env_obs(od, trends)
      fe <- left_join(fi, ode, relationship = "many-to-many", by = join_by(species, var))

      # exposure-sensitivity-response data
      esr <- esr_data(fe)

      # combine with fitted parameters and calculate predicted mortality
      fe <- fe %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop") %>%
            mutate(pred = inv_logit(pred))

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

predict_recruitment <- function(recr_draws, recr_data, trends, recr_ts, ndraws = 5){

      select <- dplyr::select
      recr_draws <- valid_draws(recr_draws)
      fi <- load_fits(recr_draws, ndraws)

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

      # combine with fitted parameters and calculate predicted
      fdp0 <- baseline_pred_rec(fi, od) %>%
            gather(stat, value, obs, pred) %>%
            unite(stat, param, stat) %>%
            spread(stat, value)

      fdp0 <- fdp0 %>%
            mutate(
                  prob_obs = as.integer(zeta_obs > 0),
                  prob_pred = ann2multi(inv_logit(zeta_pred), t),
                  mean_obs = zeta_obs,
                  mean_pred = exp(mu_pred),
                  rate_obs = prob_obs * mean_obs, # same as rate_obs = zeta_obs
                  rate_pred = prob_pred * mean_pred) %>%
            na.omit()

      fdp <- fdp0 %>%
            mutate(pred = rate_pred,
                   obs = rate_obs)


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
      fe <- left_join(fi, ode, relationship = "many-to-many", by = join_by(species, var))

      # exposure-sensitivity-response data
      esr <- esr_data(fe)

      # calculate predicted recruitment rate
      fe <- fe %>%
            group_by(i, species, plot_id, param, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop")

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
            spread(param, pred) %>%
            mutate(pred = ann2multi(inv_logit(zeta), t) * exp(mu)) %>%
            select(-mu, -zeta) %>%
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
