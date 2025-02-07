

predict_growth <- function(grow_draws, grow_data, trends, ndraws){
      # tar_load(grow_data)
      # tar_load(grow_draws)
      # trend_data <- tar_read(fia_env)$trend
      # ndraws = 5

      select <- dplyr::select

      # fitted parameter values
      # fi <- load_fits("data/draws_mortality.csv", ndraws)
      fi <- load_fits(grow_draws, ndraws)

      # observed data
      # od <- load_obs("data/model_data_mortality.csv", "mortality")
      od <- load_obs(grow_data)


      ### baseline ================================================

      # combine with fitted parameters and calculate predicted mortality
      # trans <- sqrt
      # inv_trans <- function(x)x^2
      fdp <- baseline_pred(fi, od)# %>%
            # mutate(obs = trans(obs))


      ### trends ================================

      # join env trends to obs trees
      ode <- expand_env_obs(od, trends)

      # combine with fitted parameters and calculate predicted growth
      fe <- fi %>%
            left_join(ode, relationship = "many-to-many") %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop")

      # observed trends for 2000 to 2020
      mt <- od %>%
            select(species, plot_id, tree_id, t, year, obs) %>%
            distinct() %>%
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    # obs = project(year, trans(obs), c(2000, 2020)),
                    obs = project(year, obs, c(2000, 2020)),
                    year = c(2000, 2020))

      fm <- left_join(fe, mt) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2020 - pred_2000) / 20, # option A: delta sqrt-link variable
                   drdt_obs = (obs_2020 - obs_2000) / 20) %>%
            # mutate(drdt_pred = (inv_trans(pmax(0, pred_2020)) - inv_trans(pmax(0, pred_2000))) / 20, # option B: delta linear variable
            #        drdt_obs = (inv_trans(pmax(0, obs_2020)) - inv_trans(pmax(0, obs_2000))) / 20) %>%
            filter(nreps > 1)

      ### think more about this. these are trees only remeasured once due to aging in or dying out ###
      fm <- filter(fm, drdt_pred != 0)

      # geographic bins
      fm <- fm %>%
            # left_join(trend_env %>% select(plot_id, lon, lat) %>% distinct()) %>%
            left_join(grow_data %>% select(plot_id, lon, lat) %>% distinct()) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))


      list(baseline = fdp,
           trend = fm)

}

predict_mortality <- function(mort_draws, mort_data, trends, ndraws){

      select <- dplyr::select

      # fitted parameter values
      # fi <- load_fits("data/draws_mortality.csv", ndraws)
      fi <- load_fits(mort_draws, ndraws)

      # observed data
      # od <- load_obs("data/model_data_mortality.csv", "mortality")
      od <- load_obs(mort_data)


      ### baseline ================================================

      # combine with fitted parameters and calculate predicted mortality
      fdp <- baseline_pred(fi, od) %>%
            mutate(pred = ann2multi(inv_logit(pred), t)) # convert from annual logit to multiyear probability


      ### trends ================================

      # join env trends to obs trees
      ode <- expand_env_obs(od, trends)

      # combine with fitted parameters and calculate predicted mortality
      fe <- fi %>%
            left_join(ode, relationship = "many-to-many") %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop") %>%
            mutate(pred = inv_logit(pred))
      # write_csv(fe, "data/fe_mortality.csv")

      # observed trends for 2000 to 2020
      mt <- od %>%
            select(species, plot_id, tree_id, t, year, obs) %>%
            distinct() %>%
            mutate(obs = obs / t * 2) %>% # convert from multiyear to annual ### fixme?? should it be *2 since mean mort occurrd halfway thru period?
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    obs = pmax(0, pmin(1, project(year, obs, c(2000, 2020)))),
                    year = c(2000, 2020))
      # write_csv(mt, "data/mt_mortality.csv")

      fm <- left_join(fe, mt) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2020 - pred_2000) / 20,
                   drdt_obs = (obs_2020 - obs_2000) / 20) %>%
            filter(nreps > 1)

      ### think more about this. these are trees only remeasured once due to aging in or dying out ###
      fm <- filter(fm, drdt_pred != 0)

      # geographic bins
      fm <- fm %>%
            # left_join(trend_env %>% select(plot_id, lon, lat) %>% distinct()) %>%
            left_join(mort_data %>% select(plot_id, lon, lat) %>% distinct()) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))
      # write_csv(fm, "data/eval_mortality.csv")

      list(baseline = fdp,
           trend = fm)

}

predict_recruitment <- function(recr_draws, recr_data, trends, ndraws){

      select <- dplyr::select

      # fitted parameter values
      # load_fits_rec <- function(x, i){
      #       x %>%
      #             # read_csv() %>%
      #             dplyr::select(-lp__) %>%
      #             gather(param, value, `zeta[1]`:`mu[7]`) %>%
      #             mutate(param = str_remove(param, "\\]"),
      #                    param = str_replace(param, "\\[", "_")) %>%
      #             separate(param, c("param", "var")) %>%
      #             mutate(var = recode(var, "1" = "intercept",
      #                                 "2" = "bacon", "3" = "bahet",
      #                                 "4" = "sulfur", "5" = "nitrogen",
      #                                 "6" = "temperature", "7" = "precipitation")) %>%
      #             filter(.draw %in% sample(unique(.draw), i)) %>%
      #             dplyr::select(i = .draw, species = sp, var, param, param_value = value)
      # }
      # fi <- load_fits_rec("data/draws_recruitment.csv", ndraws)
      fi <- load_fits(recr_draws, ndraws)

      # observed data
      load_obs_rec <- function(x){
            x %>%
                  mutate(intercept = 1,
                         year = (yr0 + yr1) / 2) %>%
                  rename(obs = outcome) %>%
            # x$obs <- x[[v]]
            # x %>%
                  dplyr::select(species, plot_id, #tree_id,
                                year, obs, t,
                                intercept, # dbh = dia,
                                bacon, bahet, sulfur, nitrogen, temperature = bio1, precipitation = bio12) %>%
                  gather(var, value, intercept:precipitation)
      }
      od <- load_obs_rec(recr_data) %>%
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
                  # prob_obs = as.integer(zeta_obs > 0) / t,
                  # prob_pred = inv_logit(zeta_pred),
                  mean_obs = zeta_obs,
                  mean_pred = exp(mu_pred),
                  rate_obs = prob_obs * mean_obs, # same as rate_obs = zeta_obs
                  rate_pred = prob_pred * mean_pred) %>%
            na.omit()

      fdp <- fdp0 %>%
            mutate(pred = rate_pred,
                   obs = rate_obs)


      ### trends ================================

      # # expand env trends to obs trees
      # expand_env_obs_rec <- function(od, e){
      #
      #       # # join env trends to obs trees
      #       ode <- od %>%
      #             select(species, plot_id) %>%
      #             distinct() %>%
      #             left_join(e) %>%
      #             na.omit()
      #
      #       # add intercept
      #       ode <- od %>%
      #             filter(var %in% c("dbh", "intercept")) %>%
      #             group_by(species, plot_id, var) %>%
      #             reframe(value1 = mean(value), # could project here instead of averaging
      #                     t = sum(t),
      #                     value = project(year, value, c(2000, 2020)),
      #                     year = c(2000, 2020)) %>%
      #             mutate(value = ifelse(is.finite(value), value, value1)) %>%
      #             dplyr::select(-value1) %>%
      #             bind_rows(ode)
      #
      #       ode %>% group_by(species, plot_id, year) %>% filter(length(var) == 7) %>%
      #             group_by(species, plot_id) %>%
      #             mutate(t = unique(na.omit(t))) %>% # make sure t is present for all params
      #             ungroup()
      # }

      ptrends <- trends %>%
            mutate(plot_id = str_sub(plot_id, 1, -3)) %>%
            group_by(plot_id, lon, lat, species, year, var) %>%
            summarize(value = mean(value), .groups = "drop")

      ode <- expand_env_obs(od, ptrends) %>%
            group_by(species, plot_id) %>%
            mutate(t = unique(na.omit(t))) %>% # make sure t is present for all params
            ungroup()


      # combine with fitted parameters and calculate predicted rate
      fe <- fi %>%
            # mutate(var = ifelse(var == "ba_con", "bacon2", var),
            #        var = ifelse(var == "ba_het", "bahet2", var)) %>%
            left_join(ode) %>%
            group_by(i, species, plot_id, param, year, t) %>%
            summarize(pred = sum(param_value * value),
                      .groups = "drop")

      # observed trends for 2000 to 2020
      mt <- od %>%
            select(species, plot_id, t, year, obs) %>%
            distinct() %>%
            # mutate(obs = obs / t) %>% # convert from multiyear to annual. ########### not sure this is good. may need to revisit for mort too
            # mutate(obs = as.integer(obs > 0)) %>% # multiyear zeta
            group_by(species, plot_id) %>%
            reframe(nreps = length(unique(year)),
                    tree_yrs = sum(t),
                    # obs = prj_logit(year, obs, c(2000, 2020)),
                    # obs = pmax(0, pmin(1, project(year, obs, c(2000, 2020)))),
                    obs = pmax(0, project(year, obs, c(2000, 2020))),
                    year = c(2000, 2020),
                    n = n())

      fm <- fe %>%
            spread(param, pred) %>%
            mutate(pred = ann2multi(inv_logit(zeta), t) * exp(mu)) %>%
            select(-mu, -zeta) %>%
            left_join(mt) %>%
            # left_join(fe, mt) %>%
            filter(is.finite(year)) %>%
            gather(stat, value, pred, obs) %>%
            unite(stat, stat, year) %>%
            spread(stat, value) %>%
            mutate(drdt_pred = (pred_2020 - pred_2000) / 20,
                   drdt_obs = (obs_2020 - obs_2000) / 20) %>%
            # filter(species != "Picea mariana") %>%
            filter(nreps > 1)

      ### think more about this. these are trees only remeasured once due to aging in or dying out ###
      # fm <- filter(fm, drdt_pred != 0) # not relevant for recruitment, there are no zeros

      # geographic bins
      fm <- fm %>%
            # left_join(trend_env %>% select(plot_id, lon, lat) %>% distinct()) %>%
            left_join(ptrends %>% select(plot_id, lon, lat) %>% distinct()) %>%
            mutate(bin = paste(plyr::round_any(lon, 3), plyr::round_any(lat, 3)))

      fm <- fm %>%
            filter(is.finite(drdt_pred), is.finite(drdt_obs))

      list(baseline = fdp,
           trend = fm)

}




# fitted parameter values
load_fits <- function(x, i){

      x %>%
            dplyr::select(-lp__) %>%
            gather(param, value, contains("[")) %>%
            mutate(param = str_remove(param, "\\]"),
                   param = str_replace(param, "\\[", "_")) %>%
            separate(param, c("param", "var")) %>%
            mutate(var = recode(var, "1" = "intercept",
                                "2" = "bacon", "3" = "bahet",
                                "4" = "sulfur", "5" = "nitrogen",
                                "6" = "temperature", "7" = "precipitation",
                                "8" = "dbh")) %>%
            filter(.draw %in% sample(unique(.draw), i)) %>%
            dplyr::select(i = .draw, species = sp, var, param, param_value = value)
}

# observed data
load_obs <- function(x){
      x %>%
            mutate(intercept = 1,
                   year = (year + year_next) / 2#,
                   # plot_id = str_sub(plot_id, 1, -3)
            ) %>%
            rename(obs = outcome) %>%
            dplyr::select(species, plot_id, tree_id, year, obs, t,
                          intercept, dbh = dia, bacon, bahet, sulfur, nitrogen, temperature = bio1, precipitation = bio12) %>%
            gather(var, value, intercept:precipitation)
}

baseline_pred <- function(fi, od){
      left_join(fi, od) %>%
            group_by(i, species, plot_id, tree_id, year, t) %>%
            summarize(obs = mean(obs),
                      pred = sum(param_value * value),
                      .groups = "drop")
}

baseline_pred_rec <- function(fi, od){
      left_join(fi, od) %>%
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
            left_join(e) %>%
            na.omit()

      # add trends in diameter and intercept
      ode <- od %>%
            filter(var %in% c("dbh", "intercept")) %>%
            group_by(species, plot_id, tree_id, var) %>%
            reframe(value1 = mean(value), # could project here instead of averaging
                    t = sum(t),
                    value = project(year, value, c(2000, 2020)),
                    year = c(2000, 2020)) %>%
            mutate(value = ifelse(is.finite(value), value, value1)) %>%
            dplyr::select(-value1) %>%
            bind_rows(ode)
      # ode

      ode %>%
            group_by(species, plot_id, tree_id, year) %>%
            mutate(nrows = n()) %>% ungroup() %>% #count(nrows)
            filter(nrows == nterms) %>%
            dplyr::select(-nrows)
}


