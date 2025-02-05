
evaluate <- function(recruitment_pred, growth_pred, mortality_pred){
      library(wCorr)

      add_model <- function(x, tag) lapply(x, function(y) mutate(y, model = tag))
      growth_pred <- add_model(growth_pred, "growth")
      mortality_pred <- add_model(mortality_pred, "mortality")
      recruitment_pred <- add_model(recruitment_pred, "recruitment")

      baseline <- bind_rows(growth_pred$baseline,
                            recruitment_pred$baseline,
                            mortality_pred$baseline)
      trend <- bind_rows(growth_pred$trend,
                         recruitment_pred$trend,
                         mortality_pred$trend)

      # baseline

      hl <- baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t),
                      pred = multi2ann(pred, t),
                      obs = multi2ann(obs, t))

      spm <- baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t),
                      pred = multi2ann(pred, t),
                      obs = multi2ann(obs, t))

      p <- hl %>%
            ggplot(aes(pred, obs, group = species)) +
            facet_wrap(~model, scales = "free", nrow = 2) +
            geom_abline(slope = 1, intercept = 0, color = "dodgerblue") +
            geom_line() +
            geom_point(data = spm, color = "red") +
            coord_flip() +
            theme_minimal() +
            labs(x = "predicted mortality rate (decile bins)",
                 y = "observed mortality rate within predicted bin")
      ggsave("figures/eval/baseline_scatter.pdf", p, width = 10, height = 10, units = "in")


      # trends

      p <- trend %>%
            group_by(model) %>%
            mutate(bin = decile(drdt_pred)) %>%
            group_by(model, bin) %>%
            summarize(drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      tree_yrs = sum(tree_yrs)) %>%
            ggplot(aes(drdt_pred, drdt_obs)) +
            facet_wrap(~model, scales = "free", nrow = 2) +
            geom_vline(xintercept = 0, color = "gray") +
            geom_hline(yintercept = 0, color = "gray") +
            geom_abline(slope = 1, intercept = 0, color = "dodgerblue") +
            geom_line() +
            geom_point(aes(size = tree_yrs)) +
            theme_minimal() +
            coord_flip() +
            labs(x = "predicted rate of change",
                 y = "observed rate of change",
                 size = "tree-years")
      ggsave("figures/eval/trend_scatter_overall.pdf", p, width = 10, height = 10, units = "in")

      p <- trend %>% group_by(model, species) %>% scat()
      ggsave("figures/eval/trend_scatter_species.pdf", p, width = 10, height = 10, units = "in")

      p <- trend %>% group_by(model, bin) %>% scat()
      ggsave("figures/eval/trend_scatter_geog.pdf", p, width = 10, height = 10, units = "in")

      p <- trend %>% group_by(model, species, bin) %>% scat()
      ggsave("figures/eval/trend_scatter_spgeog.pdf", p, width = 10, height = 10, units = "in")

      trend
}





scat <- function(trend){
      means <- trend %>%
            group_by(model) %>%
            summarize(drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      r = NA, n = NA)
      tr <- trend %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      n = sum(tree_yrs))
      ggplot(tr, aes(drdt_obs, drdt_pred, weight = n, color = r)) +
            facet_wrap(~model, scales = "free", nrow = 2) +
            geom_vline(xintercept = 0, color = "gray") +
            geom_hline(yintercept = 0, color = "gray") +
            geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
            geom_smooth(method = lm, color = "black") +
            geom_point(aes(size = n)) +
            geom_point(data = means, size = 5, color = "black") +
            # annotate(geom = "point", size = 5,
            #          x = weighted.mean(trend$drdt_obs, trend$tree_yrs),
            #          y = weighted.mean(trend$drdt_pred, trend$tree_yrs)) +
            scale_color_gradient2(mid = "gray") +
            theme_minimal() +
            labs(x = "observed rate of change",
                 y = "predicted rate of change",
                 size = "tree-years",
                 color = "within-group\nHL correlation")
}

scat_pct <- function(x){
      x %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      pred = mean(pred_2000 + pred_2020)/2,
                      obs = mean(obs_2000 + obs_2020)/2,
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      drdtp_pred = drdt_pred / pred,
                      drdtp_obs = drdt_obs / obs,
                      n = sum(tree_yrs)) %>%
            ggplot(aes(drdtp_obs, drdtp_pred, weight = n, color = r)) +
            geom_vline(xintercept = 0, color = "gray") +
            geom_hline(yintercept = 0, color = "gray") +
            geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
            geom_smooth(method = lm, color = "black") +
            geom_point(aes(size = n)) +
            annotate(geom = "point", size = 5,
                     x = weighted.mean(fm$drdt_obs, fm$tree_yrs) / weighted.mean(fm$obs_2000/2 + fm$obs_2020/2, fm$tree_yrs),
                     y = weighted.mean(fm$drdt_pred, fm$tree_yrs) / weighted.mean(fm$pred_2000/2 + fm$pred_2020/2, fm$tree_yrs)) +
            scale_color_gradient2(mid = "gray") +
            scale_x_continuous(labels = scales::percent) +
            scale_y_continuous(labels = scales::percent) +
            theme_minimal() +
            labs(x = "observed rate of change",
                 y = "predicted rate of change",
                 size = "tree-years",
                 color = "within-group\nHL correlation")
}
