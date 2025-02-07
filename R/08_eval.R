
evaluate <- function(recruitment_pred, growth_pred, mortality_pred){

      # combine results ==========
      add_model <- function(x, tag) lapply(x, function(y) mutate(y, model = tag))
      growth_pred <- add_model(growth_pred, "growth")
      mortality_pred <- add_model(mortality_pred, "mortality")
      recruitment_pred <- add_model(recruitment_pred, "recruitment")

      # baseline <- bind_rows(growth_pred$baseline,
      #                       recruitment_pred$baseline,
      #                       mortality_pred$baseline)

      w <- 10
      h <- 10


      # baseline ==============

      # HL
      mort <- mortality_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t),
                      pred = multi2ann(pred, t) %>% logit(),
                      obs = multi2ann(obs, t) %>% logit())
      recr <- recruitment_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin, t) %>%
            summarize(pred = mean(pred),
                      obs = mean(obs),
                      n = n()) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, n*t) %>% log10(),
                      obs = weighted.mean(obs, n*t) %>% log10())
      grow <- growth_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t))
      hl <- bind_rows(grow, recr, mort)

      # means
      mort <- mortality_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
                      t = mean(t),
                      pred = multi2ann(pred, t) %>% logit(),
                      obs = multi2ann(obs, t) %>% logit())
      grow <- growth_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
                      t = mean(t))
      recr <- recruitment_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t) %>% log10(),
                      obs = weighted.mean(obs, t) %>% log10(),
                      weight = sum(t),
                      t = mean(t))
      spm <- bind_rows(grow, recr, mort)

      p <- hl %>%
            ggplot(aes(pred, obs, group = species)) +
            facet_wrap(~model, scales = "free", ncol = 2) +
            geom_abline(slope = 1, intercept = 0, color = "dodgerblue") +
            geom_line() +
            geom_point(data = spm, color = "red", aes(size = weight)) +
            scale_size_continuous(range = c(.5, 5)) +
            coord_flip() +
            theme_minimal() +
            labs(x = "predicted rate (decile bins)",
                 y = "observed rate within predicted bin")
      ggsave("figures/eval/baseline_scatter.pdf", p, width = w, height = h, units = "in")


      # trends ===============

      trend <- bind_rows(growth_pred$trend,
                         recruitment_pred$trend,
                         mortality_pred$trend)

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
      ggsave("figures/eval/trend_scatter_overall.pdf", p, width = w, height = h, units = "in")

      p <- trend %>% group_by(model, species) %>% scat()
      ggsave("figures/eval/trend_scatter_species.pdf", p, width = w, height = h, units = "in")

      p <- trend %>% group_by(model, bin) %>% scat()
      ggsave("figures/eval/trend_scatter_geog.pdf", p, width = w, height = h, units = "in")

      p <- trend %>% group_by(model, species, bin) %>% scat()
      ggsave("figures/eval/trend_scatter_spgeog.pdf", p, width = w, height = h, units = "in")

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
            facet_wrap(~model, scales = "free", ncol = 2) +
            geom_vline(xintercept = 0, color = "gray") +
            geom_hline(yintercept = 0, color = "gray") +
            geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
            geom_smooth(method = lm, color = "black") +
            geom_point(aes(size = n)) +
            geom_point(data = means, size = 5, color = "black") +
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



z <- function(){

      fm <- bind_rows(read_csv("data/eval_growth.csv") %>% mutate(model = "growth"),
                      read_csv("data/eval_mortality.csv") %>% mutate(model = "mortality"),
                      read_csv("data/eval_recruitment.csv") %>% mutate(model = "recruitment"))


      pd <- fm %>%
            group_by(species, model) %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      n = sum(tree_yrs)) %>%
            group_by(model) %>%
            mutate(n = n / max(n),
                   group = "species")

      pd <- fm %>%
            group_by(bin, model) %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      n = sum(tree_yrs)) %>%
            group_by(model) %>%
            mutate(n = n / max(n),
                   group = "communities") %>%
            bind_rows(pd)

      pd <- fm %>%
            group_by(model) %>%
            mutate(bin = as.character(decile(drdt_pred))) %>%
            group_by(bin, model) %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      n = sum(tree_yrs)) %>%
            group_by(model) %>%
            mutate(n = n / max(n),
                   group = "prediction deciles") %>%
            bind_rows(pd)

      pd <- pd %>%
            mutate(group = factor(group, levels = c("species", "communities", "prediction deciles")))

      pdm <- pd %>%
            group_by(model, group) %>%
            summarize(drdt_obs = weighted.mean(drdt_obs, n),
                      drdt_pred = weighted.mean(drdt_pred, n))

      plt <- function(mod){
            ggplot(mapping = aes(drdt_obs, drdt_pred)) +
                  facet_grid(group~model, scales = "free") +
                  geom_vline(xintercept = 0, color = "gray") +
                  geom_hline(yintercept = 0, color = "gray") +
                  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
                  geom_smooth(data = pd %>% filter(model == mod), aes(weight = n),
                              method = lm, color = "black", alpha = .2) +
                  geom_point(data = pd %>% filter(model == mod), aes(size = n, color = r)) +
                  geom_point(data = pdm %>% filter(model == mod), size = 5) +
                  # scale_color_gradient2(mid = "gray", limits = max(abs(pd$r[is.finite(pd$r)])) * c(-1, 1)) +
                  scale_color_gradientn(colors = c("darkred","orangered", "gray", "dodgerblue", "darkblue"),
                                        limits = max(abs(pd$r[is.finite(pd$r)])) * c(-1, 1)) +
                  scale_size_continuous(limits = 0:1) +
                  theme_bw() +
                  theme(strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black")) +
                  labs(x = "observed rate of change",
                       y = "predicted rate of change",
                       size = "tree-years",
                       color = "within-group\ncorrelation")
      }

      p <- (plt("growth") +
                  coord_cartesian(xlim = c(-.01, NA)) +
                  theme(strip.background.y = element_blank(),
                        axis.title.x = element_blank())) +
            (plt("mortality") +
                   theme(strip.background.y = element_blank(),
                         axis.title.y = element_blank())) +
            (plt("recruitment") +
                   coord_cartesian(xlim = c(-2, NA),
                                   ylim = c(-4, NA)) +
                   theme(axis.title = element_blank())) +
            plot_layout(nrow = 1, guides = "collect")
      ggsave("figures/rgm/eval/combined_scatters.pdf", p, width = 10, height = 7, units = "in")


}
