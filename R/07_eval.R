
evaluate <- function(recr_pred, grow_pred, mort_pred){

      # combine results ==========
      add_model <- function(x, tag) lapply(x, function(y) mutate(y, model = tag))
      grow_pred <- add_model(grow_pred, "growth")
      mort_pred <- add_model(mort_pred, "mortality")
      recr_pred <- add_model(recr_pred, "recruitment")

      w <- 10
      h <- 10


      # baseline ==============

      # HL
      mort <- mort_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t),
                      pred = multi2ann(pred, t) %>% logit(),
                      obs = multi2ann(obs, t) %>% logit(),
                      .groups = "drop")
      recr <- recr_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin, t) %>%
            summarize(pred = mean(pred),
                      obs = mean(obs),
                      n = n(),
                      .groups = "drop") %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, n*t) %>% log10(),
                      obs = weighted.mean(obs, n*t) %>% log10(),
                      .groups = "drop")
      grow <- grow_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      t = mean(t),
                      .groups = "drop")
      hl <- bind_rows(grow, recr, mort)

      # means
      mort <- mort_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
                      t = mean(t),
                      pred = multi2ann(pred, t) %>% logit(),
                      obs = multi2ann(obs, t) %>% logit(),
                      .groups = "drop")
      grow <- grow_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
                      t = mean(t),
                      .groups = "drop")
      recr <- recr_pred$baseline %>%
            group_by(model, species) %>%
            summarize(pred = weighted.mean(pred, t) %>% log10(),
                      obs = weighted.mean(obs, t) %>% log10(),
                      weight = sum(t),
                      t = mean(t),
                      .groups = "drop")
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

      trend <- bind_rows(grow_pred$trend,
                         recr_pred$trend,
                         mort_pred$trend)

      p <- trend %>%
            group_by(model) %>%
            mutate(bin = decile(drdt_pred)) %>%
            group_by(model, bin) %>%
            summarize(drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      tree_yrs = sum(tree_yrs),
                      .groups = "drop") %>%
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

      trend <- trend %>%
            mutate(model = ifelse(model == "mortality", "survival", model),
                   model = factor(model, levels = c("recruitment", "growth", "survival")),
                   drdt_pred = ifelse(model == "survival", -drdt_pred, drdt_pred),
                   drdt_obs = ifelse(model == "survival", -drdt_obs, drdt_obs))

      trend
}




combined_scatter <- function(trend){
      # trend <- tar_read(eval)

      norm <- function(n) n / sum(n)

      bin <- function(x){
            x %>%
                  summarize(
                        r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                        drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                        drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                        n = sum(tree_yrs),
                        .groups = "drop") %>%
                  group_by(model) %>%
                  mutate(n = norm(n))
      }

      pd <- trend %>%
            group_by(species, model) %>%
            bin() %>%
            mutate(group = "species")
      pd <- trend %>%
            group_by(bin, model) %>%
            bin() %>%
            mutate(group = "landscapes") %>%
            bind_rows(pd)
      pd <- trend %>%
            group_by(model) %>%
            mutate(bin = as.character(decile(drdt_pred))) %>%
            group_by(bin, model) %>%
            bin() %>%
            mutate(group = "prediction deciles") %>%
            bind_rows(pd)

      pd <- pd %>%
            mutate(group = factor(group, levels = c("species", "landscapes", "prediction deciles")))

      pd <- filter(pd, n > .001)

      pdm <- pd %>%
            group_by(model, group) %>%
            summarize(drdt_obs = weighted.mean(drdt_obs, n),
                      drdt_pred = weighted.mean(drdt_pred, n),
                      .groups = "drop")

      labels <- pd %>%
            group_by(model, group) %>%
            summarize(

                  # within-group
                  ppos = weighted.mean(r > 0, n),

                  # among-group
                  r = weightedCorr(drdt_pred, drdt_obs,
                                   weights = n, method = "Pearson"),
                  p = summary(lm(drdt_pred ~ drdt_obs,
                                 data = bind_cols(drdt_pred, drdt_obs),
                                 weights = n))$coefficients[2,4],
                  s = weighted.mean(sign(drdt_obs) == sign(drdt_pred), n),
                  label = paste0("r = ", round(r, 2), "\n",
                                 "p = ", signif(p, 2),"\n",
                                 "CS = ", round(s*100), "%\n",
                                 "PC = ", round(ppos*100), "%"),

                  drdt_obs = min(drdt_obs),
                  drdt_pred = max(drdt_pred)) %>%
            group_by(model) %>%
            mutate(drdt_obs = min(drdt_obs),
                   .groups = "drop")

      plt <- function(mod){
            ggplot(mapping = aes(drdt_obs, drdt_pred)) +
                  facet_grid(group~model, scales = "free") +
                  geom_vline(xintercept = 0, color = "gray30") +
                  geom_hline(yintercept = 0, color = "gray30") +
                  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
                  geom_label(data = labels %>% filter(model == mod),
                             aes(label = label),
                             hjust = 0, vjust = 1, fill = "forestgreen", color = "white",
                             label.padding = unit(0.25, "lines"), lineheight = .9, size = 3) +
                  geom_smooth(data = pd %>% filter(model == mod),
                              aes(weight = n),
                              method = lm, color = "forestgreen", fill = "forestgreen", alpha = .5) +
                  geom_point(data = pdm %>% filter(model == mod),
                             size = 5, shape = 21, fill = "forestgreen", color = "black") +
                  geom_point(data = pd %>% filter(model == mod),
                             aes(size = n, fill = r), color = "black", shape = 21) +
                  geom_label(data = labels %>% filter(model == mod),
                             aes(label = label),
                             hjust = 0, vjust = 1, fill = NA, color = "white", size = 3,
                             label.padding = unit(0.25, "lines"), lineheight = .9) +
                  scale_fill_gradientn(colors = c("darkred","orangered", "gray", "dodgerblue", "darkblue"),
                                       limits = max(abs(pd$r[is.finite(pd$r)])) * c(-1, 1)) +
                  scale_size_continuous(limits = range(pd$n)) +
                  theme_bw() +
                  theme(strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black"),
                        panel.grid.minor = element_blank()) +
                  labs(x = "observed rate of change",
                       y = "predicted rate of change",
                       size = "proportion of\ntree-years",
                       fill = "within-group\ncorrelation")
      }

      p <- (plt("recruitment") +
                  theme(strip.background.y = element_blank(),
                        axis.title.x = element_blank())) +
            (plt("growth") +
                   scale_x_continuous(breaks = c(-.008, -.004, 0)) +
                   theme(strip.background.y = element_blank(),
                         axis.title.y = element_blank())) +
            (plt("survival") +
                   scale_x_continuous(breaks = c(0, .008, .016)) +
                   theme(axis.title = element_blank())) +
            plot_layout(nrow = 1, guides = "collect")

      ggsave("figures/eval/combined_scatters.pdf", p, width = 10, height = 7, units = "in")
}

# what explains variation in accuracy?
accuracy_pred <- function(eval){

      err_data <- function(x){
            x %>%
                  summarize(correlation = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                            rcorrelation = weightedCorr(drdt_pred, drdt_obs, "Spearman", weights = tree_yrs),
                            accuracy = -sqrt(weighted.mean((drdt_pred - drdt_obs)^2, tree_yrs)),
                            accuracy2 = accuracy / abs(weighted.mean(drdt_obs, tree_yrs)),
                            accuracy_pct = -sqrt(weighted.mean((drdt_pred/weighted.mean(pred_2009/2+pred_2010/2, tree_yrs) -
                                                                      drdt_obs/weighted.mean(obs_2009/2+obs_2010/2, tree_yrs))^2, tree_yrs)),
                            slope = coef(lm(drdt_pred ~ drdt_obs))[2],
                            n_sqrt = sqrt(sum(tree_yrs)),
                            obs = weighted.mean(obs_2009/2 + obs_2010/2, tree_yrs),
                            lon = weighted.mean(lon, tree_yrs),
                            lat = weighted.mean(lat, tree_yrs)) %>%
                  group_by(model) %>%
                  mutate(obs = scales::rescale(obs),
                         n_sqrt = scales::rescale(n_sqrt)) %>%
                  gather(metric, err, correlation:slope) %>%
                  gather(pred, value, n_sqrt:lat) %>%
                  group_by(metric, model) %>%
                  mutate(err = scale(err))
      }

      spp <- eval %>%
            group_by(species, model) %>%
            err_data()
      com <- eval %>%
            group_by(bin, model) %>%
            err_data()
      spp_com <- eval %>%
            group_by(species, bin, model) %>%
            err_data()

      p <- spp %>%
            ggplot(aes(value, err, color = model, fill = model)) +
            facet_grid(metric ~ pred, scales = "free") +
            geom_smooth(alpha = .25) +
            theme_bw() +
            labs(title = "species")
      ggsave("figures/eval/accuracy_predictors_species.pdf", p, width = 7, height = 6, units = "in")

      p <- com %>%
            ggplot(aes(value, err, color = model, fill = model)) +
            facet_grid(metric ~ pred, scales = "free") +
            geom_smooth(alpha = .25) +
            theme_bw() +
            labs(title = "communities")
      ggsave("figures/eval/accuracy_predictors_community.pdf", p, width = 7, height = 6, units = "in")

      p <- spp_com %>%
            ggplot(aes(value, err, color = model, fill = model)) +
            facet_grid(metric ~ pred, scales = "free") +
            geom_smooth(alpha = .25) +
            theme_bw() +
            labs(title = "species x communities")
      ggsave("figures/eval/accuracy_predictors_species-community.pdf", p, width = 7, height = 6, units = "in")



      # is there consistency in accuracy across the three demographic variables? (no, there's not)
      p <- spp %>%
            select(-pred, -value) %>%
            distinct() %>%
            filter(metric %in% c("accuracy_pct", "correlation")) %>%
            group_by(metric, model) %>% mutate(err = rank(err)/length(err)) %>%
            spread(model, err) %>% as.data.frame() %>%
            ecoclim::pairsData(xy_vars = c("growth", "survival", "recruitment"),
                               z_vars = c("species", "metric")) %>% as_tibble() %>%
            ggplot(aes(x_value, y_value)) +
            geom_point() +
            geom_smooth() +
            facet_grid(y_var ~ metric + x_var, scales = "free") +
            labs(x = "rank accuracy for one demographic rate",
                 y = "rank accuracy for another demographic rate",
                 title = "species") +
            theme_bw()
      ggsave("figures/eval/accuracy_pairs_species.pdf", p, width = 9, height = 6, units = "in")

      p <- com %>%
            select(-pred, -value) %>%
            distinct() %>%
            filter(metric %in% c("accuracy_pct", "correlation")) %>%
            group_by(metric, model) %>% mutate(err = rank(err)/length(err)) %>%
            spread(model, err) %>% as.data.frame() %>%
            ecoclim::pairsData(xy_vars = c("growth", "survival", "recruitment"),
                               z_vars = c("species", "metric")) %>% as_tibble() %>%
            ggplot(aes(x_value, y_value)) +
            geom_point() +
            geom_smooth() +
            facet_grid(y_var ~ metric + x_var, scales = "free") +
            labs(x = "rank accuracy for one demographic rate",
                 y = "rank accuracy for another demographic rate",
                 title = "communities") +
            theme_bw()
      ggsave("figures/eval/accuracy_pairs_community.pdf", p, width = 9, height = 6, units = "in")

      p <- spp_com %>%
            select(-pred, -value) %>%
            distinct() %>%
            filter(metric %in% c("accuracy_pct", "correlation")) %>%
            group_by(metric, model) %>% mutate(err = rank(err)/length(err)) %>%
            spread(model, err) %>% as.data.frame() %>%
            ecoclim::pairsData(xy_vars = c("growth", "survival", "recruitment"),
                               z_vars = c("species", "metric")) %>% as_tibble() %>%
            ggplot(aes(x_value, y_value)) +
            geom_point() +
            geom_smooth() +
            facet_grid(y_var ~ metric + x_var, scales = "free") +
            labs(x = "rank accuracy for one demographic rate",
                 y = "rank accuracy for another demographic rate",
                 title = "species x community") +
            theme_bw()
      ggsave("figures/eval/accuracy_pairs_species-community.pdf", p, width = 9, height = 6, units = "in")

}




scat <- function(trend){
      means <- trend %>%
            group_by(model) %>%
            summarize(drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      r = NA, n = NA, .groups = "drop")
      tr <- trend %>%
            summarize(r = weightedCorr(drdt_pred, drdt_obs, "Pearson", weights = tree_yrs),
                      drdt_pred = weighted.mean(drdt_pred, tree_yrs),
                      drdt_obs = weighted.mean(drdt_obs, tree_yrs),
                      n = sum(tree_yrs), .groups = "drop")
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
                      pred = mean(pred_2009 + pred_2010)/2,
                      obs = mean(obs_2009 + obs_2010)/2,
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
                     x = weighted.mean(fm$drdt_obs, fm$tree_yrs) / weighted.mean(fm$obs_2009/2 + fm$obs_2010/2, fm$tree_yrs),
                     y = weighted.mean(fm$drdt_pred, fm$tree_yrs) / weighted.mean(fm$pred_2009/2 + fm$pred_2010/2, fm$tree_yrs)) +
            scale_color_gradient2(mid = "gray") +
            scale_x_continuous(labels = scales::percent) +
            scale_y_continuous(labels = scales::percent) +
            theme_minimal() +
            labs(x = "observed rate of change",
                 y = "predicted rate of change",
                 size = "tree-years",
                 color = "within-group\nHL correlation")
}

