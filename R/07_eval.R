
evaluate <- function(recr_pred, grow_pred, mort_pred){

      # combine results ==========

      add_model <- function(x, tag) lapply(x, function(y) mutate(y, model = tag))
      grow_pred <- add_model(grow_pred, "growth")
      mort_pred <- add_model(mort_pred, "mortality")
      recr_pred <- add_model(recr_pred, "recruitment")

      w <- 10
      h <- 10


      # baseline ==============

      # HL deciles
      mort <- mort_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
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
                      weight = sum(t),
                      .groups = "drop")
      grow <- grow_pred$baseline %>%
            group_by(model, species) %>%
            mutate(bin = decile(pred)) %>%
            group_by(model, species, bin) %>%
            summarize(pred = weighted.mean(pred, t),
                      obs = weighted.mean(obs, t),
                      weight = sum(t),
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


      # weighted R2 (coefficient of determination)
      wR2 <- function(x, y, w = 1){
            tss <- sum(w * (x - mean(x))^2)
            rss <- sum(w * (y-x)^2)
            1 - rss / tss
      }

      R2 <- hl %>% group_by(model, species) %>%
            summarize(value = wR2(pred, obs),
                      weight = sum(weight),
                      stat = "R2")
      r <- bind_rows(mort_pred$baseline, recr_pred$baseline, grow_pred$baseline) %>%
            group_by(model, species) %>%
            summarize(value = cor(pred, obs),
                      weight = sum(t),
                      stat = "r")

      spm2 <- spm %>%
            group_by(model) %>%
            summarize(R2 = wR2(pred, obs),
                      r = cor(pred, obs)) %>%
            gather(stat, value, R2, r) %>%
            mutate(stat = factor(stat, levels = c("r", "R2"),
                                 labels = c("Pearson's r across trees", "R-squared across deciles")))

      p1 <- bind_rows(R2, r) %>%
            mutate(stat = factor(stat, levels = c("r", "R2"),
                                 labels = c("Pearson's r across trees", "R-squared across deciles"))) %>%
            ggplot(aes(value)) +
            facet_grid(model ~ stat, scales = "free") +
            geom_histogram(boundary = 0, binwidth = .05) +
            geom_vline(data = spm2, aes(xintercept = value), color ="red") +
            labs(x = "r                                     R^2", y = "n species")
      p2 <- hl %>%
            ggplot(aes(pred, obs, group = species)) +
            facet_wrap(~ model, ncol = 1, scales = "free", strip.position = "right") +
            geom_abline(slope = 1, intercept = 0, color = "dodgerblue") +
            geom_line() +
            geom_point(data = spm, color = "red", aes(size = weight)) +
            scale_size_continuous(range = c(.2, 2)) +
            coord_flip() +
            labs(x = "predicted rate (decile bins)",
                 y = "observed rate within predicted bin")

      p <- p2 + p1 + plot_layout(nrow = 1, widths = c(1, 2)) &
            theme_minimal() &
            theme(strip.text = element_text(color = "white"),
                  strip.background = element_rect(fill = "black", color = "black"),
                  panel.border = element_rect(color = "black", fill = NA),
                  legend.position = "none")
      ggsave(paste0("figures/eval_baseline_scat_hist.pdf"),
             p, width = 8, height = 8, units = "in")


      # trends ===============

      trend <- bind_rows(grow_pred$trend,
                         recr_pred$trend,
                         mort_pred$trend) %>%
            mutate(model = ifelse(model == "mortality", "survival", model),
                   model = factor(model, levels = c("recruitment", "growth", "survival")),
                   drdt_pred = ifelse(model == "survival", -drdt_pred, drdt_pred),
                   drdt_obs = ifelse(model == "survival", -drdt_obs, drdt_obs))
      trend
}




combined_scatter <- function(trend){

      select <- dplyr::select

      # map of grid cells
      trend %>% select(bin) %>% distinct() %>%
            separate(bin, c("x", "y"), sep = " ") %>%
            ggplot(aes(x, y)) +
            geom_tile()

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
            mutate(group = "species")# %>% filter(drdt_obs > -10)
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
                  R2 = 1 - sum((drdt_pred - drdt_obs)^2 * n) / sum((drdt_obs - mean(drdt_obs))^2 * n),
                  s = weighted.mean(sign(drdt_obs) == sign(drdt_pred), n),
                  label = paste0(#"R2 = ", round(R2, 2), "\n",
                        "r = ", round(r, 2), "\n",
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

      ggsave("figures/eval_combined_scatters.pdf", p, width = 10, height = 7, units = "in")
}
