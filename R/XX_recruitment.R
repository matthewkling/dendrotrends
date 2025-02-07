

#
# fit_recruitment_model <- function(dd, species, model_file){
#
#       spp <- species$spp
#
#       dd <- dd %>%
#             filter(species %in% spp) %>%
#             split(.$species)
#
#       model <- cmdstan_model(model_file)
#
#       fit_model <- function(ds){
#
#             ds <- ds %>% filter(is.finite(rec_ba)) %>% na.omit()
#
#             x <- bind_cols(bacon = ds$bacon,
#                            bahet = ds$bahet,
#                            sulfur = ds$sulfur,
#                            nitrogen = ds$nitrogen,
#                            bio1 = ds$bio1,
#                            bio12 = ds$bio12)
#
#             data <- list(N = nrow(ds),
#                          y = ds$rec_ba,
#                          P = ncol(x),
#                          x = x,
#                          t = ds$t)
#
#             fit <- model$sample(data = data,
#                                 iter_warmup = 500, iter_sampling = 500,
#                                 chains = 3)
#
#             fit$draws(format = "draws_df") %>%
#                   as.data.frame() %>% as_tibble() %>%
#                   mutate(sp = ds$species[1])
#       }
#
#       future::plan(multisession, workers = 9)
#       dd %>% future_map_dfr(possibly(fit_model))
# }
