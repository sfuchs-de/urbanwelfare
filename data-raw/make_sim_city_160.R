# Generate built-in 160-tract synthetic city dataset (sim_city_160)
pkgload::load_all("urbanwelfare")
sim_city_160 <- sim_urban_city(
  N = 160L,
  T = 8L,
  K = 4L,
  n_visitor_groups = 6L,
  model = "ces",
  sigma = 3.4555,
  tourism_share = 0.12,
  shock_scale = 1.0,
  noise_sd = 0.0,
  seed = 2026L
)
save(sim_city_160, file = "urbanwelfare/data/sim_city_160.rda", compress = "xz")
cat("Saved sim_city_160.rda, size:", file.info("urbanwelfare/data/sim_city_160.rda")$size, "bytes\n")
