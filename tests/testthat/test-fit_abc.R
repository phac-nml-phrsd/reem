test_that("fit_abc works", {
  date.start <- lubridate::ymd("2022-01-01")
  asof <- date.start + 90
  horizon <- 300

  prms0 <- list(
    horizon = horizon, # horizon of the simulation
    last.obs = horizon - 1, # last observation time (must be < horizon)
    B = data.frame( # Behavior change, buffer start and end dates
      # to compensate for start.delta going between -7 and 7 in fit
      date = seq(date.start - 14, date.start + horizon - 1 + 14, by = 1),
      mult = rep(1, horizon + 28)
    ),
    freq.obs.ww = 3, # average frequency of ww observation
    t.obs.cl = seq(7, 280, by = 7),
    t.obs.ww = seq(3, 200, by = 3),
    i0prop = 1e-3,
    date.start = date.start,
    start.delta = 0,
    R0 = 1.8, # Basic reproduction number
    N = 1e4, # population size
    alpha = 0.5, # transmission heterogeneity (alpha=0: homogeneous)
    I.init = c(1, 1, 3, 5), # initial incidence (overwritten in fit ABC)
    lag = 7, # Aggregation lag for clinical reports
    rho = 0.1, # mean reporting ratio
    g = get_gi(), # Generation interval distribution
    fec = get_fecalshed(), # fecal shedding kinetics
    kappa = 0.18, # decay in ww
    psi = get_psi(), # plug flow simulation,
    shed.mult = 1e-3
  )

  obj0 <- new("reem",
    name = "foo0",
    prms = prms0,
    is.fitted = FALSE
  )

  # Simulate data
  set.seed(1234)
  simepi <- obj0$simulate_epi(deterministic = FALSE)
  obs.cl <- dplyr::filter(simepi$obs.cl, date <= asof)
  obs.ww <- dplyr::filter(simepi$obs.ww, date <= asof)
  obs.ha <- dplyr::filter(simepi$obs.ha, date <= asof)

  # mess with the dates such that they
  # do not perfectly align with the simulation
  obs.cl$t <- obs.cl$t + 1
  obs.cl$date <- obs.cl$date + 1
  obs.ha$t <- obs.ha$t + 1
  obs.ha$date <- obs.ha$date + 1
  obs.ww$t <- obs.ww$t + 1
  obs.ww$date <- obs.ww$date + 1

  # Attached simulated data to new `reem` object:
  prms <- prms0
  prms$t.obs.cl <- obs.cl$t
  prms$t.obs.ww <- obs.ww$t
  obj <- new("reem",
    name = "foo",
    prms = prms,
    obs.cl = obs.cl,
    obs.ha = obs.ha,
    obs.ww = obs.ww,
    is.fitted = FALSE
  )

  # ---- Fit ----

  prm.abc <- list(
    n.abc = 1000,
    n.sim = 0, # `0` for deterministic, else`8` should be enough
    p.abc = 0.02, # 1e-2,
    n.cores = 1, # min(12, parallel::detectCores() - 1),
    use.cl = 1,
    use.ha = 1,
    use.ww = 1,
    err.type = "L2"
  )

  prms.to.fit <- list(
    R0          = list("gamma", 2, 0.251),
    alpha       = list("norm", 1, 1),
    i0prop      = list("unif", -5, -2),
    start.delta = list("unif_int", -7, 7)
  )

  thefit <- obj$fit_abc(prm.abc, prms.to.fit)

  # -- Check posterior values to known parameter values

  R0.post <- thefit$post.prms$R0
  alpha.post <- thefit$post.prms$alpha
  i0.post <- thefit$post.prms$i0prop

  if (0) { # DEBUG
    gg <- obj$plot_fit()
    gg$all
  }
  ci.check <- 0.9
  qt.check <- c(0.5 - ci.check / 2, 0.5 + ci.check / 2)

  ci.R0 <- stats::quantile(R0.post, probs = qt.check)
  ci.alpha <- stats::quantile(alpha.post, probs = qt.check)
  ci.i0 <- stats::quantile(i0.post, probs = qt.check)

  testthat::expect_true(
    ci.R0[1] <= prms0$R0 & prms0$R0 <= ci.R0[2],
    "R0 does not seem to be well fitted."
  )

  testthat::expect_true(
    ci.alpha[1] <= prms0$alpha & prms0$alpha <= ci.alpha[2],
    "`alpha` does not seem to be well fitted."
  )

  testthat::expect_true(
    ci.i0[1] <= log10(prms0$i0) & log10(prms0$R0 <= ci.i0[2]),
    "i0prop does not seem to be well fitted."
  )
})
