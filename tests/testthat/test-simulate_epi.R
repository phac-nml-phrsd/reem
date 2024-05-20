test_that("simulate_epi works", {
  
  date.start = lubridate::ymd('2022-01-01')
  asof       = lubridate::ymd('2022-03-01') 
  
  prms = list(
    horizon = 120,  # horizon of the simulation
    last.obs = 299,  # last observation time (must be < horizon)
    B       = rep(1,300), # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    t.obs.cl = seq(7,280, by = 7),
    t.obs.ww = seq(3,200, by=3),
    date.start = date.start,
    start.delta = 0, 
    R0      = 1.5, # Basic reproduction number
    N       = 1e4, # population size
    alpha   = 0.2, # transmission heterogeneity (alpha=0: homogeneous)
    I.init  = c(1,1,3,5), # initial incidence (overwritten in fit ABC)
    lag     = 7,   # Aggregation lag for clinical reports
    rho     = 0.1, # mean reporting ratio
    g       = get_gi(), # Generation interval distribution
    fec     = get_fecalshed(), # fecal shedding kinetics
    kappa   = 0.18, # decay in ww
    psi     = get_psi(), # plug flow simulation,
    shed.mult = 0.2 # deposited fecal shedding multiplier  
  )
  
  obj = new('reem', 
            name = 'forthetest', 
            prms = prms, 
            is.fitted = FALSE)
  
  expect_true(class(obj) == 'reem')
  expect_true(class(obj$obs.cl) == 'data.frame')
  expect_true(class(obj$obs.ww) == 'data.frame')
  expect_equal(nrow(obj$obs.cl), 0)   # no observations were defined
  expect_equal(nrow(obj$obs.ww), 0)   # no observations were defined
  expect_true(!obj$is.fitted)
  
  simepi  = obj$simulate_epi(deterministic = FALSE)
  expect_type(simepi, 'list')
  expect_equal(length(simepi), 3)
  check.names.simepi = all(names(simepi) %in% c('sim', 'obs.cl', 'obs.ww'))
  expect_true(check.names.simepi)
  
  
  sim = simepi$sim
  expect_true(class(sim) == 'data.frame')
  check.names.sim = all(names(sim) %in% c("t","m","I","S","A","Y","Wd","Wp","Wr","date"))
  expect_true(check.names.sim)
  expect_equal(nrow(sim), prms$horizon)
})