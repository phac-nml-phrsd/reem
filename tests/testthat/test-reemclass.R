

test_that("reem class basics",{
  
  prms = prm_model_example()
  
  obj = new('reem', 
            name = 'foo', 
            prms = prms)
  
  expect_type(obj, 'S4')
  
  n = names(obj)
  expect_true('name' %in% n)
  expect_true('prms' %in% n)
  expect_true('is.fitted' %in% n)
  expect_true('fcst.prm' %in% n)
  expect_true('fcst.obj' %in% n)
  expect_true('fit.prm' %in% n)
  expect_true('fit.obj' %in% n)
  
  expect_equal(nrow(obj$obs.cl), 0)
  expect_equal(nrow(obj$obs.ww), 0)
  expect_equal(nrow(obj$obs.ha), 0)

  })
