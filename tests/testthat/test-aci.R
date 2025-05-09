library(testthat)
library(tuneR)
library(seewave)
library(tibble)
library(CGSsounds)

# Assuming aci function is loaded from your package

# Generate a test wave object
data("birds_tanzania")

# Test 1: Basic functionality
test_that("Test ACI with default parameters", {
  result <- aci(birds_tanzania)
  expect_type(result, "list")
  expect_true(nrow(result) == 1)
  expect_true(ncol(result) >= 1) # Adjust depending on the expected number of columns
})

# Test 2: Frequency resolution impact
test_that("Test ACI with different frequency resolutions", {
  result_low_res <- aci(birds_tanzania, freq.res = 100)
  result_high_res <- aci(birds_tanzania, freq.res = 10)
  expect_true(result_low_res$value_l != result_high_res$value_l)
})

# Test 3: Noise reduction options
test_that("Test ACI with noise reduction", {
  result_no_nr <- aci(birds_tanzania, noise.red = 0)
  result_with_nr <- aci(birds_tanzania, noise.red = 1)
  expect_true(result_no_nr$value_l != result_with_nr$value_l)
})

# Test 4: Window function impact
test_that("Test ACI with different window functions", {
  result_hanning <- aci(birds_tanzania, win.fun = "hanning")
  result_blackman <- aci(birds_tanzania, win.fun = "blackman")
  expect_true(result_hanning$value_l != result_blackman$value_l)
})

# Test 5: Error handling for non-numeric min.freq
test_that("Test error handling for non-numeric frequency limits", {
  expect_error(aci(birds_tanzania, min.freq = "zero"), "min.freq is not a valid number")
})

# Run all tests
test_dir("tests/testthat")
