library(drord)
library(testthat)

context("basic tests")

test_that("basic call with wald cis", {
  set.seed(1234)
  n <- 200
  covar <- data.frame(x1 = runif(n), x2 = runif(n))
  treat <- rbinom(n, 1, 1/2)
  out <- rbinom(n, 3, plogis(-1 + covar$x1 - covar$x2 + 0.2 * treat))
  # call drord and obtain wald inference
  fit <- drord(out, treat, covar, ci = "wald")
  expect_true(class(fit) == "drord")
  expect_true(length(fit$log_odds$est) == 3)
  expect_true(all(!is.na(fit$log_odds$est)))
  expect_true(is.null(fit$log_odds$bca))

  expect_true(length(fit$weighted_mean$est$est) == 3)
  expect_true(all(!is.na(fit$weighted_mean$est$est)))
  expect_true(is.null(fit$weighted_mean$bca))

  expect_true(length(fit$mann_whitney$est) == 1)
  expect_true(all(!is.na(fit$mann_whitney$est)))
  expect_true(is.null(fit$mann_whitney$bca))
})

test_that("basic call with bca cis", {
  set.seed(1234)
  n <- 200
  covar <- data.frame(x1 = runif(n), x2 = runif(n))
  treat <- rbinom(n, 1, 1/2)
  out <- rbinom(n, 3, plogis(-1 + covar$x1 - covar$x2 + 0.2 * treat))
  # call drord and obtain bca inference
  fit <- drord(out, treat, covar, ci = "bca", nboot = 10)
  expect_true(class(fit) == "drord")
  expect_true(length(fit$log_odds$est) == 3)
  expect_true(all(!is.na(fit$log_odds$est)))
  expect_true(is.null(fit$log_odds$wald))

  expect_true(length(fit$weighted_mean$est$est) == 3)
  expect_true(all(!is.na(fit$weighted_mean$est$est)))
  expect_true(is.null(fit$weighted_mean$wald))

  expect_true(length(fit$mann_whitney$est) == 1)
  expect_true(all(!is.na(fit$mann_whitney$est)))
  expect_true(is.null(fit$mann_whitney$wald))
})

test_that("basic call to stratified estimator", {
  set.seed(1234)
  n <- 200
  covar <- data.frame(x1 = rbinom(n, 2, 0.33))
  treat <- rbinom(n, 1, 1/2)
  out <- rbinom(n, 3, plogis(-1 + covar$x1 + 0.2 * treat))
  # call drord and obtain stratified estimator
  fit <- drord(out, treat, covar, ci = "wald", stratify = TRUE)
  # call drord and obtain unadjusted stratified estimator
  fit_ua <- drord(out, treat, covar, ci = "wald", 
                  out_form = "1", stratify = TRUE)
  
  expect_true(class(fit) == "drord")
  expect_true(length(fit$log_odds$est) == 3)
  expect_true(all(!is.na(fit$log_odds$est)))
  expect_true(is.null(fit$log_odds$bca))

  expect_true(length(fit$weighted_mean$est$est) == 3)
  expect_true(all(!is.na(fit$weighted_mean$est$est)))
  expect_true(is.null(fit$weighted_mean$bca))

  expect_true(length(fit$mann_whitney$est) == 1)
  expect_true(all(!is.na(fit$mann_whitney$est)))
  expect_true(is.null(fit$mann_whitney$bca))
})

