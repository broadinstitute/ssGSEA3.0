test_that("bootstrap_es_pvalue: B=0 returns the observed ES with p=1", {
  res <- bootstrap_es_pvalue(0.42, numeric(0))
  expect_equal(res$NES, 0.42)
  expect_equal(res$p.value, 1)
})

test_that("bootstrap_es_pvalue: obs.es centered in a wide null gives a large p-value", {
  set.seed(1)
  boot.es <- rnorm(1000, mean = 0, sd = 1)
  res <- bootstrap_es_pvalue(0.01, boot.es)
  expect_gt(res$p.value, 0.5)
})

test_that("bootstrap_es_pvalue: obs.es far from a tight null gives a small p-value", {
  boot.es <- c(rep(0.01, 500), rep(-0.01, 500)) ## tight null around 0
  res <- bootstrap_es_pvalue(5, boot.es)
  expect_lt(res$p.value, 1e-6)
  expect_gt(res$NES, 1) ## positive ES normalized by mean of positive bootstrap replicates
})

test_that("bootstrap_es_pvalue matches hand-computed values for a known input", {
  boot.es <- c(-2, -1, 0, 1, 2, 3) ## mean = 0.5
  obs.es <- 4

  centered <- boot.es - mean(boot.es) ## -2.5 -1.5 -0.5 0.5 1.5 2.5
  se_boot <- sd(centered)
  z <- obs.es / se_boot
  pos.phi <- centered[centered >= 0] ## 0.5 1.5 2.5
  expected_NES <- obs.es / mean(pos.phi)
  expected_p <- 2 * pnorm(z, lower.tail = FALSE)

  res <- bootstrap_es_pvalue(obs.es, boot.es)
  expect_equal(res$NES, expected_NES)
  expect_equal(res$p.value, expected_p)
})

test_that("bootstrap_es_pvalue falls back to 0.5 when the (centered) bootstrap distribution has no same-sign replicates", {
  ## a constant bootstrap distribution centers to all-zero; ">= 0" always
  ## matches (so the positive-side fallback is unreachable), but "< 0" never
  ## does -- so a negative obs.es hits the neg.phi fallback
  boot.es <- rep(3, 10)
  res <- bootstrap_es_pvalue(-5, boot.es)
  expect_equal(res$NES, -5 / 0.5)
  expect_equal(res$p.value, 0) ## z_score = -Inf -> pnorm(-Inf) = 0
})

test_that("bootstrap_es_pvalue propagates NA if any bootstrap replicate is NA (matches original, no na.rm)", {
  boot.es <- c(1, 2, NA, -1, -2)
  res <- bootstrap_es_pvalue(3, boot.es)
  expect_true(is.na(res$NES))
  expect_true(is.na(res$p.value))
})
