sommer_fixture <- function() {
  set.seed(140)
  ids <- paste0("g", seq_len(24))
  envs <- paste0("e", seq_len(5))
  p <- expand.grid(designation = ids, environment = envs, stringsAsFactors = FALSE)
  g <- rnorm(length(ids), sd = 2)
  slope <- rnorm(length(ids), sd = 0.4)
  ei <- seq(-2, 2)
  p$predictedValue <- 20 + g[match(p$designation, ids)] +
    ei[match(p$environment, envs)] * (1 + slope[match(p$designation, ids)]) +
    rnorm(nrow(p), sd = 0.5)
  p$stdError <- rep(c(0.3, 0.5, 0.8), length.out = nrow(p))
  p$module <- "soa"
  p$analysisId <- 1
  p$trait <- "yield"
  p$reliability <- 0.8
  p$entryType <- "TEST"
  p$gid <- match(p$designation, ids)
  p$mother <- paste0("m", (p$gid - 1) %% 4 + 1)
  p$father <- paste0("f", (p$gid - 1) %% 6 + 1)
  p$pipeline <- "example"
  p$effectType <- "designation"
  list(data = list(), metadata = list(),
       status = data.frame(module = "soa", analysisId = 1, analysisIdName = "input"),
       predictions = p,
       metrics = data.frame(module = "soa", analysisId = 1, trait = "yield",
                            environment = envs, parameter = "H2", method = "test",
                            value = 0.7, stdError = NA_real_),
       modeling = data.frame(module = "soa", analysisId = 1, trait = "yield",
                             environment = envs, parameter = "designationEffectType", value = "BLUE"))
}

sommer_run <- function(o = sommer_fixture(), ...) {
  metSommer(phenoDTfile = o, analysisId = 1, trait = "yield", verbose = FALSE, ...)
}

sommer_new <- function(x, tb = "predictions") {
  x[[tb]][x[[tb]]$module == "moaSommer", , drop = FALSE]
}

test_that("weighted main effects agree with an independent REML oracle", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  z <- sommer_run(o)
  p <- sommer_new(z)
  p <- p[p$effectType == "designation", ]
  d <- o$predictions
  X <- stats::model.matrix(~environment, d)
  Z <- stats::model.matrix(~designation - 1, d)
  y <- d$predictedValue
  R <- diag(d$stdError^2)
  objective <- function(logVariance) {
    V <- R + exp(logVariance) * tcrossprod(Z)
    ViX <- solve(V, X)
    b <- solve(crossprod(X, ViX), crossprod(ViX, y))
    r <- y - X %*% b
    as.numeric(determinant(V, logarithm = TRUE)$modulus +
                 determinant(crossprod(X, ViX), logarithm = TRUE)$modulus +
                 crossprod(r, solve(V, r))) / 2
  }
  variance <- exp(stats::optimize(objective, c(-8, 8), tol = 1e-8)$minimum)
  V <- R + variance * tcrossprod(Z)
  b <- solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, y)))
  u <- variance * crossprod(Z, solve(V, y - X %*% b))
  expected <- as.numeric(mean(X %*% b) + u)
  names(expected) <- sub("^designation", "", colnames(Z))
  expect_equal(p$predictedValue, unname(expected[p$designation]), tolerance = 2e-3)
  expect_true(all(is.finite(p$stdError) & p$stdError > 0))
  expect_true(all(p$reliability >= 0 & p$reliability <= 1))
  for (column in names(d)) expect_identical(z$predictions[[column]][seq_len(nrow(d))], d[[column]])
  expect_identical(z$status$module, c("soa", "moaSommer"))
  expect_true(all(c("predictions", "metrics", "modeling", "status") %in% names(z)))
})

test_that("weights and SE controls are respected", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  weighted <- sommer_new(sommer_run(o))
  unweighted <- sommer_new(sommer_run(o, useWeights = FALSE, calculateSE = FALSE))
  weighted <- weighted[weighted$effectType == "designation", ]
  unweighted <- unweighted[unweighted$effectType == "designation", ]
  expect_gt(max(abs(weighted$predictedValue - unweighted$predictedValue)), 1e-3)
  expect_true(all(is.na(unweighted$stdError)))
  expect_true(all(is.na(unweighted$reliability)))
  o$predictions$stdError[1] <- 0
  expect_error(sommer_run(o), "positive SOA standard errors")
  expect_no_error(sommer_run(o, useWeights = FALSE))
})

test_that("custom fixed genotype effects remain in the predicted merits", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  z <- sommer_run(o, fixedTerm = list("designation"), randomTerm = list("environment"))
  p <- sommer_new(z)
  p <- p[p$effectType == "designation", ]
  expect_equal(nrow(p), 24L)
  expected <- tapply(o$predictions$predictedValue, o$predictions$designation, mean)
  expect_equal(p$predictedValue, as.numeric(expected[p$designation]), tolerance = 1e-5)
  expect_true(all(is.na(p$reliability)))
})

test_that("GxE covariance structures and FW slopes fit", {
  skip_if_not_installed("sommer", "4.5.25")
  for (cs in c("none", "Structure model_diag", "Structure model_us", "Structure model_fa")) {
    z <- suppressWarnings(sommer_run(
      randomTerm = list("designation", c("environment", "designation")),
      covMod = list("none", c(cs, "none")), nFA = list(1, 1), maxIters = 60))
    p <- sommer_new(z)
    expect_true(all(c("designation", "environment_designation") %in% p$effectType), info = cs)
    expect_true(all(is.finite(p$predictedValue)), info = cs)
    expect_equal(sum(p$effectType == "environment_designation"), 120L)
  }
  diagonal <- sommer_new(sommer_run(randomTerm = list(c("environment", "designation")),
    covMod = list(c("Structure model_diag", "none"))))
  expect_equal(sum(diagonal$effectType == "designation"), 24L)
  fw <- suppressWarnings(sommer_run(randomTerm = list("designation", c("envIndex", "designation"))))
  p <- sommer_new(fw)
  expect_equal(sum(p$effectType == "envIndex_designation"), 24L)
  expect_equal(sum(p$effectType == "environment_designation"), 120L)
  expect_true(all(is.finite(p$predictedValue)))
})

test_that("incomplete trials predict missing cells and include main-effect genetic covariance", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  o$predictions <- o$predictions[-c(1, 26, 51, 76), ]
  z <- sommer_run(o, randomTerm = list("designation", c("environment", "designation")),
                  covMod = list("none", c("Structure model_diag", "none")))
  p <- sommer_new(z)
  expect_equal(sum(p$effectType == "environment_designation"), 120L)
  correlations <- p$predictedValue[p$effectType == "GenCorrMat"]
  expect_length(correlations, 10L)
  expect_true(all(correlations > 0 & correlations <= 1))
  expect_equal(sum(p$effectType == "designation"), 24L)
  cs <- sommer_new(sommer_run(o,
    randomTerm = list("designation", c("environment", "designation"))))
  expect_equal(sum(cs$effectType == "environment_designation"), 120L)
  expect_true(all(is.finite(cs$predictedValue)))
})

test_that("trait and TPP environment filters retain the right provenance", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  o$metadata$tpp_analysis_config <- list(trait_map = c(TPPyield = "yield"),
                                       env_map = list(TPPyield = c("e1", "e2", "e3")))
  envs <- matrix(c(1, 1, 0, 1, 1), ncol = 1,
                 dimnames = list(paste0("e", 1:5), "yield"))
  z <- metSommer(o, analysisId = 1, trait = "TPPyield", envsToInclude = envs, verbose = FALSE)
  included <- sommer_new(z, "modeling")
  included <- included[included$parameter == "includedInMta" & included$value == "TRUE", ]
  expect_setequal(included$environment, c("e1", "e2"))
  expect_identical(unique(sommer_new(z)$trait), "TPPyield")
  o$metadata$tpp_analysis_config <- NULL
  expect_error(sommer_run(o, heritLB = 0.9), "Not enough|no rows")
  expect_error(sommer_run(o, meanLB = 100), "Not enough")
  expect_error(sommer_run(o, traitFamily = "system('false')"), "Unsupported trait distribution")
  expect_error(metSommer(o, 1, trait = "absent"), "not present")
})

test_that("genomic and pedigree kernels do not require ASReml", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  ids <- unique(o$predictions$designation)
  set.seed(18)
  o$data$geno <- matrix(sample(0:2, length(ids) * 80, replace = TRUE), length(ids),
                        dimnames = list(ids, paste0("s", 1:80)))
  o$status <- rbind(o$status, data.frame(module = "qaGeno", analysisId = 2, analysisIdName = "markers"))
  for (cs in c("Relationship structure_GenoA", "Relationship structure_GenoD", "Relationship structure_GenoAD")) {
    z <- suppressWarnings(sommer_run(o, analysisIdgeno = 2, covMod = list(cs)))
    expect_true(all(is.finite(sommer_new(z)$predictedValue)), info = cs)
  }
  ad <- suppressWarnings(sommer_run(o, analysisIdgeno = 2,
    fixedTerm = list("environment", "inbreeding"),
    randomTerm = list("designation", "designation"),
    covMod = list("Relationship structure_GenoA", "Relationship structure_GenoD")))
  expect_true(all(c("designation", "designationA", "designationD", "inbreeding") %in%
                    sommer_new(ad)$effectType))
  expect_error(sommer_run(o, covMod = list("Relationship structure_GenoA")), "Genotype QA/QC")
  o$data$geno[1, 1] <- NA
  expect_error(sommer_run(o, analysisIdgeno = 2, covMod = list("Relationship structure_GenoA")),
               "finite numeric matrix")
  skip_if_not_installed("pedigreemm")
  o$data$pedigree <- data.frame(id = ids, dam = NA_character_, sire = NA_character_)
  o$metadata$pedigree <- data.frame(parameter = c("designation", "mother", "father"),
                                    value = c("id", "dam", "sire"))
  z <- sommer_run(o, covMod = list("Relationship structure_Pedigree"))
  expect_true(all(is.finite(sommer_new(z)$predictedValue)))
})

test_that("Poisson, binomial and Gamma families return response scale predictions", {
  skip_if_not_installed("sommer", "4.5.25")
  for (family in c("poisson", "binomial", "Gamma")) for (weighted in c(FALSE, TRUE)) {
    o <- sommer_fixture()
    set.seed(131)
    eta <- as.numeric(factor(o$predictions$designation)) / 24
    o$predictions$predictedValue <- switch(family,
      poisson = rpois(length(eta), exp(eta)),
      binomial = rbinom(length(eta), 1, stats::plogis(eta)),
      Gamma = rgamma(length(eta), shape = 4, scale = exp(eta) / 4))
    z <- suppressWarnings(sommer_run(o, traitFamily = family, useWeights = weighted))
    p <- sommer_new(z)
    p <- p[p$effectType == "designation", ]
    expect_true(all(is.finite(p$predictedValue)), info = family)
    expect_true(all(p$predictedValue > 0), info = family)
    if (family == "binomial") expect_true(all(p$predictedValue < 1))
  }
})

test_that("parental GCA and SCA use hybrid means rather than parental SOA rows", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  p <- o$predictions[1:5, ]
  p$effectType <- NA_character_
  p$designation <- p$mother
  p$entryType <- "mother##unlabeled"
  o$predictions <- rbind(o$predictions, p)
  z <- sommer_run(o, randomTerm = list("mother", "father", "designation"))
  p <- sommer_new(z)
  expect_equal(sum(p$effectType == "designation"), 24L)
  expect_equal(sum(p$effectType == "mother_GCA"), 4L)
  expect_equal(sum(p$effectType == "father_GCA"), 6L)
  expect_true(all(is.finite(p$predictedValue)))
})

test_that("BLUP deregression is limited to the recorded environment", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  o$modeling$value[o$modeling$environment == "e1"] <- "BLUP"
  explicit <- o
  rows <- explicit$predictions$environment == "e1"
  explicit$predictions$predictedValue[rows] <-
    explicit$predictions$predictedValue[rows] / explicit$predictions$reliability[rows]
  explicit$modeling$value <- "BLUE"
  a <- sommer_new(sommer_run(o))
  b <- sommer_new(sommer_run(explicit))
  expect_equal(a$predictedValue, b$predictedValue, tolerance = 1e-8)
  o$predictions$reliability[1] <- 0
  expect_error(sommer_run(o), "reliabilities")
})

test_that("weather covariance works without markers and retains environmental connectivity", {
  skip_if_not_installed("sommer", "4.5.25")
  o <- sommer_fixture()
  w <- expand.grid(environment = paste0("e", 1:5), day = 1:3, stringsAsFactors = FALSE)
  w$temperature <- 20 + as.numeric(factor(w$environment)) + w$day / 2
  w$rainfall <- 10 - as.numeric(factor(w$environment)) + w$day
  o$data$weather <- w
  o$metadata$weather <- data.frame(parameter = c("environment", "trait", "trait"),
                                    value = c("environment", "temperature", "rainfall"))
  z <- sommer_run(o, randomTerm = list("designation", c("environment", "designation")),
                  covMod = list("none", c("weatherInfo", "none")))
  p <- sommer_new(z)
  expect_equal(sum(p$effectType == "environment_designation"), 120L)
  expect_true(all(is.finite(p$predictedValue)))
})

test_that("multiple traits remain usable by the desired-gain index", {
  skip_if_not_installed("sommer", "4.5.25")
  skip_if_not_installed("enhancer")
  o <- sommer_fixture()
  p <- o$predictions
  p$trait <- "height"
  set.seed(87)
  p$predictedValue <- p$predictedValue * 0.8 + rnorm(nrow(p), sd = 2)
  o$predictions <- rbind(o$predictions, p)
  m <- o$metrics
  m$trait <- "height"
  o$metrics <- rbind(o$metrics, m)
  m <- o$modeling
  m$trait <- "height"
  o$modeling <- rbind(o$modeling, m)
  z <- metSommer(o, analysisId = 1, trait = c("yield", "height"), verbose = FALSE)
  expect_setequal(unique(sommer_new(z)$trait), c("yield", "height"))
  expect_true("CorrTrait" %in% sommer_new(z, "metrics")$trait)
  idx <- indexDesire(z, analysisId = tail(z$status$analysisId, 1),
                     environmentToUse = "(Intercept)", effectTypeToUse = "designation",
                     trait = c("yield", "height"), desirev = c(1, 1), verbose = FALSE)
  expect_equal(tail(idx$status$module, 1), "indexD")
  expect_equal(sum(idx$predictions$module == "indexD"), 24L)
})
