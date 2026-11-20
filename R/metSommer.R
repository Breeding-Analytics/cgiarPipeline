metSommer <- function(phenoDTfile = NULL, analysisId = NULL, analysisIdgeno = NULL,
                      fixedTerm = list("environment"), randomTerm = list("designation"),
                      covMod = NULL, addG = NULL, nFA = NULL, envsToInclude = NULL,
                      trait = NULL, traitFamily = NULL, useWeights = TRUE,
                      calculateSE = TRUE, heritLB = 0.15, heritUB = 0.95,
                      meanLB = 0, meanUB = Inf, maxIters = 50, verbose = TRUE) {
  if (!requireNamespace("sommer", quietly = TRUE) ||
      utils::packageVersion("sommer") < "4.5.25") {
    stop("metSommer requires sommer >= 4.5.25.", call. = FALSE)
  }
  if (is.null(phenoDTfile) || !is.data.frame(phenoDTfile$predictions) ||
      !is.data.frame(phenoDTfile$status)) {
    stop("Please provide a bioflow phenotype object.", call. = FALSE)
  }
  soa <- phenoDTfile$status$analysisId[phenoDTfile$status$module %in% c("soa", "sta")]
  if (!length(analysisId) || any(!analysisId %in% soa)) {
    stop("Please select valid SOA analysis IDs.", call. = FALSE)
  }
  if (!length(trait) || anyNA(trait) || anyDuplicated(trait)) {
    stop("Please provide distinct traits to analyze.", call. = FALSE)
  }
  for (flag in list(useWeights, calculateSE, verbose)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("useWeights, calculateSE and verbose must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (length(maxIters) != 1L || !is.finite(maxIters) ||
      maxIters < 1 || maxIters != as.integer(maxIters)) {
    stop("maxIters must be a positive integer.", call. = FALSE)
  }
  fixedTerm <- .sommer_terms(fixedTerm, allowNone = TRUE)
  randomTerm <- .sommer_terms(randomTerm)
  if (!length(randomTerm)) stop("Please specify random effects.", call. = FALSE)
  if (is.null(covMod)) covMod <- addG
  if (is.null(covMod)) covMod <- lapply(randomTerm, function(x) rep("none", length(x)))
  if (!is.list(covMod)) covMod <- as.list(covMod)
  if (length(covMod) != length(randomTerm) ||
      any(lengths(covMod) != lengths(randomTerm))) {
    stop("Specify one covariance structure for each variable in every random term.", call. = FALSE)
  }
  covMod <- lapply(covMod, function(x) sub("^none\\.*$", "none", x))
  if (is.null(nFA)) nFA <- rep(list(1L), length(randomTerm))
  if (!is.list(nFA)) nFA <- as.list(rep(nFA, length.out = length(randomTerm)))
  if (length(nFA) != length(randomTerm)) {
    stop("nFA must have one rank per random term.", call. = FALSE)
  }
  allowed <- c("none", "Structure model_diag", "Structure model_us", "Structure model_fa",
               "Relationship structure_GenoA", "Relationship structure_GenoD",
               "Relationship structure_GenoAD", "Relationship structure_Pedigree", "weatherInfo")
  if (anyNA(unlist(covMod)) || any(!unlist(covMod) %in% allowed)) {
    stop("Unknown sommer covariance structure.", call. = FALSE)
  }
  bounds <- lapply(list(heritLB, heritUB, meanLB, meanUB), function(x) {
    if (!is.numeric(x) || !length(x) || anyNA(x) ||
        !length(x) %in% c(1L, length(trait))) {
      stop("Each filtering bound must have length one or length(trait), without NA.", call. = FALSE)
    }
    rep(x, length.out = length(trait))
  })
  if (any(bounds[[1]] > bounds[[2]]) || any(bounds[[3]] > bounds[[4]])) {
    stop("Filtering lower bounds must not exceed upper bounds.", call. = FALSE)
  }
  if (is.null(traitFamily)) traitFamily <- rep("gaussian", length(trait))
  if (length(traitFamily) != length(trait)) {
    stop("Specify one distribution for each trait.", call. = FALSE)
  }
  families <- lapply(traitFamily, .sommer_family)
  base <- phenoDTfile$predictions
  base <- base[which(base$analysisId %in% analysisId), , drop = FALSE]
  keep <- rep(TRUE, nrow(base))
  if ("effectType" %in% names(base)) {
    known <- !is.na(base$effectType) & nzchar(base$effectType)
    keep[known] <- base$effectType[known] == "designation"
  }
  if ("entryType" %in% names(base)) {
    keep <- keep & !grepl("^(mother|father)##", base$entryType, ignore.case = TRUE)
  }
  if (verbose && any(!keep)) message("Excluding ", sum(!keep), " non-designation SOA predictions.")
  base <- base[keep, , drop = FALSE]
  required <- c("trait", "environment", "designation", "predictedValue", "stdError")
  if (!all(required %in% names(base))) {
    stop("SOA predictions are missing required columns.", call. = FALSE)
  }
  config <- phenoDTfile$metadata$tpp_analysis_config
  actual <- vapply(trait, function(tt) {
    if (!is.null(config$trait_map[[tt]])) config$trait_map[[tt]] else tt
  }, character(1))
  if (any(!actual %in% base$trait)) {
    stop("Selected traits are not present in the selected SOA predictions.", call. = FALSE)
  }
  base <- .sommer_metadata(base, phenoDTfile)
  requested <- unique(c(unlist(fixedTerm), unlist(randomTerm)))
  weather <- NULL
  if ("weatherInfo" %in% unlist(covMod) ||
      length(setdiff(requested, c(names(base), "envIndex", "inbreeding")))) {
    weather <- as.data.frame(summaryWeather(phenoDTfile, wide = TRUE))
    if (!nrow(weather)) stop("No weather covariates are available.", call. = FALSE)
    names(weather) <- gsub(" ", "", names(weather), fixed = TRUE)
    weather$environment <- rownames(weather)
    extra <- setdiff(names(weather), names(base))
    base <- merge(base, weather[, c("environment", extra), drop = FALSE],
                  by = "environment", all.x = TRUE, sort = FALSE)
  }
  kernels <- .sommer_kernels(phenoDTfile, analysisIdgeno, base, randomTerm, covMod, weather)
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId), , drop = FALSE]
  outPred <- outMetrics <- outModel <- list()
  analysis <- as.numeric(Sys.time())
  while (analysis %in% phenoDTfile$status$analysisId) analysis <- analysis + 0.001
  for (it in seq_along(trait)) {
    tt <- trait[it]
    at <- actual[it]
    if (verbose) message("Fitting sommer model for ", tt)
    dat <- base[which(base$trait == at), , drop = FALSE]
    allEnvs <- unique(as.character(dat$environment))
    if (!is.null(envsToInclude)) {
      column <- if (tt %in% colnames(envsToInclude)) tt else at
      if (!column %in% colnames(envsToInclude) || is.null(rownames(envsToInclude)) ||
          anyNA(envsToInclude[, column])) {
        stop("envsToInclude must have environment row names and a nonmissing column for ", tt,
             call. = FALSE)
      }
      dat <- dat[dat$environment %in% rownames(envsToInclude)[envsToInclude[, column] > 0], , drop = FALSE]
    }
    if (!is.null(config$env_map[[tt]])) {
      dat <- dat[dat$environment %in% config$env_map[[tt]], , drop = FALSE]
    }
    h2 <- metrics[which(metrics$trait == at & metrics$parameter %in%
                          c("plotH2", "H2", "meanR2", "r2",
                            as.vector(outer(c("plotH2", "H2", "meanR2", "r2"),
                                            c("designation", "mother", "father"), paste, sep = "_")))), ]
    good <- h2$environment[which(is.finite(h2$value) &
                                  h2$value >= bounds[[1]][it] & h2$value <= bounds[[2]][it])]
    dat <- dat[dat$environment %in% good, , drop = FALSE]
    if (!nrow(dat)) stop("Not enough records remain for ", tt, ". Check H2 filters.", call. = FALSE)
    means <- metrics[which(metrics$trait == at &
                             metrics$parameter %in% c("mean", "mean_designation", "mean_mother", "mean_father")), ]
    if (!nrow(means)) {
      means <- stats::aggregate(predictedValue ~ environment, dat, mean)
      names(means)[2] <- "value"
    }
    good <- means$environment[which(is.finite(means$value) &
                                     means$value >= bounds[[3]][it] & means$value <= bounds[[4]][it])]
    dat <- dat[dat$environment %in% good, , drop = FALSE]
    for (j in seq_along(randomTerm)) {
      for (k in seq_along(randomTerm[[j]])) {
        K <- kernels[[j]][[k]]
        if (!is.null(K)) {
          variable <- randomTerm[[j]][k]
          keep <- dat[[variable]] %in% rownames(K)
          if (any(!keep)) warning(tt, ": dropping ", sum(!keep),
                                  " records without a kernel entry for ", variable, call. = FALSE)
          dat <- dat[keep, , drop = FALSE]
        }
      }
    }
    if (nrow(dat) < 3L || length(unique(dat$designation)) < 2L) {
      stop("Not enough records remain for ", tt, ". Check environment and H2/mean filters.", call. = FALSE)
    }
    key <- dat[, c("environment", "designation"), drop = FALSE]
    if (anyDuplicated(key)) {
      stop("Duplicate SOA predictions for ", tt,
           ". Select non-overlapping SOA stamps with one mean per genotype and environment.", call. = FALSE)
    }
    if (any(!is.finite(dat$predictedValue))) {
      stop("SOA predicted values must be finite for ", tt, ".", call. = FALSE)
    }
    familyName <- families[[it]]$family
    y <- dat$predictedValue
    invalid <- switch(familyName, binomial = any(y < 0 | y > 1),
                      poisson = any(y < 0 | y != round(y)),
                      Gamma = any(y <= 0), FALSE)
    if (invalid) stop("SOA values are outside the support of the ", familyName, " family.", call. = FALSE)
    model <- phenoDTfile$modeling
    blups <- model[which(model$trait == at & model$parameter == "designationEffectType" &
                            toupper(model$value) == "BLUP"), , drop = FALSE]
    deregress <- paste(dat$analysisId, dat$environment) %in% paste(blups$analysisId, blups$environment)
    if (any(deregress)) {
      r <- dat$reliability[deregress]
      if (length(r) != sum(deregress) || any(!is.finite(r) | r <= 0 | r > 1)) {
        stop("BLUP inputs require reliabilities in (0,1] for deregression.", call. = FALSE)
      }
      dat$predictedValue[deregress] <- dat$predictedValue[deregress] / r
    }
    if (stats::var(dat$predictedValue) <= 0) stop("No variation for ", tt, ".", call. = FALSE)
    dat$envIndex <- stats::ave(dat$predictedValue, dat$environment, FUN = mean) - mean(dat$predictedValue)
    if ("inbreeding" %in% requested) {
      f <- attr(kernels, "inbreeding")
      dat$inbreeding <- unname(f[as.character(dat$designation)])
    }
    missing <- setdiff(requested, names(dat))
    if (length(missing)) stop("Unknown model variables: ", paste(missing, collapse = ", "), call. = FALSE)
    categorical <- c("environment", "designation", "mother", "father", "gid", "pipeline",
                     "stage", "year", "season", "timepoint", "country", "location",
                     "trial", "study", "management", "entryType")
    for (v in requested) {
      if (anyNA(dat[[v]])) stop("Missing model values in ", v, " for ", tt, ".", call. = FALSE)
      if (v %in% categorical || is.character(dat[[v]]) || is.factor(dat[[v]])) {
        dat[[v]] <- factor(dat[[v]])
        if (any(!nzchar(as.character(dat[[v]])))) {
          stop("Empty model levels in ", v, ".", call. = FALSE)
        }
      }
    }
    e <- new.env(parent = asNamespace("sommer"))
    e$.sommer_indicator <- .sommer_indicator
    fixed <- vapply(fixedTerm, function(x) paste(.sommer_quote(x), collapse = ":"), character(1))
    fixed <- fixed[vapply(fixedTerm, function(x) {
      all(vapply(x, function(v) length(unique(dat[[v]])) > 1L, logical(1)))
    }, logical(1))]
    ff <- stats::as.formula(paste("predictedValue ~", paste(c("1", fixed), collapse = " + ")), env = e)
    random <- character(length(randomTerm))
    termK <- vector("list", length(randomTerm))
    for (j in seq_along(randomTerm)) {
      vars <- randomTerm[[j]]
      cs <- covMod[[j]]
      structured <- which(cs %in% c("Structure model_diag", "Structure model_us", "Structure model_fa"))
      # sommer orders covariance coordinates before the identity/relationship factor.
      order <- c(structured, setdiff(seq_along(vars), structured))
      parts <- character()
      identity <- character()
      Kparts <- list()
      for (k in order) {
        v <- vars[k]
        qv <- .sommer_quote(v)
        if (k %in% structured) {
          if (!is.factor(dat[[v]])) stop("Structured covariance requires a factor: ", v, call. = FALSE)
          nf <- if (is.null(nFA[[j]])) 1L else nFA[[j]]
          if (cs[k] == "Structure model_fa" &&
              (length(nf) != 1L || !is.finite(nf) || nf < 1 ||
               nf != as.integer(nf) || nf >= nlevels(dat[[v]]))) {
            stop("FA rank must be a positive integer below the number of levels.", call. = FALSE)
          }
          parts <- c(parts, switch(cs[k],
            "Structure model_diag" = paste0("dsm(", qv, ")"),
            "Structure model_us" = paste0("usm(", qv, ")"),
            "Structure model_fa" = paste0("fam(", qv, ", k=", nf, ")")))
        } else if (is.numeric(dat[[v]]) && is.null(kernels[[j]][[k]])) {
          parts <- c(parts, paste0("ism(matrix(", qv, ", ncol=1))"))
        } else {
          identity <- c(identity, v)
          K <- kernels[[j]][[k]]
          lv <- levels(factor(dat[[v]]))
          if (is.null(K)) {
            K <- diag(length(lv))
            dimnames(K) <- list(lv, lv)
          }
          else K <- K[lv, lv, drop = FALSE]
          Kparts[[length(Kparts) + 1L]] <- K
        }
      }
      if (length(identity) > 1L) {
        name <- paste0(".sommer_group_", j)
        grid <- expand.grid(lapply(identity, function(v) levels(factor(dat[[v]]))),
                            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        labs <- apply(grid, 1, paste, collapse = ":")
        if (anyDuplicated(labs)) {
          stop("Interaction labels are ambiguous; check ':' characters in model levels.", call. = FALSE)
        }
        observed <- apply(dat[, identity, drop = FALSE], 1, paste, collapse = ":")
        dat[[name]] <- factor(observed, levels = labs)
        parts <- c(parts, paste0("ism(.sommer_indicator(", name, "))"))
        K <- Reduce(function(a, b) kronecker(b, a), Kparts)
        dimnames(K) <- list(labs, labs)
      } else if (length(identity)) {
        parts <- c(parts, paste0("ism(", .sommer_quote(identity), ")"))
        K <- Kparts[[1]]
      } else {
        K <- matrix(1, 1, 1)
      }
      termK[[j]] <- K
      useKernel <- any(vapply(kernels[[j]], Negate(is.null), logical(1)))
      gu <- ""
      if (useKernel) {
        precision <- solve(K)
        attr(precision, "inverse") <- TRUE
        name <- paste0(".sommer_precision_", j)
        assign(name, precision, envir = e)
        gu <- paste0(", Gu=", name)
      }
      random[j] <- paste0("vsm(", paste(parts, collapse = ", "), gu, ")")
    }
    rf <- stats::as.formula(paste("~", paste(random, collapse = " + ")), env = e)
    gaussian <- families[[it]]$family == "gaussian"
    residual <- if (useWeights && gaussian) "~vsm(ism(units), sigma2=1, fixedSigma2=TRUE)" else "~vsm(ism(units))"
    args <- list(fixed = ff, random = rf, rcov = stats::as.formula(residual, env = e),
                 data = dat, family = families[[it]],
                 method = sommer::mmes_method(approximation = if (gaussian) "auto" else "pql"),
                 control = sommer::mmes_control(
                   convergence = list(maxIters = as.integer(maxIters)),
                   output = list(inverse = 0L),
                   execution = list(blasThreads = 1L),
                   diagnostics = list(verbose = verbose, dateWarning = FALSE)))
    if (useWeights) {
      if (any(!is.finite(dat$stdError) | dat$stdError <= 0)) {
        stop("Weighted analysis requires finite, positive SOA standard errors for ", tt, ".", call. = FALSE)
      }
      args$W <- Matrix::Diagonal(x = 1 / dat$stdError^2)
    }
    mix <- do.call(sommer::mmes, args)
    if (!isTRUE(mix$convergence)) warning("sommer did not converge for ", tt, ".", call. = FALSE)
    covariance <- lapply(seq_along(randomTerm), function(j) {
      sommer::covmatrix_mmes(mix, j, se = FALSE)$covariance
    })
    prior <- lapply(seq_along(randomTerm), function(j) kronecker(covariance[[j]], termK[[j]]))
    pred <- .sommer_predictions(mix, dat, randomTerm, covMod, prior, tt, calculateSE, families[[it]], ff, rf)
    outPred[[tt]] <- pred
    m <- list()
    addMetric <- function(parameter, value, method, environment = "across", se = NA_real_) {
      m[[length(m) + 1L]] <<- data.frame(trait = tt, environment = environment,
                                        parameter = parameter, method = method, value = value, stdError = se)
    }
    for (ef in unique(pred$effectType)) {
      pp <- pred[pred$effectType == ef, ]
      addMetric(paste0("mean_", ef), mean(pp$predictedValue), "sum(x)/n")
      addMetric(paste0("sd_", ef), stats::sd(pp$predictedValue), "sd")
      if (any(is.finite(pp$reliability))) {
        addMetric(paste0("r2_", ef), stats::median(pp$reliability, na.rm = TRUE), "1-PEV/G")
      }
    }
    vc <- summary(mix)$varcomp
    addMetric(paste0("VC_", vc$term, "_", vc$parameter), vc$estimate, "REML/PQL", se = vc$StdError)
    addMetric(c("Var_residual", "nEnv", "nEntries"),
              c(sommer::covmatrix_mmes(mix, length(randomTerm) + 1L, se = FALSE)$covariance[1, 1],
                length(unique(dat$environment)), length(unique(dat$designation))), c("REML/PQL", "n", "n"))
    envLevels <- sort(unique(as.character(dat$environment)))
    geneticCov <- matrix(0, length(envLevels), length(envLevels),
                         dimnames = list(envLevels, envLevels))
    geneticCorAvailable <- TRUE
    for (j in seq_along(covariance)) {
      C <- covariance[[j]]
      addMetric(paste0("Var_", paste(randomTerm[[j]], collapse = "_")), mean(diag(C)), "REML/PQL")
      if (any(c("designation", "mother", "father", "gid") %in% randomTerm[[j]])) {
        scale <- mean(diag(termK[[j]]))
        if ("environment" %in% randomTerm[[j]]) {
          if (nrow(C) == length(envLevels) && all(envLevels %in% rownames(C))) {
            geneticCov <- geneticCov + C[envLevels, envLevels, drop = FALSE] * scale
          } else if (nrow(C) == 1L) {
            envKernel <- kernels[[j]][[match("environment", randomTerm[[j]])]]
            if (is.null(envKernel)) {
              geneticCov <- geneticCov + diag(as.numeric(C) * scale, length(envLevels))
            } else {
              geneticCov <- geneticCov + envKernel[envLevels, envLevels, drop = FALSE] *
                as.numeric(C) * scale / mean(diag(envKernel))
            }
          } else {
            geneticCorAvailable <- FALSE
            warning("Environment correlations omitted for compound covariance coordinates in ",
                    paste(randomTerm[[j]], collapse = ":"), ".", call. = FALSE)
          }
        } else if ("envIndex" %in% randomTerm[[j]]) {
          index <- dat$envIndex[match(envLevels, as.character(dat$environment))]
          geneticCov <- geneticCov + tcrossprod(index) * as.numeric(C) * scale
        } else if (nrow(C) == 1L) {
          geneticCov <- geneticCov + as.numeric(C) * scale
        }
      }
    }
    if (geneticCorAvailable && length(envLevels) > 1L && all(diag(geneticCov) > 0) &&
        any(c("environment", "envIndex") %in% unlist(randomTerm))) {
      cor <- stats::cov2cor(geneticCov)
      ij <- which(upper.tri(cor), arr.ind = TRUE)
      gp <- data.frame(designation = rownames(cor)[ij[, 1]], predictedValue = cor[ij],
                       stdError = NA_real_, reliability = NA_real_, trait = tt,
                       effectType = "GenCorrMat", environment = colnames(cor)[ij[, 2]],
                       entryType = "unknown", gid = NA, mother = NA, father = NA, pipeline = NA)
      outPred[[tt]] <- rbind(outPred[[tt]], gp[, names(pred)])
      addMetric(paste(gp$designation, gp$environment, sep = "_"), gp$predictedValue, "GenCorrREML")
    }
    outMetrics[[tt]] <- do.call(rbind, m)
    outModel[[tt]] <- rbind(
      data.frame(trait = tt, environment = "across",
                 parameter = c("fixedFormula", "randomFormula", "family", "convergence", "engine", "predictionScale", "useWeights"),
                 value = c(paste(deparse(ff), collapse = ""), paste(deparse(rf), collapse = ""),
                           paste(families[[it]]$family, families[[it]]$link), mix$convergence,
                           paste("sommer", utils::packageVersion("sommer")),
                           "response (merits); link (coefficients)", useWeights)),
      data.frame(trait = tt, environment = allEnvs, parameter = "includedInMta",
                 value = allEnvs %in% dat$environment))
  }
  pred <- do.call(rbind, outPred)
  if (length(trait) > 1L) {
    p <- pred[pred$effectType == "designation" & pred$environment == "(Intercept)", ]
    wide <- stats::reshape(p[, c("designation", "trait", "predictedValue")],
                           idvar = "designation", timevar = "trait", direction = "wide")
    if (nrow(wide) > 1L && ncol(wide) > 2L) {
      C <- stats::cor(wide[, -1, drop = FALSE], use = "pairwise.complete.obs")
      ij <- which(upper.tri(C), arr.ind = TRUE)
      outMetrics[[length(outMetrics) + 1L]] <- data.frame(
        trait = "CorrTrait", environment = "across",
        parameter = paste(sub("predictedValue.", "", rownames(C)[ij[, 1]], fixed = TRUE),
                          sub("predictedValue.", "", colnames(C)[ij[, 2]], fixed = TRUE), sep = "_"),
        method = "PearsonCorrelation", value = C[ij], stdError = NA_real_)
    }
  }
  tables <- list(predictions = pred, metrics = do.call(rbind, outMetrics),
                 modeling = rbind(do.call(rbind, outModel),
                   data.frame(trait = "inputObject", environment = "general",
                              parameter = "analysisId", value = as.character(analysisId))))
  for (tb in names(tables)) {
    rows <- tables[[tb]]
    rows$module <- "moaSommer"
    rows$analysisId <- analysis
    phenoDTfile[[tb]] <- .sommer_append(phenoDTfile[[tb]], rows)
  }
  phenoDTfile$status <- .sommer_append(phenoDTfile$status,
    data.frame(module = "moaSommer", analysisId = analysis, analysisIdName = NA_character_))
  phenoDTfile
}

.sommer_indicator <- function(x) {
  Matrix::sparseMatrix(i = seq_along(x), j = as.integer(x), x = 1,
                       dims = c(length(x), nlevels(x)),
                       dimnames = list(NULL, levels(x)))
}

.sommer_terms <- function(x, allowNone = FALSE) {
  if (is.null(x)) return(list())
  if (!is.list(x)) x <- as.list(x)
  if (allowNone) x <- lapply(x, function(v) setdiff(v, c("none", "1", "")))
  if (allowNone) x <- x[lengths(x) > 0L]
  if (any(!vapply(x, function(v) is.character(v) && length(v) > 0L &&
                         !anyNA(v) && all(nzchar(v)) && !anyDuplicated(v), logical(1)))) {
    stop("Model terms must be lists of nonempty variable names.", call. = FALSE)
  }
  x
}

.sommer_quote <- function(x) {
  if (any(grepl("`", x, fixed = TRUE))) stop("Backticks are not supported in model variable names.", call. = FALSE)
  paste0("`", x, "`")
}

.sommer_family <- function(x) {
  # Accept the ASReml menu labels, but never evaluate user-supplied R code.
  x <- gsub("[[:space:]]", "", x)
  choices <- c("gaussian", "binomial", "Gamma", "poisson")
  for (f in choices) {
    if (identical(x, f) || identical(x, paste0("asr_", f, "(link='",
        switch(f, gaussian = "identity", binomial = "logit", Gamma = "inverse", poisson = "log"),
        "',dispersion=1.0)"))) {
      return(switch(f, gaussian = stats::gaussian(), binomial = stats::binomial(),
                    Gamma = stats::Gamma(), poisson = stats::poisson()))
    }
  }
  stop("Unsupported trait distribution. Use gaussian, binomial, Gamma or poisson.", call. = FALSE)
}

.sommer_metadata <- function(base, object) {
  meta <- object$metadata$pheno
  if (is.null(meta)) return(base)
  meta <- meta[meta$parameter %in% c("environment", "pipeline", "stage", "year", "season",
                 "timepoint", "country", "location", "trial", "study", "management") &
                 meta$value %in% names(object$data$pheno), , drop = FALSE]
  if (!"environment" %in% meta$parameter) return(base)
  other <- object$data$pheno[, meta$value, drop = FALSE]
  names(other) <- meta$parameter
  other <- unique(other[, c("environment", setdiff(names(other), names(base))), drop = FALSE])
  if (anyDuplicated(other$environment)) {
    stop("Environment metadata is not unique. Check the environment mapping.", call. = FALSE)
  }
  merge(base, other, by = "environment", all.x = TRUE, sort = FALSE)
}

.sommer_append <- function(old, rows) {
  if (is.null(old)) return(rows)
  for (n in setdiff(names(old), names(rows))) rows[[n]] <- rep(NA, nrow(rows))
  for (n in setdiff(names(rows), names(old))) old[[n]] <- rep(NA, nrow(old))
  rbind(old, rows[, names(old), drop = FALSE])
}

.sommer_kernels <- function(object, stamp, base, terms, structures, weather) {
  requested <- unlist(structures)
  M <- NULL
  if (any(grepl("^Relationship structure_Geno", requested))) {
    if (length(stamp) != 1L || is.na(stamp) ||
        !stamp %in% object$status$analysisId[object$status$module == "qaGeno"]) {
      stop("A genomic kernel requires a valid Genotype QA/QC version.", call. = FALSE)
    }
    if (inherits(object$data$geno, "genlight")) {
      hit <- cgiarBase::resolveGenoStamp(object$data$geno_imp, stamp)
      if (length(hit) != 1L) stop("The selected genotype QA/QC version is missing or ambiguous.", call. = FALSE)
      M <- as.matrix(object$data$geno_imp[[hit]])
    } else {
      M <- as.matrix(object$data$geno)
      modifications <- object$modifications$geno
      modifications <- modifications[which(modifications$analysisId %in% stamp), , drop = FALSE]
      if (!is.null(modifications) && nrow(modifications)) {
        M <- cgiarBase::applyGenoModifications(M = M, modifications = modifications)
      }
    }
    if (!is.numeric(M) || !nrow(M) || !ncol(M) || is.null(rownames(M)) ||
        anyDuplicated(rownames(M)) || any(!is.finite(M))) {
      stop("Selected QA/QC markers must be a named, finite numeric matrix.", call. = FALSE)
    }
  }
  scaleKernel <- function(K) {
    d <- mean(diag(K))
    if (!is.finite(d) || d <= 0) stop("The requested relationship kernel has zero or invalid variance.", call. = FALSE)
    K <- K / d
    diag(K) <- diag(K) + 1e-5
    K
  }
  result <- lapply(terms, function(x) vector("list", length(x)))
  inbreeding <- NULL
  for (j in seq_along(terms)) for (k in seq_along(terms[[j]])) {
    v <- terms[[j]][k]
    cs <- structures[[j]][k]
    K <- NULL
    if (grepl("^Relationship structure_Geno", cs)) {
      if (!v %in% names(base)) stop("Unknown genomic grouping variable: ", v, call. = FALSE)
      ids <- intersect(unique(as.character(base[[v]])), rownames(M))
      if (length(ids) < 2L) stop("Fewer than two model levels have marker data for ", v, ".", call. = FALSE)
      markers <- M[ids, , drop = FALSE]
      ploidy <- max(M)
      if (!is.finite(ploidy) || ploidy <= 0) stop("Invalid marker dosages.", call. = FALSE)
      if (cs == "Relationship structure_GenoA") {
        K <- sommer::A.mat(markers / (ploidy / 2) - 1)
      } else if (cs == "Relationship structure_GenoAD") {
        K <- tcrossprod(scale(log(markers + 1), center = TRUE, scale = FALSE))
      } else {
        if (ploidy == 2) {
          K <- sommer::D.mat(markers - 1)
        } else {
          p <- colMeans(markers) / ploidy
          Q <- sweep(t(markers), 1, (ploidy - 1) * p, "*")
          Q <- matrix(choose(ploidy, 2) * p^2, nrow = ncol(markers), ncol = nrow(markers)) -
            Q + t(markers) * (t(markers) - 1) / 2
          K <- crossprod(Q) / sum(choose(ploidy, 2) * p^2 * (1 - p)^2)
          dimnames(K) <- list(ids, ids)
        }
        inbreeding <- stats::setNames(rowMeans(4 * (markers / ploidy) * (1 - markers / ploidy)), ids)
      }
      K <- scaleKernel(K)
    } else if (cs == "Relationship structure_Pedigree") {
      if (!requireNamespace("pedigreemm", quietly = TRUE)) stop("Pedigree kernels require pedigreemm.", call. = FALSE)
      meta <- object$metadata$pedigree
      cols <- meta$value[match(c("designation", "mother", "father"), meta$parameter)]
      if (length(cols) != 3L || anyNA(cols) || !all(cols %in% names(object$data$pedigree))) {
        stop("Map designation, mother and father in pedigree metadata.", call. = FALSE)
      }
      p <- object$data$pedigree[, cols, drop = FALSE]
      p[] <- lapply(p, as.character)
      if (anyNA(p[[1]]) || anyDuplicated(p[[1]])) stop("Pedigree IDs must be unique and nonmissing.", call. = FALSE)
      p[[2]][p[[2]] %in% c("", "0")] <- NA
      p[[3]][p[[3]] %in% c("", "0")] <- NA
      ped <- pedigreemm::editPed(sire = p[[3]], dam = p[[2]], label = p[[1]])
      K <- as.matrix(pedigreemm::getA(pedigreemm::pedigree(sire = ped$sire, dam = ped$dam, label = ped$label)))
    } else if (cs == "weatherInfo") {
      if (is.null(weather) || !nrow(weather)) stop("No weather data for a weather kernel.", call. = FALSE)
      W <- weather[, vapply(weather, is.numeric, logical(1)), drop = FALSE]
      for (n in names(W)) {
        if (any(!is.finite(W[[n]]))) {
          warning("Imputing nonfinite weather values in ", n, " with the column mean.", call. = FALSE)
          W[[n]][!is.finite(W[[n]])] <- mean(W[[n]][is.finite(W[[n]])])
        }
      }
      W <- W[, vapply(W, function(x) is.finite(stats::sd(x)) && stats::sd(x) > 0, logical(1)), drop = FALSE]
      if (!ncol(W)) stop("No variable weather covariates for a weather kernel.", call. = FALSE)
      W <- scale(W)
      K <- tcrossprod(W)
      dimnames(K) <- list(weather$environment, weather$environment)
      K <- scaleKernel(K)
    }
    result[[j]][k] <- list(K)
  }
  attr(result, "inbreeding") <- inbreeding
  result
}

.sommer_predictions <- function(mix, dat, terms, structures, prior, trait, calculateSE, family, fixed, random) {
  fixedN <- length(mix$b)
  reference <- .sommer_prediction_grid(dat, terms, structures, all.vars(fixed)[-1L])
  if (!is.null(reference)) {
    dat <- reference
    X <- stats::model.matrix(stats::delete.response(stats::terms(fixed)), dat)
    colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"
    if (!all(rownames(mix$b) %in% colnames(X))) {
      stop("Cannot align fixed prediction columns with the fitted model.", call. = FALSE)
    }
    X <- X[, rownames(mix$b), drop = FALSE]
    expressions <- as.list(attr(stats::terms(random), "variables"))[-1L]
    Z <- lapply(expressions, function(ex) {
      structure <- eval(ex, envir = dat, enclos = environment(random))
      do.call(cbind, structure$Z)
    })
    design <- do.call(cbind, c(list(Matrix::Matrix(X, sparse = TRUE)), Z))
    if (ncol(design) != ncol(mix$W)) stop("Prediction design does not match the fitted model.", call. = FALSE)
  } else {
    design <- mix$W
  }
  n <- nrow(dat)
  # Marginal contrasts include all fitted random contributions and their joint
  # PEV, rather than averaging standard errors from individual environments.
  make <- function(groups, effect, environments = "(Intercept)") {
    ids <- unique(groups)
    A <- Matrix::sparseMatrix(i = match(groups, ids), j = seq_len(n),
                              x = 1 / as.numeric(table(factor(groups, levels = ids)))[match(groups, ids)],
                              dims = c(length(ids), n))
    D <- A %*% design
    p <- sommer::predict.mmes(mix, D = D, VarU = FALSE)
    Dr <- D
    if (fixedN) Dr[, seq_len(fixedN)] <- 0
    offset <- fixedN
    for (j in seq_along(prior)) {
      cols <- offset + seq_len(nrow(prior[[j]]))
      if (!any(c("designation", "mother", "father", "gid") %in% terms[[j]])) Dr[, cols] <- 0
      offset <- max(cols)
    }
    pr <- sommer::predict.mmes(mix, D = Dr, VarU = FALSE)
    variance <- rep(0, length(ids))
    offset <- fixedN
    for (j in seq_along(prior)) {
      cols <- offset + seq_len(nrow(prior[[j]]))
      Z <- Dr[, cols, drop = FALSE]
      variance <- variance + as.numeric(Matrix::rowSums((Z %*% prior[[j]]) * Z))
      offset <- max(cols)
    }
    reliability <- rep(NA_real_, length(ids))
    valid <- variance > 1e-12
    reliability[valid] <- pmax(0, pmin(1, 1 - diag(pr$vcov)[valid] / variance[valid]))
    value <- as.numeric(p$pvals$predicted.value)
    se <- sqrt(pmax(0, diag(p$vcov)))
    if (family$family != "gaussian") {
      se <- se * abs(family$mu.eta(value))
      value <- family$linkinv(value)
    }
    origin <- dat[match(ids, as.character(dat$designation)), , drop = FALSE]
    entry <- if ("entryType" %in% names(origin)) as.character(origin$entryType) else rep("unknown", length(ids))
    entry[is.na(entry) | !nzchar(entry)] <- "unknown"
    out <- data.frame(designation = ids, predictedValue = value,
                      stdError = if (calculateSE) se else NA_real_,
                      reliability = if (calculateSE) reliability else NA_real_,
                      trait = trait, effectType = effect, environment = environments, entryType = entry)
    for (v in c("gid", "mother", "father", "pipeline")) {
      out[[v]] <- if (v %in% names(origin)) origin[[v]] else NA_character_
    }
    out
  }
  result <- list(make(as.character(dat$designation), "designation"))
  # Export fixed coefficients separately for reports and cross-prediction users.
  if (fixedN) {
    D <- Matrix::Matrix(0, fixedN, ncol(design), sparse = TRUE)
    D[, seq_len(fixedN)] <- Matrix::Diagonal(fixedN)
    pp <- sommer::predict.mmes(mix, D = D, VarU = FALSE)
    p <- result[[1]][rep(1L, fixedN), , drop = FALSE]
    p$designation <- rownames(mix$b)
    p$predictedValue <- as.numeric(mix$b)
    p$stdError <- if (calculateSE) sqrt(pmax(0, diag(pp$vcov))) else NA_real_
    p$reliability <- NA_real_
    p$effectType <- vapply(p$designation, function(x) {
      if (x %in% c("(Intercept)", "Intercept")) return("(Intercept)")
      variables <- all.vars(fixed)[-1L]
      hit <- variables[startsWith(x, variables)]
      label <- if (length(hit)) hit[which.max(nchar(hit))] else x
      if (label == "designation") "designation_BLUE" else label
    }, character(1))
    p$entryType <- "unknown"
    p[, c("gid", "mother", "father", "pipeline")] <- NA
    result[[length(result) + 1L]] <- p
  }
  if (any(vapply(terms, function(x) any(c("environment", "envIndex") %in% x) &&
                                 "designation" %in% x, logical(1)))) {
    groups <- paste(dat$environment, dat$designation, sep = ":")
    p <- make(groups, "environment_designation",
              as.character(dat$environment[match(unique(groups), groups)]))
    origin <- dat[match(unique(groups), groups), , drop = FALSE]
    p$entryType <- if ("entryType" %in% names(origin)) as.character(origin$entryType) else "unknown"
    for (v in c("gid", "mother", "father", "pipeline")) {
      if (v %in% names(origin)) p[[v]] <- origin[[v]]
    }
    result[[length(result) + 1L]] <- p
  }
  # Keep component solutions distinct (A, D, parental GCA and FW slopes).
  for (j in seq_along(terms)) if (length(terms[[j]]) == 1L || "envIndex" %in% terms[[j]]) {
    u <- mix$uList[[j]]
    slope <- as.numeric(u)
    if (ncol(u) != 1L) next
    p <- result[[1]][rep(1L, nrow(u)), , drop = FALSE]
    p$designation <- rownames(u)
    origin <- dat[match(p$designation, as.character(dat$designation)), , drop = FALSE]
    for (v in c("gid", "mother", "father", "pipeline")) {
      p[[v]] <- if (v %in% names(origin)) origin[[v]] else NA
    }
    p$entryType <- if ("entryType" %in% names(origin)) as.character(origin$entryType) else "unknown"
    p$entryType[is.na(p$entryType) | !nzchar(p$entryType)] <- "unknown"
    p$predictedValue <- slope
    if ("envIndex" %in% terms[[j]]) {
      p$effectType <- "envIndex_designation"
    } else {
      v <- terms[[j]][1]
      suffix <- switch(structures[[j]][1],
        "Relationship structure_GenoA" = "A", "Relationship structure_GenoD" = "D",
        "Relationship structure_GenoAD" = "AD", "Relationship structure_Pedigree" = "Ped",
        "weatherInfo" = "Weather", "Idv")
      p$effectType <- if (v %in% c("mother", "father")) paste0(v, "_GCA") else paste0(v, suffix)
    }
    D <- Matrix::Matrix(0, nrow(u) * ncol(u), ncol(design), sparse = TRUE)
    offset <- fixedN + sum(vapply(prior[seq_len(j - 1L)], nrow, integer(1)))
    D[, offset + seq_along(slope)] <- Matrix::Diagonal(length(slope))
    prediction <- sommer::predict.mmes(mix, D = D, VarU = FALSE)
    p$stdError <- if (calculateSE) sqrt(pmax(0, diag(prediction$vcov))) else NA_real_
    p$reliability <- NA_real_
    valid <- diag(prior[[j]]) > 1e-12
    if (calculateSE) {
      p$reliability[valid] <- pmax(0, pmin(1, 1 - diag(prediction$vcov)[valid] / diag(prior[[j]])[valid]))
    }
    result[[length(result) + 1L]] <- p
  }
  do.call(rbind, result)
}

.sommer_prediction_grid <- function(dat, terms, structures, fixedVars) {
  ids <- unique(as.character(dat$designation))
  envs <- unique(as.character(dat$environment))
  grid <- expand.grid(designation = ids, environment = envs, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  needed <- unique(c(unlist(terms), fixedVars))
  for (v in setdiff(names(dat), c("designation", "environment"))) {
    values <- dat[[v]]
    byEnv <- !anyDuplicated(unique(dat[, c("environment", v), drop = FALSE])$environment)
    byId <- !anyDuplicated(unique(dat[, c("designation", v), drop = FALSE])$designation)
    if (byEnv) {
      grid[[v]] <- values[match(grid$environment, as.character(dat$environment))]
    } else if (byId) {
      grid[[v]] <- values[match(grid$designation, as.character(dat$designation))]
    } else if (v %in% needed) {
      warning("Using observed-record marginal predictions: ", v,
              " cannot be mapped uniquely to entries or environments.", call. = FALSE)
      return(NULL)
    }
  }
  for (v in intersect(names(dat), names(grid))) {
    if (is.factor(dat[[v]])) grid[[v]] <- factor(grid[[v]], levels = levels(dat[[v]]))
  }
  for (j in seq_along(terms)) {
    name <- paste0(".sommer_group_", j)
    if (name %in% names(dat)) {
      identity <- terms[[j]][!structures[[j]] %in%
                              c("Structure model_diag", "Structure model_us", "Structure model_fa")]
      identity <- identity[!vapply(dat[identity], is.numeric, logical(1))]
      grid[[name]] <- factor(apply(grid[, identity, drop = FALSE], 1, paste, collapse = ":"),
                             levels = levels(dat[[name]]))
    }
  }
  grid
}
