.getMode <- function(beta, n) {
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta,
                     method = "meanshift", kernel = "gaussian"
                   ) / sqrt(n))

  mode
}


.estimateNBDispersion <- function(counts,
                                  design,
                                  normalization = "TMM") {
  require(edgeR)
  d <- edgeR::DGEList(counts)
  if (normalization == "TMM") d <- edgeR::calcNormFactors(d)
  d <- edgeR::estimateDisp(d, design)
  return(d$tagwise.dispersion)
}


### NON-PARAMETRIC BOOTSTRAP
### for vector of coefs
.calcBias <- function(beta, ids, n) {
  .getMode(beta[ids], n = n)
}
.nonParametricBootBeta <- function(beta, n, R = 4e3) {
  library(boot)
  bootRes <- boot(beta, # we are resampling rows in 'boot'
    .calcBias,
    R = 4e3,
    n = n
  )
  varMode <- var(bootRes$t[, 1])
  return(varMode)
}

### PARAMETRIC BOOTSTRAP
.parametricBootstrap <- function(beta,
                                 design,
                                 sigma2,
                                 weights,
                                 n,
                                 R = 4e3,
                                 L = NULL) {

  X <- design

  ## loop over populations to get variance-covariance
  covBetaList <- list()
  for (pp in seq_len(nrow(beta))) {
    curW <- diag(weights[pp, ])
    curCovBeta <- sigma2[pp] * solve(t(X) %*% curW %*% X) # agrees with topTable
    covBetaList[[pp]] <- curCovBeta
  }

  ## simulate for each population
  if (is.null(L)) {
    simBetaList <- list()
    for (pp in seq_len(nrow(beta))) {
      curSimBeta <- mixtools::rmvnorm(n = R, mu = beta[pp, ],
                                      sigma = covBetaList[[pp]])
      simBetaList[[pp]] <- curSimBeta
    }
  } else {
    simBetaList <- list()
    for (pp in seq_len(nrow(beta))) {
      curSimBeta <- mixtools::rmvnorm(n = R, mu = beta[pp, ],
                                      sigma = t(L) %*% covBetaList[[pp]] %*% L)
      simBetaList[[pp]] <- curSimBeta
    }
  }

  ## calculate mode for each bootstrap iteration
  modeMat <- matrix(NA, nrow = R, ncol = ncol(beta))
  for (bb in 1:R) {
    curBeta <- do.call(rbind, lapply(simBetaList, function(x) x[bb, ]))
    curModes <- apply(curBeta, 2, function(x) .getMode(x, n))
    modeMat[bb, ] <- curModes
  }

  ## calculate covariance
  covMat <- matrix(NA,
    nrow = nrow(beta), ncol = ncol(beta),
    dimnames = dimnames(beta)
  )
  for (pp in seq_along(length(simBetaList))) {
    for (cc in seq_len(ncol(modeMat))) {
      covMat[pp, cc] <- cov(simBetaList[[pp]][, cc], modeMat[, cc])
    }
  }

  ## variance of mode
  varMode <- apply(modeMat, 2, var)
  names(varMode) <- colnames(beta)

  return(list(
    "varMode" = varMode,
    "covMode" = covMat
  ))
}
