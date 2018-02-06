#' CRISPhieRmix
#'
#' A tool for identifying interesting genes in large pooled CRISPRi and CRISPRa screen
#' using hierarchical mixture models
#'
#' @import sn
#' @import statmod
#'
#' @name CRISPhieRmix
NULL
# > NULL

logSumLogVec <- function(logVec){
  log_max = max(logVec)
  logVec = logVec[-which.max(logVec)]
  return(log_max + log(1 + sum(exp(logVec - log_max))))
}

skewtMixExpectationStep2comp <- function(x, skewtFit, mu, sigma, pq){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dst(x, dp = skewtFit$dp, log = TRUE)
  log_denom = apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)
  return(exp(log_pos_prob - log_denom))
}

skewNormalMixExpectationStep2comp <- function(x, skewNormalFit, mu, sigma, pq){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  log_denom = apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)
  return(exp(log_pos_prob - log_denom))
}

skewtMixMaxStep2comp <- function(x, posProbs, skewtFit = NULL){
  mu = 0
  sigma = 0
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = mean(posProbs*(x - mu)^2)/mean(posProbs)
  if(is.null(skewtFit)){
    skewtFit = sn::st.mple(y = x, w = 1 - posProbs)
  }
  return(list(pq = pq, mu = mu, sigma = sigma, skewtFit = skewtFit))
}

skewNormalMixMaxStep2comp <- function(x, posProbs, skewNormalFit = NULL){
  mu = 0
  sigma = 0
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = mean(posProbs*(x - mu)^2)/mean(posProbs)
  if(is.null(skewNormalFit)){
    skewNormalFit = sn::sn.mple(y = x, w = 1 - posProbs)
  }
  return(list(pq = pq, mu = mu, sigma = sigma, skewNormalFit = skewNormalFit))
}

skewtLogLike <- function(x, pq, skewtFit, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dst(x, dp = skewtFit$dp, log = TRUE)
  return(sum(log(exp(log_pos_prob) + exp(log_null_prob))))
}

skewNormalLogLike <- function(x, pq, skewNormalFit, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  return(sum(log(exp(log_pos_prob) + exp(log_null_prob))))
}

skewtEM <- function(x, skewtFit = NULL, max_iter = 1000, tol = 1e-10,
                    pq = 0.1, mu = 4, sigma = 1, VERBOSE = FALSE){
  n_obs = length(x)
  providedFit = !is.null(skewtFit)
  loglike = -1e100;
  iter = 0;
  posProbs = rep(pq, times = n_obs)
  repeat{
    prevloglike = loglike
    # max step
    if(providedFit){
      updated_params = skewtMixMaxStep2comp(x, posProbs, skewtFit)
    }
    else{
      updated_params = skewtMixMaxStep2comp(x, posProbs, skewtFit = NULL)
    }
    pq = updated_params$pq
    mu = updated_params$mu
    sigma = updated_params$sigma
    skewtFit = updated_params$skewtFit
    # expectation step
    posProbs = skewtMixExpectationStep2comp(x, skewtFit, mu, sigma, pq)

    loglike = skewtLogLike(x, pq, skewtFit, mu, sigma)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if((loglike - prevloglike)/n_obs < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(guide_posterior = posProbs, skewtFit = skewtFit,
              pq = pq, mu = mu, sigma = sigma))
}

skewNormalEM <- function(x, skewNormalFit = NULL, max_iter = 1000, tol = 1e-10,
                    pq = 0.1, mu = 4, sigma = 1, VERBOSE = FALSE){
  n_obs = length(x)
  providedFit = !is.null(skewNormalFit)
  loglike = -1e100;
  iter = 0;
  posProbs = rep(pq, times = n_obs)
  repeat{
    prevloglike = loglike
    # max step
    if(providedFit){
      updated_params = skewNormalMixMaxStep2comp(x, posProbs, skewNormalFit)
    }
    else{
      updated_params = skewNormalMixMaxStep2comp(x, posProbs, skewNormalFit = NULL)
    }
    pq = updated_params$pq
    mu = updated_params$mu
    sigma = updated_params$sigma
    skewNormalFit = updated_params$skewNormalFit
    # expectation step
    posProbs = skewNormalMixExpectationStep2comp(x, skewNormalFit, mu, sigma, pq)

    loglike = skewNormalLogLike(x, pq, skewNormalFit, mu, sigma)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if((loglike - prevloglike)/n_obs < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(guide_posterior = posProbs, skewNormalFit = skewNormalFit,
              pq = pq, mu = mu, sigma = sigma))
}

geneExpectationSkewtMix <- function(x, geneIds, q, p, skewtFit, mu = 5, sigma = 1){
  log_pos_probs = sapply(unique(geneIds),
                         function(i) sum(sapply(x[which(geneIds == i)],
                                                function(y) logSumLogVec(c(log(q) + dnorm(y, mu, sigma, log = TRUE),
                                                                           log(1 - q) + sn::dst(y,
                                                                                                dp = skewtFit$dp,
                                                                                                log = TRUE))))))
  log_null_probs = sapply(unique(geneIds), function(i) sum(sn::dst(x[which(geneIds == i)],
                                                                  dp = skewtFit$dp, log = TRUE)))
  log_denom = sapply(1:length(unique(geneIds)),
                    function(i)
                      logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])))
  return(exp(log(p) + log_pos_probs - log_denom))
}

geneExpectationSkewNormalMix <- function(x, geneIds, q, p, skewNormalFit, mu = 5, sigma = 1){
  log_pos_probs = sapply(unique(geneIds),
                         function(i) sum(sapply(x[which(geneIds == i)],
                                                function(y) logSumLogVec(c(log(q) + dnorm(y, mu, sigma, log = TRUE),
                                                                           log(1 - q) + sn::dsn(y,
                                                                                                dp = skewNormalFit$dp,
                                                                                                log = TRUE))))))
  log_null_probs = sapply(unique(geneIds), function(i) sum(sn::dst(x[which(geneIds == i)],
                                                                  dp = skewNormalFit$dp, log = TRUE)))
  log_denom = sapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])))
  return(exp(log(p) + log_pos_probs - log_denom))
}



integratedGeneExpectationSkewtMix <- function(x, geneIds, skewtFit, mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationSkewtMix(x, geneIds, q.grid[i], pq/q.grid[i],
                                                        skewtFit, mu, sigma),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, mean))
}

integratedGeneExpectationSkewNormalMix <- function(x, geneIds, skewNormalFit, mu = 5,
                                                   sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationSkewNormalMix(x, geneIds, q.grid[i], pq/q.grid[i],
                                                             skewNormalFit, mu, sigma),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, mean))
}

geneExpectationNormalMix <- function(x, geneIds, q = 0.5, p = 0.1,
                                     mu0 = 0, sigma0 = 1, mu = 5, sigma = 1){
  log_null_probs = sapply(unique(geneIds),
                         function(g) sum(dnorm(x[which(geneIds == g)], mu0, sigma0, log = TRUE)))
  log_pos_probs = sapply(unique(geneIds),
                         function(g) sum(sapply(x[which(geneIds == g)],
                                                function(y) logSumLogVec(c(log(q) + dnorm(y, mu, sigma, log = TRUE),
                                                                           log(1 - q) + dnorm(y, mu0, sigma0, log = TRUE))))))
  log_denom = sapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])))
  return(exp(log(p) + log_pos_probs - log_denom))
}


integratedGeneExpectationNormalMix <- function(x, geneIds, mu0 = 0, sigma0 = 1,
                                              mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationNormalMix(x, geneIds, q.grid[i], pq/q.grid[i],
                                                         mu0, sigma0, mu, sigma),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, mean))
}


averageReplicatesSkewtMix <- function(x, geneIds, neg.ctrl = NULL,
                                      skewtFits = NULL, max_iter = 100, tol = 1e-10,
                                      pq = 0.1, mu0 = 4, sigma0 = 1, nMesh = 1000,
                                      VERBOSE = FALSE){
  nReps = dim(x)[2]
  if(!is.null(neg.ctrl)){
    stopifnot(dim(neg.ctrl)[2] == nReps)
    skewtFits = lapply(1:nReps,
                       function(i) sn::st.mple(y = neg.ctrl[,i]))
  }
  skewtMixFits = lapply(1:nReps,
                        function(i) skewtEM(x[,i], skewtFit = skewtFits[[i]],
                                            max_iter = max_iter, tol = tol,
                                            pq = pq, mu = mu0, sigma = sigma0,
                                            VERBOSE = VERBOSE))
  names(skewtMixFits) = paste0("rep", 1:nReps)
  genePosteriors = sapply(1:nReps,
                          function(i) integratedGeneExpectationSkewtMix(x = x[,i], geneIds = geneIds,
                                                                        skewtFit = skewtMixFits[[i]][["skewtFit"]],
                                                                        mu = skewtMixFits[[i]][["mu"]],
                                                                        sigma = skewtMixFits[[i]][["sigma"]],
                                                                        pq = skewtMixFits[[i]][["pq"]],
                                                                        nMesh = nMesh),
                          simplify = TRUE)
  return(list(genes = unique(geneIds),
              score = apply(genePosteriors, 1, mean),
              skewtMixFits = skewtMixFits))
}

averageReplicatesSkewNormalMix <- function(x, geneIds, neg.ctrl = NULL,
                                           skewNormalFits = NULL, max_iter = 100, tol = 1e-10,
                                           pq = 0.1, mu0 = 4, sigma0 = 1, nMesh = 1000,
                                           VERBOSE = FALSE){
  nReps = dim(x)[2]
  if(!is.null(neg.ctrl)){
    stopifnot(dim(neg.ctrl)[2] == nReps)
    skewNormalFits = lapply(1:nReps,
                       function(i) sn::sn.mple(y = neg.ctrl[,i]))
  }
  skewNormalMixFits = lapply(1:nReps,
                             function(i) skewNormalEM(x[,i], skewNormalFit = skewNormalFits[[i]],
                                                      max_iter = max_iter, tol = tol,
                                                      pq = pq, mu = mu0, sigma = sigma0,
                                                      VERBOSE = VERBOSE))
  names(skewNormalMixFits) = paste0("rep", 1:nReps)
  genePosteriors = sapply(1:nReps,
                          function(i) integratedGeneExpectationSkewNormalMix(x = x[,i], geneIds = geneIds,
                                                                             skewNormalFit = skewNormalMixFits[[i]][["skewtFit"]],
                                                                             mu = skewNormalMixFits[[i]][["mu"]],
                                                                             sigma = skewNormalMixFits[[i]][["sigma"]],
                                                                             pq = skewNormalMixFits[[i]][["pq"]],
                                                                             nMesh = nMesh),
                          simplify = TRUE)
  return(list(genes = unique(geneIds),
              score = apply(genePosteriors, 1, mean),
              skewNormalMixFits = skewNormalMixFits))
}

#' normalMix
#'
#' fit a normal hierarchical 2 groups model to log2fc changes from a CRISPRi/a screen
#'
#' @param x log2fc of gene targetting guides
#' @param geneIds gene ids corresponding to x
#' @param pq initial value of p*q
#' @param mu0 initial value of mean for the null distribution, default = 0
#' note that label switching can occur and this should be set with prior knowledge
#' @param sigma0 initial value of sd for the null distribution, default = 1
#' @param mu initial value of the mean for the alternative distribution, default = 5
#' @param sigma initial value of the sd for the alternative distribution, default = 1
#' @param nMesh number of points for computing the gene posterior probabilities, default = 100
#' @param VERBOSE default = FALSE
#' @param PLOT default = FALSE
#'
#' @return a list containing genes, the corresponding posterior probabilities,
#' and the normal mixture fit.
#'
#' @examples
#' Rosenbluh2017CRISPRi.sim.DESeq.log2fc.NormalMix =
#' normalMix(x = Rosenbluh2017CRISPRiSim$x, geneIds = Rosenbluh2017CRISPRiSim$geneIds,
#' mu0 = 0, sigma0 = 0.5, mu = -2, sigma = 1, nMesh = 200)
#'
normalMix <- function(x, geneIds, pq = 0.1, mu0 = 0, sigma0 = 1,
                      mu = 5, sigma = 1, nMesh = 100,
                      VERBOSE = FALSE, PLOT = FALSE){
  require(mixtools)
  normalMixFit = mixtools::normalmixEM(x, k = 2, mu = c(mu0, mu),
                                        sigma = c(sigma0, sigma))
  if(PLOT){
    plot(normalMixFit, density = TRUE)[2]
  }
  genePosteriors = integratedGeneExpectationNormalMix(x, geneIds = geneIds,
                                                      mu0 = normalMixFit[["mu"]][[1]],
                                                      mu = normalMixFit[["mu"]][[2]],
                                                      sigma0 = normalMixFit[["sigma"]][[1]],
                                                      sigma = normalMixFit[["sigma"]][[2]],
                                                      pq = normalMixFit[["lambda"]][[2]],
                                                      nMesh = nMesh)
  return(list(genes = unique(geneIds),
              score = genePosteriors,
              normalMixFit = normalMixFit))
}

empiricalExpectationStep <- function(x, gene_info, null_coefficients,
                                     null_log_norm_factor, mu, sigma, pq){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_neg_prob = log(1 - pq) +
                 apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                 null_coefficients[1] - null_log_norm_factor
  log_denom = apply(cbind(log_pos_prob, log_neg_prob), 1, logSumLogVec)
  return(exp(log_pos_prob - log_denom))
}

empiricalMaxStep <- function(x, gene_info, posProbs){
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = mean(posProbs*(x - mu)^2)/mean(posProbs)
  return(list(pq = pq, mu = mu, sigma = sigma))
}


emprirical2GroupsMixtureLogLike <- function(x, pq, null_coefficients,
                                     null_log_norm_factor, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_neg_prob = log(1 - pq) +
                 apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                 null_coefficients[1] - null_log_norm_factor
  return(sum(log(exp(log_pos_prob) + exp(log_null_prob))))
}


emprical2GroupEMmix <- function(x, gene_info, null_coefficients, null_log_norm_factor,
                                max_iter = 1000, tol = 1e-10, pq = 0.1,
                                mu = 4, sigma = 1, VERBOSE = FALSE){
  n_obs = length(x)
  loglike = emprirical2GroupsMixtureLogLike(x, pq, null_coefficients, null_log_norm_factor, mu, sigma)
  prevloglike = -1e100;
  iter = 0;
  repeat{
    prevloglike = loglike
    # expectation step
    posProbs = empiricalExpectationStep(x, gene_info, null_coefficients,
                                        null_log_norm_factor, mu, sigma, pq)
    # max step
    updated_params = empiricalMaxStep(x, gene_info, posProbs)
    pq = updated_params$pq
    mu = updated_params$mu
    sigma = updated_params$sigma
    loglike = emprirical2GroupsMixtureLogLike(x, pq, null_coefficients, null_log_norm_factor, mu, sigma)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if((loglike - prevloglike) < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(gene = unique(gene_info),
              posProbs = posProbs,
              pq = pq, mu = mu, sigma = sigma,
              null_coefficients = null_coefficients,
              null_log_norm_factor = null_log_norm_factor,
              loglike = loglike))
}

empirical3GroupsExpectationStep <- function(x, gene_info, null_coefficients,
                                            null_log_norm_factor, muNeg, muPos,
                                            sigmaNeg, sigmaPos, qpNeg, qpPos){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigmaPos, log = TRUE)
  log_neg_prob = log(qpNeg) + dnorm(x, mean = muNeg, sd = sigmaNeg, log = TRUE)

  log_null_prob = log(1 - qpPos - qpNeg) +
    apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
    null_coefficients[1] - null_log_norm_factor
  log_denom = apply(cbind(log_pos_prob, log_neg_prob, log_null_prob), 1, logSumLogVec)
  return(list(posProbs = exp(log_pos_prob - log_denom),
              negProbs = exp(log_neg_prob - log_denom)))
}

empirical3GroupsMaxStep <- function(x, posProbs, negProbs){
  qpPos = mean(posProbs)
  qpNeg = mean(negProbs)
  muPos = mean(posProbs*x)/mean(posProbs)
  sigmaPos = mean(posProbs*(x - muPos)^2)/mean(posProbs)
  muNeg = mean(negProbs*x)/mean(negProbs)
  sigmaNeg = mean(negProbs*(x - muPos)^2)/mean(posProbs)  
  return(list(qpPos = qpPos, qpNeg = qpNeg, 
              muPos = muPos, muNeg = muNeg,
              sigmaPos = sigmaPos, sigmaNeg = sigmaNeg))
}


emprirical3GroupsMixtureLogLike <- function(x, null_coefficients, null_log_norm_factor, 
                                            qpPos, qpNeg, muPos, muNeg,
                                            sigmaPos, sigmaNeg){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigmaPos, log = TRUE)
  log_neg_prob = log(qpNeg) + dnorm(x, mean = muNeg, sd = sigmaNeg, log = TRUE)
  
  log_null_prob = log(1 - qpPos - qpNeg) +
    apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
    null_coefficients[1] - null_log_norm_factor

  return(sum(log(exp(log_pos_prob) + exp(log_neg_prob) + exp(log_null_prob))))
}


emprical3GroupEMmix <- function(x, null_coefficients, null_log_norm_factor,
                                max_iter = 1000, tol = 1e-10, qpPos = 0.05, qpNeg = 0.05,
                                muPos = 5, muNeg = -5, sigmaPos = 1, sigmahNeg = 1,
                                VERBOSE = FALSE){
  n_obs = length(x)
  loglike = emprirical3GroupsMixtureLogLike(x, null_coefficients, null_log_norm_factor, 
                                            qpPos, qpNeg, muPos, muNeg,
                                            sigmaPos, sigmaNeg)
  prevloglike = -1e100;
  iter = 0;
  repeat{
    prevloglike = loglike
    # expectation step
    guideExpectations = empirical3GroupsExpectationStep(x, gene_info, null_coefficients,
                                                        null_log_norm_factor, muNeg, muPos,
                                                        sigmaNeg, sigmaPos, qpNeg, qpPos)
    posProbs = guideExpectations$posProbs
    negProbs = guideExpectations$negProbs
    
    # max step
    updated_params = empirical3GroupsMaxStep(x, posProbs, negProbs)
    qpPos = updated_params$qpPos 
    qpNeg = updated_params$qpNeg 
    muPos = updated_params$muPos
    muNeg = updated_params$muNeg
    sigmaPos = updated_params$sigmaPos
    sigmaNeg = updated_params$sigmaNeg

    loglike = emprirical3GroupsMixtureLogLike(x, null_coefficients, null_log_norm_factor, 
                                              qpPos, qpNeg, muPos, muNeg,
                                              sigmaPos, sigmaNeg)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if((loglike - prevloglike) < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(gene = unique(gene_info),
              posProbs = posProbs,
              negProbs = negProbs,
              qpPos = qpPos,
              qpNeg = qpNeg,
              muPos = muPos,
              muNeg = muNeg,
              sigmaPos = sigmaPos, 
              sigmaNeg = sigmaNeg,
              null_coefficients = null_coefficients,
              null_log_norm_factor = null_log_norm_factor,
              loglike = loglike))
}

geneExpectationEmpiricalMix <- function(x, geneIds, q = 0.5, p = 0.1, mu = 5, sigma = 1,
                                        null_coefficients, null_log_norm_factor){
  log_null_probs = sapply(unique(geneIds),
                         function(g) sum(apply(t(null_coefficients[-1]*t(poly(x[which(geneIds == g)], degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                                           null_coefficients[1] - null_log_norm_factor))
  log_pos_probs = sapply(unique(geneIds),
                         function(g) sum(sapply(x[which(geneIds == g)],
                                                function(y) logSumLogVec(c(log(q) + dnorm(y, mu, sigma, log = TRUE),
                                                                           log(1 - q) + apply(t(null_coefficients[-1]*t(poly(y, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                                                                           null_coefficients[1] - null_log_norm_factor)))))
  log_denom = sapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])))
  return(exp(log(p) + log_pos_probs - log_denom))
}

gaussQuadGeneExpectationEmpiricalMix <- function(x, geneIds, null_coefficients, null_log_norm_factor,
                                                 mu = 5, sigma = 1, lowerLim = 0.1, upperLim = 1, 
                                                 nMesh = 100){
  quad.points.weights = statmod::gauss.quad.prob(nMesh, dist = "uniform", l = lowerLim, u = upperLim)
  #q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  #q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:nMesh,
                    function(i) quad.points.weights$weights[i]*geneExpectationEmpiricalMix(x, geneIds, q = quad.points.weights$nodes[i],
                                                                                           p = pq/quad.points.weights$nodes[i],
                                                                                           mu = mu, sigma = sigma, null_coefficients,
                                                                                           null_log_norm_factor),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, sum))
}


integratedGeneExpectationEmpiricalMix <- function(x, geneIds, null_coefficients, null_log_norm_factor,
                                                  mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationEmpiricalMix(x, geneIds, q = q.grid[i], p = pq/q.grid[i],
                                                            mu = mu, sigma = sigma, null_coefficients,
                                                            null_log_norm_factor),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, mean))
}

fitNegCtrl <- function(neg.ctrl, maxDegree = 20, minDegree = 4,
                       breaks = NULL, VERBOSE = FALSE, PLOT = FALSE){
  stopifnot(length(neg.ctrl) > maxDegree)
  if(VERBOSE){
    cat("negCtrl = ", head(neg.ctrl), "\n")
  }
  if(length(breaks) > 1){
    h = hist(neg.ctrl, breaks = breaks, plot = FALSE)
    if(VERBOSE){
      cat("using supplied breaks \n")
    }
  }
  else if(length(breaks) == 1){
    # I find extending the breaks out helps to force the largest term to be negative
    breaks = seq(from = min(neg.ctrl) - 0.5, to = max(neg.ctrl) + 0.5, length = breaks)
    h = hist(neg.ctrl, breaks = breaks, plot = FALSE)
    if(VERBOSE){
      cat("using supplied break length")
    }
  }
  else{
    breaks = seq(from = min(neg.ctrl) - 0.5, to = max(neg.ctrl) + 0.5, length = 101)
    h = hist(neg.ctrl, breaks = breaks, plot = FALSE)
    if(VERBOSE){

    }
  }
  midpoints = (breaks[-length(breaks)] + breaks[-1])/2
  d = breaks[2] - breaks[1]
  N = length(neg.ctrl)
  degree = minDegree
  repeat{
    lindsey.fit = glm(h$counts ~ poly(midpoints, degree = degree, raw = TRUE),
                    family = poisson(link = "log"))
    if(VERBOSE){
      cat("fit of degree ", degree, ":", "\n")
      cat(lindsey.fit$coefficients)
      cat("\n")
    }
    if(lindsey.fit$coefficients[length(lindsey.fit$coefficients)] < 0){
      break
    }
    if(degree == maxDegree){
      break
      cat("problem converging, increase max degree \n")
    }
    degree = degree + 2
  }
  if(PLOT){
    hist(neg.ctrl, breaks = breaks, probability = TRUE)
    lines(breaks, exp(apply(t(lindsey.fit$coefficients[-1]*t(poly(breaks, degree = degree, raw = TRUE))), 1, sum)
                      + lindsey.fit$coefficients[1] - (log(N) + log(d))), col = "red", lwd  = 2, lty = 2)
  }
  return(list(coefficients = lindsey.fit$coefficients,
              log_norm_factor = log(N) + log(d)))
}

setBimodalParams <- function(mu, sigma, pq){
  if(length(mu) == 1){
    muPos = abs(mu)
    muNeg = -abs(mu)
  }
  else{
    muPos = max(mu)
    stopifnot(muPos > 0)
    muNeg = min(mu)
    stopifnot(muNeg < 0)
  }
  if(length(sigma == 1)){
    sigmaPos = sigma
    sigmaNeg = sigma
  }
  else{
    sigmaPos = sigma[which.max(mu)]
    sigmaNeg = sigma[which.min(mu)]
  }
  if(length(pq) == 1){
    qpPos = pq/2
    qpNeg = pq/2
  }
  else{
    qpPos = pq[which.max(mu)]
    qpNeg = pq[which.min(mu)]
  }
  return(list(qpPos = qpPos,
              qpNeg = qpNeg,
              muPos = muPos,
              muNeg = muNeg,
              sigmaPos = sigmaPos, 
              sigmaNeg = sigmaNeg))
}

CRISPhieRmix2Groups <- function(x, geneIds, negCtrlFit,
                                max_iter = 100, tol = 1e-10, 
                                pq = 0.1, mu = -4, sigma = 1,
                                nMesh = 100, maxDegree = 20, 
                                minDegree = 4,breaks = 101,
                                VERBOSE = FALSE, PLOT = FALSE){
  mixFit = emprical2GroupEMmix(x, gene_info = geneIds,
                               null_coefficients = negCtrlFit[["coefficients"]],
                               null_log_norm_factor = negCtrlFit[["log_norm_factor"]],
                               max_iter = max_iter, tol = tol, pq = pq,
                               mu = mu, sigma = sigma, VERBOSE = VERBOSE)
  if(VERBOSE){
    cat("fit mixture \n")
    cat("mu = ", mixFit$mu,
        ", sigma = ", mixFit$sigma,
        ", pq = ", mixFit$pq, "\n")
  }
  if(PLOT){
    b = seq(from = min(x) - 0.1, to = max(x) + 0.1, length = 81)
    hist(x, breaks = b, probability = TRUE, main = "mixture fit to observations")
    lines(b, mixFit$pq*dnorm(b, mixFit$mu, mixFit$sigma), lwd = 2, col = "darkgreen")
    lines(b, (1 - mixFit$pq)*exp(apply(t(mixFit$null_coefficients[-1]*t(poly(b, degree = length(mixFit$null_coefficients) - 1, raw = TRUE))), 1, sum) + mixFit$null_coefficients[1] - mixFit$null_log_norm_factor), col = "red", lwd  = 2)
  }
  
  genePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds,
                                                        negCtrlFit[["coefficients"]],
                                                        negCtrlFit[["log_norm_factor"]],
                                                        mu = mixFit[["mu"]],
                                                        sigma = mixFit[["sigma"]],
                                                        lowerLim = mixFit[["pq"]],
                                                        nMesh = nMesh)
  
  return(list(locfdr = 1 - genePosteriors,
              genePosteriors = genePosteriors,
              mixFit = mixFit))
}

CRISPhieRmix3Groups <- function(x, geneIds, negCtrlFit,
                                max_iter = 100, tol = 1e-10, pq = c(0.05, 0.05), 
                                mu = c(-4, 4), sigma = c(1, 1),
                                nMesh = 100, breaks = 101,
                                VERBOSE = FALSE, PLOT = FALSE){
  params = setBimodalParams(mu, sigma, pq)  
  mixFit = emprical3GroupEMmix(x, null_coefficients = negCtrlFit[["coefficients"]], 
                               null_log_norm_factor = negCtrlFit[["log_norm_factor"]],
                               max_iter = max_iter, tol = tol, 
                               qpPos = params$qpPos, qpNeg = params$qpNeg,
                               muPos = params$muPos, muNeg = params$muNeg, 
                               sigmaPos = params$sigmaPos, sigmahNeg = params$sigmaNeg,
                               VERBOSE = VERBOSE)
  if(VERBOSE){
    cat("fit mixture \n")
    cat("mu = ", mixFit$muPos, "\t", mixFit$muNeg, "\n",
        "sigma = ", mixFit$sigmaPos, "\t", mixFit$sigmaNeg, "\n",
        "pq = ", mixFit$qpPos, "\t", mixFit$qpNeg, "\n")
  }
  if(PLOT){
    b = seq(from = min(x) - 0.1, to = max(x) + 0.1, length = 81)
    hist(x, breaks = b, probability = TRUE, main = "mixture fit to observations")
    lines(b, mixFit$qpPos*dnorm(b, mixFit$muPos, mixFit$sigmaPos) + mixFit$qpNeg*dnorm(b, mixFit$muNeg, mixFit$sigmaNeg), lwd = 2, col = "darkgreen")
    lines(b, (1 - mixFit$qpPos - mixFit$qpNeg)*exp(apply(t(mixFit$null_coefficients[-1]*t(poly(b, degree = length(mixFit$null_coefficients) - 1, raw = TRUE))), 1, sum) + mixFit$null_coefficients[1] - mixFit$null_log_norm_factor), col = "red", lwd  = 2)
  }
  
  # can we compute the following independently? 
  posGenePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds,
                                                           negCtrlFit[["coefficients"]],
                                                           negCtrlFit[["log_norm_factor"]],
                                                           mu = mixFit[["muPos"]],
                                                           sigma = mixFit[["sigmaPos"]],
                                                           lowerLim = mixFit[["qpPos"]],
                                                           upperLim = 1 - mixFit[["qpNeg"]],
                                                           nMesh = nMesh)
  negGenePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds,
                                                           negCtrlFit[["coefficients"]],
                                                           negCtrlFit[["log_norm_factor"]],
                                                           mu = mixFit[["muNeg"]],
                                                           sigma = mixFit[["sigmaNeg"]],
                                                           lowerLim = mixFit[["qpNeg"]],
                                                           upperLim = mixFit[["qpPos"]]
                                                           nMesh = nMesh)
  
  return(list(locfdr = 1 - posGenePosteriors - negGenePosteriors,
              posGenePosteriors = posGenePosteriors,
              negGenePosteriors = negGenePosteriors,
              mixFit = mixFit))
  
}

#' CRISPhieRmix
#'
#' a hierarchical mixture model for analysing large-scale CRISPRi/a pooled screen
#'
#' @param x log2 fold changes of guides targeting genes
#' @param geneIds gene ids corresponding to x
#' @param negCtrl log2 fold changes of negative control guides
#' @param max_iter maximum number of iterations for EM algorithm, default = 100
#' @param tol tolerance for convergence of EM algorithm, default = 1e-10
#' @param pq initial value of p*q, default = 0.1
#' @param mu initial value of mu for the interesting genes, default = 4
#' @param sigma initial value of sigma for the interesting genes, default = 1
#' @param nMesh the number of points to use in numerical integration of posterior probabilities, default = 100
#' @param maxDegree the maximum degree of the Lindsey fit, default = 20
#' @param minDegree the minimum degree of the Lindsey fit, default = 4
#' @param breaks either a number indicating the number of breaks or the sequence of breaks, default = 101
#' @param VERBOSE boolean variable for VERBOSE mode, default = FALSE
#' @param PLOT boolean variable to produce plots, default = FALSE
#' @return a list containing genes, the corresponding posterior probabilities of being non-null,
#' and the mixture fit
#'
#' @return a list containing genes, the corresponding posterior probabilities,
#' and the mixture fit.
#'
#' @author Timothy Daley, \email{tdaley@stanford.edu}
#'
#' @examples
#' Rosenbluh2017CRISPRi.sim.DESeq.log2fc.CRISPhieRmix =
#' CRISPhieRmix(x = Rosenbluh2017CRISPRiSim$x, geneIds = Rosenbluh2017CRISPRiSim$geneIds,
#' negCtrl = Rosenbluh2017CRISPRiSim$negCtrl, mu = -2, sigma = 0.5, nMesh = 200)
#'
#' @export
CRISPhieRmix <- function(x, geneIds, negCtrl = NULL,
                         max_iter = 100, tol = 1e-10, pq = 0.1, mu = -4, sigma = 1,
                         nMesh = 100, maxDegree = 20, minDegree = 4,
                         breaks = 101,
                         BIMODAL = FALSE,
                         VERBOSE = FALSE, PLOT = FALSE){}
  #stopifnot(!is.null(neg.ctrl) | !is.null(empiricalNegCtrlFits))
  if(!is.null(negCtrl)){
    negCtrlFit = fitNegCtrl(as.numeric(negCtrl), maxDegree = maxDegree, minDegree = minDegree,
                            breaks = 101, VERBOSE = VERBOSE, PLOT = PLOT)
    if(VERBOSE){
      cat("fit negative control distributions \n")
    }
    if(BIMODAL){
      if(VERBOSE){
        cat("3 groups \n")
      }
      CRISPhieRmixFit = CRISPhieRmix3Groups(x = x, geneIds = geneIds, negCtrlFit = negCtrlFit,
                                            max_iter = max_iter, tol = tol, pq = pq, 
                                            mu = mu, sigma = sigma, nMesh = nMesh,
                                            breaks = 101, VERBOSE = FALSE, PLOT = FALSE)
      
    } else{
      stopifnot(length(mu) == 1)
      if(VERBOSE){
        cat("2 groups \n")
      }
      CRISPhieRmixFit = CRISPhieRmix2Groups(x, geneIds = geneIds, negCtrlFit = negCtrlFit, 
                                            max_iter = max_iter, tol = tol, pq = pq, 
                                            mu = mu, sigma = sigma, nMesh = nMesh, breaks = breaks,
                                            VERBOSE = FALSE, PLOT = FALSE)
    }
  } else{
    require(mixtools)
    normalMixFit = mixtools::normalmixEM(x, k = 2, mu = c(0, mu),
                                         sigma = c(1, sigma))
    if(PLOT){
      plot(normalMixFit, density = TRUE)[2]
    }
    genePosteriors = integratedGeneExpectationNormalMix(x, geneIds = geneIds,
                                                        mu0 = normalMixFit[["mu"]][[1]],
                                                        mu = normalMixFit[["mu"]][[2]],
                                                        sigma0 = normalMixFit[["sigma"]][[1]],
                                                        sigma = normalMixFit[["sigma"]][[2]],
                                                        pq = normalMixFit[["lambda"]][[2]],
                                                        nMesh = nMesh)
    mixFit = normalMixFit
    CRISPhieRmixFit = list(genes = unique(geneIds), 
                           locfdr = 1 - genePosteriors,
                           geneScores = genePosteriors,
                           mixFit = normalMixFit)
  }
  return(CRISPhieRmixFit)
}




