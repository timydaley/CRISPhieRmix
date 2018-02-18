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

LogSumWeightedAverage <- function(logprobs, q){
  stopifnot(dim(logprobs)[2] == 2)
  r = sum(vapply(1:dim(logprobs)[1], function(i) logSumLogVec(c(log(q) + logprobs[i, 1],
                                                                log(1 - q) + logprobs[i, 2])),
                 FUN.VALUE = double(1)))
  return(r)
}

geneExpectations2Group <-function(x, geneIds, q, p,
                                  log_alt_guide_probs,
                                  log_null_guide_probs){
  log_null_gene_probs = by(log_null_guide_probs, geneIds, sum)
  logprobs = data.frame(log_alt_guide_probs = log_alt_guide_probs,
                        log_null_guide_probs = log_null_guide_probs)
  log_pos_gene_probs = by(d, geneIds, function(y) LogSumWeightedAverage(y, q))
  log_denom = vapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_gene_probs[i], log(1 - p) + log_null_gene_probs[i])),
                     FUN.VALUE = double(1))
  return(exp(log(p) + log_pos_probs - log_denom))
}

geneExpectations3Groups <-function(x, geneIds, q, p,
                                  log_pos_guide_probs,
                                  log_neg_guide_probs,
                                  log_null_guide_probs){
  log_null_gene_probs = by(log_null_guide_probs, geneIds, sum)
  logprobs = data.frame(log_alt_guide_probs = log_alt_guide_probs,
                        log_null_guide_probs = log_null_guide_probs)
  log_pos_gene_probs = by(d, geneIds, function(y) LogSumWeightedAverage(y, q))
  log_denom = vapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_gene_probs[i], log(1 - p) + log_null_gene_probs[i])),
                     FUN.VALUE = double(1))
  return(exp(log(p) + log_pos_probs - log_denom))
}

gaussQuadGeneExpectation2Groups <- function(x, geneIds, 
                                            log_alt_guide_probs,
                                            log_null_guide_probs,
                                            lowerLim = 0.1, upperLim = 1, 
                                            nMesh = 100){

  quad.points.weights = statmod::gauss.quad.prob(nMesh, dist = "uniform", l = lowerLim, u = upperLim)
  nGenes = length(unique(geneIds))
  EZ_g.mat = vapply(1:nMesh,
                    function(i) quad.points.weights$weights[i]*geneExpectations2Group(x, geneIds,
                                                                                      q = quad.points.weights$nodes[i],
                                                                                      p = lowerLim/quad.points.weights$nodes[i],
                                                                                      log_alt_guide_probs = log_alt_guide_probs,
                                                                                      log_null_guide_probs = log_null_guide_probs),
                    FUN.VALUE = double(nGenes))
  EZ_g.mat = t(EZ_g.mat)
  return(apply(EZ_g.mat, 2, sum))
}

integratedGeneExpectation2Groups <- function(x, geneIds, 
                                             log_alt_guide_probs,
                                             log_null_guide_probs,
                                             lowerLim = 0.1, upperLim = 1, 
                                             nMesh = 100){
  q.grid = seq(from = lowerLim, to = upperLim, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  nGenes = length(unique(geneIds))
  EZ_g.mat = vapply(1:length(q.grid),
                    function(i) geneExpectations2Group(x, geneIds, q = q.grid[i], p = lowerLim/q.grid[i],
                                                       log_alt_guide_probs = log_alt_guide_probs,
                                                       log_null_guide_probs = log_null_guide_probs),
                    FUN.VALUE = double(nGenes))
  
  return(apply(EZ_g.mat, 1, mean))
}

geneExpectationNormalMix <- function(x, geneIds, q = 0.5, p = 0.1,
                                     mu0 = 0, sigma0 = 1, mu = 5, sigma = 1){
  log_null_guide_probs = dnorm(x, mu0, sigma0, log = TRUE)
  log_pos_guide_probs = dnorm(x, mu, sigma, log = TRUE)
  
  log_null_gene_probs = by(log_null_guide_probs, geneIds, sum)
  d = data.frame(log_pos_guide_probs = log_pos_guide_probs, 
                 log_null_guide_probs = log_null_guide_probs)
  log_pos_gene_probs = by(d, geneIds, function(y) LogSumWeightedAverage(y$log_pos_guide_probs, y$log_null_guide_probs, 2))
  #log_null_probs = vapply(unique(geneIds), function(g) sum(dnorm(x[which(geneIds == g)], mu0, sigma0, log = TRUE)))
  #log_pos_probs = vapply(unique(geneIds),
  #                      function(g) sum(vapply(x[which(geneIds == g)],
  #                                             function(y) logSumLogVec(c(log(q) + dnorm(y, mu, sigma, log = TRUE),
  #                                                                        log(1 - q) + dnorm(y, mu0, sigma0, log = TRUE))))))
  
  d = data.frame(log_pos_gene_probs = log_pos_gene_probs, 
                 log_null_gene_probs = log_null_gene_probs)
  log_denom = vapply(1:length(log_pos_gene_probs),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_gene_probs[i], log(1 - p) + log_null_gene_probs[i])), 
                     FUN.VALUE = double(1))
  return(exp(log(p) + log_pos_probs - log_denom))
}


integratedGeneExpectationNormalMix <- function(x, geneIds, mu0 = 0, sigma0 = 1,
                                              mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  log_null_guide_probs = dnorm(x, mu0, sigma0, log = TRUE)
  log_alt_guide_probs = dnorm(x, mu, sigma, log = TRUE)
  
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectations2Group(x, geneIds, q = q.grid[i], p = pq/q.grid[i],
                                                       log_alt_guide_probs = log_alt_guide_probs,
                                                       log_null_guide_probs = log_null_guide_probs),
                    simplify = TRUE)

  return(apply(EZ_g.mat, 1, mean))
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

empiricalExpectationStep <- function(x, null_coefficients,
                                     null_log_norm_factor, mu, sigma, pq){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_neg_prob = log(1 - pq) +
                 apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                 null_coefficients[1] - null_log_norm_factor
  log_denom = apply(cbind(log_pos_prob, log_neg_prob), 1, logSumLogVec)
  return(exp(log_pos_prob - log_denom))
}

empiricalMaxStep <- function(x, posProbs){
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = sqrt(mean(posProbs*(x - mu)^2)/mean(posProbs))
  return(list(pq = pq, mu = mu, sigma = sigma))
}


empirical2GroupsMixtureLogLike <- function(x, pq, null_coefficients,
                                     null_log_norm_factor, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) +
                 apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
                 null_coefficients[1] - null_log_norm_factor
  return(sum(apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)))
}


empirical2GroupEMmix <- function(x, null_coefficients, null_log_norm_factor,
                                max_iter = 1000, tol = 1e-10, pq = 0.1,
                                mu = 4, sigma = 1, VERBOSE = FALSE){
  n_obs = length(x)
  loglike = empirical2GroupsMixtureLogLike(x, pq, null_coefficients, null_log_norm_factor, mu, sigma)
  prevloglike = -1e100;
  iter = 0;
  repeat{
    prevloglike = loglike
    # expectation step
    posProbs = empiricalExpectationStep(x, null_coefficients,
                                        null_log_norm_factor, mu, sigma, pq)
    # max step
    updated_params = empiricalMaxStep(x, posProbs)
    pq = updated_params$pq
    mu = updated_params$mu
    sigma = updated_params$sigma
    loglike = empirical2GroupsMixtureLogLike(x, pq, null_coefficients, null_log_norm_factor, mu, sigma)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if(abs(loglike - prevloglike) < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(posProbs = posProbs,
              pq = pq, mu = mu, sigma = sigma,
              null_coefficients = null_coefficients,
              null_log_norm_factor = null_log_norm_factor,
              loglike = loglike))
}

empirical3GroupsExpectationStep <- function(x, null_coefficients,
                                            null_log_norm_factor, 
                                            qpPos, qpNeg, 
                                            muPos, muNeg,
                                            sigmaPos, sigmaNeg){
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
  muPos = mean(posProbs*x)/qpPos
  sigmaPos = sqrt(mean(posProbs*(x - muPos)^2)/qpPos)
  muNeg = mean(negProbs*x)/qpNeg
  sigmaNeg = sqrt(mean(negProbs*(x - muNeg)^2)/qpNeg)
  return(list(qpPos = qpPos, qpNeg = qpNeg, 
              muPos = muPos, muNeg = muNeg,
              sigmaPos = sigmaPos, sigmaNeg = sigmaNeg))
}


empirical3GroupsMixtureLogLike <- function(x, null_coefficients, null_log_norm_factor, 
                                            qpPos, qpNeg, 
                                            muPos, muNeg,
                                            sigmaPos, sigmaNeg){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigmaPos, log = TRUE)
  log_neg_prob = log(qpNeg) + dnorm(x, mean = muNeg, sd = sigmaNeg, log = TRUE)
  
  log_null_prob = log(1 - qpPos - qpNeg) +
    apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
    null_coefficients[1] - null_log_norm_factor

  return(sum(log(exp(log_pos_prob) + exp(log_neg_prob) + exp(log_null_prob))))
}


empirical3GroupEMmix <- function(x, null_coefficients, null_log_norm_factor,
                                max_iter = 1000, tol = 1e-10, qpPos = 0.05, qpNeg = 0.05,
                                muPos = 5, muNeg = -5, sigmaPos = 1, sigmaNeg = 1,
                                VERBOSE = FALSE){
  n_obs = length(x)
  loglike = empirical3GroupsMixtureLogLike(x, null_coefficients, null_log_norm_factor, 
                                            qpPos, qpNeg, 
                                            muPos, muNeg,
                                            sigmaPos, sigmaNeg)
  prevloglike = -1e100;
  iter = 0;
  repeat{
    prevloglike = loglike
    # expectation step
    guideExpectations = empirical3GroupsExpectationStep(x, null_coefficients, null_log_norm_factor, 
                                                        qpPos, qpNeg, 
                                                        muPos, muNeg,
                                                        sigmaPos, sigmaNeg)
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

    loglike = empirical3GroupsMixtureLogLike(x, null_coefficients, null_log_norm_factor, 
                                              qpPos, qpNeg, muPos, muNeg,
                                              sigmaPos, sigmaNeg)
    if(VERBOSE){
      cat("iter = ", iter, "\n")
      cat("loglike = ", loglike, "\n")
      cat("muPos = ", muPos, "\n",
          "muNeg = ", muNeg, "\n",
          "sigmaPos = ", sigmaPos, "\n",
          "sigmaNeg = ", sigmaNeg, "\n",
          "qpPos = ", qpPos, "\n", 
          "qpNeg = ", qpNeg, "\n")
    }
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if(abs(loglike - prevloglike) < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(posProbs = posProbs,
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

geneExpectationEmpiricalMix <- function(x, geneIds, 
                                        q, p,
                                        log_norm_probs,
                                        log_null_guide_probs){
  log_null_probs = by(log_null_guide_probs, geneIds, sum)
  d = data.frame()
  log_pos_probs = 
  log_pos_probs = vapply(unique(geneIds),
                         function(g) sum(vapply(which(geneIds == g),
                                                function(i) logSumLogVec(c(log(q) + log_norm_probs[i],
                                                                           log(1 - q) + log_null_guide_probs[i])),
                                                FUN.VALUE = double(1))),
                         FUN.VALUE = double(1))
  log_denom = vapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])), 
                     FUN.VALUE = double(1))
  return(exp(log(p) + log_pos_probs - log_denom))
}

gaussQuadGeneExpectationEmpiricalMix <- function(x, geneIds, null_coefficients, null_log_norm_factor,
                                                 mu = 5, sigma = 1, lowerLim = 0.1, upperLim = 1, 
                                                 nMesh = 100){
  log_null_guide_probs = apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
    null_coefficients[1] - null_log_norm_factor
  log_norm_probs = dnorm(x, mu, sigma, log = TRUE)
  quad.points.weights = statmod::gauss.quad.prob(nMesh, dist = "uniform", l = lowerLim, u = upperLim)
  #q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  #q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = vapply(1:nMesh,
                    function(i) quad.points.weights$weights[i]*geneExpectationEmpiricalMix(x, geneIds, q = quad.points.weights$nodes[i],
                                                                                           p = lowerLim/quad.points.weights$nodes[i],
                                                                                           log_norm_probs = log_norm_probs,
                                                                                           log_null_guide_probs = log_null_guide_probs),
                    FUN.VALUE = double(length(unique(geneIds))))

  return(apply(EZ_g.mat, 1, sum))
}


integratedGeneExpectationEmpiricalMix <- function(x, geneIds, null_coefficients, null_log_norm_factor,
                                                  mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  log_null_guide_probs = apply(t(null_coefficients[-1]*t(poly(x, degree = length(null_coefficients) - 1, raw = TRUE))), 1, sum) +
    null_coefficients[1] - null_log_norm_factor
  log_norm_probs = dnorm(x, mu, sigma, log = TRUE)
  nGenes = length(unique(geneIds))
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = vapply(1:length(q.grid),
                    function(i) geneExpectationEmpiricalMix(x, geneIds, q = q.grid[i], p = pq/q.grid[i],
                                                            log_norm_probs = log_norm_probs,
                                                            log_null_guide_probs = log_null_guide_probs),
                    FUN.VALUE = double(nGenes))

  return(apply(EZ_g.mat, 1, mean))
}

fitNegCtrl <- function(neg.ctrl, maxDegree = 20, minDegree = 4,
                       breaks = NULL, VERBOSE = FALSE){
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
  mixFit = empirical2GroupEMmix(x, 
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
    lines(b, mixFit$pq*dnorm(b, mixFit$mu, mixFit$sigma) + 
            (1 - mixFit$pq)*exp(apply(t(mixFit$null_coefficients[-1]*t(poly(b, degree = length(mixFit$null_coefficients) - 1, raw = TRUE))), 1, sum) + mixFit$null_coefficients[1] - mixFit$null_log_norm_factor), col = "darkviolet", lty = 2, lwd = 2)
  }
  
  genePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds,
                                                        negCtrlFit[["coefficients"]],
                                                        negCtrlFit[["log_norm_factor"]],
                                                        mu = mixFit[["mu"]],
                                                        sigma = mixFit[["sigma"]],
                                                        lowerLim = mixFit[["pq"]],
                                                        nMesh = nMesh)
  
  return(list(genes = unique(geneIds),
              locfdr = 1 - genePosteriors,
              genePosteriors = genePosteriors,
              mixFit = mixFit))
}

CRISPhieRmix3Groups <- function(x, geneIds, negCtrlFit,
                                max_iter = 100, tol = 1e-10, pq = c(0.05, 0.05), 
                                mu = c(-4, 4), sigma = c(1, 1),
                                nMesh = 100, breaks = 101,
                                VERBOSE = FALSE, PLOT = FALSE){
  params = setBimodalParams(mu, sigma, pq)  
  if(VERBOSE){
    cat("muPos = ", params$muPos, "\n",
        "muNeg = ", params$muNeg, "\n",
        "sigmaPos = ", params$sigmaPos, "\n",
        "sigmaNeg = ", params$sigmaNeg, "\n",
        "qpPos = ", params$qpPos, "\n", 
        "qpNeg = ", params$qpNeg, "\n")
  }
  mixFit = empirical3GroupEMmix(x, null_coefficients = negCtrlFit[["coefficients"]], 
                               null_log_norm_factor = negCtrlFit[["log_norm_factor"]],
                               max_iter = max_iter, tol = tol, 
                               qpPos = params$qpPos, qpNeg = params$qpNeg,
                               muPos = params$muPos, muNeg = params$muNeg, 
                               sigmaPos = params$sigmaPos, sigmaNeg = params$sigmaNeg,
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
  log_null_guide_probs = apply(t(negCtrlFit[["coefficients"]][-1]*t(poly(x, degree = length(negCtrlFit[["coefficients"]]) - 1, raw = TRUE))), 1, sum) +
    negCtrlFit[["coefficients"]][1] - negCtrlFit[["log_norm_factor"]]
  log_pos_guide_probs = dnorm(x, mean = mixFit[["muPos"]], sd = mixFit[["sigmaPos"]], log = TRUE)
  log_neg_guide_probs = dnorm(x, mean = mixFit[["muPos"]], sd = mixFit[["sigmaPos"]], log = TRUE)
  
  posGenePosteriors = gaussQuadGeneExpectation2Groups(x, geneIds = geneIds,
                                                      log_alt_guide_probs = log_pos_guide_probs,
                                                      log_null_guide_probs = log_null_guide_probs,
                                                      lowerLim = mixFit[["qpPos"]],
                                                      upperLim = 1 - mixFit[["qpNeg"]], 
                                                      nMesh = 100)
  
#  posGenePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds,
#                                                           negCtrlFit[["coefficients"]],
#                                                           negCtrlFit[["log_norm_factor"]],
#                                                           mu = mixFit[["muPos"]],
#                                                           sigma = mixFit[["sigmaPos"]],
#                                                           lowerLim = mixFit[["qpPos"]],
#                                                           upperLim = 1 - mixFit[["qpNeg"]],
#                                                           nMesh = nMesh)
  negGenePosteriors = gaussQuadGeneExpectation2Groups(x, geneIds = geneIds,
                                                      log_alt_guide_probs = log_neg_guide_probs,
                                                      log_null_guide_probs = log_null_guide_probs,
                                                      lowerLim = mixFit[["qpPos"]],
                                                      upperLim = 1 - mixFit[["qpNeg"]], 
                                                      nMesh = 100)
#    negGenePosteriors = gaussQuadGeneExpectationEmpiricalMix(x, geneIds, negCtrlFit[["coefficients"]], negCtrlFit[["log_norm_factor"]],
#                                                             mu = mixFit[["muNeg"]], sigma = mixFit[["sigmaNeg"]],
#                                                             lowerLim = mixFit[["qpNeg"]], upperLim = 1 - mixFit[["qpPos"]],
#                                                             nMesh = nMesh)
  
  return(list(genes = unique(geneIds),
              locfdr = 1 - posGenePosteriors - negGenePosteriors,
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
                         VERBOSE = FALSE, PLOT = FALSE){
  #stopifnot(!is.null(neg.ctrl) | !is.null(empiricalNegCtrlFits))
  if(!is.null(negCtrl)){
    negCtrlFit = fitNegCtrl(as.numeric(negCtrl), maxDegree = maxDegree, minDegree = minDegree,
                            breaks = 101, VERBOSE = VERBOSE)
    if(PLOT){
      s = seq(from = min(x), to = max(x), by = 0.1)
      hist(negCtrl, breaks = 80, probability = TRUE, xlim = c(min(x), max(x)), main = "negative control fit")
      lines(s, exp(apply(t(negCtrlFit$coefficients[-1]*t(poly(s, degree = length(negCtrlFit$coefficients) - 1, raw = TRUE))), 1, sum)
                   + negCtrlFit$coefficients[1] - negCtrlFit$log_norm_factor), col = "red", lwd  = 2, lty = 2)
    }
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
                                            breaks = 101, VERBOSE = VERBOSE, PLOT = PLOT)
      
    } else{
      stopifnot(length(mu) == 1)
      if(VERBOSE){
        cat("2 groups \n")
      }
      CRISPhieRmixFit = CRISPhieRmix2Groups(x, geneIds = geneIds, negCtrlFit = negCtrlFit, 
                                            max_iter = max_iter, tol = tol, pq = pq, 
                                            mu = mu, sigma = sigma, nMesh = nMesh, breaks = breaks,
                                            VERBOSE = VERBOSE, PLOT = PLOT)
      
    }
  } else{
    if(VERBOSE){
      cat("no negative controls provided, fitting hierarchical normal model \n")
    }
    require(mixtools)
    normalMixFit = mixtools::normalmixEM(x, k = 2, mu = c(0, mu),
                                         sigma = c(1, sigma))
    if(PLOT){
      plot(normalMixFit, density = TRUE)[2]
    }
    log_alt_guide_probs = dnorm(x, mean = normalMixFit[["mu"]][[2]], sd = normalMixFit[["sigma"]][[2]], log = TRUE)
    log_null_guide_probs = dnorm(x, mean = normalMixFit[["mu"]][[1]], sd = normalMixFit[["sigma"]][[1]], log = TRUE)
    
    genePosteriors = gaussQuadGeneExpectation2Groups(x, geneIds = geneIds,
                                                      log_alt_guide_probs = log_alt_guide_probs,
                                                      log_null_guide_probs = log_null_guide_probs,
                                                      lowerLim = normalMixFit[["lambda"]][[2]], 
                                                      upperLim = 1, nMesh = 100)
#    genePosteriors = integratedGeneExpectationNormalMix(x, geneIds = geneIds, mu0 = normalMixFit[["mu"]][[1]],
#                                                        mu = normalMixFit[["mu"]][[2]], sigma0 = normalMixFit[["sigma"]][[1]],
#                                                        sigma = normalMixFit[["sigma"]][[2]],  pq = normalMixFit[["lambda"]][[2]], nMesh = nMesh)
    CRISPhieRmixFit = list(genes = unique(geneIds), 
                           locfdr = 1 - genePosteriors,
                           geneScores = genePosteriors,
                           mixFit = normalMixFit)
  }
  CRISPhieRmixFit$locfdr = vapply(CRISPhieRmixFit$locfdr, function(x) max(x, 0),
                                  FUN.VALUE = double(1))
  return(CRISPhieRmixFit)
}




