logSumLogVec <- function(logVec){
  log_max = max(logVec)
  logVec = logVec[-which.max(logVec)]
  return(log_max + log(1 + sum(exp(logVec - log_max))))
}

skewNormalMixExpectationStep2comp <- function(x, skewNormalFit, mu, sigma, pq){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  log_denom = apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)
  return(exp(log_pos_prob - log_denom))
}

skewNormalMixExpectationStep3comp <- function(x, skewNormalFit, 
                                              qpPos, qpNeg,
                                              muPos, muNeg,
                                              sigmaPos, sigmaNeg){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigmaPos, log = TRUE)
  log_neg_prob = log(qpNos) + dnorm(x, mean = muNeg, sd = sigmaNeg, log = TRUE)
  log_null_prob = log(1 - qpPos - qpNeg) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  log_denom = apply(cbind(log_pos_prob, log_neg_prob, log_null_prob), 1, logSumLogVec)
  return(list(posProbs = exp(log_pos_prob - log_denom),
              negProbs = exp(log_neg_prob - log_denom)))
}

skewNormalMixMaxStep2comp <- function(x, posProbs, skewNormalFit = NULL){
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = mean(posProbs*(x - mu)^2)/mean(posProbs)
  if(is.null(skewNormalFit)){
    skewNormalFit = sn::sn.mple(y = x, w = 1 - posProbs)
  }
  return(list(pq = pq, mu = mu, sigma = sigma, skewNormalFit = skewNormalFit))
}

skewNormalMixMaxStep3comp <- function(x, posProbs, negProbs, skewNormalFit = NULL){
  qpPos = mean(posProbs)
  qpNeg = mean(negProbs)
  muPos = mean(posProbs*x)/qpPos
  muPos = mean(negProbs*x)/qpNeg
  sigmaPos = mean(posProbs*(x - muPos)^2)/qpPos
  sigmaNeg = mean(negProbs*(x - muNeg)^2)/qpNeg
  if(is.null(skewNormalFit)){
    skewNormalFit = sn::sn.mple(y = x, w = 1 - posProbs)
  }
  return(list(qpPos = qpPos, qpNeg = qpNeg,
              muPos = muPos, muNeg = muNeg,
              sigmaPos = sigmaPos, sigmaNeg = sigmaNeg, 
              skewNormalFit = skewNormalFit))
}

skew2compNormalLogLike <- function(x, pq, skewNormalFit, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  return(sum(apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)))
}

skew3compNormalLogLike <- function(x, skewNormalFit, 
                                   qpPos, qpNeg,
                                   muPos, muNeg,
                                   sigmaPos, sigmaNeg){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigmaPos, log = TRUE)
  log_neg_prob = log(qpNos) + dnorm(x, mean = muNeg, sd = sigmaNeg, log = TRUE)
  log_null_prob = log(1 - qpPos - qpNeg) + sn::dsn(x, dp = skewNormalFit$dp, log = TRUE)
  
  return(sum(apply(cbind(log_pos_prob, log_neg_prob, log_null_prob), 1, logSumLogVec)))
}

skewNormalEM2comp <- function(x, skewNormalFit = NULL, max_iter = 1000, tol = 1e-10,
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
    
    loglike = skew2compNormalLogLike(x, pq, skewNormalFit, mu, sigma)
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
  return(list(posProbs = posProbs, skewNormalFit = skewNormalFit,
              pq = pq, mu = mu, sigma = sigma))
}

skewNormalEM3comp <- function(x, skewNormalFit = NULL, 
                              max_iter = 1000, tol = 1e-10,
                              qpPos = 0.05, qpNeg = 0.05,
                              muPos = 5, muNeg = -5, 
                              sigmaPos = 1, sigmaNeg = 1, 
                              VERBOSE = FALSE){
  n_obs = length(x)
  providedFit = !is.null(skewNormalFit)
  loglike = -1e100;
  iter = 0;
  posProbs = rep(qpPos, times = n_obs)
  negProbs = rep(qpNeg, times = n_obs)
  repeat{
    prevloglike = loglike
    # max step
    if(providedFit){
      updated_params = skewNormalMixMaxStep3comp(x, posProbs, negProbs, skewNormalFit = skewNormalFit)
    }
    else{
      updated_params = skewNormalMixMaxStep3comp(x, posProbs, negProbs, skewNormalFit = NULL)
    }
    qpPos = updated_params$qpPos 
    qpNeg = updated_params$qpNeg
    muPos = updated_params$muPos
    muNeg = updated_params$muNeg
    sigmaPos = updated_params$sigmaPos
    sigmaNeg = updated_params$sigmaNeg 
    skewNormalFit = updated_params$skewNormalFit
    # expectation step
    updated_expectations = skewNormalMixExpectationStep3comp(x, skewNormalFit, qpPos, qpNeg,
                                                             muPos, muNeg, sigmaPos, sigmaNeg)
    
    loglike = skew3compNormalLogLike(x, skewNormalFit, qpPos, qpNeg,
                                     muPos, muNeg, sigmaPos, sigmaNeg)
    iter = iter + 1
    if(VERBOSE & (iter %% 50 == 0)){
      cat("iter: ", iter, "\n")
    }
    if(abs(loglike - prevloglike)/n_obs < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(guide_posterior = posProbs, skewNormalFit = skewNormalFit,
              pq = pq, mu = mu, sigma = sigma))
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

integratedGeneExpectationSkewNormalMix <- function(x, geneIds, skewNormalFit, mu = 5,
                                                   sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationSkewNormalMix(x, geneIds, q.grid[i], pq/q.grid[i],
                                                             skewNormalFit, mu, sigma),
                    simplify = TRUE)
  
  return(apply(t(EZ_g.mat), 2, mean))
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

