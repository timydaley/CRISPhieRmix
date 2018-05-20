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

skewtMixExpectationStep3comp <- function(x, skewtFit, 
                                         qpPos, qpNeg,
                                         muPos, muNeg, 
                                         sigma){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigma, log = TRUE)
  log_neg_prob = log(qpNeg) + dnorm(x, mean = muNeg, sd = sigma, log = TRUE)
  log_null_prob = log(1 - qpPos - qpNeg) + sn::dst(x, dp = skewtFit$dp, log = TRUE)
  
  log_denom = apply(cbind(log_pos_prob, log_neg_prob, log_null_prob), 1, logSumLogVec)
  return(list(posProbs = exp(log_pos_prob - log_denom),
              negProbs = exp(log_neg_prob - log_denom)))
}

skewtMixMaxStep2comp <- function(x, posProbs, skewtFit = NULL){
  pq = mean(posProbs)
  mu = mean(posProbs*x)/mean(posProbs)
  sigma = sqrt(mean(posProbs*(x - mu)^2)/mean(posProbs))
  if(is.null(skewtFit)){
    skewtFit = sn::st.mple(y = x, w = 1 - posProbs)
  }
  return(list(pq = pq, mu = mu, sigma = sigma, skewtFit = skewtFit))
}

skewtMixMaxStep3comp <- function(x, posProbs, negProbs, skewtFit = NULL){
  qpPos = mean(posProbs)
  qpNeg = mean(negProbs)
  muPos = mean(posProbs*x)/qpPos
  muNeg = mean(negProbs*x)/qpNeg
  sigma = sqrt(mean(posProbs*(x - muPos)^2 + negProbs*(x - muNeg)^2)/(qpPos + qpNeg))

  if(is.null(skewtFit)){
    cat("refitting skewtFit \n")
    skewtFit = sn::st.mple(y = x, w = 1 - posProbs - negProbs)
  }
  return(list(qpPos = qpPos, qpNeg = qpNeg,
              muPos = muPos, muNeg = muNeg,
              sigma = sigma,
              skewtFit = skewtFit))
}

skewt2compLogLike <- function(x, pq, skewtFit, mu, sigma){
  log_pos_prob = log(pq) + dnorm(x, mean = mu, sd = sigma, log = TRUE)
  log_null_prob = log(1 - pq) + sn::dst(x, dp = skewtFit$dp, log = TRUE)
  return(sum(apply(cbind(log_pos_prob, log_null_prob), 1, logSumLogVec)))
}

skewt3compLogLike <- function(x, skewtFit, 
                              qpPos, qpNeg,
                              muPos, muNeg, 
                              sigma){
  log_pos_prob = log(qpPos) + dnorm(x, mean = muPos, sd = sigma, log = TRUE)
  log_neg_prob = log(qpNeg) + dnorm(x, mean = muNeg, sd = sigma, log = TRUE)
  log_null_prob = log(1 - qpPos - qpNeg) + sn::dst(x, dp = skewtFit$dp, log = TRUE)
  
  return(sum(apply(cbind(log_pos_prob, log_neg_prob, log_null_prob), 1, logSumLogVec)))
}


skewtEM2comp <- function(x, skewtFit = NULL, max_iter = 1000, tol = 1e-10,
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
    
    loglike = skewt2compLogLike(x, pq, skewtFit, mu, sigma)
    iter = iter + 1
    if(VERBOSE){
      cat("iter: ", iter, "\n")
      cat("loglike = ", loglike, "\n")
      cat("prevloglike = ", prevloglike, "\n")
      cat("mu = ", mu, "\n")
      cat("sigma = ", sigma, "\n")
      cat("pq = ", pq, "\n")
    }
    if(abs(loglike - prevloglike)/n_obs < tol | iter > max_iter){
      if(VERBOSE){
        cat("stop after iteration ", iter, "\n")
      }
      break
    }
  }
  return(list(posProbs = posProbs, skewtFit = skewtFit,
              pq = pq, mu = mu, sigma = sigma))
}

skewtEM3comp <- function(x, skewtFit = NULL, 
                         max_iter = 1000, tol = 1e-10,
                         qpPos = 0.05, qpNeg = 0.05,
                         muPos = 5, muNeg = -5,
                         sigma = 1,
                         VERBOSE = FALSE){
  n_obs = length(x)
  providedFit = !is.null(skewtFit)
  prevloglike = -1e100;
  loglike = skewt3compLogLike(x, skewtFit, qpPos, qpNeg,
                              muPos, muNeg, sigma)
  iter = 0;
  initial_expectations = skewtMixExpectationStep3comp(x, skewtFit, qpPos, qpNeg,
                                                      muPos, muNeg, sigma)
  posProbs = initial_expectations$posProbs
  negProbs = initial_expectations$negProbs
  
  repeat{
    prevloglike = loglike
    # max step
    if(providedFit){
      updated_params = skewtMixMaxStep3comp(x, posProbs, negProbs, skewtFit)
    }
    else{
      updated_params = skewtMixMaxStep3comp(x, posProbs, negProbs, skewtFit = NULL)
    }
    qpPos = updated_params$qpPos
    qpNeg = updated_params$qpNeg
    muPos = updated_params$muPos
    muNeg = updated_params$muNeg
    sigma = updated_params$sigma
    skewtFit = updated_params$skewtFit
    # expectation step
    updated_expectations = skewtMixExpectationStep3comp(x, skewtFit, qpPos, qpNeg,
                                                        muPos, muNeg, sigma)
    posProbs = updated_expectations$posProbs
    negProbs = updated_expectations$negProbs
    
    loglike = skewt3compLogLike(x, skewtFit, qpPos, qpNeg,
                                muPos, muNeg, sigma)
    iter = iter + 1
    if(VERBOSE){
      cat("iter = ", iter, "\n",
          "prevloglike = ", prevloglike, "\n",
          "loglike = ", loglike, "\n",
          "qpPos = ", qpPos, "\n",
          "qpNeg = ", qpNeg, "\n", 
          "muPos = ", muPos, "\n",
          "muNeg = ", muNeg, "\n",
          "sigma = ", sigma, "\n")
    }
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
  return(list(posProbs = posProbs, 
              negProbs = negProbs,
              skewtFit = skewtFit,
              qpPos = qpPos, qpNeg = qpNeg,
              muPos = muPos, muNeg = muNeg, 
              sigma = sigma))
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


integratedGeneExpectationSkewtMix <- function(x, geneIds, skewtFit, mu = 5, sigma = 1, pq = 0.1, nMesh = 1000){
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationSkewtMix(x, geneIds, skewtFit = skewtFit, 
                                                        q = q.grid[i], p = pq/q.grid[i],
                                                        mu = mu, sigma = sigma),
                    simplify = TRUE)
  EZ_g.mat = t(EZ_g.mat)
  return(apply(EZ_g.mat, 2, mean))
}


setBimodalParams <- function(mu, pq){
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
              muNeg = muNeg))
}

geneExpectationSkewtMix <- function(x, geneIds, 
                                    p, q,
                                    log_skewt_probs,
                                    log_norm_probs){
  log_null_probs = sapply(unique(geneIds),
                          function(g) sum(log_skewt_probs[which(geneIds == g)]))
  log_pos_probs = sapply(unique(geneIds),
                         function(g) sum(sapply(which(geneIds == g),
                                                function(i) logSumLogVec(c(log(q) + log_norm_probs[i],
                                                                           log(1 - q) + log_skewt_probs[i])))))
  log_denom = sapply(1:length(unique(geneIds)),
                     function(i)
                       logSumLogVec(c(log(p) + log_pos_probs[i], log(1 - p) + log_null_probs[i])))
  return(exp(log(p) + log_pos_probs - log_denom))
}

gaussQuadGeneExpectationSkewtMix <- function(x, geneIds, skewtFit,
                                             mu = 5, sigma = 1, 
                                             lowerLim = 0.1, upperLim = 1, 
                                             nMesh = 100){
  log_skewt_probs = sn::dst(x, dp = skewtFit$dp, log = TRUE)
  log_norm_probs = dnorm(x, mu, sigma, log = TRUE)
  quad.points.weights = statmod::gauss.quad.prob(nMesh, dist = "uniform", l = lowerLim, u = upperLim)
  #q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  #q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:nMesh,
                    function(i) quad.points.weights$weights[i]*geneExpectationSkewtMix(x, geneIds, 
                                                                                       q = quad.points.weights$nodes[i],
                                                                                       p = lowerLim/quad.points.weights$nodes[i],
                                                                                       log_skewt_probs = log_skewt_probs,
                                                                                       log_norm_probs = log_norm_probs),
                    simplify = TRUE)
  EZ_g.mat = t(EZ_g.mat)
  return(apply(EZ_g.mat, 2, sum))
}


integratedGeneExpectationSkewt2groupsMix <- function(x, geneIds, skewtFit,
                                                  mu = 5, sigma = 1, pq = 0.1, nMesh = 100){
  log_skewt_probs = sn::dst(x, dp = skewtFit$dp, log = TRUE)
  log_norm_probs = dnorm(x, mu, sigma, log = TRUE)
  q.grid = seq(from = pq, to = 1, length = nMesh + 2)
  q.grid = q.grid[-c(1, length(q.grid))]
  EZ_g.mat = sapply(1:length(q.grid),
                    function(i) geneExpectationSkewtMix(x, geneIds, skewtFit = skewtFit,
                                                        q = q.grid[i], p = pq/q.grid[i],
                                                        log_skewt_probs = log_skewt_probs,
                                                        log_norm_probs = log_norm_probs),
                    simplify = TRUE)
  EZ_g.mat = t(EZ_g.mat)
  return(apply(EZ_g.mat, 2, mean))
}




CRISPhieRmixSkewt <- function(x, geneIds, negCtrl,
                              max_iter = 100, tol = 1e-10, pq = 0.1, mu = -4, sigma = 1,
                              nMesh = 100, breaks = 101, BIMODAL = FALSE,
                              VERBOSE = FALSE, PLOT = FALSE){
  #stopifnot(!is.null(neg.ctrl) | !is.null(empiricalNegCtrlFits))
  negCtrlFit = sn::st.mple(y = negCtrl)
  if(PLOT){
    s = seq(from = min(x), to = max(x), by = 0.1)
    hist(negCtrl, breaks = 80, probability = TRUE, xlim = c(min(x), max(x)), main = "negative control fit")
    lines(s, sn::dst(s, dp = negCtrlFit$dp), col = "red", lwd  = 2, lty = 2)
  }
  if(VERBOSE){
    cat("fit negative control distributions \n")
  }
  if(BIMODAL){
    if(VERBOSE){
      cat("3 groups \n")
    }
    params = setBimodalParams(mu, pq)
    skewtMix = skewtEM3comp(x, skewtFit = negCtrlFit, max_iter = max_iter, tol = tol,
                            qpPos = params$qpPos, qpNeg = params$qpNeg, 
                            muPos = params$muPos, muNeg = params$muNeg, sigma = sigma,
                            VERBOSE = VERBOSE)
    if(PLOT){
      b = seq(from = min(x) - 0.1, to = max(x) + 0.1, length = 81)
      hist(x, breaks = b, probability = TRUE, main = "mixture fit to observations")
      lines(b, skewtMix$qpPos*dnorm(b, skewtMix$muPos, skewtMix$sigma) + skewtMix$qpNeg*dnorm(b, skewtMix$muNeg, skewtMix$sigma), 
            lwd = 2, col = "darkgreen")
      lines(b, (1 - skewtMix$qpPos - skewtMix$qpNeg)*sn::dst(b, dp = negCtrlFit$dp), col = "red", lwd  = 2)
      lines(b, skewtMix$qpPos*dnorm(b, skewtMix$muPos, skewtMix$sigma) + skewtMix$qpNeg*dnorm(b, skewtMix$muNeg, skewtMix$sigma) +
              (1 - skewtMix$qpPos - skewtMix$qpNeg)*sn::dst(b, dp = negCtrlFit$dp), col = "darkviolet", lty = 2, lwd = 2)
    }
    posGenePosteriors = gaussQuadGeneExpectationSkewtMix(x, geneIds,
                                                         negCtrlFit,
                                                         mu = skewtMix$muPos,
                                                         sigma = skewtMix$sigma,
                                                         lowerLim = skewtMix$qpPos,
                                                         upperLim = 1 - skewtMix$qpNeg,
                                                         nMesh = nMesh)
    negGenePosteriors = gaussQuadGeneExpectationSkewtMix(x, geneIds,
                                                         negCtrlFit,
                                                         mu = skewtMix$muNeg,
                                                         sigma = skewtMix$sigma,
                                                         lowerLim = skewtMix$qpNeg,
                                                         upperLim = 1 - skewtMix$qpPos,
                                                         nMesh = nMesh)
    skewtMixFit = list(genes = unique(geneIds),
                       locfdr = sapply(1 - posGenePosteriors - negGenePosteriors, function(y) max(0, y)),
                       posGenePosteriors = posGenePosteriors,
                       negGenePosteriors = negGenePosteriors,
                       skewtMixFit = skewtMix)
      
  } else{
    stopifnot(length(mu) == 1)
    if(VERBOSE){
      cat("2 groups \n")
    }
    skewtMix = skewtEM2comp(x, skewtFit = negCtrlFit, max_iter = max_iter, tol = tol,
                            pq = pq, mu = mu, sigma = sigma, VERBOSE = VERBOSE)
    if(VERBOSE){
      cat("EM converged \n")
      cat("mu = ", skewtMix$mu, "\n")
      cat("sigma = ", skewtMix$sigma, "\n")
      cat("pq = ", skewtMix$pq, "\n")
    }
    if(PLOT){
      b = seq(from = min(x) - 0.1, to = max(x) + 0.1, length = 81)
      hist(x, breaks = b, probability = TRUE, main = "mixture fit to observations")
      lines(b, skewtMix$pq*dnorm(b, skewtMix$mu, skewtMix$sigma), 
            lwd = 2, col = "darkgreen")
      lines(b, (1 - skewtMix$pq)*sn::dst(b, dp = negCtrlFit$dp), col = "red", lwd  = 2)
      lines(b, skewtMix$pq*dnorm(b, skewtMix$mu, skewtMix$sigma) + (1 - skewtMix$pq)*sn::dst(b, dp = negCtrlFit$dp), col = "darkviolet", lty = 2, lwd = 2)
    }
    genePosteriors = gaussQuadGeneExpectationSkewtMix(x, geneIds, skewtFit = negCtrlFit,
                                                      mu = skewtMix$mu, sigma = skewtMix$sigma, 
                                                      lowerLim = skewtMix$pq, upperLim = 1, 
                                                      nMesh = nMesh)   
    skewtMixFit = list(genes = unique(geneIds),
                       locfdr = sapply(1 - genePosteriors, function(y) max(0, y)),
                       genePosteriors = genePosteriors,
                       skewtMixFit = skewtMix)
  }
  return(skewtMixFit)
}





