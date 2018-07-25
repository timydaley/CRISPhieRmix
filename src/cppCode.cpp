
#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// from smithlab_cpp
double logSumLogVec(NumericVector vals) {
  const NumericVector::iterator x = 
    std::max_element(vals.begin(), vals.begin() + vals.size());
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < vals.size(); ++i)
    if (i != max_idx)
      sum += exp(vals[i] - max_val);
  return max_val + log(sum);
}

double logSumWeightedAverage(double log_x, 
                             double log_y,
                             double q){
  NumericVector vals(2);
  vals(0) = log(q) + log_x;
  vals(1) = log(1 - q) + log_y;
  return logSumLogVec(vals);
}

// [[Rcpp::export]]
NumericVector integratedExpectation(NumericVector geneIds,
                                    NumericVector log_alt_guide_probs,
                                    NumericVector log_null_guide_probs,
                                    NumericVector q,
                                    NumericVector p,
                                    NumericVector weights){
  Rcerr << "length of geneIds = " << geneIds.size() << std::endl;
  Rcerr << "length of log_alt_guide_probs = " << log_alt_guide_probs.size() << std::endl;
  Rcerr << "length of log_null_guide_probs = " << log_null_guide_probs.size() << std::endl;
  Rcerr << "length of q = " << q.size() << std::endl;
  Rcerr << "length of p = " << p.size() << std::endl;
  Rcerr << "length of weights = " << weights.size() << std::endl;
  
  assert(q.size() == p.size());
  assert(q.size() == weights.size());
  int nGenes = max(geneIds);
  Rcerr << "nGenes = " << nGenes << std::endl;
  
  
  NumericVector genePosteriors(nGenes);
  for(size_t i = 0; i < q.size(); i++){
    NumericVector logPosGeneProbs(nGenes);
    NumericVector logNullGeneProbs(nGenes);
    for(size_t j = 0; j < logPosGeneProbs.size(); j++){
      logPosGeneProbs(j) = log(p(i));
      logNullGeneProbs(j) = log(1 - p(i));
    }
    for(size_t j = 0; j < geneIds.size(); j++){
      logPosGeneProbs(geneIds(j) - 1) += logSumWeightedAverage(log_alt_guide_probs(j),
                                                               log_null_guide_probs(j), q(i));
      logNullGeneProbs(geneIds(j) - 1) += log_null_guide_probs(j);
    }
    for(size_t j = 0; j < genePosteriors.size(); j++){
      NumericVector y(2);
      y(0) = logPosGeneProbs(j);
      y(1) = logNullGeneProbs(j);
      double logDenom = logSumLogVec(y);
      genePosteriors(j) += weights(i)*exp(logPosGeneProbs(j) - logDenom);
    }
  }

  return genePosteriors;
}

