
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
NumericVector integratedExpectation2groups(NumericVector geneIds,
                                           NumericVector log_alt_guide_probs,
                                           NumericVector log_null_guide_probs,
                                           NumericVector q,
                                           double tau,
                                           NumericVector weights){
// make sure the size of the weights agree weight the size of the mesh
  assert(q.size() == weights.size());
  int nGenes = max(geneIds);  
  
  NumericVector marginalized_alt_probs(nGenes);
  NumericVector marginalized_null_probs(nGenes);
  for(size_t i = 0; i < q.size(); i++){
    NumericVector logAltGeneProbs(nGenes);
    NumericVector logNullGeneProbs(nGenes);
    for(size_t j = 0; j < geneIds.size(); j++){
      logAltGeneProbs(geneIds(j) - 1) += log(q(i)*exp(log_alt_guide_probs(j))
                                             + (1 - q(i))*exp(log_null_guide_probs(j)));
      logNullGeneProbs(geneIds(j) - 1) += log_null_guide_probs(j);
    }
    for(size_t g = 0; g < nGenes; g++){
      marginalized_alt_probs(g) += exp(log(weights(i)) + log(tau) - log(q(i))
                                       + logAltGeneProbs(g));
      marginalized_null_probs(g) += exp(log(weights(i)) + log(1 - tau/q(i))
                                        + logNullGeneProbs(g));
    }
  }
  NumericVector genePosteriors(nGenes);
  for(size_t g = 0; g < nGenes; g++){
    genePosteriors(g) = marginalized_alt_probs(g)/(marginalized_alt_probs(g) + marginalized_null_probs(g));
  }
  
  return genePosteriors;

}

// [[Rcpp::export]]
NumericVector integratedExpectation3groups(NumericVector geneIds,
                                           NumericVector log_pos_guide_probs,
                                           NumericVector log_neg_guide_probs,
                                           NumericVector log_null_guide_probs,
                                           NumericVector q,
                                           double tau_pos,
                                           double tau_neg,
                                           NumericVector weights){
  
  assert(q.size() == weights.size());
  int nGenes = max(geneIds);
  
  NumericVector marginalized_neg_probs(nGenes);
  NumericVector marginalized_pos_probs(nGenes);
  NumericVector marginalized_null_probs(nGenes);
  for(size_t i = 0; i < q.size(); i++){
    NumericVector logPosGeneProbs(nGenes);
    NumericVector logNegGeneProbs(nGenes);
    NumericVector logNullGeneProbs(nGenes);
    for(size_t j = 0; j < geneIds.size(); j++){
      logPosGeneProbs(geneIds(j) - 1) += log(q(i)*exp(log_pos_guide_probs(j))
                                             + (1 - q(i))*exp(log_null_guide_probs(j)));
      logNegGeneProbs(geneIds(j) - 1) += log(q(i)*exp(log_neg_guide_probs(j))
                                             + (1 - q(i))*exp(log_null_guide_probs(j)));
      logNullGeneProbs(geneIds(j) - 1) += log_null_guide_probs(j);
    }
    for(size_t g = 0; g < nGenes; g++){
      marginalized_pos_probs(g) += exp(log(weights(i)) + log(tau_pos) - log(q(i))
                                       + logPosGeneProbs(g));
      marginalized_neg_probs(g) += exp(log(weights(i)) + log(tau_neg) - log(q(i))
                                       + logNegGeneProbs(g));
      marginalized_null_probs(g) += exp(log(weights(i)) + log(1 - (tau_pos + tau_neg)/q(i))
                                        + logNullGeneProbs(g));
    }
  }
  NumericVector genePosteriors(nGenes);
  for(size_t g = 0; g < nGenes; g++){
    genePosteriors(g) = marginalized_pos_probs(g)/(marginalized_pos_probs(g) + marginalized_neg_probs(g) + marginalized_null_probs(g));
  }
  
  return genePosteriors;
}

