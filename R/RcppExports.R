# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

integratedExpectation2groups <- function(geneIds, log_alt_guide_probs, log_null_guide_probs, q, tau, weights) {
    .Call('_CRISPhieRmix_integratedExpectation2groups', PACKAGE = 'CRISPhieRmix', geneIds, log_alt_guide_probs, log_null_guide_probs, q, tau, weights)
}

integratedExpectation3groups <- function(geneIds, log_pos_guide_probs, log_neg_guide_probs, log_null_guide_probs, q, tau_pos, tau_neg, weights) {
    .Call('_CRISPhieRmix_integratedExpectation3groups', PACKAGE = 'CRISPhieRmix', geneIds, log_pos_guide_probs, log_neg_guide_probs, log_null_guide_probs, q, tau_pos, tau_neg, weights)
}

