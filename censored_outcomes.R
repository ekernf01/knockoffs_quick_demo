setwd("~/Desktop/jhu/research/projects/knockoffs/applications/package_demo/")
library("magrittr")
library("dplyr")
library("ggplot2")
ggplot2::theme_update(text = element_text(family = "ArialMT"))

# Suppose we want to judge calibration on a dataset where we have positive
# and negative examples of associated variables, but not all pairs are 
# observed and the overall balance may be off.
# Can we still obtain "honest" results?

# Simulation is based on DREAM5 e coli data
withr::with_dir(
  "~/Desktop/jhu/research/datalake/dream5/DREAM5_network_inference_challenge",
  {
    ecoli_expression = read.table("Network3/input data/net3_expression_data.tsv", header = T) 
    ecoli_tf = read.table("Network3/input data/net3_transcription_factors.tsv") 
    ecoli_network = read.table("Network3/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv")
  }
)
head(ecoli_expression)
dim(ecoli_expression)
ecoli_expression %<>% sweep(2, colMeans(ecoli_expression), FUN = "-")
ecoli_tf_expression = ecoli_expression[ecoli_tf[[1]]] %>% as.matrix
ecoli_tf_covariance = corpcor::cor.shrink(ecoli_tf_expression)
ecoli_tf_covariance = matrix(ecoli_tf_covariance, dimnames = dimnames(ecoli_tf_covariance), nrow = nrow(ecoli_tf_covariance))
# But, all data are simulated and Gaussian with 0 mean and known covariance.
set.seed(0)
ecoli_tf_covariance_sqrt = chol(ecoli_tf_covariance)

runKnockoffFilter = function(){
  X = matrix(rnorm(prod(dim(ecoli_tf_expression))),
             nrow = nrow(ecoli_tf_expression)) %*% 
    ecoli_tf_covariance_sqrt
  X_k = knockoff::create.gaussian(X, 0, ecoli_tf_covariance)
  # Simple outcome, active sets of size 1
  Y = X
  # Knockoff association statistics -- difference in OLS coeffs
  W = coef(lm(Y ~ X + 0)) -  coef(lm(Y ~ X_k + 0)) 
  W %<>% reshape2::melt(value.name = "knockoff_stat")
  colnames(W)[1:2] = c("feature", "target")
  W %<>% mutate(feature = gsub("X", "", feature))
  W %<>% mutate(is_confirmed = feature==target)
  W
}

#' Take a network and make it "unbalanced" by removing
#' different fractions of positive or negative examples.
#'
makeUnbalancedNetwork = function(network, prop_positives, prop_negatives, seed = 0){
  network %<>% subset(!is.na(is_confirmed))
  if(all(network$is_confirmed, na.rm = T)){
    stop("Cannot make it unbalanced with no negative examples.\n")
  }
  set.seed(seed)
  include_if_positive = rbinom(prob = prop_positives, size = 1, n = nrow(network))
  include_if_negative = rbinom(prob = prop_negatives, size = 1, n = nrow(network))
  network$include = ifelse( network$is_confirmed, include_if_positive, include_if_negative) %>% as.logical
  network %<>% subset(include)
  network
}

# Check the calibration: if q<0.3, is it about 30% wrong?
na2zero = function(X, filler = 0) {
  X[is.na(X)] = filler
  X
}
checkCalibration = function (W, targeted_fdrs = (1:10)/10) {
  fdr = rep(NA,10)
  for (j in seq_along(targeted_fdrs)) {
    fdr[[j]] = 
      W %>% 
      subset(q < targeted_fdrs[[j]]) %>% 
      extract2("is_confirmed") %>% 
      mean %>% 
      subtract(1, .) %>%
      na2zero
  }
  fdr
}

# Try it with different types of bias in the reference data.
{
  results = list()
  i = 0
  for(replicate in 1:10){
    for(bias in c("positive", "negative", "no bias")){
      i = i + 1
      W = runKnockoffFilter()
      # compute q on all hypotheses
      W$q = W$knockoff_stat %>% rlookc::knockoffQvals() 
      # Omit hypotheses with no availability in the gold standard
      W %<>% makeUnbalancedNetwork(
        prop_positives = ifelse(bias=="positive", 1, 0.2), 
        prop_negatives = ifelse(bias=="negative", 1, 0.2), 
        seed = replicate
      )
      results[[i]] = 
        data.frame(expected_fdr = (1:10)/10, 
                   observed_fdr = checkCalibration(W), 
                   bias = bias, 
                   method = "FDR on all\nhypotheses",
                   replicate = replicate)
      # recompute q only for hypotheses with available data
      i = i + 1
      W$q = W$knockoff_stat %>% rlookc::knockoffQvals() 
      results[[i]] = 
        data.frame(expected_fdr = (1:10)/10, 
                   observed_fdr = checkCalibration(W), 
                   bias = bias, 
                   method = "FDR on testable\nhypotheses",
                   replicate = replicate)
    }
  }
  results %<>% data.table::rbindlist()
}

ggplot(results) + 
  geom_boxplot(aes(x = expected_fdr, y = observed_fdr, group = expected_fdr)) + 
  facet_grid(method~bias) + 
  geom_abline(aes(slope = 1, intercept = 0)) + 
  ggtitle("Calibration with biased gold standards", "Bias towards which type of example:") + 
  scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
  scale_y_continuous(breaks = (0:2)/2, limits = 0:1) + 
  coord_fixed()
ggsave("censored_gold_standards.pdf", width = 6, height = 4)
ggsave("censored_gold_standards.svg", width = 6, height = 4)

