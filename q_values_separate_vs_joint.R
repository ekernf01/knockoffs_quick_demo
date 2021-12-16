setwd("~/Desktop/jhu/research/projects/knockoffs/applications/package_demo/")
library("magrittr")
library("dplyr")
library("ggplot2")
ggplot2::theme_update(text = element_text(family = "ArialMT"))

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
# But, all data are simulated and Gaussian with0 mean and known covariance.
set.seed(0)
ecoli_tf_covariance_sqrt = chol(ecoli_tf_covariance)
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
write.csv(W, "threshold_demo.csv")
# Demonstration of different pooling methods
ggplot(W %>% subset(target %in% paste0("G", 1:5))) + 
  geom_histogram(aes(x = knockoff_stat, fill = target)) + 
  ggtitle("Combining multiple knockoff runs")
ggsave("merge_illustration.png", width = 5, height = 5)
W$q_merged = W$knockoff_stat %>% rlookc::knockoffQvals()
W %<>% 
  dplyr::group_by(target) %>%
  dplyr::mutate(q_separate = rlookc::knockoffQvals(knockoff_stat)) %>%
  dplyr::ungroup()
W %<>% mutate(is_correct = feature==target)
ggplot(W) + 
  geom_point(aes(x = q_merged, y = q_separate, color = is_correct)) + 
  xlim(0:1) + ylim(0:1) + 
  ggtitle("Merging is a more powerful approach")
ggsave("merge_power.pdf", width = 5, height = 5)
W %<>% dplyr::arrange(q_merged)
W$empirical_fdr_merged   = not(W$is_correct) %>% cumsum %>% divide_by(1:nrow(W))
W %<>% dplyr::arrange(q_separate)
W$empirical_fdr_separate = not(W$is_correct) %>% cumsum %>% divide_by(1:nrow(W))
W_long = W %>% 
  tidyr::pivot_longer(cols = starts_with("q_"),
                      names_prefix = "q_",
                      names_to = "q_value_strategy",
                      values_to = "q_value") %>%  
  tidyr::pivot_longer(cols = starts_with("empirical_fdr_"), 
                      names_prefix = "empirical_fdr_",
                      names_to = "empirical_fdr_strategy",
                      values_to = "empirical_fdr") %>%
  subset(empirical_fdr_strategy==q_value_strategy) 
stratified_sample = function(x){
  to_return = c()
  for (bin in levels(cut(0:1, 10))){
    bin = bin %>% gsub("\\[|\\]|\\(|\\)", "", .) %>% strsplit(",") %>% extract2(1) %>% as.numeric
    lo = bin[[1]]
    hi = bin[[2]]
    stratum = which(lo < x & x < hi)
    to_return %<>% c(sample(stratum, min(100, length(stratum))))
  }
  to_return
}  
# Reduce the number of points or inkscape will go         v  e  r  y    s  l  o  w  l  y 
W_long = W_long[union( stratified_sample(W_long$q_value), stratified_sample(W_long$empirical_fdr) ), ]
ggplot(W_long) +  
  geom_point(aes(q_value, empirical_fdr, colour = q_value_strategy)) +  
  xlim(0:1) + ylim(0:1) + 
  geom_abline(aes(slope = 1, intercept = 0)) + 
  ggtitle("Calibration")
ggsave("merge_calibration.pdf", width = 5, height = 5)


