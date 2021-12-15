setwd("~/Desktop/jhu/research/projects/knockoffs/applications/q_values")

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
ecoli_expression %<>% sweep(2, colMeans(ecoli_tf_expression), FUN = "-")
ecoli_tf_expression = ecoli_expression[ecoli_tf[[1]]] %>% as.matrix
ecoli_tf_covariance = cov(ecoli_tf_expression)

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
# Demons of different pooling methods
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
ggsave("merge_power.png", width = 5, height = 5)
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
ggplot(W_long) +  
  geom_point(aes(q_value, empirical_fdr, colour = q_value_strategy)) +  
  xlim(0:1) + ylim(0:1) + 
  geom_abline(aes(slope = 1, intercept = 0)) + 
  ggtitle("Merging is better calibrated")
ggsave("merge_calibration.png", width = 5, height = 5)

W_long %>% 
  group_by(q_value_strategy, target) %>%
  summarize(
    false = sum(q_value<0.5 & !is_correct), 
    true = sum(q_value<0.5 & is_correct), 
  ) %>% 
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols = c("false", "true"),
                      names_to = "type", 
                      values_to = "discoveries") %>%
  ggplot() + 
  geom_histogram(aes(x = discoveries, 
                     fill = type)) + 
  facet_grid(~q_value_strategy) + 
  xlab("Number of discoveries") + 
  ylab("Number of targets")
ggsave("discovery_count_by_method.png", width = 5, height = 5)

