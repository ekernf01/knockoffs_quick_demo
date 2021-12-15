setwd("~/Desktop/jhu/research/projects/knockoffs/applications/package_demo/")
library("ggplot2")
library("magrittr")
results = data.frame(dimension = c(10, 20, 50, 100, 200),
                     time_rlookc =NA, time_knockoff = NA,
                     memory_rlookc = NA, memory_knockoff = NA)

for(i in seq(nrow(results))){
  X = rnorm(results[i, "dimension"]*1000) %>% matrix(ncol = results[i, "dimension"], nrow = 1000)
  sigma = cov(X)
  leave_one_out_naive = function(X){
    lapply(seq(ncol(X)), function(i) knockoff::create.gaussian(X[,-i], mu = 0, Sigma = sigma[-i, -i] ))
  }
  results[i, "time_rlookc"] = microbenchmark::microbenchmark(
    {knockoffs = rlookc::generateLooks(X, mu = 0, Sigma = sigma, output_type = "knockoffs_compact")},
    times = 1
  )$time %>% mean
  results[i, "time_knockoff"] = microbenchmark::microbenchmark(
    {knockoffs = leave_one_out_naive(X)},
    times = 1
  )$time %>% mean
  results[i, "memory_rlookc"] = peakRAM::peakRAM(
    {knockoffs = rlookc::generateLooks(X, mu = 0, Sigma = sigma, output_type = "knockoffs_compact")}
  )$Peak_RAM_Used_MiB
  results[i, "memory_knockoff"] = peakRAM::peakRAM(
    {knockoffs = leave_one_out_naive(X)}
  )$Peak_RAM_Used_MiB
}
write.csv(results, "leave_one_out.csv")
results = read.csv("leave_one_out.csv", row.names = 1)
results %>%
  tidyr::pivot_longer(cols = 2:5, names_to = "name", values_to = "value") %>%
  tidyr::separate("name", into = c("resource", "package")) %>%
  dplyr::mutate(resource = gsub("memory", "memory (MiB)", resource)) %>%
  dplyr::mutate(resource = gsub("time", "time (ms)", resource)) %>%
  ggplot() +
  geom_point(aes(y = value, x = dimension, color = package)) +
  facet_wrap(~resource, scales = "free_y") +
  scale_y_log10() +
  scale_x_log10() +
  ylab("Resource consumption") +
  ggtitle("Leave-one-out knockoff construction")
ggsave("leave_one_out.pdf", width = 4, height = 2)
ggsave("leave_one_out.svg", width = 4, height = 2)
