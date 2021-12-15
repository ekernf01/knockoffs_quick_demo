setwd("~/Desktop/jhu/research/projects/knockoffs/applications/package_demo/")
library("ggplot2")
library("magrittr")
results = data.frame(dimension = c(2, 5, 10, 20, 50, 100)*10,
                     time_rlookc =NA, time_knockoff = NA,
                     memory_rlookc = NA, memory_knockoff = NA)
for(i in seq(nrow(results))){
  X = rnorm(results[i, "dimension"]*10) %>% matrix(ncol = results[i, "dimension"], nrow = 10)
  results[i, "time_rlookc"] = microbenchmark::microbenchmark(
    {knockoffs = rlookc::createHighDimensionalKnockoffs(X, output_type = "knockoffs", lambda = 0.1)},
    times = 10
  )$time %>% mean
  results[i, "time_knockoff"] = microbenchmark::microbenchmark(
    {knockoffs = knockoff::create.second_order(X)},
    times = 10
  )$time %>% mean
  results[i, "memory_rlookc"] = peakRAM::peakRAM(
    {knockoffs = rlookc::createHighDimensionalKnockoffs(X, output_type = "knockoffs", lambda = 0.1)}
  )$Peak_RAM_Used_MiB
  results[i, "memory_knockoff"] = peakRAM::peakRAM(
    {knockoffs = knockoff::create.second_order(X)}
  )$Peak_RAM_Used_MiB
}
write.csv(results, "high_dimensional.csv")
results = read.csv("high_dimensional.csv", row.names = 1)
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
  ggtitle("High dimensional knockoff construction")
ggsave("high_dimensional.pdf", width = 4, height = 2)
ggsave("high_dimensional.svg", width = 4, height = 2)
