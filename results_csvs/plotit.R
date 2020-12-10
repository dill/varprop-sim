# plot results from the simulations

library(ggplot2)

# load all data
results_files <- list.files(pattern=".csv$", full.names=TRUE)
results <- c()
for(fn in results_files){
  results <- rbind(results,
                   read.csv(fn, header=FALSE))
}
# give correct names
names(results) <- c("ind", "model", "Nbirds", "Nhat", "var", "n_obs", "betaset", "id")
# make numeric
results$Nhat <- as.numeric(results$Nhat)
results$var <- as.numeric(results$var)
results$n_obs <- as.numeric(results$n_obs)

# what quantile of the estimate does truth lie in?
get_N_quantile <- function(N, Nhat, cvN){

  meantrans <- function(m, cv) log(m)-0.5*log(cv^2+1)
  setrans <- function(cv) sqrt(log(cv^2+1))

  plnorm(N, meantrans(Nhat, cvN), setrans(cvN))
}

results$quantile <- get_N_quantile(results$Nbirds,
                                   results$Nhat,
                                   sqrt(results$var)/results$Nhat)

# do we have 200 results for each model/df combo?
table(results[,c("model","betaset")])

# plot n_obs
n_obs_plot <- ggplot(subset(results, model==" varprop")) +
  geom_histogram(aes(n_obs)) +
  facet_wrap(~betaset, scale="free_x") +
  labs(x="Number of observations", y="Count") +
  theme_minimal()


# plot results
ggplot(results) +
  geom_histogram(aes(quantile), binwidth=0.15) +
  geom_vline(xintercept=0.5, colour="red") +
  facet_grid(model~betaset) +
  labs(x="Quantile", y="Count") +
  theme_minimal()


results$CV <- sqrt(results$var)/results$Nhat

ggplot(results) +
  geom_histogram(aes(CV)) +
  facet_grid(model~betaset, scale="free_x") +
  labs(x="Quantile", y="Count") +
  theme_minimal()


ggplot(results) +
  geom_histogram(aes(Nbirds-Nhat),bins=20 ) +
  geom_vline(xintercept=0, colour="red") +
  facet_grid(model~betaset) +
  labs(x="Bias", y="Count") +
  theme_minimal()

