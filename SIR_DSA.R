# setwd("~/OneDrive - The University of Nottingham/Scripts/SIRapproxSEIR")

library(ggplot2)
library(cmdstanr)

# Remove all variables and clear the workspace
rm(list = ls())

# Confirm that all variables are removed
ls()


# SEIR Model Parameters
# alpha <- 0.25  # Incubation rate
# beta <- 0.8   # Infection rate
# gamma <- 0.4  # Recovery rate
# rho <- 0.01

# SEIR Model Parameters
alpha <- 2.25  # Incubation rate
beta <- 1.5   # Infection rate
gamma <- 0.5  # Recovery rate
rho <- 0.01


R0 <- beta/gamma
print(R0)

# Initial population values
n <- 10000     # Total population size
init_susceptible <- n 

init_exposed <- 0
init_infected <- ceiling(init_susceptible*rho)
init_removed <- 0

filename <- paste0("data/SEIR_sims_n_", n, "_alpha_", alpha, "_beta_", beta, "_gamma_", gamma, "_rho_", rho, ".csv")
seir_data <- read.csv(filename)
####### Extracting event times 
maxT <- 15
infection_times <- seir_data$Time[seir_data$Event == 1]
infection_times <- infection_times[infection_times <maxT]
hist(infection_times)


I2R_times <- seir_data$Time[seir_data$Event == 3]
I2R_times <- I2R_times[I2R_times < maxT]
hist(I2R_times)

######################
#####Stan parameters
N <- 2500
K <- 2500

### Take a random sample of infection and recovery times
if (length(infection_times) < N){
  infection_times_sample <- sort(unique(infection_times), decreasing = FALSE)
  N <- length(infection_times_sample)
} else {
  infection_times_sample <- sort(sample(unique(infection_times), size = N, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(infection_times_sample)
}


hist(infection_times_sample)


if (length(I2R_times) < K){
  I2R_times_sample <- sort(unique(I2R_times), decreasing = FALSE)
  K <- length(I2R_times_sample)
} else {
  I2R_times_sample <- sort(sample(unique(I2R_times), size = K, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  K <- length(I2R_times_sample)
}

hist(I2R_times_sample)

###################################################
############ Inference with only infection times
######### Create Stan data 
stan_data <- list(N = N, 
                  infection_times = infection_times_sample)

nChains <- 2
nIter <- 7500

#file <- file.path(getwd(),"infection_times_DSA.stan")
file <- file.path(getwd(),"SIR_DSA_infection_times.stan")
mod_1 <- cmdstan_model(file)

print(mod_1$print())

fit_1 <- mod_1$sample(
  data = stan_data, 
  seed = 123, 
  chains = nChains, 
  parallel_chains = nChains,
  iter_sampling = nIter, 
  adapt_delta = 0.9999,
  refresh = 500 # print update every 500 iters
)

fit_1$summary()

print(fit_1$summary())

fit_1$cmdstan_diagnose()


# fit_vb <- mod_1$variational(data = stan_data, seed = 123)
# print(fit_vb$summary())
# hist(fit_vb$draws("beta"))
# hist(fit_vb$draws("alpha"))
# hist(fit_vb$draws("gamma"))
# hist(fit_vb$draws("rho"))
# hist(fit_vb$draws("R0"))






posterior_samples_1 <- fit_1$draws(format = "df")

# write.csv(posterior_samples_1, file = "plots/SIR/SIR_DSA_posterior_samples.csv")

#posterior_samples_1 <- read.csv("plots/SIR/SIR_DSA_posterior_samples.csv")

R0_plot_1<-ggplot(posterior_samples_1, aes(x=R0)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() 
R0_plot_1
# ggsave("plots/SIR/SIR_R0_histogram_uninformative.pdf", plot=R0_plot_1, device="pdf", width = 6, height = 4)

alpha_plot_1<-ggplot(posterior_samples_1, aes(x=alpha)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Latent period", y = "Density")+
  geom_point(aes(x = mean(alpha), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
alpha_plot_1
# ggsave("plots/SIR/SIR_alpha_histogram_uninformative.pdf", plot=alpha_plot_1, device="pdf", width = 6, height = 4)

beta_plot_1<-ggplot(posterior_samples_1, aes(x=beta)) + 
  geom_density(aes(y=..density..), alpha = 0.5,  fill = "salmon") +
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(beta), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
beta_plot_1
# ggsave("plots/SIR/SIR_beta_histogram_uninformative.pdf", plot=beta_plot_1, device="pdf", width = 6, height = 4)


gamma_plot_1<-ggplot(posterior_samples_1, aes(x=gamma)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Infectious period", y = "Density")+
  geom_point(aes(x = mean(gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot_1
# ggsave("plots/SIR/SIR_gamma_histogram_uninformative.pdf", plot=gamma_plot_1, device="pdf", width = 6, height = 4)


rho_plot_1<-ggplot(posterior_samples_1, aes(x=rho)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot_1
# ggsave("plots/SIR/SIR_rho_histogram_uninformative.pdf", plot=rho_plot_1, device="pdf", width = 6, height = 4)



#########################################
############ Inference using both infection and recovery times
#########################################

#### Create Stan data
stan_data <- list(N = N, 
                  infection_times = infection_times_sample,
                  K = K, 
                  I2R_times = I2R_times_sample)

#stan_data <- list(N = N, infection_times = infection_times_sample,K = K, I2R_times = I2R_times_sample,alpha_0 = 0.2,beta_0 = 0.7,gamma_0 = 0.34,rho_0 = 0.009)

nChains <- 1
nIter <- 2500

#file <- file.path(getwd(),"infection_times_DSA.stan")
file <- file.path(getwd(),"SIR_full_DSA.stan")
mod <- cmdstan_model(file)

print(mod$print())

fit <- mod$sample(
  data = stan_data, 
  seed = 123, 
  chains = nChains, 
  parallel_chains = nChains,
  iter_sampling = nIter, 
  adapt_delta = 0.99,
  refresh = 500 # print update every 500 iters
)

fit$summary()

print(fit$summary())

posterior_samples <- fit$draws(format = "df")

# write.csv(posterior_samples, file = "plots/SIR/SIR_full_DSA_posterior_samples.csv")

#posterior_samples <- read.csv("plots/SIR/SIR_full_DSA_posterior_samples.csv")
# posterior_samples <- read.csv("plots/FMD/FMD_DSA_posterior_samples.csv")
# summary_table <- summary(posterior_samples)
# write.csv(summary_table, file = "plots/FMD/posterior_summary.csv")


R0_plot<-ggplot(posterior_samples, aes(x=R0)) + 
  geom_density(aes(x=R0, y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() 
R0_plot

ggsave("plots/SIR/SIR_full_R0_histogram_uninformative.pdf", plot=R0_plot, device="pdf", width = 6, height = 4)

alpha_plot<-ggplot(posterior_samples, aes(x=alpha)) + 
  geom_density(aes(y=..density..),  alpha = 0.5, fill = "salmon") +
  labs(x="Latent period", y = "Density")+
  geom_point(aes(x = mean(alpha), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
alpha_plot

ggsave("plots/SIR/SIR_full_alpha_histogram_uninformative.pdf", plot=alpha_plot, device="pdf", width = 6, height = 4)

beta_plot<-ggplot(posterior_samples, aes(x=beta)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(beta), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
beta_plot

ggsave("plots/SIR/SIR_full_beta_histogram_uninformative.pdf", plot=beta_plot, device="pdf", width = 6, height = 4)


gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_density(aes(y=..density..), alpha = 0.5,  fill = "salmon") +
  labs(x="Infectious period", y = "Density")+
  geom_point(aes(x = mean(gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot

ggsave("plots/SIR/SIR_full_gamma_histogram_uninformative.pdf", plot=gamma_plot, device="pdf", width = 6, height = 4)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_density(aes(y=..density..), alpha = 0.5, fill = "salmon") +
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(rho), y = 0), pch = 17, size = 10, col = 2) +
  geom_point(aes(x = rho, y = 0), pch = 13, size = 10, col = 2)+
  theme_classic()
rho_plot

ggsave("plots/SIR/SIR_full_rho_histogram_uninformative.pdf", plot=rho_plot, device="pdf", width = 6, height = 4)



