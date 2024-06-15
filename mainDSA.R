# setwd("~/OneDrive - The University of Nottingham/Scripts/SIRapproxSEIR")

library(ggplot2)
library(cmdstanr)
library(deSolve)
library(survival)
# Remove all variables and clear the workspace
rm(list = ls())

# Confirm that all variables are removed
ls()

source("helper_functions.R")


# SEIR Model Parameters
# alpha <- 0.25  # Incubation rate
# beta <- 0.8   # Infection rate
# gamma <- 0.4  # Recovery rate
# rho <- 0.01

# SEIR Model Parameters
alpha <- 2  # Incubation rate
beta <- 0.8   # Infection rate
gamma <- 0.4  # Recovery rate
rho <- 0.01

R0 <- beta/gamma
print(R0)


# Set the maximum simulation time
max_time <- 50

# Initial population values
n <- 1000    # Total population size
init_susceptible <- n 

init_exposed <- 0
init_infected <- ceiling(init_susceptible*rho)
init_removed <- 0


comparison_plots(alpha = alpha, beta = beta, gamma = gamma, rho = rho,
                 n = n, init_exposed = init_exposed, init_infected = init_infected, init_removed = init_removed, max_time = max_time)

# Initial population values
n <- 5000    # Total population size
init_susceptible <- n 

init_exposed <- 0
init_infected <- ceiling(init_susceptible*rho)
init_removed <- 0


comparison_plots(alpha = alpha, beta = beta, gamma = gamma, rho = rho,
                 n = n, init_exposed = init_exposed, init_infected = init_infected, init_removed = init_removed, max_time = max_time)


# Initial population values
n <- 10000    # Total population size
init_susceptible <- n 

init_exposed <- 0
init_infected <- ceiling(init_susceptible*rho)
init_removed <- 0


comparison_plots(alpha = alpha, beta = beta, gamma = gamma, rho = rho,
                 n = n, init_exposed = init_exposed, init_infected = init_infected, init_removed = init_removed, max_time = max_time)



# Initial population values
n <- 50000    # Total population size
init_susceptible <- n 

init_exposed <- 0
init_infected <- ceiling(init_susceptible*rho)
init_removed <- 0


comparison_plots(alpha = alpha, beta = beta, gamma = gamma, rho = rho,
                 n = n, init_exposed = init_exposed, init_infected = init_infected, init_removed = init_removed, max_time = max_time)




max_time <- 50
# Simulate the SEIR model using Gillespie's algorithm
seir_data <- simulate_SEIR_gillespie(alpha, beta, gamma, n, init_exposed, init_infected, init_removed, max_time)

print(seir_data)

# Simulate the SIR model using Gillespie's algorithm
#sir_data <- simulate_SIR_gillespie(beta, gamma, n, init_infected, init_removed, max_time)

#print(sir_data)

# Simulate the SVR model using Gillespie's algorithm
svr_data <- simulate_SVR_gillespie(alpha, beta, gamma, rho, n, init_infected, init_removed, max_time)

print(svr_data)


filename <- paste0("data/SEIR_sims_n_", n, "_alpha_", alpha, "_beta_", beta, "_gamma_", gamma, "_rho_", rho, ".csv")
write.csv(seir_data, file = filename)

####################################################
############ Comparison of the models 
####################################################
####### Extracting event times 

seir_infection_times <- seir_data$Time[seir_data$Event == 1]
hist(seir_infection_times)

svr_infection_times <- svr_data$Time[svr_data$Event == 1]
hist(svr_infection_times)


infection_times_hist <- ggplot() +
  geom_density(data = data.frame(x = seir_infection_times), aes(x = x, y = ..density..), fill = "steelblue", alpha = 0.5) +
  geom_density(color="steelblue",alpha=.8, lwd = 1.2, adjust = 1.75)+
  geom_density(data = data.frame(x = svr_infection_times), aes(x = x, y = ..density..), fill = "salmon", alpha = 0.5) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+ 
  labs(x = "Infection times", y = "Density") +
  scale_fill_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
  theme_classic() +
  theme(legend.position = "top")

infection_times_hist
fname <- paste0("plots/infection_times_hist_comparison_n_", n, ".pdf")

ggsave(fname, plot=infection_times_hist, device="pdf", width = 6, height = 4)

####### Extracting recovery times 

seir_recovery_times <- seir_data$Time[seir_data$Event == 3]

svr_recovery_times <- svr_data$Time[svr_data$Event == 2]

recovery_times_hist <- ggplot() +
  geom_histogram(data = data.frame(x = seir_recovery_times), aes(x = x, y = ..density..), fill = "steelblue", alpha = 0.5, binwidth = 0.5) +
  geom_histogram(data = data.frame(x = svr_recovery_times), aes(x = x, y = ..density..), fill = "salmon", alpha = 0.5, binwidth = 0.5) +
  labs(x = "Recovery times", y = "Density") +
  scale_fill_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
  theme_classic() +
  theme(legend.position = "top")

recovery_times_hist
fname <- paste0("plots/recovery_times_hist_comparison_n_", n, ".pdf")
ggsave(fname, plot=recovery_times_hist, device="pdf", width = 6, height = 4)




################################################
############# Compare survival functions 
################################################
# Create failure data (time-to-event) and event indicator
seir_failure_data <- data.frame(time = seir_infection_times,
                           event = rep(1, times = length(seir_infection_times)))

svr_failure_data <- data.frame(time = svr_infection_times,
                                event = rep(1, times = length(svr_infection_times)))

# Create a Surv object
seir_surv_obj <- with(seir_failure_data, Surv(time, event))
svr_surv_obj <- with(svr_failure_data, Surv(time, event))


# Estimate the Kaplan-Meier survival curve (empirical survival function)
seir_survival_fit <- survfit(seir_surv_obj ~ 1)
svr_survival_fit <- survfit(svr_surv_obj ~ 1)

df1 <- data.frame(Time = seir_survival_fit$time, 
                          Survival_Probability = seir_survival_fit$surv,
                          Model = "SEIR model")
df2 <- data.frame(Time = svr_survival_fit$time, 
                  Survival_Probability = svr_survival_fit$surv,
                  Model = "SVR model")

combined_df <- rbind(df1, df2)

# Create the plot using ggplot2
emp_surv <- ggplot(combined_df, aes(x = Time, y = Survival_Probability, color = Model)) +
  geom_step() +
  labs(x = "Time", y = "Empirical Survival Function") +
  scale_color_manual(values = c("SEIR model" = "steelblue", "SVR model" = "salmon")) +
  theme_minimal()
emp_surv
fname <- paste0("plots/empirical_survival_comparison_n_", n, ".pdf")
ggsave(fname, plot=emp_surv, device="pdf", width = 6, height = 4)







#seir_data <- read.csv(filename)

ggplot(seir_data, aes(x = Time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), linewidth = 1) +
  geom_line(aes(y = Exposed, color = "Exposed"), linewidth = 1) +
  geom_line(aes(y = Infected, color = "Infected"), linewidth = 1) +
  geom_line(aes(y = Removed, color = "Removed"), linewidth = 1) +
  labs(x = "Time", y = "Population", color = "Compartment") +
  scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange",
                                "Infected" = "red", "Removed" = "green")) +
  theme_minimal()


ggplot(sir_data, aes(x = Time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), linewidth = 1) +
  geom_line(aes(y = Infected, color = "Infected"), linewidth = 1) +
  geom_line(aes(y = Removed, color = "Removed"), linewidth = 1) +
  labs(x = "Time", y = "Population", color = "Compartment") +
  scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange",
                                "Infected" = "red", "Removed" = "green")) +
  theme_minimal()



ggplot(svr_data, aes(x = Time)) +
  geom_line(aes(y = Susceptible, color = "Susceptible"), linewidth = 1) +
  geom_line(aes(y = Total_Infected, color = "Infected"), linewidth = 1) +
  geom_line(aes(y = Removed, color = "Removed"), linewidth = 1) +
  labs(x = "Time", y = "Population", color = "Compartment") +
  scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange",
                                "Infected" = "red", "Removed" = "green")) +
  theme_minimal()



