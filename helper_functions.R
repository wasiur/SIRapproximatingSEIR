
# Function to simulate the SEIR model using Gillespie's algorithm
simulate_SEIR_gillespie <- function(alpha, beta, gamma, n, init_exposed, init_infected, init_removed, max_time) {
  time <- 0
  S <- n
  E <- init_exposed
  I <- init_infected
  R <- init_removed
  
  event_times <- numeric(0)
  susceptible <- numeric(0)
  exposed <- numeric(0)
  infected <- numeric(0)
  removed <- numeric(0)
  event_types <- numeric(0)
  
  while (time < max_time & I >0 & S>0) {
    rates <- c(beta * S * I / n, alpha * E, gamma * I)
    total_rate <- sum(rates)
    if (total_rate <= 0) break
    
    time_step <- rexp(1, rate = total_rate)
    time <- time + time_step
    
    event <- sample(1:3, size = 1, prob = rates / total_rate)
    if (event == 1) {
      S <- S - 1
      E <- E + 1
    } else if (event == 2) {
      E <- E - 1
      I <- I + 1
    } else if (event == 3) {
      I <- I - 1
      R <- R + 1
    }
    
    event_types <- c(event_types, event)
    event_times <- c(event_times, time)
    susceptible <- c(susceptible, S)
    exposed <- c(exposed, E)
    infected <- c(infected, I)
    removed <- c(removed, R)
  }
  
  return(data.frame(Time = c(0, event_times), Susceptible = c(init_susceptible, susceptible),
                    Exposed = c(init_exposed, exposed), Infected = c(init_infected, infected),
                    Removed = c(init_removed, removed), Event = c(0, event_types)))
}


# Function to simulate the SIR model using Gillespie's algorithm
simulate_SIR_gillespie <- function(beta, gamma, n, init_infected, init_removed, max_time) {
  time <- 0
  S <- n
  I <- init_infected
  R <- init_removed
  
  event_times <- numeric(0)
  susceptible <- numeric(0)
  infected <- numeric(0)
  removed <- numeric(0)
  event_types <- numeric(0)
  
  while (time < max_time & I >0 & S>0) {
    rates <- c(beta * S * I / n, gamma * I)
    total_rate <- sum(rates)
    if (total_rate <= 0) break
    
    time_step <- rexp(1, rate = total_rate)
    time <- time + time_step
    
    event <- sample(1:2, size = 1, prob = rates / total_rate)
    if (event == 1) {
      S <- S - 1
      I <- I + 1
    } else if (event == 2) {
      I <- I - 1
      R <- R + 1
    } 
    
    event_types <- c(event_types, event)
    event_times <- c(event_times, time)
    susceptible <- c(susceptible, S)
    infected <- c(infected, I)
    removed <- c(removed, R)
  }
  
  return(data.frame(Time = c(0, event_times), Susceptible = c(init_susceptible, susceptible),
                    Infected = c(init_infected, infected),
                    Removed = c(init_removed, removed), Event = c(0, event_types)))
}

# SEIR ODE System
seir_ode <- function (t, x, params) {
  infec <- params["beta"]*x[1]*x[3]
  symp <- params["alpha"]*x[2]
  recov <- params["gamma"]*x[3]
  list(c(-infec,infec-symp,symp-recov,recov, infec-recov))
}


# Function to simulate the SIR model using Gillespie's algorithm
simulate_SVR_gillespie <- function(alpha, beta, gamma, rho, n, init_infected, init_removed, max_time) {
  time <- 0
  S <- n
  V <- init_infected
  R <- init_removed
  
  # Initial population values
  S_0 <- 1
  E_0 <- 0
  I_0 <- rho
  R_0 <- 0
  V_0 <- rho
  
  params <- c(alpha=alpha, beta=beta,gamma=gamma)
  init_state <- c(S_0, E_0, I_0, R_0, V_0)
  time_points <- seq(0, max_time, by = 0.01)
  
  # Solve the SEIR ODEs
  seir_solution <- ode(y = init_state, times = time_points, func = seir_ode, parms = params)
  df_seir_sol <- data.frame(seir_solution)
  
  event_times <- numeric(0)
  susceptible <- numeric(0)
  infected <- numeric(0)
  removed <- numeric(0)
  event_types <- numeric(0)
  
  while (time < max_time & V >0 & S>0) {
    
    approx_i <- approx(df_seir_sol$time, df_seir_sol$X3, xout = c(time))
    approx_v <- approx(df_seir_sol$time, df_seir_sol$X5, xout = c(time))
    factor <- approx_i$y/approx_v$y
    
    rates <- c(beta * factor * S * V / n, gamma * factor * V)
    total_rate <- sum(rates)
    if (total_rate <= 0) break
    
    time_step <- rexp(1, rate = total_rate)
    time <- time + time_step
    
    event <- sample(1:2, size = 1, prob = rates / total_rate)
    if (event == 1) {
      S <- S - 1
      V <- V + 1
    } else if (event == 2) {
      V <- V - 1
      R <- R + 1
    } 
    
    event_types <- c(event_types, event)
    event_times <- c(event_times, time)
    susceptible <- c(susceptible, S)
    infected <- c(infected, V)
    removed <- c(removed, R)
  }
  
  return(data.frame(Time = c(0, event_times), Susceptible = c(init_susceptible, susceptible),
                    Total_Infected = c(init_infected, infected),
                    Removed = c(init_removed, removed), Event = c(0, event_types)))
}



comparison_plots <- function(alpha, beta, gamma, rho, n, init_exposed, init_infected, init_removed, max_time){
  seir_data <- simulate_SEIR_gillespie(alpha, beta, gamma, n, init_exposed, init_infected, init_removed, max_time)
  svr_data <- simulate_SVR_gillespie(alpha, beta, gamma, rho, n, init_infected, init_removed, max_time)
  seir_infection_times <- seir_data$Time[seir_data$Event == 1]
  svr_infection_times <- svr_data$Time[svr_data$Event == 1]
  infection_times_hist <- ggplot() +
    geom_density(data = data.frame(x = seir_infection_times), aes(x = x, y = ..density..), fill = "steelblue", alpha = 0.5) +
    geom_density(color="steelblue",alpha=.8, lwd = 1.2, adjust = 1.75)+
    geom_density(data = data.frame(x = svr_infection_times), aes(x = x, y = ..density..), fill = "salmon", alpha = 0.5) +
    geom_density(color="salmon", alpha=.8, lwd = 1.2, adjust = 1.75)+
    labs(x = "Infection times", y = "Density") +
    scale_colour_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
    scale_fill_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
    theme_classic() +
    theme(legend.position = "top")
  
  infection_times_hist
  fname <- paste0("plots/infection_times_hist_comparison_n_", n, ".pdf")
  ggsave(fname, plot=infection_times_hist, device="pdf", width = 4, height = 3)
  
  ####### Extracting recovery times 
  seir_recovery_times <- seir_data$Time[seir_data$Event == 3]
  svr_recovery_times <- svr_data$Time[svr_data$Event == 2]
  
  recovery_times_hist <- ggplot() +
    geom_density(data = data.frame(x = seir_recovery_times), aes(x = x, y = after_stat(density)), fill = "steelblue", alpha = 0.5) +
    geom_density(color="steelblue",alpha=.8, lwd = 1.2, adjust = 1.75)+
    geom_density(data = data.frame(x = svr_recovery_times), aes(x = x, y = after_stat(density)), fill = "salmon", alpha = 0.5) +
    geom_density(color="salmon",alpha=.8, lwd = 1.2, adjust = 1.75)+
    labs(x = "Recovery times", y = "Density") +
    scale_colour_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
    scale_fill_manual(values = c("steelblue", "salmon"), labels = c("SEIR model", "SVR model")) +
    theme_classic() +
    theme(legend.position = "top")
  
  recovery_times_hist
  fname <- paste0("plots/recovery_times_hist_comparison_n_", n, ".pdf")
  ggsave(fname, plot=recovery_times_hist, device="pdf", width = 4, height = 3)
  
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
                    Model = "SIR model")
  
  combined_df <- rbind(df1, df2)
  
  # Create the plot using ggplot2
  emp_surv <- ggplot(combined_df, aes(x = Time, y = Survival_Probability, color = Model)) +
    geom_step() +
    labs(x = "Time", y = "Empirical Survival Function") +
    scale_color_manual(values = c("SEIR model" = "steelblue", "SIR model" = "salmon")) +
    theme_classic()
  emp_surv
  fname <- paste0("plots/empirical_survival_comparison_n_", n, ".pdf")
  ggsave(fname, plot=emp_surv, device="pdf", width = 4, height = 3)
  
}


