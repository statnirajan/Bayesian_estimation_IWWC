library(tidyverse)
sample_data <- function(ssize, theta, gamma) {
  u <- runif(ssize)
  x <- rep(NA, ssize)
  
  for (j in 1:ssize) {
    if (u[j] <= 0.5) {    
      x[j] <- theta * (log(exp(1) / (2 * u[j])))^(-1 / gamma)          
    } else {
      x[j] <- theta * (log(exp(1) / (2 * (1 - u[j]))))^(1 / gamma)
    }
  }
  
  return(x)
}

num_simulations <- 250

est_gamma_ls <- numeric(num_simulations)
est_theta_ls <- numeric(num_simulations)
est_gamma <- numeric(num_simulations)
est_theta <- numeric(num_simulations)


all_mse_gamma <- numeric(num_simulations)
all_bias_gamma <- numeric(num_simulations)
all_mse_theta <- numeric(num_simulations)
all_bias_theta <- numeric(num_simulations)
k_p <- numeric(num_simulations)

interval_gamma_l = matrix(NA,num_simulations, 2)
interval_theta_l = matrix(NA, num_simulations, 2)
interval_gamma_b = matrix(NA, num_simulations, 2)
interval_theta_b = matrix(NA, num_simulations, 2)

gamma_coverage_l = numeric(num_simulations)
theta_coverage_l = numeric(num_simulations)
gamma_coverage_b = numeric(num_simulations)
theta_coverage_b = numeric(num_simulations)

acceptance_rate_eta_2 <- numeric(num_simulations)
acceptance_rate_eta_1 <- numeric(num_simulations)

all_gamma_samples <- list()
all_theta_samples <- list()



## True values
theta_true = 25
gamma_true = 0.5

set.seed(3234668)

for (s in 1:num_simulations) {
  
  ##############################################################################
  repeat{
    ##############################################################################
    # Step 1: Sample a new data set
    #set.seed(1994 + s) #Remove
    ssize = 50                                                    #n=10
    xdata <- sample_data(ssize=ssize,theta=theta_true,gamma=gamma_true)  
    Xa = sort(xdata)
    a=1:length(Xa)
    
    x1 <-as.data.frame(cbind(a, Xa))
    
    n <- length(Xa)
    gammahat <-rep(NA,n)
    thetahat <- rep(NA,n)
    yi <- log(Xa)
    for (m in 1:(n-1)){
      
      x1<-x1 %>%
        mutate(X1i=-log(log(((ssize+0.4)*exp(1))/(2*(a-0.3)))))  
      x1<-x1 %>%
        mutate(X2i=log(log(((ssize+0.4)*exp(1))/(2*(ssize-a+0.7)))))
      x1<-x1 %>%
        mutate(yi)
      
      c1=subset(x1,a<=m,select=c(X1i,yi))
      
      sum(c1$X1i)
      
      c2=subset(x1,a>m,select=c(X2i,yi))
      
      sum(c2$X2i)
      xbar=(sum(c1$X1i)+sum(c2$X2i))/ssize
      xbar
      ybar=sum(yi)/ssize
      ybar
      c1<-c1 %>%
        mutate(up1=(X1i-xbar)*(yi-ybar))
      c2<-c2 %>%
        mutate(up2=(X2i-xbar)*(yi-ybar))
      up=sum(c1$up1)+sum(c2$up2)
      c1<-c1 %>%
        mutate(up11=(X1i-xbar)^2)
      c2<-c2 %>%
        mutate(up22=(X2i-xbar)^2)
      lw=sum(c1$up11)+sum(c2$up22)
      amhat=up/lw
      gammahat[m]=1/amhat
      
      
      bmhat=ybar-amhat*xbar
      thetahat[m]=exp(bmhat)
      
    }
    cbind(x1[-n],thetahat)
    K=rep(NA,n-1)
    Z<-as.data.frame(x1)
    Xa
    for (p in 1:(n-1)){
      if (Xa[p] < thetahat[p] & thetahat[p] < Xa[p+1]){
        K[p]<-1
      }else{
        K[p]<-0
      }
    }
    K
    m_position=which(K %in% 1)
    if (length(m_position) == 0) {
      # Handle the case where m_position is empty.
      # For example, you can skip the current iteration of the simulation.
      next
    }
    m_position
    
    ######### Least Squares Estimate #######
    
    theta_ls <- thetahat[m_position]
    gamma_ls <- gammahat[m_position]
    
    ## Standard Errors
    amhat = 1/gamma_ls
    bmhat = log(theta_ls)
    pred.y1 = amhat*x1$X1i[1:m_position] + bmhat
    pred.y2 = amhat*x1$X2i[(m_position + 1):nrow(x1)] + bmhat
    pred = c(pred.y1,pred.y2)
    pred
    res = x1$yi - pred
    sig = sqrt((sum(res[1:m_position]^2) + sum (res[(m_position + 1):nrow(x1)]^2))/(nrow(x1)-2))
    
    l.x1=x1$X1i[1:m_position]-mean(c(x1$X1i[1:m_position],x1$X2i[(m_position + 1):nrow(x1)]))
    s1 = sum(l.x1^2)
    l.x2=x1$X2i[(m_position + 1):nrow(x1)]-mean(c(x1$X1i[1:m_position],x1$X2i[(m_position + 1):nrow(x1)]))
    s2 = sum(l.x2^2)
    sqrt(sum(s1,s2))
    se_amhat = sig/sqrt(sum(s1,s2))  # This is SE(amhat)
    se_gammahat = se_amhat/(amhat^2) #This is SE(gammamhat) 
    
    se_bmhat = sqrt((sum((x1$X1i[1:m_position])^2) + sum((x1$X2i[(m_position + 1):nrow(x1)])^2)) / (n * sum(s1,s2)))*sig # This is SE(bmhat)
    se_thetahat = exp(bmhat) * se_bmhat   #This is SE(thetahat) 
    
    interval_gamma_l[s,] = c(gamma_ls-1.96*se_gammahat,gamma_ls+1.96*se_gammahat)
    interval_theta_l[s,] = c(theta_ls-1.96*se_thetahat,theta_ls+1.96*se_thetahat)
    
    gamma_coverage_l[s] = if (interval_gamma_l[s,1] < gamma_true & interval_gamma_l[s,2] > gamma_true) 1 else 0
    theta_coverage_l[s] = if (interval_theta_l[s,1] < theta_true & interval_theta_l[s,2] > theta_true) 1 else 0
    
    # estimates
    est_gamma_ls[s] = gamma_ls
    est_theta_ls[s] = theta_ls
    
    # Data and Initialization
    initial_eta_1 <- log(thetahat[m_position])
    initial_eta_2 <- log(gammahat[m_position])
    N <- 10000
    m <- m_position
    n <- ssize
    
    # Log of full conditional/posterior for eta.1
    log_posterior_eta_1 <- function(eta_1, xdata, n, m, eta_2) {
      (eta_1 * exp(eta_2) * (2 * m - n)) -
        sum((exp(eta_1) / xdata[1:m])^(exp(eta_2))) -
        sum((xdata[(m+1):n] / exp(eta_1))^(exp(eta_2)))
    }
    
    # Log of the full conditional/posterior for eta.2
    log_posterior_eta_2 <- function(eta_2, xdata, n, m, eta_1) {
      (eta_1 * exp(eta_2) * (2 * m - n)) +
        ((n + 1) * eta_2) +
        (exp(eta_2) * sum(log(xdata[(m+1):n]))) -
        (exp(eta_2) * sum(log(xdata[1:m]))) -
        sum((exp(eta_1) / xdata[1:m])^(exp(eta_2))) -
        sum((xdata[(m+1):n] / exp(eta_1))^(exp(eta_2)))
    }
    
    # MCMC Initialization
    current_eta_1 <- initial_eta_1
    current_eta_2 <- initial_eta_2
    acceptance_eta_1 <- acceptance_eta_2 <- numeric(N)
    eta_1_samples <- eta_2_samples <- numeric(N)
    
    # MCMC Loop
    set.seed(2994 + s)
    for (j in 1:N) {
      
      # Update for eta.2
      candidate_eta_2 <- rnorm(1, mean = current_eta_2, sd = sqrt(0.04))
      acceptance_prob_eta_2 <- exp(log_posterior_eta_2(candidate_eta_2, xdata, n, m, current_eta_1) -
                                     log_posterior_eta_2(current_eta_2, xdata, n, m, current_eta_1))
      if (runif(1) < min(acceptance_prob_eta_2, 1)) {
        current_eta_2 <- candidate_eta_2
        acceptance_eta_2[j] <- 1  # Record the acceptance
      } else {
        acceptance_eta_2[j] <- 0  # Record the rejection
      }
      eta_2_samples[j] <- current_eta_2
      
      # Update for eta.1
      candidate_eta_1 <- rnorm(1, mean = current_eta_1, sd = sqrt(0.22))
      acceptance_prob_eta_1 <- exp(log_posterior_eta_1(candidate_eta_1, xdata, n, m, current_eta_2) -
                                     log_posterior_eta_1(current_eta_1, xdata, n, m, current_eta_2))
      if (runif(1) < min(acceptance_prob_eta_1, 1)) {
        current_eta_1 <- candidate_eta_1
        acceptance_eta_1[j] <- 1  # Record the acceptance
      } else {
        acceptance_eta_1[j] <- 0  # Record the rejection
      }
      eta_1_samples[j] <- current_eta_1
    }
    
    # Post-Processing
    burn_in <- 2001:N
    acceptance_rate_eta_2[s] <- mean(acceptance_eta_2[burn_in])
    acceptance_rate_eta_1[s] <- mean(acceptance_eta_1[burn_in])
    
    ################################################
    if ((acceptance_rate_eta_1[s] > 0.35 & acceptance_rate_eta_1[s] < 0.50) &
        (acceptance_rate_eta_2[s] > 0.35 & acceptance_rate_eta_2[s] < 0.50)) 
      break
  }
  #########################################################
  gamma_samples <- exp(eta_2_samples)
  theta_samples <- exp(eta_1_samples)
  
  k_p[s] = m
  all_gamma_samples[[s]] <- gamma_samples
  all_theta_samples[[s]] <- theta_samples
  
  # Point estimates
  est_gamma[s] = mean(gamma_samples[burn_in])
  est_theta[s] = mean(theta_samples[burn_in])
  
  
  interval_gamma_b[s,] = quantile(gamma_samples[burn_in],c(0.025,0.975), name=FALSE)
  interval_theta_b[s,] = quantile(theta_samples[burn_in],c(0.025,0.975), name=FALSE)
  
  
  gamma_coverage_b[s] = if (interval_gamma_b[s,1] < gamma_true & interval_gamma_b[s,2] > gamma_true) 1 else 0
  theta_coverage_b[s] = if (interval_theta_b[s,1] < theta_true & interval_theta_b[s,2] > theta_true) 1 else 0
  
  # Calculate MSE and Bias for this simulation
  all_mse_gamma[s] <- (mean(gamma_samples[burn_in]) - gamma_true)^2
  all_bias_gamma[s] <- mean(gamma_samples[burn_in]) - gamma_true
  all_mse_theta[s] <- (mean(theta_samples[burn_in]) - theta_true)^2
  all_bias_theta[s] <- mean(theta_samples[burn_in]) - theta_true

}



###################### ALL LS ESTIMATES ###################################
(ls_cp_gamma = mean(gamma_coverage_l))
(ls_bias_gamma = mean(est_gamma_ls-gamma_true))
(ls_mse_gamma = mean((est_gamma_ls-gamma_true)^2))

(ls_cp_theta = mean(theta_coverage_l))
(ls_bias_theta = mean(est_theta_ls-theta_true))
(ls_mse_theta = mean((est_theta_ls-theta_true)^2))


###################### ALL Bayes ESTIMATES_Jeffrey's Prior #################
(bayes_cp_gamma = mean(gamma_coverage_b))
(bayes_bias_gamma <- mean(all_bias_gamma))
(bayes_mse_gamma <- mean(all_mse_gamma))

(bayes_cp_theta = mean(theta_coverage_b))
(bayes_bias_theta <- mean(all_bias_theta))
(bayes_mse_theta <- mean(all_mse_theta))

###########################################################################

acceptance_rate_eta_2
acceptance_rate_eta_1
summary(acceptance_rate_eta_2)
summary(acceptance_rate_eta_1)
mean(acceptance_rate_eta_2)
mean(acceptance_rate_eta_1)







iteration_with_lowest_acceptance_rate_eta_2 <- which.min(acceptance_rate_eta_2)
iteration_with_lowest_acceptance_rate_eta_2

iteration_with_lowest_acceptance_rate_eta_1 <- which.min(acceptance_rate_eta_1)
iteration_with_lowest_acceptance_rate_eta_1


iteration_with_highest_acceptance_rate_eta_2 <- which.max(acceptance_rate_eta_2)
iteration_with_highest_acceptance_rate_eta_2

iteration_with_highest_acceptance_rate_eta_1 <- which.max(acceptance_rate_eta_1)
iteration_with_highest_acceptance_rate_eta_1
# Trace plots after the MCMC loop
par(mfrow=c(2,1))  # Set up the plotting area for two plots

# Trace plot for gamma

plot(all_gamma_samples[[53]], type="l", col="blue", main="Trace Plot for Gamma", ylab="Gamma", xlab="Iteration")

# Trace plot for theta

plot(all_theta_samples[[74]], type="l", col="red", main="Trace Plot for Theta", ylab="Theta", xlab="Iteration")


# Trace plot for gamma - highest

plot(all_gamma_samples[[175]], type="l", col="blue", main="Trace Plot for Gamma", ylab="Gamma", xlab="Iteration")

# Trace plot for theta

plot(all_theta_samples[[29]], type="l", col="red", main="Trace Plot for Theta", ylab="Theta", xlab="Iteration")


# Average over all simulations

(avg_est_gamma <- mean(est_gamma))
quantile(est_gamma,c(0.025,0.975))
(avg_est_theta <- mean(est_theta))
quantile(est_theta,c(0.025,0.975))

avg_mse_gamma <- mean(all_mse_gamma)
avg_bias_gamma <- mean(all_bias_gamma)
avg_mse_theta <- mean(all_mse_theta)
avg_bias_theta <- mean(all_bias_theta)

# Print the results
cat("Average MSE for Gamma:", avg_mse_gamma, "\n")
cat("Average Bias for Gamma:", avg_bias_gamma, "\n")
cat("Average MSE for Theta:", avg_mse_theta, "\n")
cat("Average Bias for Theta:", avg_bias_theta, "\n")
