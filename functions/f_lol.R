
# Calculate negative log-likelihood given the data and the model.
f_lol <- function(setup, estims, fin_noise, macro_noise, data, quarters, daysQ){
  
  ### Initialization ###
  
  c_1 <- estims[1]
  c_2 <- estims[2]
  h <- estims[3]
  delta_und <- 0 # undersmoothing parameter (0 for no undersmoothing)
  sim <- setup$estim$sim # number of simulations
  
  ## Macroeconomic Part ##
  
  macro_obs <- setup$macro_model$obs
  macro_params <- list(tau = c(setup$macro_model$tau), 
                       kappa = c(setup$macro_model$kappa), 
                       phi_y = c(setup$macro_model$phi_y), 
                       phi_pi = c(setup$macro_model$phi_pi), 
                       eta = c(setup$macro_model$eta), 
                       iota = c(setup$macro_model$iota), 
                       mu = c(setup$macro_model$mu), 
                       gamma = c(setup$macro_model$gamma),
                       nu = c(setup$macro_model$nu),
                       rho = c(setup$macro_model$rho),
                       sig_y = c(setup$macro_model$sig_y),
                       sig_p = c(setup$macro_model$sig_pi),
                       sig_r = c(setup$macro_model$sig_r))
  
  # macro data
  ypr <- cbind(data$y, data$pi, data$r)
  yps <- cbind(data$y, data$pi, data$delta_s - data$pi)
  
  # initialize data structure for current series
  msim <- matrix(0, nrow = macro_obs, 1)
  ADA <- matrix(0, macro_obs, 3) 
  U <- matrix(0, macro_obs, 9) 
  fra <- matrix(0, macro_obs, 9) 
  
  # set initial values
  fra[1:3, ] <- 1/3
  
  
  ## Financial Part ##
  
  fin_obs <- setup$fin_model$obs 
  fin_params <- list(psi = c(setup$fin_model$psi), 
                     chi = c(setup$fin_model$chi), 
                     sig_f = c(setup$fin_model$sig_f), 
                     sig_c = c(setup$fin_model$sig_c), 
                     delta_0 = c(setup$fin_model$delta_0), 
                     delta_h = c(setup$fin_model$delta_h), 
                     delta_m = c(setup$fin_model$delta_m), 
                     ups = c(setup$fin_model$ups), 
                     beta = c(setup$fin_model$beta), 
                     p_1 = c(setup$fin_model$p_1))
  
  # fin data
  p <- matrix(data$p, ncol = 1)
  s_q <- data$s_q
  delta_s <- data$delta_s
  
  # initialize data structures for current series
  fsim <- matrix(0, nrow = fin_obs, 1)
  n_f <- matrix(0, fin_obs, 1)
  d_f <- matrix(0, fin_obs, 1)
  d_c <- matrix(0, fin_obs, 1)
  a <- matrix(0, fin_obs, 1)
  
  # set initial values
  a[1, ] <- 0
  n_f[1, ] <- 0.5
  d_f[1, ] <- 0
  d_c[1, ] <- 0
  
  # # initialize financial data for first three quarters
  fin_3m_out <- f_EstFinData_3q(fin_params, a, p, n_f, d_f, d_c, sim, fin_noise, ypr, quarters, daysQ)
  n_f <- fin_3m_out$n_f
  d_f <- fin_3m_out$d_f
  d_c <- fin_3m_out$d_c
  a <- fin_3m_out$a
  
  
  ### Generation of the Datasets for Financial and Macroeconomic Series ###
  for (q in 4:macro_obs){
    
    p_star <- h * yps[q-1, 1]
    
    # generate fin data for q-th quarter
    Q <- quarters[q]
    dQ <- which(daysQ == Q)
    for (t in dQ){ #! (63*(q-1)+1):(63*q)

      fin_out <- f_FinModel(fin_params,
                            p[t, , drop = FALSE],
                            a[t - 1, , drop = FALSE],
                            p[t - 1, , drop = FALSE],
                            n_f[t - 1, , drop = FALSE],
                            d_f[t - 1, , drop = FALSE],
                            d_c[t - 1, , drop = FALSE],
                            p_star, 
                            sim,
                            fin_noise[, , t, drop = FALSE])
      n_f[t, ] <- fin_out$n_f_t
      a[t, ] <- fin_out$a_t
      d_f[t, ] <- fin_out$d_f_t
      d_c[t, ] <- fin_out$d_c_t
      p_sim <- fin_out$p_t
      
      # bandwidth
      hp <- (4 / (5 * sim^(1 + delta_und)))^(1 / 7) * sd(p_sim)
      # if (q == length(quarters) & t == length(daysQ)){ P <-  p[t] } else {P <- p[t+1]} #!# ugly solution -> redo
      K <- dnorm(p_sim, mean = p[t+1], sd = hp) 
      fsim[t] <- max(1e-32, mean(prod(K)))
    }

    # recalculate average mean values
    yps_avg_l <- matrix(c(mean(yps[1:(q-1), 1]), mean(yps[1:(q-1), 2]), mean(yps[1:(q-1), 3])))
    yps_avg_ll <- matrix(c(mean(yps[1:(q-2), 1]), mean(yps[1:(q-2), 2]), mean(yps[1:(q-2), 3])))
    
    # calculate next macro observation
    
    yps_l <-  t(yps[q-1, , drop = FALSE])
    yps_ll <- t(yps[q-2, , drop = FALSE])
    yps_lll <- t(yps[q-3, , drop = FALSE])
    yps_avg_l
    yps_avg_ll
    ADA_ll <- t(ADA[q-2, , drop = FALSE])
    U_ll <- t(U[q-2, , drop = FALSE])
    fra_ll <- t(fra[q-2, , drop = FALSE])
    noise <- macro_noise[, , q, drop = FALSE]
    s_q_ <- s_q[q]
    d_s <- delta_s[q]
    
    macro_out <- f_MacroModel(macro_params, 
                              t(yps[q-1, , drop = FALSE]), 
                              t(yps[q-2, , drop = FALSE]), 
                              t(yps[q-3, , drop = FALSE]), 
                              yps_avg_l, 
                              yps_avg_ll, 
                              t(ADA[q-2, , drop = FALSE]), 
                              t(U[q-2, , drop = FALSE]), 
                              t(fra[q-2, , drop = FALSE]), 
                              sim, 
                              macro_noise[, , q, drop = FALSE], 
                              s_q[q], 
                              delta_s[q], 
                              c_1, 
                              c_2)
    ypr_sim <- macro_out$ypr_q
    ADA[q, ] <- macro_out$ADA_q
    U[q, ] <- macro_out$U_q
    fra[q, ] <- macro_out$fra_q
    
    # bandwidth
      # Bandwidth computed following Wand & Jones (1995, p. 98) + possible undersmoothing by delta_und, see https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation#Rule_of_thumb
      # h = (4/(3 * sim^(1+delta_und)))^(1/5) * std(psim)
    hy <- (4 / (5 * sim^(1 + delta_und)))^(1 / 7) * sd(ypr_sim[, 1])
    hpi <- (4 / (5 * sim^(1 + delta_und)))^(1 / 7) * sd(ypr_sim[, 2])
    hr <- (4 / (5 * sim^(1 + delta_und)))^(1 / 7) * sd(ypr_sim[, 3])
    
    K <- matrix(0, 3, sim)
    K[1, ] <- dnorm(ypr_sim[, 1], mean = ypr[q, 1], sd = hy)
    K[2, ] <- dnorm(ypr_sim[, 2], mean = ypr[q, 2], sd = hpi)
    K[3, ] <- dnorm(ypr_sim[, 3], mean = ypr[q, 3], sd = hr)
    
    msim[q] <- max(1e-32, mean(apply(K, 2, prod)))
  }
  
  macro_logL <- mean(log(msim[4:macro_obs]))
  fin_logL <- mean(log(fsim[(63*(4-1)+1):fin_obs]))
  logL <- mean(macro_logL, fin_logL)
  
  return(-logL)
}

# Generate financial data for the first three quarters
f_EstFinData_3q <- function(fin_params, a, p, n_f, d_f, d_c, sim, fin_noise, ypr, quarters, daysQ){
  
  for (q in 1:3){
    
    p_star <- ypr[q, 1] # set p_star as a constant for first three months
    Q <- quarters[q]
    dQ <- which(daysQ == Q)
    
    if (q == 1) {dQ <- dQ[2:length(dQ)]}
    
    # iteration to calculate first three months of fin data
    for (t in dQ) {
      
      # calculate next observation
      fin_out <- f_FinModel(fin_params,
                            p[t, , drop = FALSE],
                            a[t - 1, , drop = FALSE],
                            p[t - 1, , drop = FALSE],
                            n_f[t - 1, , drop = FALSE],
                            d_f[t - 1, , drop = FALSE],
                            d_c[t - 1, , drop = FALSE],
                            p_star, 
                            1, 
                            fin_noise[, , t, drop = FALSE])
      
      n_f[t, ] <- fin_out$n_f_t
      a[t, ] <- fin_out$a_t
      d_f[t, ] <- fin_out$d_f_t
      d_c[t, ] <- fin_out$d_c_t
      if (sim == 1) {
        p[t + 1, ] <- fin_out$p_t_sim
      }
    }
  }
  
  return(list(a = a, p = p, n_f = n_f, d_f = d_f, d_c = d_c))
}

