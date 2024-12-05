

f_GenData <- function(setup, c_1, c_2, h, quarters){
  
  ### Initialization ###
  sim <- 1
  
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

  
  # initialize data structure for current series
  yps <- matrix(0, macro_obs, 3) 
  ypr <- matrix(0, macro_obs, 3) 
  ADA <- matrix(0, macro_obs, 3) 
  U <- matrix(0, macro_obs, 9) 
  fra <- matrix(0, macro_obs, 9) 
  
  # set initial values
  fra[1:3, ] <- 1/3
  
  # generate random noise (normally distributed)
  # macro_noise <- rbind(rnorm(macro_obs, mean = 0, sd = macro_params$sig_y), 
  #                      rnorm(macro_obs, mean = 0, sd = macro_params$sig_p), 
  #                      rnorm(macro_obs, mean = 0, sd = macro_params$sig_r))
  macro_noise <- array(rnorm(3 * sim * setup$macro_model$obs, mean = 0, sd = 1), 
                       dim = c(3, sim, setup$macro_model$obs))
  
  
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
  
  # initialize data structures for current series
  p <- matrix(0, fin_obs+1, 1)
  n_f <- matrix(0, fin_obs, 1)
  d_f <- matrix(0, fin_obs, 1)
  d_c <- matrix(0, fin_obs, 1)
  a <- matrix(0, fin_obs, 1)
  
  # set initial values
  p[1:2, ] <- fin_params$p_1
  a[1, ] <- 0
  n_f[1, ] <- 0.5
  d_f[1, ] <- 0
  d_c[1, ] <- 0
  
  # generate random noise (normally distributed)
  # fin_noise <- rbind(rnorm(fin_obs, mean = 0, sd = fin_params$sig_f), 
  #                    rnorm(fin_obs, mean = 0, sd = fin_params$sig_c))
  fin_noise <- array(rnorm(3 * sim * setup$fin_model$obs, mean = 0, sd = 1), 
                     dim = c(2, sim, setup$fin_model$obs))
  
  
  # initialize financial data for first three quarters
  fin_3m_out <- f_GenFinData_3q(fin_params, a, p, n_f, d_f, d_c, sim, fin_noise)
  p <- fin_3m_out$p
  n_f <- fin_3m_out$n_f
  d_f <- fin_3m_out$d_f
  d_c <- fin_3m_out$d_c
  a <- fin_3m_out$a
  
  # initialize s_q and delta_s for first three quarters
  sq_out <- f_sq(macro_obs, p)
  s_q <- sq_out$s_q
  delta_s <- sq_out$delta_s
  
  #s_q <- fin_data$s_q
  #delta_s <- fin_data$delta_s
  #p <- fin_data$p
  
  ### Generation of the Datasets for Financial and Macroeconomic Series ###
  
  # generate macro data for q-th quarter
  for (q in 4:macro_obs){
    
    p_star <- h * ypr[q-1, 1]
    
    # generate fin data for q-th quarter
    for (t in (63*(q-1)+1):(63*q)){
      
      # sim <- 1
      # p_t <- p[t, , drop = FALSE]
      # a_l <- a[t - 1, , drop = FALSE]
      # p_l <- p[t - 1, , drop = FALSE]
      # n_f_l <- n_f[t - 1, , drop = FALSE]
      # d_f_l <- d_f[t - 1, , drop = FALSE]
      # d_c_l <- d_c[t - 1, , drop = FALSE]
      # p_star <- 0
      # noise <- fin_noise[, , t, drop = FALSE]
      
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
      p[t+1, ] <- fin_out$p_t_sim
    }
    
    # update s_q and delta_s
    s_q[q, 1] <- sum(p[(63*(q-1)+1):(63*q)])/63
    delta_s[q] <- diff(s_q[(q-1):q, 1])
    
    # recalculate average mean values
    yps_avg_l <- matrix(c(mean(yps[1:(q-1), 1]), mean(yps[1:(q-1), 2]), mean(yps[1:(q-1), 3])))
    yps_avg_ll <- matrix(c(mean(yps[1:(q-2), 1]), mean(yps[1:(q-2), 2]), mean(yps[1:(q-2), 3])))
    
    ###
    # setup$estim$sim <- 1
    # 
    # yps_l <-  t(yps[q-1, , drop = FALSE])
    # yps_ll <- t(yps[q-2, , drop = FALSE])
    # yps_lll <- t(yps[q-3, , drop = FALSE])
    # yps_avg_l
    # yps_avg_ll
    # ADA_ll <- t(ADA[q-2, , drop = FALSE])
    # U_ll <- t(U[q-2, , drop = FALSE])
    # fra_ll <- t(fra[q-2, , drop = FALSE])
    # sim <- 1
    # noise <- macro_noise[, , q, drop = FALSE]
    # s_q_ <- s_q[q]
    # d_s <- delta_s[q]
    # c_1 <- 0.3
    # c_2 <- 16
    ###
    
    # calculate next macro observation
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
    ypr[q, ] <- macro_out$ypr_q
    ADA[q, ] <- macro_out$ADA_q
    U[q, ] <- macro_out$U_q
    fra[q, ] <- macro_out$fra_q
    yps[q, ] <- macro_out$yps_q
  }
  
  # get quarter names
  dQ <- rep(0, length(quarters)*63)
  for (q in seq_along(quarters)){
    dQ[(63*(q-1)+1):(63*q)] <- quarters[q]
  }
  
  data <- list(p = p, s_q = s_q, delta_s = delta_s, y = ypr[, 1], pi = ypr[, 2], r = ypr[, 3], dQ = dQ)
  return(data)
}
  

# Generate financial data for the first three quarters
f_GenFinData_3q <- function(fin_params, a, p, n_f, d_f, d_c, sim, fin_noise){
  
  p_star <- 0 # set p_star as a constant for first three months
  
  # iteration to calculate first three months of fin data
  for (t in 2:(63*3)) {
    
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
  
  return(list(a = a, p = p, n_f = n_f, d_f = d_f, d_c = d_c))
}

# compute s_q and delta_s for first three quarters of fin data
f_sq <- function(macro_obs, p){
  
  s_q <- matrix(0, macro_obs, 1)
  for (q in 1:3){
    t_start <- (63*(q-1)+1)
    t_end <- (63*q)
    s_q[q, 1] <- sum(p[t_start:t_end])/63
  }
  delta_s_3m <- diff(s_q[1:3, 1])
  delta_s <- matrix(0, macro_obs, 1)
  delta_s[2:3, ] <- delta_s_3m
  
  return(list(s_q = s_q, delta_s = delta_s))
}