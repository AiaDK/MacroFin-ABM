
f_FinModel <- function(fin_params, p_t, a_l, p_l, n_f_l, d_f_l, d_c_l, p_star, sim, noise) {
  
  # parameters
  psi <- fin_params$psi
  chi <- fin_params$chi
  sig_f <- fin_params$sig_f
  sig_c <- fin_params$sig_c
  delta_0 <- fin_params$delta_0
  delta_h <- fin_params$delta_h
  delta_m <- fin_params$delta_m
  ups <- fin_params$ups
  beta <- fin_params$beta
  
  
  # Update the Agent Population Variable
  n_f_t <- 1 / (1 + exp(-beta * a_l))
  
  # Update the Relative Attractiveness
  a_t <- delta_0 + delta_h * (2 * n_f_t - 1) + delta_m * (p_t - p_star) ^ 2
  
  # shocks
  noise[1, , ] <- noise[1, , ] * sig_f
  noise[2, , ] <- noise[2, , ] * sig_c
  
  df_t_sim <- matrix(0, nrow = 1, ncol = sim) 
  dc_t_sim <- matrix(0, nrow = 1, ncol = sim)
  p_t_sim <- matrix(0, nrow = sim, ncol = 1)
  
  for (i in 1:sim) {
    # Update the Net Demands
    df_t <- psi * (p_star - p_t) + noise[1, i, ]
    df_t_sim[, i] <- df_t
    dc_t <- chi * (p_t - p_l) + noise[2, i, ]
    dc_t_sim[, i] <- dc_t
    
    # Update the price
    p_tt <- p_t + ups* (n_f_t * df_t + (1 - n_f_t) * dc_t)
    p_t_sim[i, 1] <- p_tt
  }
  d_f_t <- mean(df_t_sim)
  d_c_t <- mean(dc_t_sim)
  
  return(list(a_t = a_t, p_t_sim = p_t_sim, n_f_t = n_f_t, d_f_t = d_f_t, d_c_t = d_c_t))
}

