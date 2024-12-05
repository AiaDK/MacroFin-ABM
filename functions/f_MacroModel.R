
# Generate output of the New Keynesian model
f_MacroModel <- function(macro_params, yps_l, yps_ll, yps_lll, yps_avg_l, yps_avg_ll, ADA_ll, U_ll, fra_ll, sim, noise, s_q_, d_s, c_1, c_2){ # 's_q_, delta_s1'
  
  # parameters
  nu <- macro_params$nu 
  tau <- macro_params$tau
  kappa <- macro_params$kappa
  phi_y <- macro_params$phi_y
  phi_p <- macro_params$phi_p 
  sig_y <- macro_params$sig_y 
  sig_p <- macro_params$sig_p 
  sig_r <- macro_params$sig_r 
  eta <- macro_params$eta 
  iota <- macro_params$iota 
  mu <- macro_params$mu 
  rho <- macro_params$rho 
  gamma <- macro_params$gamma 
  
  # state-space matrices
  A <- matrix(c(1.0, 0.0, tau,
                -kappa, 1.0, 0.0,
                -phi_y, -phi_p, 1.0), nrow = 3, byrow = TRUE)
  
  B <- matrix(c(1.0, tau, c_1,                                        # c_1 here!!
                0.0, nu, 0.0,
                0.0, 0.0, 0.0), nrow = 3, byrow = TRUE)
  
  C <- matrix(c(0.0, 0.0, 0.0,
                0.0, -c_2, 0.0,                                       # c_2 here!!
                0.0, 0.0, 0.0), nrow = 3, byrow = TRUE)
  
  D <- matrix(c(-1.0, 0.0, 0.0,
                0.0, -1.0, 0.0,
                0.0, 0.0, -1.0), nrow = 3, byrow = TRUE)
  
  
  # PREVIOUS period
  
  # initialisation
  U_l <- matrix(0, nrow = nrow(U_ll), ncol = ncol(U_ll))
  fra_l <- matrix(0, nrow = nrow(fra_ll), ncol = ncol(fra_ll))
  
  #v#
  ADA_y_l <- eta * yps_ll[1] + (1.0 - eta) * ADA_ll[1]
  TF_y_l <- yps_ll[1] + iota * (yps_ll[1] - yps_lll[1])
  LAA_y_l <- iota * (yps_avg_ll[1] + yps_ll[1]) + (yps_ll[1] - yps_lll[1])
  ADA_p_l <- eta * yps_ll[2] + (1.0 - eta) * ADA_ll[2]
  TF_p_l <- yps_ll[2] + iota * (yps_ll[2] - yps_lll[2])
  LAA_p_l <- iota * (yps_avg_ll[2] + yps_ll[2]) + (yps_ll[2] - yps_lll[2])
  ADA_dsp_l <- eta * yps_ll[3] + (1.0 - eta) * ADA_ll[3] 
  TF_dsp_l <- yps_ll[3] + iota * (yps_ll[3] - yps_lll[3]) 
  LAA_dsp_l <- iota * (yps_avg_ll[3] + yps_ll[3]) + (yps_ll[3] - yps_lll[3]) 
  
  ADA_l <- matrix(c(ADA_y_l, ADA_p_l, ADA_dsp_l)) 
  
  # forecast performance measures 
  U_l[1] <- rho * U_ll[1] - (ADA_y_l - yps_l[1])^2
  U_l[2] <- rho * U_ll[2] - (TF_y_l - yps_l[1])^2
  U_l[3] <- rho * U_ll[3] - (LAA_y_l - yps_l[1])^2
  U_l[4] <- rho * U_ll[4] - (ADA_p_l - yps_l[2])^2
  U_l[5] <- rho * U_ll[5] - (TF_p_l - yps_l[2])^2
  U_l[6] <- rho * U_ll[6] - (LAA_p_l - yps_l[2])^2
  U_l[7] <- rho * U_ll[7] - (ADA_dsp_l - yps_l[3])^2 
  U_l[8] <- rho * U_ll[8] - (TF_dsp_l - yps_l[3])^2 
  U_l[9] <- rho * U_ll[9] - (LAA_dsp_l - yps_l[3])^2 
  
  # updating of fractions 
  norm123pre <- (exp(gamma * U_l[1]) + exp(gamma * U_l[2]) + exp(gamma * U_l[3]))
  norm456pre <- (exp(gamma * U_l[4]) + exp(gamma * U_l[5]) + exp(gamma * U_l[6]))
  norm789pre <- (exp(gamma * U_l[7]) + exp(gamma * U_l[8]) + exp(gamma * U_l[9])) 
  fra_l[1] <- exp(gamma * U_l[1]) / norm123pre
  fra_l[2] <- exp(gamma * U_l[2]) / norm123pre
  fra_l[3] <- exp(gamma * U_l[3]) / norm123pre
  fra_l[4] <- exp(gamma * U_l[4]) / norm456pre
  fra_l[5] <- exp(gamma * U_l[5]) / norm456pre
  fra_l[6] <- exp(gamma * U_l[6]) / norm456pre
  fra_l[7] <- exp(gamma * U_l[7]) / norm789pre 
  fra_l[8] <- exp(gamma * U_l[8]) / norm789pre 
  fra_l[9] <- exp(gamma * U_l[9]) / norm789pre 
  
  
  # NEW period
  
  # initialization
  ypr_q <- matrix(0, nrow = sim, ncol = 3)  # variable to store our X output
  U_q <- matrix(0, nrow = nrow(U_ll), ncol = ncol(U_ll))
  fra_q <- matrix(0, nrow = nrow(fra_ll), ncol = ncol(fra_ll))
  
  yps_q <- matrix(0, nrow = sim, ncol = 3)  # variable to store output for further calculation
  
  # shocks
  noise[1, , ] <- noise[1, , ] * sig_y
  noise[2, , ] <- noise[2, , ] * sig_p
  noise[3, , ] <- noise[3, , ] * sig_r
  
  # expectations for BRF #v#
  ADA_y <- eta * yps_l[1] + (1.0 - eta) * ADA_l[1]
  TF_y <- yps_l[1] + iota * (yps_l[1] - yps_ll[1])
  LAA_y <- iota * (yps_avg_l[1] + yps_l[1]) + (yps_l[1] - yps_ll[1])
  ADA_p <- eta * yps_l[2] + (1.0 - eta) * ADA_l[2]
  TF_p <- yps_l[2] + iota * (yps_l[2] - yps_ll[2])
  LAA_p <- iota * (yps_avg_l[2] + yps_l[2]) + (yps_l[2] - yps_ll[2])
  ADA_dsp <- eta * yps_l[3] + (1.0 - eta) * ADA_l[3] 
  TF_dsp <- yps_l[3] + iota * (yps_l[3] - yps_ll[3]) 
  LAA_dsp <- iota * (yps_avg_l[3] + yps_l[3]) + (yps_l[3] - yps_ll[3]) 
  
  #!#
  E_y <- fra_l[1] * ADA_y + fra_l[2] * TF_y + fra_l[3] * LAA_y
  E_p <- fra_l[4] * ADA_p + fra_l[5] * TF_p + fra_l[6] * LAA_p
  E_dsp <- fra_l[7] * ADA_dsp + fra_l[8] * TF_dsp + fra_l[9] * LAA_dsp 
  forecasts <- matrix(c(E_y, E_p, E_dsp)) 
  ADA_q <- matrix(c(ADA_y, ADA_p, ADA_dsp)) 
  
  s_q_ <- matrix(rep(s_q_,3))
  
  # final computation
  for (i in 1:sim) {
    
    X <- -solve(A) %*% (B %*% forecasts + C %*% s_q_ + D %*% noise[, i, , drop = FALSE]) 
    
    # model output (y, pi, r)  
    ypr_q[i, 1] <- X[1]
    ypr_q[i, 2] <- X[2]
    ypr_q[i, 3] <- X[3]
    # output for calculation (y, pi, delta_s - pi) 
    yps_q[i, 1] <- X[1]
    yps_q[i, 2] <- X[2]
    yps_q[i, 3] <- d_s - X[2] 
  }
  

  # for SIMULATION and BRF
  #if (sim == 1) {
  
  # forecast performance measures
  U_q[1] <- rho * U_l[1] - (ADA_y - yps_q[1])^2
  U_q[2] <- rho * U_l[2] - (TF_y - yps_q[1])^2
  U_q[3] <- rho * U_l[3] - (LAA_y - yps_q[1])^2
  U_q[4] <- rho * U_l[4] - (ADA_p - yps_q[2])^2
  U_q[5] <- rho * U_l[5] - (TF_p - yps_q[2])^2
  U_q[6] <- rho * U_l[6] - (LAA_p - yps_q[2])^2
  U_q[7] <- rho * U_l[7] - (ADA_dsp - yps_q[3])^2 
  U_q[8] <- rho * U_l[8] - (TF_dsp - yps_q[3])^2 
  U_q[9] <- rho * U_l[9] - (LAA_dsp - yps_q[3])^2 
  
  
  # updating of fractions
  norm123fin <- (exp(gamma * U_q[1]) + exp(gamma * U_q[2]) + exp(gamma * U_q[3]))
  norm456fin <- (exp(gamma * U_q[4]) + exp(gamma * U_q[5]) + exp(gamma * U_q[6]))
  norm789fin <- (exp(gamma * U_q[7]) + exp(gamma * U_q[8]) + exp(gamma * U_q[9])) 
  fra_q[1] <- (exp(gamma * U_q[1]) / norm123fin)
  fra_q[2] <- (exp(gamma * U_q[2]) / norm123fin)
  fra_q[3] <- (exp(gamma * U_q[3]) / norm123fin)
  fra_q[4] <- (exp(gamma * U_q[4]) / norm456fin)
  fra_q[5] <- (exp(gamma * U_q[5]) / norm456fin)
  fra_q[6] <- (exp(gamma * U_q[6]) / norm456fin)
  fra_q[7] <- (exp(gamma * U_q[7]) / norm789fin) 
  fra_q[8] <- (exp(gamma * U_q[8]) / norm789fin) 
  fra_q[9] <- (exp(gamma * U_q[9]) / norm789fin) 
  #}
  
  
  return(list(ypr_q = ypr_q, yps_q = yps_q, ADA_q = ADA_q, U_q = U_q, fra_q = fra_q))
}

