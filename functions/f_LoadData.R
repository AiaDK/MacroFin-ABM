

# Load empirical data 
f_LoadData <- function(setup){
  
  ### Load data from the file specified in the fin_setup ###
  
  # load matrix of empirical observations from csv file 
  fin_emp_series <- data.frame(read.csv(paste0("data/", setup$estim$fin_data)))
  fin_emp_series <- fin_emp_series[order(fin_emp_series$Date),]
  
  # adjustment of the period to fit along with macroeconomic data  
  start_date <- as.Date("1954-07-01")
  end_date <- as.Date("2022-06-30")
  fin_data <- fin_emp_series[-which(fin_emp_series$Date < start_date | fin_emp_series$Date > end_date), ]
  SP <- fin_data$.SP500.Close 
  
  # compute log prices
  log_p <- log(fin_data$.SP500.Close)
  #fin_data$.SP500.Close <- log_p
  
  # compute s and delta_s
  quarters <- quarter(fin_data$Date, with_year = TRUE)
  fin_data <- cbind(fin_data, quarters)
  quarters_names <- unique(quarters)
  s_q <- matrix(0, length(quarters_names), 1)
  for (i in seq_along(quarters_names)){
    s_q[i] <- mean(fin_data[fin_data$quarters == quarters_names[i], 2])
  }
  delta_s <- matrix(c(0, diff(s_q)))
  
  
  ### Load data from the file specified in the macro_setup ###
  
  # load matrix of empirical observations from csv file 
  macro_emp_series <- read.csv(paste0("data/", setup$estim$macro_data))
  macro_data <- macro_emp_series
  
  # separate output
  y <- macro_data$Y_US_all_obs
  pi <- macro_data$Pi_US_all_obs
  r <- macro_data$R_US_all_obs
  
  return(list(SP = SP, p = log_p, s_q = s_q, delta_s = delta_s, y = y, pi = pi, r = r, quarters = quarters, quarters_names = quarters_names, p_dates = as_date(fin_data$Date)))
}






