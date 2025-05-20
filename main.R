
### Initialization ###
library(yaml) 
library(miceadds)
library(abind)
library(lubridate)
library(philentropy)

experiment <- "test_4"
# Source functions
source.all("functions")
# Load setup
setup <- yaml.load_file("setup_test.yaml")


### Get Data ###

if (is.null(setup$estim$macro_data) & is.null(setup$estim$fin_data)){
  
  # set estimates values for pseudo-empirical data
  c_1 <- 0.309 # Table 2, p. 18
  c_2 <- 16.3 # Table 2, p. 18
  h <- 0.589 # Table 2, p. 18
  
  # generate pseudo-empirical data
  quarters <- read.csv("data/quarters.csv", header = FALSE)
  quarters <- quarters$V1
  set.seed(setup$estim$seed)
  data <- f_GenData(setup, c_1, c_2, h, quarters)
  
  # plot(data$p,
  #       main = "S&P",
  #       xlab = "daily",
  #       ylab = "log price",
  #       type = "l",
  #       col = "blue")
  
  daysQ <- data$dQ
  
  cat("##- Pseudo-empirical data generated! -## \n")
  
} else {
  
  #! DOES NOT WORK WITH REAL DATA -> fix
  
  # load empirical data 
  data <- f_LoadData(setup)
  cat("##- Empirical data loaded! -## \n")
  
  # extra for plots (temp)
  setup$macro_model$obs <- length(data$y)
  setup$fin_model$obs <- length(data$p)
  
  quarters <- data$quarters_names
  daysQ <- data$quarters
}


### Estimation ###
cat(paste0("##- Executing experiment '", experiment, "' -## \n"))
start_time <- Sys.time()

# GA Parameters
bounds <- list(c(-1, 1.2), # c_1
               c(0, 30.0), # c_2
               c(0, 1.5)) # h
n_pop <- 50 # population size
n_bits <- 16 # bits per variable
n_iter <- 100 # number of iterations
r_cross <- 0.9 # crossover rate
r_mut <- 1.0 / (n_bits * length(bounds)) # mutation rate
r_selpr <- 1  # Selection pressure
objective <- f_lol

# Run the genetic algorithm to estimate parameters
set.seed(setup$estim$seed)
result <- f_GA(setup, data, objective, bounds, n_bits, n_iter, n_pop, r_cross, r_mut, r_selpr, quarters, daysQ)
best <- result[[1]]
score <- result[[2]]
decoded <- decode(bounds, n_bits, best)
cat(sprintf("Final Score: f(%s) = %f\n", toString(round(decoded, 5)), score))

# Save results
# write.csv(results, paste0("results/", setup$foldername, "/results.csv")) 
end_time <- Sys.time()
cat(paste0("##- Done estimating model and generating results! Experiment: ", experiment, ". Time taken: ", round(end_time - start_time, 4), " -## \n"))
cat(paste0("##- Results saved to ", "results/", setup$foldername, "/results.csv -## \n"))






