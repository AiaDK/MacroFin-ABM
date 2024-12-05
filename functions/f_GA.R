
# Genetic Algorithm
f_GA <- function(setup, data, objective, bounds, n_bits, n_iter, n_pop, r_cross, r_mut, r_selpr, quarters, daysQ){
  
  # Initialize population
  population <- list()
  for (i in 1:n_pop) {
    population[[i]] <- sample(0:1, n_bits * length(bounds), replace = TRUE)
  }
  
  # noise
  #set.seed(setup$estim$seed)
  macro_noise <- array(rnorm(3 * setup$estim$sim * setup$macro_model$obs, mean = 0, sd = 1), 
                       dim = c(3, setup$estim$sim, setup$macro_model$obs))
  #set.seed(setup$estim$seed)
  fin_noise <- array(rnorm(3 * setup$estim$sim * setup$fin_model$obs, mean = 0, sd = 1), 
                     dim = c(2, setup$estim$sim, setup$fin_model$obs))

  
  # Keep track of the best solution
  best <- population[[1]]
  best_eval <- objective(setup, decode(bounds, n_bits, best), fin_noise, macro_noise, data, quarters, daysQ)
  
  for (gen in 1:n_iter) {
    
    # Decode the population to continuous values
    decoded <- list()
    for (i in 1:n_pop) {
      decoded[[i]] <- decode(bounds, n_bits, population[[i]])
    }
    
    # Evaluate all candidates in the population
    scores <- c()
    for (i in seq_along(decoded)) {
      scores[i] <- objective(setup, decoded[[i]], fin_noise, macro_noise, data, quarters, daysQ)
    }
    
    # Check for a new best solution
    for (i in 1:n_pop) {
      if (scores[i] < best_eval) {
        best <- population[[i]]
        best_eval <- scores[i]
        cat(sprintf("Generation %d> New Best: f(%s) = %f\n", gen, toString(round(decoded[[i]], 4)), scores[i]))
      }
    }
    
    # Calculate probabilities for roulette wheel selection
    avg_score <- mean(scores)
    if (avg_score != 0) {
      scores <- scores / avg_score
    }
    probs <- exp(-r_selpr * scores)
    
    # Select parents
    selection <- list()
    for (i in 1:n_pop) {
      selection[[i]] <- roulette_wheel_selection(probs, population)
    }
    
    # Create the next generation
    children <- list()
    for (i in seq(1, n_pop, by = 2)) {
      # Get selected parents in pairs
      p1 <- selection[[i]]
      p2 <- selection[[i + 1]]
      # Crossover and mutation
      offspring <- crossover(p1, p2, r_cross)
      for (c in offspring) {
        c <- mutate(c, r_mut)
        # Store for the next generation
        children <- append(children, list(c))
      }
    }
    # Replace the population with the new generation
    population <- children
  }
  
  return(list(best, best_eval))
}


decode <- function(bounds, n_bits, bitstring) {
  decoded <- c()
  largest <- 2^n_bits
  for (i in seq_along(bounds)) {
    # Extract the substring for the current variable
    start <- (i - 1) * n_bits + 1
    end <- i * n_bits
    substring <- bitstring[start:end]
    # Convert bitstring to an integer
    chars <- paste0(substring, collapse = "")
    integer <- strtoi(chars, base = 2)
    # Scale integer to the desired range
    value <- bounds[[i]][1] + (integer / largest) * (bounds[[i]][2] - bounds[[i]][1])
    # Store the decoded value
    decoded <- c(decoded, value)
  }
  return(decoded)
}

# Roulette wheel selection for choosing parents based on their probabilities
roulette_wheel_selection <- function(p, population) {
  # Cumulative sum of probabilities
  c <- cumsum(p)  
  # Random number in the range of the sum of probabilities
  r <- sum(p) * runif(1)  
  # Find the index where the random number fits
  ind <- which(r <= c)[1]  
  # Select individual
  selected <- population[[ind]]
  return(selected)  
}

# Crossover two parents to create two children
crossover <- function(p1, p2, r_cross) {
  # Children are copies of parents by default
  c1 <- p1
  c2 <- p2
  # Check for recombination
  if (runif(1) < r_cross) {
    # Select crossover point that is not at the end of the string
    pt <- sample(2:(length(p1)-1), 1)
    # Perform crossover
    c1 <- c(p1[1:pt], p2[(pt+1):length(p2)])
    c2 <- c(p2[1:pt], p1[(pt+1):length(p1)])
  }
  return(list(c1, c2))
}

# Mutation operator
mutate <- function(bitstring, r_mut) {
  for (i in seq_along(bitstring)) {
    # Check for a mutation
    if (runif(1) < r_mut) {
      # Flip the bit
      bitstring[i] <- 1 - bitstring[i]
    }
  }
  return(bitstring)
}
