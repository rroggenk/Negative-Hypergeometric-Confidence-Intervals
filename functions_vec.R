###################################
#---------------------------------#
# Preliminary Functions           #
#---------------------------------#
###################################

ngh_pmf = function(x, N, M, m) {
  numerator = choose(m + x - 1, m - 1) * choose(N - m - x, M - m)
  denominator = choose(N, M)
  result = numerator / denominator
  return(result)
}


ngh_cdf = function(x, N, M, m, lower_tail = TRUE) {
  # m = # total successes (unknown) - (our notation: M)
  m_pmf = M  
  
  # n = # total failures - (our notation: X = N - M)
  n_pmf = N - m_pmf 
  
  # r = # fixed successes (our notation: m)
  r_pmf = m   
  
  # x = # balls being drawn (our notation: n = m + x) 
  x_pmf = r_pmf + x
  
  return(pnhyper(q = x_pmf, m = m_pmf, n = n_pmf, r = r_pmf, lower.tail = lower_tail))
}


sum_ngh_pmf <- function(N, M, m, min_x, max_x) {
  sum_pmf = 0
  for (x in min_x:max_x) {
    sum_pmf = sum_pmf + ngh_pmf(x, N, M, m)
  }
  return(sum_pmf)
}



###################################
#---------------------------------#
# M Unknown                       #
#---------------------------------#
###################################

###################################
#---------------------------------#
# Normal Approximation (MLE)      #
#---------------------------------#
###################################

M_MLE_function <- function(m, x, N) {
  M_MLE = (m / (m + x)) * N
  M_MLE = ceiling(M_MLE)
  return(M_MLE)
}


SE_MLE_function <- function(M_MLE, N, m) {
  part1 = (N * (M_MLE + 1)) / (N + 1)
  part2 = ((N - M_MLE) * (M_MLE - m + 1)) / (m * (N + 1) * (M_MLE + 2))
  part2 = sqrt(part2)
  return(part1 * part2)
}


CI_cov_prob_MLE <- function(N, m, conf_level = 0.95) {
  results = data.frame(x = 0:(N), lower_bound = NA, upper_bound = NA)
  
  for (xi in 0:(N)) {
    # Speical End Cases: When x in N-m to N
    if ((xi >= (N - m)) & (xi <= N)) {
      CI_lb = N - xi
      CI_ub = N - xi
    }
    
    else {
      M_MLE = M_MLE_function(m, xi, N)
      SE_MLE = SE_MLE_function(M_MLE, N, m)
      
      z = qnorm(1 - ((1 - conf_level)/2))
      
      CI_lb = M_MLE - (z * SE_MLE)
      CI_ub = M_MLE + (z * SE_MLE)
    }
    # Store the results
    results[xi + 1, "lower_bound"] = CI_lb
    results[xi + 1, "upper_bound"] = CI_ub
  }
  
  return(results)
}



###################################
#---------------------------------#
# Normal Approximation (Unbiased) #
#---------------------------------#
###################################

M_unbiased_function <- function(m, x, N) {
  M_unbiased = ((m - 1) / (m + x - 1)) * N
  return(M_unbiased)
}


SE_unbiased_function <- function(M_unbiased, N, m) {
  part1 = (N * (m - 1) * (M_unbiased + 1)) / ((m * N) - M_unbiased + m - 1)^2
  part2 = (m * (N - M_unbiased) * (M_unbiased - m + 1) * (N + 1)) / (M_unbiased + 2)
  part2 = sqrt(part2)
  return(part1 * part2)
}


CI_cov_prob_unbiased <- function(N, m, conf_level = 0.95) {
  results = data.frame(x = 0:(N), lower_bound = NA, upper_bound = NA)
  
  for (xi in 0:(N)) {
    # Speical End Cases: When x in N-m to N
    if ((xi >= (N - m)) & (xi <= N)) {
      CI_lb = N - xi
      CI_ub = N - xi
    }
    
    else {
      M_unbiased = M_unbiased_function(m, xi, N)
      SE_unbiased = SE_unbiased_function(M_unbiased, N, m)
      
      z = qnorm(1 - ((1 - conf_level)/2))
      
      CI_lb = M_unbiased - (z * SE_unbiased)
      CI_ub = M_unbiased + (z * SE_unbiased)
    }
    
    # Store the results
    results[xi + 1, "lower_bound"] = CI_lb
    results[xi + 1, "upper_bound"] = CI_ub
  }
  
  return(results)
}



###################################
#---------------------------------#
# Analog to Clopper Pearson       #
#---------------------------------#
###################################

CI_cov_prob <- function(N, m, conf_level = 0.95) {
  target_probability = (1 - conf_level) / 2
  
  #results = data.frame(x = 0:(N - m), lower_bound = NA, upper_bound = NA)
  results = data.frame(x = 0:(N), lower_bound = NA, upper_bound = NA)
  
  #for (xi in 0:(N - m)) {
  for (xi in 0:(N)) {
    # Speical End Cases: When x in N-m to N
    if ((xi >= (N - m)) & (xi <= N)) {
      lower_bound = N - xi
      upper_bound = N - xi
    }
    
    else {
      lower_bound = m
      upper_bound = N
      
      # Find lower bound
      for (M in m:N) {
        area_left = ngh_cdf(x = xi, N = N, M = M, m = m, lower_tail = TRUE)
        
        #if (area_left == target_probability) {
        if (isTRUE(all.equal(area_left, target_probability))) {
          lower_bound = M
          break
        }
        
        else if (area_left > target_probability) {
          lower_bound = M
          break
        }
      }
      
      # Find upper bound
      for (M in N:m) {
        area_right = ngh_cdf(x = xi - 1, N = N, M = M, m = m, lower_tail = FALSE)
        
        #if (area_right == target_probability) {
        if (isTRUE(all.equal(area_right, target_probability))) {
          upper_bound = M
          break
        }
        
        else if (area_right > target_probability) {
          upper_bound = M
          break
        }
      }
    }
    
    # Store the results
    results[xi + 1, "lower_bound"] = lower_bound
    results[xi + 1, "upper_bound"] = upper_bound
  }
  
  return(results)
}



###################################
#---------------------------------#
# Minimal Cardinality             #
#---------------------------------#
###################################

all_mc_ac_vec <- function(N, m, conf_level = 0.95) {
  # Initialize the final results
  results <- data.frame(
    M              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )
  
  # Start with these constraints
  min_a <- 0
  min_b <- 0
  
  # Loop over M from N down to 0
  for (M_val in seq(N, 0, by = -1)) {
    
    # We'll build temp_results for the current M_val
    temp_results <- data.frame(
      M              = integer(),
      a              = integer(),
      b              = integer(),
      cardinality    = integer(),
      coverage_prob  = numeric()
    )
    
    # Special End Case: When M < m
    if (M_val < m) {
      a_val <- N - M_val
      b_val <- N - M_val
      
      coverage_prob <- sum_ngh_pmf(N, M_val, m, a_val, b_val)
      cardinality   <- b_val - a_val + 1
      
      temp_results <- rbind(
        temp_results,
        data.frame(
          M              = M_val,
          a              = a_val,
          b              = b_val,
          cardinality    = cardinality,
          coverage_prob  = coverage_prob
        )
      )
      
    } else {
      # If M >= m, we build all (a,b) pairs with:
      #   a in [min_a, N - M_val]
      #   b in [min_b, N - M_val]
      #   b >= a
      # Then compute coverage in a vectorized fashion.
      
      ab_grid <- expand.grid(
        a = seq.int(min_a, N - M_val),
        b = seq.int(min_b, N - M_val)
      ) %>%
        filter(b >= a)
      
      # If there are no (a,b) pairs, skip
      if (nrow(ab_grid) > 0) {
        coverage_vec <- mapply(
          FUN = function(a_val, b_val) {
            sum_ngh_pmf(N, M_val, m, a_val, b_val)
          },
          ab_grid$a,
          ab_grid$b
        )
        
        temp_results <- data.frame(
          M              = M_val,
          a              = ab_grid$a,
          b              = ab_grid$b,
          cardinality    = ab_grid$b - ab_grid$a + 1,
          coverage_prob  = coverage_vec
        )
      }
    }
    
    # If M_val >= m, we apply the coverage/confidence filter
    if (M_val >= m) {
      temp_results <- temp_results %>%
        filter(coverage_prob >= conf_level) %>%
        group_by(M) %>%
        slice_min(order_by = cardinality, with_ties = TRUE) %>%
        ungroup()
      
      # Update min_a, min_b if we found valid intervals
      if (nrow(temp_results) > 0) {
        min_a <- max(min_a, min(temp_results$a))
        min_b <- max(min_b, min(temp_results$b))
      }
    }
    
    # Append current iteration’s data to the final results
    results <- rbind(results, temp_results)
  }
  
  # Separate results into M >= m and M < m
  results_M_ge_m <- results %>% filter(M >= m)
  results_M_lt_m <- results %>% filter(M < m)
  
  # Combine and order
  filtered_results <- bind_rows(results_M_ge_m, results_M_lt_m) %>%
    arrange(desc(M))
  
  # Add x_set column
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}


minimal_cardinality_ci_vec <- function(N, m, conf_level = 0.95, procedure = "MST") {
  # Chooses which minimal cardinality procedure 
  if (procedure == "MST") {
    results = mst_ac_vec(N, m, conf_level)
  } 
  else if (procedure == "CG") {
    results = cg_ac_vec(N, m, conf_level)
  } 
  else if (procedure == "BK") {
    results = bk_ac_vec(N, m, conf_level)
  } 
  else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Loops through each x 
  for (x in 0:N) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the M of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub = first_occurrence$M
      ci_lb = last_occurrence$M
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}



###################################
#---------------------------------#
# MST                             #
#---------------------------------#
###################################

mst_ac_vec <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves 
  results = all_mc_ac_vec(N, m, conf_level)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(M = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each M 
  for (current_M in N:0) {
    # Only looks at acceptance curves for the current M 
    subset_results = results %>% 
      filter(M == current_M)
    
    # If only one acceptance curve, then that is the acceptance curve  
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply MST procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the highest 
    # coverage probability 
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(desc(coverage_prob), desc(a), desc(b)) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(desc(M))
  
  return(final_results)
}



###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_vec <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac_vec(N, m, conf_level)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(M = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each M
  for (current_M in N:0) {
    # Only looks at acceptance curves for the current M 
    subset_results = results %>% 
      filter(M == current_M)
    
    # If only one acceptance curve, then that is the acceptance curve 
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply CG procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the largest 
    # possible a and b
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(a, b) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(desc(M))
  
  return(final_results)
}


###################################
#---------------------------------#
# BK                              #
#---------------------------------#
###################################

bk_ac_vec <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac_vec(N, m, conf_level)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(M = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each M
  for (current_M in N:0) {
    # Only looks at acceptance curves for the current M 
    subset_results = results %>% 
      filter(M == current_M)
    
    # If only one acceptance curve, then that is the acceptance curve 
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply MST procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the smallest 
    # possible a and b
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(desc(a), desc(b)) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(desc(M))
  
  return(final_results)
}



###################################
#---------------------------------#
# Blaker                          #
#---------------------------------#
###################################

blaker_ac_vec <- function(N, m, conf_level = 0.95) {
  alpha <- 1 - conf_level
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1) FIRST PASS: Build 'results' with (M, x, min_tail_prob, pmf_x, acceptance_set)
  #    in a more vectorized way.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  results <- data.frame(
    M              = integer(),
    x              = integer(),
    min_tail_prob  = numeric(),
    pmf_x          = numeric(),   # we'll store PMF here to avoid repeated calls later
    acceptance_set = character()
  )
  
  for (M_val in seq.int(0, N)) {
    
    # If M_val < m, there's exactly one x = N - M_val, with no min_tail_prob needed
    if (M_val < m) {
      x_val <- N - M_val
      # Append a single row
      results <- rbind(
        results,
        data.frame(
          M              = M_val,
          x              = x_val,
          min_tail_prob  = NA,   # as in your original code
          pmf_x          = NA,   # won't be used if M_val < m
          acceptance_set = as.character(x_val)
        )
      )
      
    } else {
      # M_val >= m => we consider x from 0 to (N - M_val)
      x_seq <- seq.int(0, N - M_val)
      
      # Vectorized calls to ngh_cdf() for "left" and "right" tails:
      area_left_vec <- sapply(
        x_seq,
        function(xx) ngh_cdf(x = xx,    N = N, M = M_val, m = m, lower_tail = TRUE )
      )
      area_right_vec <- sapply(
        x_seq,
        function(xx) ngh_cdf(x = xx - 1, N = N, M = M_val, m = m, lower_tail = FALSE)
      )
      
      # Vector of min_tail_prob
      min_tail_prob_vec <- pmin(area_left_vec, area_right_vec)
      
      # Precompute PMF for each x as well
      pmf_vec <- sapply(
        x_seq,
        function(xx) ngh_pmf(x = xx, N = N, M = M_val, m = m)
      )
      
      # Build a temporary data frame for all x in [0, N - M_val]
      tmp_df <- data.frame(
        M              = M_val,
        x              = x_seq,
        min_tail_prob  = min_tail_prob_vec,
        pmf_x          = pmf_vec,
        acceptance_set = NA_character_
      )
      
      # Bind once per M_val
      results <- rbind(results, tmp_df)
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2) SECOND PASS: Build final_results.
  #    We do exactly what your original code did, but we skip repeated PMF calls,
  #    because 'results' already has 'pmf_x'.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  final_results <- data.frame(
    M               = integer(),
    acceptance_set  = character(),
    a               = integer(),
    b               = integer(),
    cardinality     = integer(),
    coverage_prob   = numeric(),
    x_set           = character(),
    gap             = logical()
  )
  
  # We iterate from M=0 up to M=N, same as your code
  for (fixed_M in seq.int(0, N)) {
    
    acceptance_set <- c()
    x_set_str      <- ""
    
    if (fixed_M < m) {
      # The acceptance set is exactly { N - M }
      acceptance_set_num <- N - fixed_M
      
      a_val         <- acceptance_set_num
      b_val         <- acceptance_set_num
      cardinality   <- 1
      coverage_prob <- 0  # same as your code
      gap_val       <- FALSE
      acceptance_set_str <- as.character(acceptance_set_num)
      x_set_str     <- paste(a_val, b_val, sep = "-")
      
    } else {
      # M >= m
      # Subset 'results' to the rows for this M
      sub_df <- results %>% filter(M == fixed_M)
      
      # For each x in [0, N - M], run the original logic
      # "for (fixed_x in 0:(N - fixed_M)) { ... }"
      x_seq <- seq.int(0, N - fixed_M)
      
      for (fixed_x in x_seq) {
        # min_tail_prob for this x
        min_tail_prob_fixed_x <- sub_df$min_tail_prob[sub_df$x == fixed_x]
        
        # All x' with min_tail_prob <= min_tail_prob_fixed_x
        # (the original code calls 'results %>% filter(...)', we do it here in memory)
        possible_xs_df <- sub_df[sub_df$min_tail_prob <= min_tail_prob_fixed_x, ]
        
        # Sum pmf_x for these x'
        prob_sum <- sum(possible_xs_df$pmf_x, na.rm = TRUE)
        
        # If prob_sum > alpha, we add 'fixed_x' to acceptance_set
        if (prob_sum > alpha) {
          acceptance_set <- c(acceptance_set, fixed_x)
        }
      }
      
      acceptance_set <- unique(as.numeric(acceptance_set))  # just in case
      acceptance_set_str <- paste(acceptance_set, collapse = ",")
      
      # Calculate coverage_prob by summing precomputed pmf_x
      # Only for those x in acceptance_set
      coverage_prob <- 0
      if (length(acceptance_set) > 0) {
        # We'll match again on sub_df$x
        coverage_prob <- sum(
          sub_df$pmf_x[sub_df$x %in% acceptance_set],
          na.rm = TRUE
        )
      }
      
      # Build a, b, cardinality
      a_val       <- min(acceptance_set)
      b_val       <- max(acceptance_set)
      cardinality <- length(acceptance_set)
      
      # Check for a "gap"
      gap_val <- any(diff(sort(acceptance_set)) > 1)
      
      # Construct x_set
      if (!gap_val) {
        # No gap => single interval
        x_set_str <- paste(a_val, b_val, sep = "-")
      } else {
        # There is a gap => multiple intervals
        acceptance_set_sorted <- sort(acceptance_set)
        intervals <- c()
        start_int <- acceptance_set_sorted[1]
        
        for (i in seq.int(2, length(acceptance_set_sorted))) {
          if (acceptance_set_sorted[i] != acceptance_set_sorted[i - 1] + 1) {
            intervals <- c(intervals, paste(start_int, acceptance_set_sorted[i - 1], sep = "-"))
            start_int <- acceptance_set_sorted[i]
          }
        }
        # Add the last interval
        intervals <- c(intervals, paste(start_int, acceptance_set_sorted[length(acceptance_set_sorted)], sep = "-"))
        x_set_str <- paste(intervals, collapse = ", ")
      }
    }
    
    # Append to final_results
    final_results <- rbind(
      final_results,
      data.frame(
        M               = fixed_M,
        acceptance_set  = acceptance_set_str,
        a               = a_val,
        b               = b_val,
        cardinality     = cardinality,
        coverage_prob   = coverage_prob,
        x_set           = x_set_str,
        gap             = gap_val
      )
    )
  }
  
  # Finally, arrange in descending order of M
  final_results <- final_results %>%
    arrange(desc(M))
  
  return(final_results)
}


blaker_ci_vec <- function(N, m, conf_level = 0.95) {
  results = blaker_ac_vec(N, m, conf_level)
  
  # Check if there are any gaps
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Loops through each x 
  for (x in 0:N) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the M of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub = first_occurrence$M
      ci_lb = last_occurrence$M
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}



###################################
#---------------------------------#
# CMC                             #
#---------------------------------#
###################################

cmc_ac_vec <- function(N, m, conf_level = 0.95) {
  # Initialize the final results
  results <- data.frame(
    M              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )
  
  # Start with these constraints
  min_a <- 0
  min_b <- 0
  
  # Loop over M from N down to 0
  for (M_val in seq.int(N, 0, by = -1)) {
    # A temporary data frame for the current M_val
    temp_results <- data.frame(
      M              = integer(),
      a              = integer(),
      b              = integer(),
      cardinality    = integer(),
      coverage_prob  = numeric()
    )
    
    # Special End Case: M_val < m => we have exactly (a, b) = (N - M_val, N - M_val)
    if (M_val < m) {
      a_val <- N - M_val
      b_val <- a_val
      
      coverage_val  <- sum_ngh_pmf(N, M_val, m, a_val, b_val)
      cardinality   <- b_val - a_val + 1
      
      temp_results <- rbind(
        temp_results,
        data.frame(
          M              = M_val,
          a              = a_val,
          b              = b_val,
          cardinality    = cardinality,
          coverage_prob  = coverage_val
        )
      )
      
    } else {
      # M_val >= m => create a grid of all valid (a, b) pairs with:
      #   a in [min_a, N - M_val],
      #   b in [min_b, N - M_val],
      #   b >= a (to ensure non-decreasing).
      
      ab_grid <- expand.grid(
        a = seq.int(min_a, N - M_val),
        b = seq.int(min_b, N - M_val)
      ) %>%
        filter(b >= a)
      
      # If there's no valid (a, b), skip
      if (nrow(ab_grid) > 0) {
        # Vectorized coverage probability computation using mapply
        coverage_vec <- mapply(
          FUN = function(a_val, b_val) {
            sum_ngh_pmf(N, M_val, m, a_val, b_val)
          },
          ab_grid$a,
          ab_grid$b
        )
        
        # Build temp_results for all (a, b)
        temp_results <- data.frame(
          M              = M_val,
          a              = ab_grid$a,
          b              = ab_grid$b,
          cardinality    = ab_grid$b - ab_grid$a + 1,
          coverage_prob  = coverage_vec
        )
      }
    }
    
    # Filter out the sets with coverage probability >= conf_level
    # Among those, choose the acceptance curve with the highest a and the lowest b
    # (the same as "the one with the lowest coverage_prob that is still >= conf_level"
    # in your description).
    if (M_val >= m && nrow(temp_results) > 0) {
      temp_results <- temp_results %>%
        filter(coverage_prob >= conf_level) %>%
        group_by(M) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()
      
      # Update min_a and min_b if we found any valid intervals
      if (nrow(temp_results) > 0) {
        min_a <- max(min_a, min(temp_results$a))
        min_b <- max(min_b, min(temp_results$b))
      }
    }
    
    # Append the current iteration’s data to the final results
    results <- rbind(results, temp_results)
  }
  
  # Arrange in descending order of M
  filtered_results <- results %>%
    arrange(desc(M))
  
  # Add a column "x_set" to show "a-b"
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}


cmc_ci_vec <- function(N, m, conf_level = 0.95) {
  results = cmc_ac_vec(N, m, conf_level)
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Loops through each x 
  for (x in 0:N) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the M of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub = first_occurrence$M
      ci_lb = last_occurrence$M
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}



###################################
#---------------------------------#
# N Unknown                       #
#---------------------------------#
###################################

###################################
#---------------------------------#
# Analog to Clopper Pearson       #
#---------------------------------#
###################################

CI_Analog_CP_N_Unknown <- function(M, m, conf_level = 0.95, max_N = 1000) {
  target_probability = (1 - conf_level) / 2
  max_x = max_N - M
  
  results = data.frame(x = 0:max_x, lower_bound = NA, upper_bound = NA)
  previous_upper_bound = 0
  
  for (xi in 0:max_x) {
    
    lower_bound = M + xi
    # Find lower bound
    for (N in (M + xi):(M + max_x)) {
      area_right = ngh_cdf(x = xi - 1, N = N, M = M, m = m, lower_tail = FALSE)
      
      if (isTRUE(all.equal(area_right, target_probability))) {
        lower_bound = N
        break
      } 
      else if (area_right > target_probability) {
        lower_bound = N
        break
      }
    }
    
    upper_bound = lower_bound
    # Find upper bound
    for (N in (lower_bound):(M + max_x)) {  
      area_left = ngh_cdf(x = xi, N = N, M = M, m = m, lower_tail = TRUE)
      
      if (isTRUE(all.equal(area_left, target_probability))) {
        upper_bound = N
        break
      } 
      else if (area_left < target_probability) {
        upper_bound = N - 1  
        break
      }
    }
    
    # Stop the iteration if the upper bound starts decreasing (because of max_N)
    if (upper_bound < previous_upper_bound) {
      break
    }
    
    # Update the previous upper bound for the next iteration
    previous_upper_bound = upper_bound
    
    # Store the results
    results[xi + 1, "lower_bound"] = lower_bound
    results[xi + 1, "upper_bound"] = upper_bound
  }
  
  # Filter out any rows where the upper bound is invalid
  results <- results[!is.na(results$upper_bound), ]
  
  return(results)
}


coverage_prob_ACP_N_unknown <- function(M, N, m, conf_level = 0.95, max_N = 1000) {
  found_N_in_last_CI <- TRUE
  
  while (found_N_in_last_CI) {
    # Calculates all confidence intervals 
    ci_results <- CI_Analog_CP_N_Unknown(M, m, conf_level, max_N)
    
    # Get the confidence interval of the last x
    last_x_ci <- ci_results[nrow(ci_results), ]
    
    # Check if N is within the last x's confidence interval
    if (N >= last_x_ci$lower_bound & N <= last_x_ci$upper_bound) {
      # If N is still in the CI, increase max_N and try again
      max_N <- max_N + 100
    } else {
      # Stop increasing if N is not in the last x's CI anymore
      found_N_in_last_CI <- FALSE
    }
  }
  
  # Once max_N is large enough, continue with the original calculation
  ci_results <- CI_Analog_CP_N_Unknown(M, m, conf_level, max_N)
  
  # Finds all x's where N is in the confidence interval 
  covered_x <- ci_results %>%
    filter(lower_bound <= N & upper_bound >= N) %>%
    pull(x)
  
  if (length(covered_x) == 0) {
    return(data.frame(N = N, coverage_prob = NA, min_x = NA, max_x = NA))
  }
  
  # Finds the min and max of covered x's to know which lines to connect in plot
  min_x <- min(covered_x, na.rm = TRUE)
  max_x <- max(covered_x, na.rm = TRUE)
  
  # Sums the probabilities (pmf's) of all x's where N is in the CI 
  total_prob <- sum(unlist(lapply(covered_x, function(x) ngh_pmf(x, N, M, m))))
  
  return(data.frame(N = N, coverage_prob = total_prob, min_x = min_x, max_x = max_x))
}



###################################
#---------------------------------#
# Minimal Cardinality             #
#---------------------------------#
###################################

all_mc_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Initialize the output data frame
  results <- data.frame(
    N              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )
  
  # Start points for a and b
  min_a <- 0
  min_b <- 0
  
  # We'll not allow a or b to exceed max_x
  max_x <- max_N - M
  
  # Iterate over N
  for (N in M:max_N) {
    
    # Create all possible (a, b) pairs subject to
    #   a >= min_a, b >= min_b, and b >= a, up to max_x
    ab_grid <- expand.grid(
      a = seq.int(min_a, max_x),
      b = seq.int(min_b, max_x)
    ) %>%
      dplyr::filter(b >= a)
    
    # If there's nothing to evaluate, skip
    if (nrow(ab_grid) == 0) {
      next
    }
    
    # Compute coverage probabilities in a vectorized manner
    # (mapply effectively loops internally, but is often faster 
    #  than writing an explicit R loop. If sum_ngh_pmf() 
    #  supports vectorization, you can replace mapply with a 
    #  direct call to sum_ngh_pmf(N, M, m, ab_grid$a, ab_grid$b).)
    coverage <- mapply(
      function(a_val, b_val) {
        sum_ngh_pmf(N, M, m, a_val, b_val)
      },
      ab_grid$a,
      ab_grid$b
    )
    
    # Add cardinality (b - a + 1) and coverage to ab_grid
    temp_results <- data.frame(
      N             = N,
      a             = ab_grid$a,
      b             = ab_grid$b,
      cardinality   = ab_grid$b - ab_grid$a + 1,
      coverage_prob = coverage
    )
    
    # Keep only those with coverage >= conf_level
    # and pick the row(s) with the smallest cardinality.
    temp_results <- temp_results %>%
      dplyr::filter(coverage_prob >= conf_level) %>%
      dplyr::group_by(N) %>%
      dplyr::slice_min(order_by = cardinality, with_ties = TRUE) %>%
      dplyr::ungroup()
    
    # If we found any valid (a, b) intervals, update min_a and min_b
    # and append to results
    if (nrow(temp_results) > 0) {
      min_a <- max(min_a, min(temp_results$a))
      min_b <- max(min_b, min(temp_results$b))
      
      results <- dplyr::bind_rows(results, temp_results)
    }
  }
  
  # Add a column "x_set" that shows "a-b"
  filtered_results <- results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}


minimal_cardinality_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 1000, 
                                            procedure = "MST") {
  # Chooses which minimal cardinality procedure 
  if (procedure == "MST") {
    results = mst_ac_N_unknown_vec(M, m, conf_level, max_N)
  } 
  else if (procedure == "CG") {
    results = cg_ac_N_unknown_vec(M, m, conf_level, max_N)
  } 
  else if (procedure == "BK") {
    results = bk_ac_N_unknown_vec(M, m, conf_level, max_N)
  } 
  else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Find the second-highest "a"
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  if (length(unique_a_values) > 1) {
    max_x <- unique_a_values[2]  # Second-highest value
  } else {
    max_x <- unique_a_values[1]  # Fallback if only one value exists
  }
  
  # Loops through each x 
  for (x in 0:max_x) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the N of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb = first_occurrence$N
      ci_ub = last_occurrence$N
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}


###################################
#---------------------------------#
# MST                             #
#---------------------------------#
###################################

mst_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves 
  results = all_mc_ac_N_unknown_vec(M, m, conf_level, max_N)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(N = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each N 
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current N
    subset_results = results %>% 
      filter(N == current_N)
    
    # If only one acceptance curve, then that is the acceptance curve  
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply MST procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the highest 
    # coverage probability 
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(desc(coverage_prob), desc(a), desc(b)) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(N)
  
  return(final_results)
}


###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac_N_unknown_vec(M, m, conf_level, max_N)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(N = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each N
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current M 
    subset_results = results %>% 
      filter(N == current_N)
    
    # If only one acceptance curve, then that is the acceptance curve 
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply CG procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the smallest 
    # possible a and b
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(a, b) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(N)
  
  return(final_results)
}


###################################
#---------------------------------#
# BK                              #
#---------------------------------#
###################################

bk_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac_N_unknown_vec(M, m, conf_level, max_N)
  
  # Initializes data frame that will be outputted 
  final_results = data.frame(N = integer(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character())
  
  # Initializes min_a and min_b, starting at 0-0
  min_a = 0
  min_b = 0
  
  # Loops through with each N
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current M 
    subset_results = results %>% 
      filter(N == current_N)
    
    # If only one acceptance curve, then that is the acceptance curve 
    if (nrow(subset_results) == 1) {
      chosen_row = subset_results
    } 
    
    # If has more than one option, apply BK procedure
    # Filters so a and b are non-decreasing, and then chooses the row with the largest 
    # possible a and b
    else {
      chosen_row = subset_results %>% 
        filter(a >= min_a, b >= min_b) %>% 
        arrange(desc(a), desc(b)) %>% 
        slice(1)
    }
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    if (nrow(chosen_row) > 0) {
      min_a = max(min_a, chosen_row$a)
      min_b = max(min_b, chosen_row$b)
      final_results = rbind(final_results, chosen_row)
    }
  }
  
  final_results = final_results %>%
    arrange(N)
  
  return(final_results)
}



###################################
#---------------------------------#
# Blaker                          #
#---------------------------------#
###################################

blaker_ac_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  alpha <- 1 - conf_level
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1) FIRST PASS: Build a data frame 'results' with:
  #    (N, x, min_tail_prob, pmf_x)
  #    Instead of looping repeatedly with rbind, we:
  #    - Create combinations for N in [M, max_N] and x in [0, N - M].
  #    - Compute min_tail_prob and pmf_x in a vectorized fashion.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  N_values <- seq.int(M, max_N)
  
  # Build all (N, x) pairs in one go
  combos <- do.call(
    rbind,
    lapply(N_values, function(N_val) {
      data.frame(
        N = N_val,
        x = seq.int(0, N_val - M)
      )
    })
  )
  
  # Compute min_tail_prob for each (N, x) pair
  # and also store pmf_x to avoid repeated calls later.
  # We'll do this in two steps for clarity.
  combos <- combos %>%
    rowwise() %>%
    mutate(
      area_left  = ngh_cdf(x, N, M, m, lower_tail = TRUE),
      area_right = ngh_cdf(x - 1, N, M, m, lower_tail = FALSE),
      min_tail_prob = min(area_left, area_right),
      pmf_x         = ngh_pmf(x, N, M, m)
    ) %>%
    ungroup()
  
  # We'll keep the "acceptance_set" column as NA here (to match your original structure)
  # though it's not used until the second pass.
  results <- combos %>%
    mutate(acceptance_set = NA_character_) %>%
    select(N, x, min_tail_prob, pmf_x, acceptance_set)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2) SECOND PASS: Build 'final_results' by scanning each fixed N
  #    and determining acceptance sets. The logic is unchanged.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  final_results <- data.frame(
    N              = integer(),
    acceptance_set = character(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric(),
    x_set          = character(),
    gap            = logical()
  )
  
  for (fixed_N in seq.int(M, max_N)) {
    # Acceptance set for this N
    acceptance_set <- integer(0)
    
    # Subset results for this N just once
    sub_df <- results %>% filter(N == fixed_N)
    
    # For each x in [0, fixed_N - M], check the condition
    x_seq <- seq.int(0, fixed_N - M)
    for (fixed_x in x_seq) {
      # min_tail_prob for x
      min_tp_fx <- sub_df$min_tail_prob[sub_df$x == fixed_x]
      
      # All x's with min_tail_prob <= min_tail_prob_fixed_x
      possible_xs <- sub_df$x[sub_df$min_tail_prob <= min_tp_fx]
      
      # Sum up the pmf of all these x's
      pmf_sum <- sum(sub_df$pmf_x[sub_df$x %in% possible_xs])
      
      # If sum of pmf is > alpha, include fixed_x in acceptance set
      if (pmf_sum > alpha) {
        acceptance_set <- c(acceptance_set, fixed_x)
      }
    }
    
    # Convert to numeric (unique), build acceptance_set string
    acceptance_set <- sort(unique(as.numeric(acceptance_set)))
    acceptance_set_str <- paste(acceptance_set, collapse = ",")
    
    # a, b, cardinality
    a_val         <- min(acceptance_set)
    b_val         <- max(acceptance_set)
    cardinality   <- length(acceptance_set)
    
    # Coverage probability: sum pmf_x for x in acceptance_set
    coverage_prob <- sum(sub_df$pmf_x[sub_df$x %in% acceptance_set])
    
    # Check if there is a gap
    gap_val <- any(diff(acceptance_set) > 1)
    
    # Build x_set (interval notation)
    x_set_str <- ""
    if (!gap_val) {
      # No gap => single interval
      x_set_str <- paste(a_val, b_val, sep = "-")
    } else {
      # Gap => multiple intervals
      intervals <- c()
      start_int <- acceptance_set[1]
      for (i in seq.int(2, length(acceptance_set))) {
        if (acceptance_set[i] != acceptance_set[i - 1] + 1) {
          intervals <- c(intervals, paste(start_int, acceptance_set[i - 1], sep = "-"))
          start_int <- acceptance_set[i]
        }
      }
      # Add last interval
      intervals <- c(intervals, paste(start_int, acceptance_set[length(acceptance_set)], sep = "-"))
      x_set_str <- paste(intervals, collapse = ", ")
    }
    
    # Add row to final_results
    final_results <- rbind(
      final_results,
      data.frame(
        N               = fixed_N,
        acceptance_set  = acceptance_set_str,
        a               = a_val,
        b               = b_val,
        cardinality     = cardinality,
        coverage_prob   = coverage_prob,
        x_set           = x_set_str,
        gap             = gap_val
      )
    )
  }
  
  final_results <- final_results %>%
    arrange(N)
  
  return(final_results)
}


blaker_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250, procedure = "MST") {
  results = blaker_ac_N_unkown_vec(M, m, conf_level)
  
  # Check if there are any gaps
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Find the second-highest "a"
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  if (length(unique_a_values) > 1) {
    max_x <- unique_a_values[2]  # Second-highest value
  } else {
    max_x <- unique_a_values[1]  # Fallback if only one value exists
  }
  
  # Loops through each x 
  for (x in 0:max_x) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the N of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb = first_occurrence$N
      ci_ub = last_occurrence$N
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}



###################################
#---------------------------------#
# CMC                             #
#---------------------------------#
###################################

cmc_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  # Initialize the final results
  results <- data.frame(
    N              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )
  
  # Start with these constraints
  min_a <- 0
  min_b <- 0
  max_x <- max_N - M
  
  # Loop over N from M up to max_N
  for (N_val in seq.int(M, max_N)) {
    # Build all (a, b) pairs (with b >= a) in a single step
    ab_grid <- expand.grid(
      a = seq.int(min_a, max_x),
      b = seq.int(min_b, max_x)
    ) %>%
      filter(b >= a)
    
    # Compute coverage probability in a vectorized manner
    if (nrow(ab_grid) > 0) {
      coverage_vec <- mapply(
        FUN = function(a_val, b_val) {
          sum_ngh_pmf(N_val, M, m, a_val, b_val)
        },
        ab_grid$a,
        ab_grid$b
      )
      
      temp_results <- data.frame(
        N             = N_val,
        a             = ab_grid$a,
        b             = ab_grid$b,
        cardinality   = ab_grid$b - ab_grid$a + 1,
        coverage_prob = coverage_vec
      )
    } else {
      # If there's no valid (a,b) pair, skip
      temp_results <- data.frame(
        N             = integer(),
        a             = integer(),
        b             = integer(),
        cardinality   = integer(),
        coverage_prob = numeric()
      )
    }
    
    # Filter out sets with coverage_prob >= conf_level
    # Then pick the acceptance curve with the highest 'a' and the lowest 'b'
    temp_results <- temp_results %>%
      filter(coverage_prob >= conf_level) %>%
      group_by(N) %>%
      filter(a == max(a)) %>%
      filter(b == min(b)) %>%
      ungroup()
    
    # Update min_a, min_b if we found any valid intervals
    if (nrow(temp_results) > 0) {
      min_a <- max(min_a, min(temp_results$a))
      min_b <- max(min_b, min(temp_results$b))
    }
    
    # Append to main results
    results <- rbind(results, temp_results)
  }
  
  # Arrange by N in ascending order
  filtered_results <- results %>%
    arrange(N)
  
  # Add x_set column "a-b"
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}


cmc_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250, procedure = "MST") {
  results = cmc_ac_N_unknown_vec(M, m, conf_level, max_N)
  
  # Initializes data frame that will be outputted 
  ci_results = data.frame(x = integer(), 
                          ci_lb = integer(), 
                          ci_ub = integer(), 
                          ci = character(),
                          stringsAsFactors = FALSE)
  
  # Find the second-highest "a"
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  if (length(unique_a_values) > 1) {
    max_x <- unique_a_values[2]  # Second-highest value
  } else {
    max_x <- unique_a_values[1]  # Fallback if only one value exists
  }
  
  # Loops through each x 
  for (x in 0:max_x) {
    # Finds first interval where x appears 
    first_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(1)
    
    # Finds last interval where x appears 
    last_occurrence = results %>% 
      filter(a <= x, x <= b) %>% 
      slice(n())
    
    # Finds the N of the corresponding above intervals and saves those as the CI bounds
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb = first_occurrence$N
      ci_ub = last_occurrence$N
      ci = paste0("[", ci_lb, ", ", ci_ub, "]")
      
      ci_results = rbind(ci_results, data.frame(x = x, 
                                                ci_lb = ci_lb, 
                                                ci_ub = ci_ub, 
                                                ci = ci))
    }
  }
  
  return(ci_results)
}





