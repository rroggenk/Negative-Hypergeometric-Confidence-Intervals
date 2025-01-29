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
    if (x >= 0 && x <= (N - M)) {
      sum_pmf = sum_pmf + ngh_pmf(x, N, M, m)
    }
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


CI_cov_prob_MLE_vec <- function(N, m, conf_level = 0.95) {
  
  # Create a vector of all possible x in [0, N].
  x_seq <- 0:N
  
  # Pre-allocate vectors for lower and upper bounds.
  lower_bounds <- numeric(length(x_seq))
  upper_bounds <- numeric(length(x_seq))
  
  # Precompute z outside the loop (same for all x).
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # Loop once over x_seq, filling in bounds.
  for (i in seq_along(x_seq)) {
    xi <- x_seq[i]
    
    # Special End Cases
    if (xi >= (N - m) && xi <= N) {
      # If x is in [N - m, N], set lb = ub = N - x
      lower_bounds[i] <- N - xi
      upper_bounds[i] <- N - xi
    } else {
      M_MLE_val <- M_MLE_function(m, xi, N)
      SE_MLE_val <- SE_MLE_function(M_MLE_val, N, m)
      
      lower_bounds[i] <- M_MLE_val - (z * SE_MLE_val)
      upper_bounds[i] <- M_MLE_val + (z * SE_MLE_val)
    }
  }
  
  # Build and return the result data frame
  results <- data.frame(
    x            = x_seq,
    lower_bound  = lower_bounds,
    upper_bound  = upper_bounds
  )
  
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


CI_cov_prob_unbiased_vec <- function(N, m, conf_level = 0.95) {
  
  # Create a vector of x in [0, N].
  x_seq <- 0:N
  
  # Pre-allocate vectors for lower and upper bounds.
  lower_bounds <- numeric(length(x_seq))
  upper_bounds <- numeric(length(x_seq))
  
  # Precompute z once.
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # Loop once over x_seq
  for (i in seq_along(x_seq)) {
    xi <- x_seq[i]
    
    # Special End Cases
    if (xi >= (N - m) && xi <= N) {
      lower_bounds[i] <- N - xi
      upper_bounds[i] <- N - xi
    } else {
      M_unb_val <- M_unbiased_function(m, xi, N)
      SE_unb_val <- SE_unbiased_function(M_unb_val, N, m)
      
      lower_bounds[i] <- M_unb_val - (z * SE_unb_val)
      upper_bounds[i] <- M_unb_val + (z * SE_unb_val)
    }
  }
  
  # Build the final data frame
  results <- data.frame(
    x            = x_seq,
    lower_bound  = lower_bounds,
    upper_bound  = upper_bounds
  )
  
  return(results)
}



###################################
#---------------------------------#
# Analog to Clopper Pearson       #
#---------------------------------#
###################################

CI_cov_prob_vec <- function(N, m, conf_level = 0.95) {
  target_probability <- (1 - conf_level) / 2
  
  # Pre-allocate storage
  x_values      <- seq.int(0, N)
  lower_bounds  <- numeric(length(x_values))
  upper_bounds  <- numeric(length(x_values))
  
  # Loop once over all x
  for (i in seq_along(x_values)) {
    xi <- x_values[i]
    
    if (xi >= (N - m) && xi <= N) {
      # Special case
      lower_bounds[i] <- N - xi
      upper_bounds[i] <- N - xi
    } else {
      # Defaults
      lb <- m
      ub <- N
      
      # Find lower bound
      for (M_val in seq.int(m, N)) {
        area_left <- ngh_cdf(x = xi, N = N, M = M_val, m = m, lower_tail = TRUE)
        if (isTRUE(all.equal(area_left, target_probability)) ||
            (area_left > target_probability)) {
          lb <- M_val
          break
        }
      }
      
      # Find upper bound
      for (M_val in seq.int(N, m, by = -1)) {
        area_right <- ngh_cdf(x = xi - 1, N = N, M = M_val, m = m, lower_tail = FALSE)
        if (isTRUE(all.equal(area_right, target_probability)) ||
            (area_right > target_probability)) {
          ub <- M_val
          break
        }
      }
      
      lower_bounds[i] <- lb
      upper_bounds[i] <- ub
    }
  }
  
  # Build the final data frame
  results <- data.frame(
    x           = x_values,
    lower_bound = lower_bounds,
    upper_bound = upper_bounds
  )
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
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
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
  # Choose which minimal cardinality procedure
  if (procedure == "MST") {
    results <- mst_ac_vec(N, m, conf_level)
  } else if (procedure == "CG") {
    results <- cg_ac_vec(N, m, conf_level)
  } else if (procedure == "BK") {
    results <- bk_ac_vec(N, m, conf_level)
  } else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }
  
  # We'll store rows in a list, then bind them once at the end
  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all non-NULL rows
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}




###################################
#---------------------------------#
# MST                             #
#---------------------------------#
###################################

mst_ac_vec <- function(N, m, conf_level = 0.95) {
  results <- all_mc_ac_vec(N, m, conf_level)
  
  # We'll collect rows in a list
  row_list <- vector("list", length = N + 1)
  idx <- 1
  
  min_a <- 0
  min_b <- 0
  
  for (current_M in seq.int(N, 0, by = -1)) {
    subset_results <- results %>% filter(M == current_M)
    
    if (nrow(subset_results) == 1) {
      chosen_row <- subset_results
    } else {
      chosen_row <- subset_results %>%
        filter(a >= min_a, b >= min_b) %>%
        arrange(desc(coverage_prob), desc(a), desc(b)) %>%
        slice(1)
    }
    
    if (nrow(chosen_row) > 0) {
      min_a <- max(min_a, chosen_row$a)
      min_b <- max(min_b, chosen_row$b)
      
      row_list[[idx]] <- chosen_row
      idx <- idx + 1
    }
  }
  
  # Combine chosen rows
  final_results <- do.call(rbind, row_list[1:(idx - 1)])
  final_results <- final_results %>% arrange(desc(M))
  
  return(final_results)
}




###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_vec <- function(N, m, conf_level = 0.95) {
  results <- all_mc_ac_vec(N, m, conf_level)
  
  row_list <- vector("list", length = N + 1)
  idx <- 1
  
  min_a <- 0
  min_b <- 0
  
  for (current_M in seq.int(N, 0, by = -1)) {
    subset_results <- results %>% filter(M == current_M)
    
    if (nrow(subset_results) == 1) {
      chosen_row <- subset_results
    } else {
      chosen_row <- subset_results %>%
        filter(a >= min_a, b >= min_b) %>%
        arrange(a, b) %>%
        slice(1)
    }
    
    if (nrow(chosen_row) > 0) {
      min_a <- max(min_a, chosen_row$a)
      min_b <- max(min_b, chosen_row$b)
      
      row_list[[idx]] <- chosen_row
      idx <- idx + 1
    }
  }
  
  final_results <- do.call(rbind, row_list[1:(idx - 1)])
  final_results <- final_results %>% arrange(desc(M))
  return(final_results)
}



###################################
#---------------------------------#
# BK                              #
#---------------------------------#
###################################

bk_ac_vec <- function(N, m, conf_level = 0.95) {
  results <- all_mc_ac_vec(N, m, conf_level)
  
  row_list <- vector("list", length = N + 1)
  idx <- 1
  
  min_a <- 0
  min_b <- 0
  
  for (current_M in seq.int(N, 0, by = -1)) {
    subset_results <- results %>% filter(M == current_M)
    
    if (nrow(subset_results) == 1) {
      chosen_row <- subset_results
    } else {
      chosen_row <- subset_results %>%
        filter(a >= min_a, b >= min_b) %>%
        arrange(desc(a), desc(b)) %>%
        slice(1)
    }
    
    if (nrow(chosen_row) > 0) {
      min_a <- max(min_a, chosen_row$a)
      min_b <- max(min_b, chosen_row$b)
      
      row_list[[idx]] <- chosen_row
      idx <- idx + 1
    }
  }
  
  final_results <- do.call(rbind, row_list[1:(idx - 1)])
  final_results <- final_results %>% arrange(desc(M))
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
  results <- blaker_ac_vec(N, m, conf_level)
  
  # If there's any gap, just return the message
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }
  
  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
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
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
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
  results <- cmc_ac_vec(N, m, conf_level)
  
  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
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

CI_Analog_CP_N_Unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  target_probability <- (1 - conf_level) / 2
  max_x <- max_N - M
  
  x_values <- seq.int(0, max_x)
  lower_bounds <- rep(NA, length(x_values))
  upper_bounds <- rep(NA, length(x_values))
  
  previous_upper_bound <- 0
  
  for (i in seq_along(x_values)) {
    xi <- x_values[i]
    
    # Initialize lower_bound
    lb <- M + xi
    # Find lower bound
    for (N_val in seq.int(M + xi, M + max_x)) {
      area_right <- ngh_cdf(x = xi - 1, N = N_val, M = M, m = m, lower_tail = FALSE)
      
      if (isTRUE(all.equal(area_right, target_probability)) || (area_right > target_probability)) {
        lb <- N_val
        break
      }
    }
    
    # Initialize upper_bound
    ub <- lb
    # Find upper bound
    for (N_val in seq.int(lb, M + max_x)) {
      area_left <- ngh_cdf(x = xi, N = N_val, M = M, m = m, lower_tail = TRUE)
      
      if (isTRUE(all.equal(area_left, target_probability))) {
        ub <- N_val
        break
      } else if (area_left < target_probability) {
        ub <- N_val - 1
        break
      }
    }
    
    # Stop if the upper bound starts decreasing
    if (ub < previous_upper_bound) {
      # Mark the rest as NA and break
      break
    }
    previous_upper_bound <- ub
    
    lower_bounds[i] <- lb
    upper_bounds[i] <- ub
  }
  
  # Build final results
  results <- data.frame(
    x = x_values,
    lower_bound = lower_bounds,
    upper_bound = upper_bounds
  )
  
  # Filter out rows where upper_bound is NA
  results <- results[!is.na(results$upper_bound), ]
  return(results)
}



coverage_prob_ACP_N_unknown_vec <- function(M, N, m, conf_level = 0.95, max_N = 1000) {
  found_N_in_last_CI <- TRUE
  
  # Keep increasing max_N until N is no longer in the last CI
  while (found_N_in_last_CI) {
    ci_results <- CI_Analog_CP_N_Unknown_vec(M, m, conf_level, max_N)
    last_x_ci <- ci_results[nrow(ci_results), ]
    
    if (N >= last_x_ci$lower_bound && N <= last_x_ci$upper_bound) {
      max_N <- max_N + 100
    } else {
      found_N_in_last_CI <- FALSE
    }
  }
  
  # Then compute final coverage probability
  ci_results <- CI_Analog_CP_N_Unknown_vec(M, m, conf_level, max_N)
  
  covered_x <- ci_results %>%
    dplyr::filter(lower_bound <= N & upper_bound >= N) %>%
    dplyr::pull(x)
  
  if (length(covered_x) == 0) {
    return(data.frame(N = N, coverage_prob = NA, min_x = NA, max_x = NA))
  }
  
  min_x <- min(covered_x)
  max_x <- max(covered_x)
  
  total_prob <- sum(
    sapply(covered_x, function(x_val) ngh_pmf(x_val, N, M, m))
  )
  
  return(data.frame(N = N, coverage_prob = total_prob, min_x = min_x, max_x = max_x))
}




###################################
#---------------------------------#
# Minimal Cardinality             #
#---------------------------------#
###################################

all_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # We'll collect final data for each N in a list
  all_results_list <- vector("list", length = (max_N - M + 1))
  
  # Index to store in our list
  idx <- 1
  
  # Loop from N = M to N = max_N
  for (N_val in seq.int(M, max_N)) {
    # The largest possible b = N_val - M
    max_x <- N_val - M
    
    # Build all (a, b) pairs with 0 <= a <= max_x, 0 <= b <= max_x, and b >= a
    ab_grid <- expand.grid(
      a = seq.int(0, max_x),
      b = seq.int(0, max_x)
    ) %>%
      filter(b >= a)
    
    # Calculate coverage_prob and cardinality
    coverage_vec <- mapply(
      FUN = function(a_val, b_val) {
        sum_ngh_pmf(N_val, M, m, a_val, b_val)
      },
      ab_grid$a,
      ab_grid$b
    )
    
    cardinality_vec <- ab_grid$b - ab_grid$a + 1
    
    # Combine columns into a temporary data frame
    temp_results <- data.frame(
      N             = N_val,
      a             = ab_grid$a,
      b             = ab_grid$b,
      cardinality   = cardinality_vec,
      coverage_prob = coverage_vec
    )
    
    # Filter out sets based on coverage
    temp_results <- temp_results %>%
      filter(
        coverage_prob >= conf_level,
        coverage_prob <= 1,
        coverage_prob >= 0,
        a <= b
      )
    
    # Store results if we got any valid rows
    if (nrow(temp_results) > 0) {
      all_results_list[[idx]] <- temp_results
      idx <- idx + 1
    }
  }
  
  # Combine everything into one data frame
  # (Only up to idx - 1, since we might not have used the entire list)
  if (idx > 1) {
    combined_results <- do.call(rbind, all_results_list[1:(idx - 1)])
  } else {
    # If no results found at all, return an empty data frame with the same structure
    combined_results <- data.frame(
      N = integer(),
      a = integer(),
      b = integer(),
      cardinality = integer(),
      coverage_prob = numeric()
    )
  }
  
  # Add a "x_set" column
  filtered_results <- combined_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}



minimal_cardinality_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 1000, 
                                                procedure = "MST") {
  # Choose the procedure
  if (procedure == "MST") {
    results <- mst_ac_N_unknown(M, m, conf_level, max_N)
  } else if (procedure == "CG") {
    results <- cg_ac_N_unknown(M, m, conf_level, max_N)
  } else if (procedure == "BK") {
    results <- bk_ac_N_unknown(M, m, conf_level, max_N)
  } else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }
  
  # Determine max_x from the second-highest 'a'
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  max_x <- if (length(unique_a_values) > 1) unique_a_values[2] else unique_a_values[1]
  
  x_values <- seq.int(0, max_x)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb <- first_occurrence$N
      ci_ub <- last_occurrence$N
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}



###################################
#---------------------------------#
# MST                             #
#---------------------------------#
###################################

mst_ac_N_unknown <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves 
  results = all_ac_N_unknown_vec(M, m, conf_level, max_N)
  
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
  previous_cardinality = 0
  
  # Loops through with each N 
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current N
    subset_results = results %>% 
      filter(N == current_N)
    
    # print(subset_results)
    # print(min(subset_results$cardinality))
    
    chosen_row <- subset_results %>%
      # filter(cardinality >= previous_cardinality) %>%
      filter(a >= min_a, b >= min_b) %>% 
      filter(cardinality == min(cardinality)) %>%  
      arrange(desc(coverage_prob), desc(a), desc(b)) %>% 
      slice(1)
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    min_a = chosen_row$a
    min_b = chosen_row$b
    previous_cardinality = chosen_row$cardinality
    final_results = rbind(final_results, chosen_row)
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

cg_ac_N_unknown <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves 
  results = all_ac_N_unknown_vec(M, m, conf_level, max_N)
  
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
  previous_cardinality = 0
  
  # Loops through with each N 
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current N
    subset_results = results %>% 
      filter(N == current_N)
    
    # print(subset_results)
    # print(min(subset_results$cardinality))
    
    chosen_row <- subset_results %>%
      # filter(cardinality >= previous_cardinality) %>%
      filter(a >= min_a, b >= min_b) %>% 
      filter(cardinality == min(cardinality)) %>%  
      arrange(a, b) %>%
      slice(1)
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    min_a = chosen_row$a
    min_b = chosen_row$b
    previous_cardinality = chosen_row$cardinality
    final_results = rbind(final_results, chosen_row)
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

bk_ac_N_unknown <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # Gets all minimal cardinality acceptance curves 
  results = all_ac_N_unknown_vec(M, m, conf_level, max_N)
  
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
  previous_cardinality = 0
  
  # Loops through with each N 
  for (current_N in M:max_N) {
    # Only looks at acceptance curves for the current N
    subset_results = results %>% 
      filter(N == current_N)
    
    # print(subset_results)
    # print(min(subset_results$cardinality))
    
    chosen_row <- subset_results %>%
      # filter(cardinality >= previous_cardinality) %>%
      filter(a >= min_a, b >= min_b) %>% 
      filter(cardinality == min(cardinality)) %>%  
      arrange(desc(a), desc(b)) %>% 
      slice(1)
    
    # Updates min_a and min_b and then adds row to final outputted data frame 
    min_a = chosen_row$a
    min_b = chosen_row$b
    previous_cardinality = chosen_row$cardinality
    final_results = rbind(final_results, chosen_row)
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


blaker_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  results <- blaker_ac_N_unkown_vec(M, m, conf_level, max_N)
  
  # If there's any gap, just return
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }
  
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  max_x <- if (length(unique_a_values) > 1) unique_a_values[2] else unique_a_values[1]
  
  x_values <- seq.int(0, max_x)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb <- first_occurrence$N
      ci_ub <- last_occurrence$N
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
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
  
  # Loop over N from M up to max_N
  for (N_val in seq.int(M, max_N)) {
    max_x <- N_val - M
    # Build all (a, b) pairs (with b >= a) in a single step
    ab_grid <- expand.grid(
      a = seq.int(min_a, max_x),
      # a = seq.int(min_a, min_a+1),
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
      filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1)
    
    # If no rows left, skip the group_by part
    if (nrow(temp_results) > 0) {
      temp_results <- temp_results %>%
        group_by(N) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()
    }
    
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



cmc_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  results <- cmc_ac_N_unknown_vec(M, m, conf_level, max_N)
  
  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  max_x <- if (length(unique_a_values) > 1) unique_a_values[2] else unique_a_values[1]
  
  x_values <- seq.int(0, max_x)
  row_list <- vector("list", length(x_values))
  
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    
    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())
    
    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb <- first_occurrence$N
      ci_ub <- last_occurrence$N
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")
      
      row_list[[i]] <- data.frame(
        x = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}






