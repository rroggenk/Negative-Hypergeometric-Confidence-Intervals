###################################
#---------------------------------#
# Preliminary Functions           #
#---------------------------------#
###################################

ngh_pmf <- function(x, N, M, m) {
  # Conditions: check supports
  if (x > (N - M) || N < M || m > M) {
    return(0)
  }
  
  numerator <- choose(m + x - 1, m - 1) * choose(N - m - x, M - m)
  denominator <- choose(N, M)
  result <- numerator / denominator
  
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

all_mc_ac <- function(N, m, conf_level = 0.95) {
  # Initializing the data frame that the function will output 
  results = data.frame(M = integer(), 
                       a = integer(), 
                       b = integer(), 
                       cardinality = integer(), 
                       coverage_prob = numeric())
  
  # Setting the initial min_a and min_b, always starts with 0-0
  min_a = 0
  min_b = 0
  
  for (M in N:0) {
    # Initializing a data frame that temporarily stores results
    temp_results = data.frame(M = integer(), 
                              a = integer(), 
                              b = integer(), 
                              cardinality = integer(), 
                              coverage_prob = numeric())
    
    # Find all possible acceptance curves that have non-decreasing a and b
    
    # Speical End Cases: When M < m
    if (M < m) {
      a = N - M
      b = N - M
      coverage_prob = sum_ngh_pmf(N, M, m, a, b)
      cardinality = b - a + 1
      temp_results = rbind(temp_results, data.frame(M = M, 
                                                    a = a, 
                                                    b = b, 
                                                    cardinality = cardinality, 
                                                    coverage_prob = coverage_prob))
    } 
    else {
      # Loops through the a first, making sure it only starts at min_a so that a is 
      # non-decreasing, stops at a N-M
      for (a in min_a:(N - M)) {
        # Loops through b: starting at the max of a and min_b to make sure that b >= a and 
        # so that b is non-decreasing, stops at N-M
        for (b in max(a, min_b):(N - M)) {
          # Calculates coverage probability, cardinality, and stores it in data frame
          coverage_prob = sum_ngh_pmf(N, M, m, a, b)
          cardinality = b - a + 1
          temp_results = rbind(temp_results, data.frame(M = M, 
                                                        a = a, 
                                                        b = b, 
                                                        cardinality = cardinality, 
                                                        coverage_prob = coverage_prob))
        }
      }
    }
    
    # Filter out the sets with the smallest cardinality and coverage probability >= conf_level 
    if (M >= m) {
      temp_results = temp_results %>%
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
        group_by(M) %>%
        slice_min(order_by = cardinality, with_ties = TRUE) %>%
        ungroup()
      
      # Updates min_a and min_b for each iteration 
      if (nrow(temp_results) > 0) {
        min_a = max(min_a, min(temp_results$a))
        min_b = max(min_b, min(temp_results$b))
      }
    }
    
    results = rbind(results, temp_results)
  }
  
  # Separate the results into two data frames
  results_M_ge_m = results %>% 
    filter(M >= m)
  
  results_M_lt_m = results %>% 
    filter(M < m)
  
  # Combine the two data frames
  filtered_results = bind_rows(results_M_ge_m, results_M_lt_m) %>%
    arrange(desc(M))
  
  # Adds a column of the x set
  filtered_results = filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}


minimal_cardinality_ci <- function(N, m, conf_level = 0.95, procedure = "MST") {
  # Chooses which minimal cardinality procedure 
  if (procedure == "MST") {
    results = mst_ac(N, m, conf_level)
  } 
  else if (procedure == "CG") {
    results = cg_ac(N, m, conf_level)
  } 
  else if (procedure == "BK") {
    results = bk_ac(N, m, conf_level)
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

mst_ac <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves 
  results = all_mc_ac(N, m, conf_level)
  
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

cg_ac <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac(N, m, conf_level)
  
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

bk_ac <- function(N, m, conf_level = 0.95) {
  # Gets all minimal cardinality acceptance curves
  results = all_mc_ac(N, m, conf_level)
  
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

blaker_ac <- function(N, m, conf_level = 0.95) {
  alpha = 1 - conf_level 
  
  # Initializing data frame to store info through iterations
  results = data.frame(M = integer(), 
                       x = integer(), 
                       min_tail_prob = numeric(), 
                       acceptance_set = character())
  
  # Calculates min_tail_prob for each M and x
  # Iterating through each M
  for (M in 0:N) {
    # Special case if M < m, x bounds are known so so min_tail_prob is needed 
    if (M < m) {
      x = N - M
      results = rbind(results, data.frame(M = M, 
                                          x = x,
                                          min_tail_prob = NA,
                                          acceptance_set = as.character(x)))
    } 
    else {
      # Iterating through each x 
      for (x in 0:(N - M)) {
        # Calculates min_tail_prob
        area_left = ngh_cdf(x = x, N = N, M = M, m = m, lower_tail = TRUE)
        area_right = ngh_cdf(x = x - 1, N = N, M = M, m = m, lower_tail = FALSE)
        min_tail_prob = min(area_left, area_right)
        
        # Stores min_tail_prob in a dataframe along with coresponding M and x
        results = rbind(results, data.frame(M = M, 
                                            x = x,
                                            min_tail_prob = min_tail_prob,
                                            acceptance_set = NA))
      }
    }
  }
  
  # Initializes data frame for output 
  final_results = data.frame(M = integer(), 
                             acceptance_set = character(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character(), 
                             gap = logical())
  
  # Calculate the acceptance sets based on min_tail_prob
  # Iterating through each M 
  for (fixed_M in 0:N) {
    acceptance_set = c()
    x_set = ""
    
    # Special Case: if M < m, acceptance set is already known 
    if (fixed_M < m) {
      # Directly include x = N - M for M < m
      acceptance_set = as.character(N - fixed_M)
      acceptance_set_num = as.numeric(acceptance_set)
      
      a = min(acceptance_set_num)
      b = max(acceptance_set_num)
      cardinality = length(acceptance_set_num)
      coverage_prob = 0
      gap = FALSE
      acceptance_set_str = acceptance_set
      x_set = paste(a, b, sep = "-")
    } 
    
    # When M >= m
    else {
      # Iterating through each x 
      for (fixed_x in 0:(N - fixed_M)) {
        # Finds the min_tail_prob for the fixed_M and fixed_x
        min_tail_prob_fixed_x = results %>% 
          filter(M == fixed_M & x == fixed_x) %>% 
          select(min_tail_prob) %>% 
          pull()
        
        # Finds all x's with min_tail_prob as small as the fixed_x (includes fixed_x)
        possible_xs = results %>% 
          #filter(M == fixed_M & x != fixed_x & min_tail_prob <= min_tail_prob_fixed_x) %>% 
          # if i use above line, would need a condition when fixed_M == N
          filter(M == fixed_M & min_tail_prob <= min_tail_prob_fixed_x) %>% 
          select(x) %>% 
          pull()
        
        # Sums up the pmf of all x's with min-tail prob as small as fixed_x
        prob_sum = sum(unlist(lapply(possible_xs, function(px) ngh_pmf(x = px, 
                                                                       N = N, 
                                                                       M = fixed_M, 
                                                                       m = m))))
        
        # If the sum of pmf is greater than alpha, add the fixed_x to the acceptance set 
        if (prob_sum > alpha) {
          acceptance_set = c(acceptance_set, fixed_x)
        }
      }
      
      acceptance_set = as.numeric(acceptance_set)
      acceptance_set_str = paste(unique(acceptance_set), collapse = ",")
      
      # Calculating cardinality and coverage probability 
      a = min(acceptance_set)
      b = max(acceptance_set)
      cardinality = length(acceptance_set)
      coverage_prob = sum(unlist(lapply(acceptance_set, function(px) ngh_pmf(x = px, 
                                                                             N = N, 
                                                                             M = fixed_M, 
                                                                             m = m))))
      # gap is a boolean (TRUE / FALSE) of whether there is a gap present in acceptance set
      gap = any(diff(sort(acceptance_set)) > 1)
      
      # If there is no gap, formatting x_set
      if (!gap) {
        x_set = paste(a, b, sep = "-")
      } 
      # When there is a gap, formatting x_set 
      else {
        intervals = c()
        start = acceptance_set[1]
        for (i in 2:length(acceptance_set)) {
          if (acceptance_set[i] != acceptance_set[i - 1] + 1) {
            intervals = c(intervals, paste(start, acceptance_set[i - 1], sep = "-"))
            start = acceptance_set[i]
          }
        }
        intervals = c(intervals, 
                      paste(start, acceptance_set[length(acceptance_set)], sep = "-"))
        x_set = paste(intervals, collapse = ", ")
      }
    }
    
    # Adding final_results to the data frame to be outputted 
    final_results = rbind(final_results, 
                          data.frame(M = fixed_M, 
                                     acceptance_set = acceptance_set_str, 
                                     a = a, 
                                     b = b, 
                                     cardinality = cardinality, 
                                     coverage_prob = coverage_prob, 
                                     x_set = x_set, 
                                     gap = gap))
  }
  final_results = final_results %>%
    arrange(desc(M))
  
  return(final_results)
}


blaker_ci <- function(N, m, conf_level = 0.95) {
  results = blaker_ac(N, m, conf_level)
  
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

cmc_ac <- function(N, m, conf_level = 0.95) {
  # Initializing the data frame that the function will output 
  results = data.frame(M = integer(), 
                       a = integer(), 
                       b = integer(), 
                       cardinality = integer(), 
                       coverage_prob = numeric())
  
  # Setting the initial min_a and min_b, always starts with 0-0
  min_a = 0
  min_b = 0
  
  for (M in N:0) {
    # Initializing a data frame that temporarily stores results
    temp_results = data.frame(M = integer(), 
                              a = integer(), 
                              b = integer(), 
                              cardinality = integer(), 
                              coverage_prob = numeric())
    
    # Find all possible acceptance curves that have non-decreasing a and b
    
    # Speical End Cases: When M < m
    if (M < m) {
      a = N - M
      b = N - M
      coverage_prob = sum_ngh_pmf(N, M, m, a, b)
      cardinality = b - a + 1
      temp_results = rbind(temp_results, data.frame(M = M, 
                                                    a = a, 
                                                    b = b, 
                                                    cardinality = cardinality, 
                                                    coverage_prob = coverage_prob))
    } 
    else {
      # Loops through the a first, making sure it only starts at min_a so that a is 
      # non-decreasing, stops at a N-M
      for (a in min_a:(N - M)) {
        # Loops through b: starting at the max of a and min_b to make sure that b >= a and 
        # so that b is non-decreasing, stops at N-M
        for (b in max(a, min_b):(N - M)) {
          # Calculates coverage probability, cardinality, and stores it in data frame
          coverage_prob = sum_ngh_pmf(N, M, m, a, b)
          cardinality = b - a + 1
          temp_results = rbind(temp_results, data.frame(M = M, 
                                                        a = a, 
                                                        b = b, 
                                                        cardinality = cardinality, 
                                                        coverage_prob = coverage_prob))
        }
      }
    }
    
    
    # Filter out the sets with the coverage probability >= conf_level
    # For each M, filters out to choose the acceptance curve with the highest a and the 
    # lowest b (which is the same as the one with the lowest coverage prob that is still
    # above the confidence level)
    if (M >= m) {
      temp_results = temp_results %>%
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
        group_by(M) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()
      
      # Updates min_a and min_b for each iteration
      if (nrow(temp_results) > 0) {
        min_a = max(min_a, min(temp_results$a))
        min_b = max(min_b, min(temp_results$b))
      }
    }
    
    results = rbind(results, temp_results)
  }
  
  # Arranges data by descreasing M
  filtered_results = results %>%
    arrange(desc(M))
  
  # Adds a column of the x set
  filtered_results = filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
  
}


cmc_ci <- function(N, m, conf_level = 0.95) {
  results = cmc_ac(N, m, conf_level)
  
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

minimal_cardinality_ci_N_unkown <- function(M, m, conf_level = 0.95, max_N = 1000, 
                                            procedure = "MST") {
  # Chooses which minimal cardinality procedure 
  if (procedure == "MST") {
    results = mst_ac_N_unknown_direct(M, m, conf_level, max_N)
  } 
  else if (procedure == "CG") {
    results = cg_ac_N_unknown_direct(M, m, conf_level, max_N)
  } 
  else if (procedure == "BK") {
    results = bk_ac_N_unknown_direct(M, m, conf_level, max_N)
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

mst_ac_N_unknown_direct <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # This function loops through N from M to max_N. For each N it searches 
  # for a viable acceptance curve (a, b) having the smallest possible cardinality.
  # a and b are constrained to be non-decreasing across N.
  # The candidate (a, b) is chosen by scanning over (a, b) pairs with fixed
  # cardinality (b - a + 1) and then selecting the one with highest coverage probability,
  # then highest a and b if ties exist.
  
  results <- list()
  
  # Initialize state: a and b from the previous acceptance curve.
  min_a <- 0
  min_b <- 0
  previous_cardinality <- NULL  # for N = M, we start with candidate_card = 1
  
  for (current_N in M:max_N) {
    max_x <- current_N - M  # a and b are chosen from 0:max_x
    # (Normally, max_x should be >= 0 because current_N >= M.)
    
    # For N == M, the only possibility is a = b = 0 (cardinality 1).
    candidate_card <- if (!is.null(previous_cardinality)) previous_cardinality else 1
    candidate_found <- FALSE
    best_candidate <- NULL
    
    # Loop over candidate cardinalities (the maximum possible is max_x+1)
    while (candidate_card <= (max_x + 1) && !candidate_found) {
      # a is chosen from a_min to a_max such that b = a + candidate_card - 1 <= max_x.
      # Also, we require a >= min_a and b >= min_b.
      a_min <- max(min_a, ifelse(candidate_card > 1, min_b - candidate_card + 1, 0))
      a_max <- max_x - (candidate_card - 1)
      
      # If the lower bound for a is greater than the upper bound, no candidate exists for this cardinality.
      if (a_min > a_max) {
        candidate_card <- candidate_card + 1
        next
      }
      
      # Initialize temporary storage for the best candidate for this candidate_card.
      best_cov <- -Inf
      best_a <- NA
      best_b <- NA
      
      # Loop through all possible a values
      for (a_val in a_min:a_max) {
        b_val <- a_val + candidate_card - 1
        if (b_val < min_b) next  # ensure b does not drop below previous min_b
        
        # Calculate the coverage probability for this (a, b)
        cov_prob <- sum_ngh_pmf(current_N, M, m, a_val, b_val)
        
        # Check if the (a, b) meets the coverage requirements.
        if (cov_prob >= conf_level && cov_prob <= 1 && cov_prob >= 0) {
          # Among candidates of the same cardinality, we select the one with:
          # 1. Maximum coverage probability,
          # 2. then maximum a, and
          # 3. then maximum b.
          if ((cov_prob > best_cov) ||
              (cov_prob == best_cov && a_val > best_a) ||
              (cov_prob == best_cov && a_val == best_a && b_val > best_b)) {
            best_cov <- cov_prob
            best_a <- a_val
            best_b <- b_val
          }
        }
      }  # end loop over a
      
      # If we found at least one candidate for this cardinality, record the best one.
      if (!is.na(best_a)) {
        best_candidate <- list(
          N = current_N,
          a = best_a,
          b = best_b,
          cardinality = candidate_card,
          coverage_prob = best_cov
        )
        candidate_found <- TRUE
      } else {
        # No candidate found at this cardinality; increase the cardinality and try again.
        candidate_card <- candidate_card + 1
      }
    }  # end while over candidate_card
    
    # If a candidate was found for this N, update the state and add the candidate to the results.
    if (!is.null(best_candidate)) {
      results[[length(results) + 1]] <- best_candidate
      min_a <- best_candidate$a
      min_b <- best_candidate$b
      previous_cardinality <- best_candidate$cardinality
    }
    # (If no candidate is found, this N is simply skipped.)
  }  # end loop over N
  
  # If no acceptance curves were found at all, return an empty data.frame with the expected columns.
  if (length(results) == 0) {
    return(data.frame(
      N = integer(),
      a = integer(),
      b = integer(),
      cardinality = integer(),
      coverage_prob = numeric(),
      x_set = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Convert the list of accepted curves to a data.table and add the x_set column.
  dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
  dt[, x_set := paste(a, b, sep = "-")]
  setcolorder(dt, c("N", "a", "b", "cardinality", "coverage_prob", "x_set"))
  
  # Return as a data.frame.
  return(as.data.frame(dt))
}


###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_N_unknown_direct <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # This function finds minimal‐cardinality acceptance curves following the CG procedure.
  # For each N from M to max_N, it searches for an (a, b) (with b = a + k - 1, where k is the candidate cardinality)
  # such that the coverage probability (computed by sum_ngh_pmf) meets the conf_level.
  # Among all candidates at the minimal viable cardinality, it chooses the one with the largest a, then largest b.
  # In addition, a and b are not allowed to decrease across N.
  
  results <- list()
  
  # Initialize state variables
  min_a <- 0
  min_b <- 0
  previous_cardinality <- NA  # for N=M, we start with candidate_card = 1
  
  for (current_N in M:max_N) {
    max_x <- current_N - M  # possible a, b range is 0:max_x
    candidate_card <- if (!is.na(previous_cardinality)) previous_cardinality else 1
    candidate_found <- FALSE
    best_candidate <- NULL
    
    # Loop over candidate cardinalities until a viable acceptance curve is found.
    while (candidate_card <= (max_x + 1) && !candidate_found) {
      # For candidate cardinality k, ensure a >= min_a and also a >= min_b - (k - 1) so that b >= min_b.
      a_min <- max(min_a, min_b - candidate_card + 1)
      a_max <- max_x - (candidate_card - 1)
      
      if (a_min > a_max) {
        candidate_card <- candidate_card + 1
        next
      }
      
      # For CG, choose the candidate with the largest a; if tied, choose the one with largest b.
      for (a_val in a_min:a_max) {
        b_val <- a_val + candidate_card - 1
        if (b_val < min_b) next
        
        cov_prob <- sum_ngh_pmf(current_N, M, m, a_val, b_val)
        if (cov_prob >= conf_level && cov_prob <= 1 && cov_prob >= 0) {
          if (is.null(best_candidate)) {
            best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                   cardinality = candidate_card, coverage_prob = cov_prob)
          } else {
            # For CG, prefer larger a; if a is equal then prefer larger b.
            if (a_val > best_candidate$a ||
                (a_val == best_candidate$a && b_val > best_candidate$b)) {
              best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                     cardinality = candidate_card, coverage_prob = cov_prob)
            }
          }
        }
      }
      
      if (!is.null(best_candidate)) {
        candidate_found <- TRUE
      } else {
        candidate_card <- candidate_card + 1
      }
    }  # end while over candidate_card
    
    if (!is.null(best_candidate)) {
      results[[length(results) + 1]] <- best_candidate
      min_a <- best_candidate$a
      min_b <- best_candidate$b
      previous_cardinality <- candidate_card
    }
  }  # end loop over N
  
  if (length(results) == 0) {
    return(data.frame(
      N = integer(),
      a = integer(),
      b = integer(),
      cardinality = integer(),
      coverage_prob = numeric(),
      x_set = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
  dt[, x_set := paste(a, b, sep = "-")]
  setcolorder(dt, c("N", "a", "b", "cardinality", "coverage_prob", "x_set"))
  
  return(as.data.frame(dt))
}


###################################
#---------------------------------#
# BK                              #
#---------------------------------#
###################################

bk_ac_N_unknown_direct <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # This function finds minimal‐cardinality acceptance curves following the BK procedure.
  # For each N from M to max_N, it searches for an (a, b) (with b = a + k - 1, where k is the candidate cardinality)
  # such that the coverage probability (computed by sum_ngh_pmf) meets the conf_level.
  # Among all candidates at the minimal viable cardinality, it chooses the one with the smallest a, then smallest b.
  # Additionally, a and b are not allowed to decrease across N.
  
  results <- list()
  
  # Initialize state variables
  min_a <- 0
  min_b <- 0
  previous_cardinality <- NA  # for N=M, we start with candidate_card = 1
  
  for (current_N in M:max_N) {
    max_x <- current_N - M  # possible values for a and b lie in 0:max_x
    candidate_card <- if (!is.na(previous_cardinality)) previous_cardinality else 1
    candidate_found <- FALSE
    best_candidate <- NULL
    
    # Loop over candidate cardinalities until a viable acceptance curve is found.
    while (candidate_card <= (max_x + 1) && !candidate_found) {
      # For a candidate cardinality k, we must have a >= min_a and also
      # b = a + k - 1 >= min_b, so a >= min_b - (k - 1).
      a_min <- max(min_a, min_b - candidate_card + 1)
      a_max <- max_x - (candidate_card - 1)
      
      if (a_min > a_max) {
        candidate_card <- candidate_card + 1
        next
      }
      
      # For BK, we select the candidate with the smallest a then smallest b.
      for (a_val in a_min:a_max) {
        b_val <- a_val + candidate_card - 1
        if (b_val < min_b) next  # enforce non-decreasing b
        
        cov_prob <- sum_ngh_pmf(current_N, M, m, a_val, b_val)
        if (cov_prob >= conf_level && cov_prob <= 1 && cov_prob >= 0) {
          # Viable candidate found at this cardinality.
          if (is.null(best_candidate)) {
            best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                   cardinality = candidate_card, coverage_prob = cov_prob)
          } else {
            # Choose the candidate with smaller a; if tied, choose smaller b.
            if (a_val < best_candidate$a ||
                (a_val == best_candidate$a && b_val < best_candidate$b)) {
              best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                     cardinality = candidate_card, coverage_prob = cov_prob)
            }
          }
        }
      }
      
      if (!is.null(best_candidate)) {
        candidate_found <- TRUE
      } else {
        candidate_card <- candidate_card + 1
      }
    }  # end while over candidate_card
    
    # If a viable candidate was found for this N, update state and record the candidate.
    if (!is.null(best_candidate)) {
      results[[length(results) + 1]] <- best_candidate
      min_a <- best_candidate$a
      min_b <- best_candidate$b
      previous_cardinality <- candidate_card
    }
    # (If no candidate is found for current_N, that N is skipped.)
  }  # end loop over N
  
  # If no curves were found, return an empty data.frame with the expected columns.
  if (length(results) == 0) {
    return(data.frame(
      N = integer(),
      a = integer(),
      b = integer(),
      cardinality = integer(),
      coverage_prob = numeric(),
      x_set = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
  dt[, x_set := paste(a, b, sep = "-")]
  setcolorder(dt, c("N", "a", "b", "cardinality", "coverage_prob", "x_set"))
  
  return(as.data.frame(dt))
}



###################################
#---------------------------------#
# Blaker                          #
#---------------------------------#
###################################

# When finding all x's with min tail prob as small, including the fixed_x 
blaker_ac_N_unkown <- function(M, m, conf_level = 0.95, max_N = 250) {
  alpha = 1 - conf_level 
  
  # Initializing data frame to store info through iterations
  results = data.frame(N = integer(), 
                       x = integer(), 
                       min_tail_prob = numeric(), 
                       acceptance_set = character())
  
  # Calculates min_tail_prob for each N and x
  # Iterating through each N
  for (N in M:max_N) {
    # Iterating through each x 
    for (x in 0:(N - M)) {
      # Calculates min_tail_prob
      area_left = ngh_cdf(x = x, N = N, M = M, m = m, lower_tail = TRUE)
      area_right = ngh_cdf(x = x - 1, N = N, M = M, m = m, lower_tail = FALSE)
      min_tail_prob = min(area_left, area_right)
      
      # Stores min_tail_prob in a dataframe along with coresponding M and x
      results = rbind(results, data.frame(N = N, 
                                          x = x,
                                          min_tail_prob = min_tail_prob,
                                          acceptance_set = NA))
    }
  }
  
  # Initializes data frame for output 
  final_results = data.frame(N = integer(), 
                             acceptance_set = character(), 
                             a = integer(), 
                             b = integer(), 
                             cardinality = integer(), 
                             coverage_prob = numeric(), 
                             x_set = character(), 
                             gap = logical())
  
  # Calculate the acceptance sets based on min_tail_prob
  # Iterating through each N
  for (fixed_N in M:max_N) {
    acceptance_set = c()
    x_set = ""
    
    # Iterating through each x 
    for (fixed_x in 0:(fixed_N - M)) {
      # Finds the min_tail_prob for the fixed_N and fixed_x
      min_tail_prob_fixed_x = results %>% 
        filter(N == fixed_N & x == fixed_x) %>% 
        select(min_tail_prob) %>% 
        pull()
      
      # Finds all x's with min_tail_prob as small as the fixed_x (includes fixed_x)
      possible_xs = results %>% 
        filter(N == fixed_N & min_tail_prob <= min_tail_prob_fixed_x) %>% 
        select(x) %>% 
        pull()
      
      # Sums up the pmf of all x's with min-tail prob as small as fixed_x
      prob_sum = sum(unlist(lapply(possible_xs, function(px) ngh_pmf(x = px, 
                                                                     N = fixed_N, 
                                                                     M = M, 
                                                                     m = m))))
      
      # If the sum of pmf is greater than alpha, add the fixed_x to the acceptance set 
      if (prob_sum > alpha) {
        acceptance_set = c(acceptance_set, fixed_x)
      }
    }
    
    acceptance_set = as.numeric(acceptance_set)
    acceptance_set_str = paste(unique(acceptance_set), collapse = ",")
    
    # Calculating cardinality and coverage probability 
    a = min(acceptance_set)
    b = max(acceptance_set)
    cardinality = length(acceptance_set)
    coverage_prob = sum(unlist(lapply(acceptance_set, function(px) ngh_pmf(x = px, 
                                                                           N = fixed_N, 
                                                                           M = M, 
                                                                           m = m))))
    # gap is a boolean (TRUE / FALSE) of whether there is a gap present in acceptance set
    gap = any(diff(sort(acceptance_set)) > 1)
    
    # If there is no gap, formatting x_set
    if (!gap) {
      x_set = paste(a, b, sep = "-")
    } 
    # When there is a gap, formatting x_set 
    else {
      intervals = c()
      start = acceptance_set[1]
      for (i in 2:length(acceptance_set)) {
        if (acceptance_set[i] != acceptance_set[i - 1] + 1) {
          intervals = c(intervals, paste(start, acceptance_set[i - 1], sep = "-"))
          start = acceptance_set[i]
        }
      }
      intervals = c(intervals, 
                    paste(start, acceptance_set[length(acceptance_set)], sep = "-"))
      x_set = paste(intervals, collapse = ", ")
    }
    
    # Adding final_results to the data frame to be outputted 
    final_results = rbind(final_results, 
                          data.frame(N = fixed_N, 
                                     acceptance_set = acceptance_set_str, 
                                     a = a, 
                                     b = b, 
                                     cardinality = cardinality, 
                                     coverage_prob = coverage_prob, 
                                     x_set = x_set, 
                                     gap = gap))
  }
  final_results = final_results %>%
    arrange(N)
  
  return(final_results)
}


blaker_ci_N_unkown <- function(M, m, conf_level = 0.95, max_N = 250) {
  results = blaker_ac_N_unkown(M, m, conf_level, max_N)
  
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

cmc_ac_N_unknown <- function(M, m, conf_level = 0.95, max_N = 250) {
  # Initializing the data frame that the function will output 
  results = data.frame(N = integer(), 
                       a = integer(), 
                       b = integer(), 
                       cardinality = integer(), 
                       coverage_prob = numeric())
  
  # Setting the initial min_a and min_b, always starts with 0-0
  min_a = 0
  min_b = 0
  
  for (N in M:max_N) {
    max_x = N - M
    # cat("N: ", N, "\n")
    
    # Initializing a data frame that temporarily stores results
    temp_results = data.frame(N = integer(), 
                              a = integer(), 
                              b = integer(), 
                              cardinality = integer(), 
                              coverage_prob = numeric())
    
    # Loop through possible 'a' values (ensuring non-decreasing a)
    for (a in min_a:min(min_a + 5, max_x)) {
      # Loop through b's (starting at max(a, min_b) to ensure b >= a and non-decreasing b)
      for (b in max(a, min_b):max_x) {
        coverage_prob = sum_ngh_pmf(N, M, m, a, b)
        cardinality = b - a + 1
        temp_results = rbind(temp_results, data.frame(N = N, 
                                                      a = a, 
                                                      b = b, 
                                                      cardinality = cardinality, 
                                                      coverage_prob = coverage_prob))
        # Added break condition: once coverage meets or exceeds the conf_level, stop further b's
        if (coverage_prob >= conf_level) {
          break
        }
      }
    }
    
    # Filter out the sets with the coverage probability >= conf_level
    temp_results <- temp_results %>%
      filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1)
    
    # If any valid rows exist, choose the acceptance curve with max a and min b
    if (nrow(temp_results) > 0) {
      temp_results <- temp_results %>%
        group_by(N) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()
    }
    
    # Update min_a and min_b for the next iteration
    if (nrow(temp_results) > 0) {
      min_a = max(min_a, min(temp_results$a))
      min_b = max(min_b, min(temp_results$b))
    }
    results = rbind(results, temp_results)
  }
  
  # Arrange results by N and add x_set column
  filtered_results = results %>%
    arrange(N) %>%
    mutate(x_set = paste(a, b, sep = "-"))
  
  return(filtered_results)
}



cmc_ci_N_unkown <- function(M, m, conf_level = 0.95, max_N = 250) {
  results = cmc_ac_N_unknown(M, m, conf_level, max_N)
  
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





