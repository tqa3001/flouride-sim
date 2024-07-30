# Far TODOs
# - cost/DALY, DALY/1000, DALY diff
# - Maybe add more threads

library(tidyverse)
setwd(this.path::here())
source("helper.R")

# Input parameters (plenty embedded in sample functions too, consider a refactor)
current_coverage <- "urban_only"
current_scenario <- "baseline"
n_times <- 3
n_cycles <- 15
death_age <- 110
rate_with_symptom <- 0.324
disability_weight <- 0.057
discount_rate <- 0.03

# Enums
v_coverage <- c("urban_only", "full_coverage")
v_scenarios <- c("baseline", "reduction_all_ages", "double_mean_cost", "include_xray")
v_pop_types <- c("rural", "urban")
v_states <- c("healthy", "aged_care", "caries", "edentulism", "dead")
v_outputs <- c("prevalence", "daly", "cost")
n_states <- length(v_states)

# Util (double check correctness)
sample_gamma <- function(mean, sd) {
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  rgamma(1, shape = shape, scale = scale)
}

sample_beta <- function(mean, sd) {
  alpha <- ((1 - mean)/sd^2 - 1/mean) * mean^2
  beta <- alpha * (1/mean - 1)
  rbeta(1, shape1 = alpha, shape2 = beta)
}

sample_normal <- function(mean, sd) {
  rnorm(1, mean = mean, sd = sd)
}

# ref DARTH paper
get_wcc_vector <- function(len) {
  is_even <- function(x) x %% 2 == 0
  is_odd <- function(x) x %% 2 != 0
  v_cycles <- seq(1, len + 1)
  v_wcc <- is_even(v_cycles)*(2/3) + is_odd(v_cycles)*(4/3)
  v_wcc[1] <- v_wcc[len + 1] <- 1/3
  return(v_wcc)
}

# before national oral health plan 69%, with health plan 89%,
# full coverage 100% pop
get_population <- function(population_type, is_indigenous) {
  get_random_cohort_population <- function(age) {
    sample(1:10000, n_states, replace = TRUE)
  }
  result <- matrix(NA, death_age, n_states, dimnames=list(1:death_age,v_states))
  for (age in 1:nrow(init_states)) {
    result[age,] <- get_cohort_population(age)
  }
  return(result)
}

# Caries
caries_inc_rate_by_age <- function(age, rates) {
  case_when(
    age >= 0 & age <= 14 ~ rates[1],
    age >= 15 & age <= 64 ~ rates[2],
    age >= 65 ~ rates[3],
    TRUE ~ NA_real_
  )
}

sample_caries_rates <- function() {
  c(
    sample_normal(0.26, 0.16),
    sample_normal(0.44, 0.01),
    sample_normal(0.10, 0.02)
  )
}

get_caries_symptom_duration <- function(age) {
  case_when(
    age >= 0 & age <= 16 ~ 28 / 365,
    age >= 17 ~ 55 / 365
  )
}

sample_residential_aged_care_rate <- function() {
  c(
    sample_beta(0.013, 0.0001),
    sample_beta(0.16, 0.0004)
  )
}

get_residential_aged_care_rate <- function(age, rates) {
  case_when(
    age < 60 ~ 0,
    age >= 60 & age <= 79 ~ rates[1],
    age >= 80 ~ rates[2]
  )
}

# Cohorts
# 5 years age, sex, 2003 start

# Edentulism
# Used in original paper to adjust caries rate
# Here we use as a prob transition
# How to "include a downward trend in edentulism to the year 2041"?

edentulism_rate_by_age <- function(age, rates) {
  case_when(
    age < 15 ~ 0, # double check
    age >= 15 & age <= 34 ~ rates[1],
    age >= 35 & age <= 54 ~ rates[2],
    age >= 55 & age <= 74 ~ rates[3],
    age >= 75 ~ 0.016,
    TRUE ~ NA_real_
  )
}

sample_edentulism_rates <- function() {
  c(
    0, # sample_beta(0, 0.0003), # beta mean 0?
    sample_beta(0.017, 0.0023),
    sample_beta(0.14, 0.0064),
    sample_beta(0.36, 0.016)
  )
}

apply_risk_ratio <- function(rate, rr) { rate * rr }

indigenous_RR_by_age <- function(age) {
  case_when(
    age >= 1 & age <= 9 ~ 2.4,
    age >= 10 & age <= 14 ~ 1.8,
    age >= 15 & age <= 34 ~ 2.7,
    age >= 35 & age <= 54 ~ 2.2,
    age >= 55 ~ 2.0,
    TRUE ~ NA_real_
  )
}

remote_RR_by_age <- function(age) {
  case_when(
    age >= 1 & age <= 9 ~ 0.9,
    age >= 10 & age <= 14 ~ 1.1,
    age >= 15 & age <= 24 ~ 0.3,
    age >= 25 & age <= 44 ~ 0.5,
    age >= 45 ~ 0.5,
    TRUE ~ NA_real_
  )
}

# Cost sampling
sample_costs <- function(scenario) {
  multiplier <- if (scenario == "double_mean_cost") 2 else 1
  c(
    flouridation_urban = sample_gamma(0.26 * multiplier, 0.05),
    flouridation_rural = sample_gamma(26 * multiplier, 5),
    oral_exam = sample_gamma(45, 13),
    metallic_restoration = sample_gamma(109, 21),
    opg_exposure = sample_gamma(82, 15)
  )
}

get_one_time_cost <- function(costs, population_type, population_count, scenario) {
  cost_per_person <- case_when(
    population_type == "urban" ~ costs[["flouridation_urban"]],
    population_type == "rural" ~ costs[["flouridation_rural"]],
    TRUE ~ NA_real_
  )
  return(cost_per_person * population_count)
}

get_treatment_cost <- function(costs, scenario) {
  cost <- costs[["oral_exam"]] + costs[["metallic_restoration"]]
  if (scenario == "include_xray") {
    cost <- cost + costs[[opg_exposure]]
  }
  return(cost)
}

# Mortality: Perhaps follow a distribution, higher rate when older?
get_mortality_rate <- function(age) {
  return(0.005 * age)
}

# Generate transitions for each cycle/age group
generate_transitions_once <- function(coverage, scenario, population) {
  
  risk_difference_caries_free_prop_change <- sample_normal(0.154, 0.0237)
  caries_rates <- sample_caries_rates()
  edentulism_rates <- sample_edentulism_rates()
  aged_care_rates <- sample_residential_aged_care_rate()
  root_caries_rate_aged <- sample_normal(0.29, 0.69)
  coronal_caries_rate_aged <- sample_normal(0.71, 1.09)
  
  default_P <- matrix(0, nrow = n_states, ncol = n_states,
                      dimnames = list(v_states, v_states))
  
  get_transitions_for_age <- function(age) {
    
    P <- default_P
    p_D <- get_mortality_rate(age) %>% rate_to_prob
    r_caries <- caries_inc_rate_by_age(age, caries_rates)
    r_edent <- edentulism_rate_by_age(age, edentulism_rates)
    r_caries_aged_care <- root_caries_rate_aged + coronal_caries_rate_aged
    rr_healthy_caries <- 1
    
    if (population$population_type == "rural") {
      rr_healthy_caries <- rr_healthy_caries * remote_RR_by_age(age)
      if (population$is_indigenous) {
        rr_healthy_caries <- rr_healthy_caries * indigenous_RR_by_age(age)
      }
    }
    
    is_fluoridated <- (
      (coverage == "urban_only" && population$population_type == "urban") ||
      (coverage == "full_coverage")
    )
    is_effective <- (
      (current_scenario != "reduction_all_ages" && age < 16) || 
      (current_scenario == "reduction_all_ages")
    )
    # should this be subtract or multiply?
    if (is_fluoridated && is_effective) {
        rr_healthy_caries <- rr_healthy_caries * 
          risk_difference_caries_free_prop_change
    }
    
    P["healthy", "edentulism"] <- r_edent %>% rate_to_prob() * (1 - p_D)
    P["healthy", "aged_care"] <- get_residential_aged_care_rate(age, aged_care_rates) %>% 
      rate_to_prob() * (1 - p_D)
    P["healthy", "caries"] <- r_caries %>% 
      apply_risk_ratio(rr_healthy_caries) %>%
      rate_to_prob() * (1 - p_D)
    P["healthy", "healthy"] <- (1 - sum(P["healthy", ])) * (1 - p_D)
    
    P["caries", "edentulism"] <- r_edent %>%
      rate_to_prob() * (1 - p_D)
    P["caries", "caries"] <- (1 - sum(P["caries", ])) * (1 - p_D)
    
    P["aged_care", "caries"] <- r_caries_aged_care %>% rate_to_prob() * (1 - p_D)
    P["aged_care", "aged_care"] <- (1 - sum(P["aged_care", ])) * (1 - p_D)
    
    P["healthy", "dead"] <- p_D
    P["caries", "dead"] <- p_D
    P["aged_care", "dead"] <- p_D
    P["edentulism", "dead"] <- p_D
    # add later: edent to edent
    P["dead", "dead"] <- 1
    
    return(P)
  }
  
  transitions_all_ages <- vapply(1:death_age, get_transitions_for_age, default_P)
  
  return(transitions_all_ages)
}

simulate_once <- function(coverage, scenario, population) {
  
  traces <- vector("list", length = death_age)
  P <- generate_transitions_once(coverage, scenario, population)
  costs <- sample_costs(scenario)
  
  # Each age forms a cohort
  for (start_age in 1:(death_age-1)) {
    M <- matrix(nrow = death_age, ncol = n_states)
    colnames(M) <- v_states
    M[start_age,] <- population$init_states[start_age,]
    for (age in start_age:min(death_age - 1, start_age + n_cycles)) {
      M[age + 1,] <- M[age,] %*% P[,,age]
    }
    traces[[start_age]] <- M
  }
  
  labels <- append(v_outputs, c("cycle"))
  result <- data.frame(matrix(ncol = length(labels), nrow = n_cycles))
  colnames(result) <- labels
  wcc_coefficients <- get_wcc_vector(n_cycles)
  
  total_pop <- sum(population$init_states)
  cumulative_cost <- get_one_time_cost(
    costs, population$population_type, total_pop, current_scenario)
  cumulative_daly <- 0
  
  for (i in 1:n_cycles) {
    
    current_prevalence <- 0
    new_cases <- 0
    
    for (start_age in 1:(death_age - i)) {
      trace <- traces[[start_age]]
      
    # Calculate prevalence at the i-th cycle
      current_prevalence <- current_prevalence + trace[[start_age + i, "caries"]]
    
    # Calculate cost in AUD after i cycles
      new_cases_in_cohort <- trace[[start_age + i - 1, "healthy"]] * 
        P["healthy", "caries", age]
      new_cases <- new_cases + new_cases_in_cohort
      cumulative_cost <- cumulative_cost + 
        get_treatment_cost(costs, current_scenario) * 
          new_cases * (1 + discount_rate)^i * wcc_coefficients[i]
      
    # DALY
      caries_cohort_pop <- trace[[start_age + i, "caries"]]
      symptomatic_caries_cohort_pop <- caries_cohort_pop * rate_with_symptom
      duration <- get_caries_symptom_duration(start_age + i)
      daly_for_cohort <- symptomatic_caries_cohort_pop * disability_weight * duration
      cumulative_daly <- daly_for_cohort * wcc_coefficients[i]
    }
    
    data <- c(
      cycle = i, 
      prevalence = as.numeric(current_prevalence),
      daly = cumulative_daly,
      cost = cumulative_cost
    )
    
    result[i, names(data)] <- data
  }
  
  print(population$population_type)
  print(population$is_indigenous)
  print(result)
  return(result)
}

simulate_monte_carlo <- function(coverage, scenario) {
  
  result <- vector("list", n_times)
  
  for (cycle in 1:n_times) {
    # Calculate results for different populations
    for (population_type in v_pop_types) {
      for (is_indigenous in c(T, F)) {
        population <- list(
          init_states = get_population(population_type, is_indigenous),
          population_type = population_type,
          is_indigenous = is_indigenous
        )
        result_for_population <- simulate_once(coverage, scenario, population)
        if (is.null(result[[cycle]])) {
          result[[cycle]] <- result_for_population
        }
        else {
          result[[cycle]] <- result[[cycle]] + result_for_population
        }
      }
    }
  }
  
  combined_df <- bind_rows(result, .id = "replicate")
  result <- combined_df %>%
    group_by(cycle) %>%
    summarise(
      prevalence_mean = mean(prevalence),
      daly_mean = mean(daly),
      cost_mean = mean(cost)
    )
  
  return(result)
}


result <- simulate_monte_carlo(current_coverage, current_scenario)

# Visualization
print("final result table:")
result

# Visualizing prevalence over time
# Baseline year on the paper is 2003
baseline_year = 2003
result |>
  ggplot(mapping = aes(x = baseline_year + cycle, y = prevalence_mean)) +
  geom_line() +
  labs(title = "Prevalence over the years",
       x = "year",
       y = "mean prevalence")

# - DALY over time (before after and 100%)

# - Cost over time (before after and 100%)

# Further visualisation
# - tornado stuff ?
  