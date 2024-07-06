# TODO:
# - Add half-cycle correction
# - Fill in sampling functions
#
# - Add modelling for the 15% thingy (follow paper)
# - Make DALY correct
# - Add visualisation
# BIG:
# - modify functions to handle multiple scenarios

# Far TODOs
# - cost/DALY, DALY/1000, DALY diff
# - Maybe add more threads

library(tidyverse)
setwd(this.path::here())
source("helper.R")

# Input parameters (plenty embedded in sample functions too, consider a refactor)
scenario <- "standard"
n_times <- 3
n_cycles <- 15
death_age <- 110
discount_rate <- 0.03
cost_per_person_urban <- 0.26
cost_per_person_rural <- 26
prop_with_symptom <- 0.324
disability_weight <- 0.057

# Enums
v_states <- c("healthy", "caries", "edentulism", "dead")
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

# analysing
# different scenarios: 
# - before national oral health plan 69%, with health plan 89%, full coverage 100% pop
# - baseline: only <16yo gets 15% reduction in proportion of caries-free
# - extend: >= 16yo receive same reduction
# - extend: x2 urban estimate cost from 0.26 AUD to 0.52 AUD per person
# todo: can we use enum or something smh

# v_coverage <- c("pre_health_plan", "urban", "full")
# v_effective_group <- c("children", "all")
# v_cost_type <- 

get_cohort_population <- function(age) {
  sample(1:10000, 4, replace = TRUE) # put real population here
}

init_states <- matrix(NA, death_age, n_states, 
                      dimnames=list(1:death_age,v_states))

for (age in 1:nrow(init_states)) {
  init_states[age,] = get_cohort_population(age)
}

print("initial populations:")
print(init_states)

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
    sample_beta(0, 0.0003),
    sample_beta(0.017, 0.0023),
    sample_beta(0.14, 0.0064),
    sample_beta(0.36, 0.016)
  )
}

# Risk factors
# - Caries rates: non-Indigenous population, the nonremote Indigenous population, and remote Indigenous population 
# Why: to reflect difference in communities <1000 and >=1000
# How exactly? 

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
sample_treatment_costs <- function() {
  c(
    oral_exam = sample_gamma(45, 13),
    metallic_restoration = sample_gamma(109, 21),
    opg_exposure = sample_gamma(82, 15)
  )
}

get_one_time_cost <- function(total_population, urban_population) {
  cost_urban <- cost_per_person_urban * total_population
  cost_all <- cost_per_person_urban * urban_population + 
    cost_per_person_rural * (total_population - urban_population)
  
  # todo: add case switch base on scenario
  
  return(cost_urban)
}

# Mortality: Perhaps follow a distribution, higher rate when older?
get_mortality_rate <- function(age) {
  return(0.005 * age)
}

# Generate transitions for each cycle/age group
# Missing: RR for caries -> eden, RRs based on scenario
generate_transitions_once <- function() {
  
  caries_rates <- sample_caries_rates()
  edentulism_rates <- sample_edentulism_rates()
  default_P <- matrix(0, nrow = n_states, ncol = n_states,
                      dimnames = list(v_states, v_states))
  
  get_transitions_for_age <- function(age) {
    P <- default_P
    p_D <- get_mortality_rate(age) %>% rate_to_prob
    r_edent <- edentulism_rate_by_age(age, edentulism_rates)
    
    P["healthy", "edentulism"] <- r_edent %>% rate_to_prob() * (1 - p_D)
    P["healthy", "caries"] <- caries_inc_rate_by_age(age, caries_rates) %>% 
      apply_risk_ratio(1) %>%
      rate_to_prob() * (1 - p_D)
    P["healthy", "healthy"] <- (1 - sum(P["healthy", ])) * (1 - p_D)
    
    P["caries", "edentulism"] <- r_edent %>% 
      apply_risk_ratio(1) %>%
      rate_to_prob() * (1 - p_D)
    P["caries", "caries"] <- (1 - P["caries", "edentulism"]) * (1 - p_D)
    
    P["healthy", "dead"] <- p_D
    P["caries", "dead"] <- p_D
    P["edentulism", "dead"] <- p_D
    P["dead", "dead"] <- 1
    
    return(P)
  }
  
  transitions_all_ages <- vapply(1:death_age, get_transitions_for_age, default_P)
  
  return(transitions_all_ages)
}

simulate_once <- function() {
  
  traces = vector("list", length = death_age)
  
  P <- generate_transitions_once()
  costs <- sample_treatment_costs()
  
  # split like paper? cohort = each 5 years
  for (start_age in 1:(death_age-1)) {
    M <- matrix(nrow = death_age, ncol = n_states)
    colnames(M) <- v_states
    M[start_age,] <- init_states[start_age,]
    for (age in start_age:min(death_age - 1, start_age + n_cycles)) {
      M[age + 1,] <- M[age,] %*% P[,,age]
    }
    traces[[start_age]] <- M
  }
  
  labels <- append(v_outputs, c("cycle"))
  result <- data.frame(matrix(ncol = length(labels), nrow = n_cycles))
  colnames(result) <- labels
  
  # Variable for keeping track of the total cost
  # Using random values for args
  cumulative_cost <- get_one_time_cost(10000, 8000)
  cumulative_daly <- 0
  
  for (i in 1:n_cycles) {
    
    current_prevalence <- 0
    new_cases <- 0
    
    for (start_age in 1:(death_age - i)) {
      trace <- traces[[start_age]]
      
    # Calculate prevalence at the i-th cycle
      current_prevalence <- current_prevalence + trace[start_age + i, "caries"]
    
    # Calculate cost in AUD after i cycles
      new_cases_in_cohort <- trace[start_age + i - 1, "healthy"] * 
        P["healthy", "caries", age]
      new_cases <- new_cases + new_cases_in_cohort
      
    # Calculate DALY after i cycles
    # DALY
    # 25% of dental caries patients required an emergency visit gets 1h pain per day for 1.5 years
    # We also have time spent with caries symptoms from clinic data for children and adolescents (34) and adults (35)
    # Disability weight = 0.057 (0 = perfect, 1 = death) FOR SYMPTOMATIC CARIES 
      cumulative_daly <- cumulative_daly + 
        disability_weight * (trace[start_age + i, "caries"])
    }
    
    cumulative_cost <- cumulative_cost + 
      (costs["oral_exam"] + costs["metallic_restoration"]) * 
        new_cases * (1 + discount_rate)^i
    
    data <- c(
      cycle = i, 
      prevalence = as.numeric(current_prevalence),
      daly = cumulative_daly,
      cost = cumulative_cost
    )
    
    result[i, names(data)] <- data
  }
  
  return(result)
}

simulate_monte_carlo <- function() {
  
  result <- vector("list", n_times)
  
  for (i in 1:n_times) {
    result[[i]] <- simulate_once()
  }
  
  print(result[[1]])
  print("dsfsdfsld")
  
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


result <- simulate_monte_carlo()

# Visualization
result

# Visualizing prevalence over time
# Baseline year on the paper is 2003
baseline_year = 2003
result |>
  ggplot(mapping = aes(x = baseline_year + cycle, y = prevalence_mean)) +
  geom_line() +
  labs(title = "Prevalence over the years")

# - DALY over time (before after and 100%)

# - Cost over time (before after and 100%)

# Further visualisation
# - tornado stuff ?
  