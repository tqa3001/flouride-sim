# File Loading ------------------------------------------------------------

load_paths <- \(.) list.files(path = "data", pattern = "data", full.names = TRUE)

# Loads the files "low_ses_data, high_ses_data and whole_data and cleans
# them up into a single tibble
load_data <- function(paths, death_age)
{
  # Obsession with pipes
  paths |>
    set_names(basename) |>
    map(read_excel) |>
    list_rbind(names_to = "ses") |>
    mutate(ses = sub("^([^_]+).*", "\\1", ses)) |>
    mutate(across(c(ses, sex), factor)) |>
    clean_names() |>
    select(!c(caries_prevalence, caries_incidence)) |>
    rename(
      caries_incidence = caries_diff_inc, 
      caries_prevalence = caries_diff_prev,
      mortality = background_mortality
    ) |>
    mutate(max_age = max(age)) |>
    group_by(sex, ses) |>
    slice(c(1:n(), rep(n(), death_age - max(age)))) |>
    (\(x) mutate(x, age = -1 + min(age) + 1:(dim(x)[1]/6)))() |>
    mutate(starting_age = age <= max_age) |>
    select(!max_age) |>
    ungroup()
}


# Cohorts -----------------------------------------------------------------

generate_cohorts <- function()
{
  expand.grid(
    c("male", "female"), 
    c("high", "low", "overall"), 
    c("soc", "sug", "ssb")
  ) |>
    tibble() |>
    rename(sex = Var1, ses = Var2, trt = Var3) |>
    # Add MI but only for low ses
    rbind(tibble(
      sex = c("male", "female"),
      ses = c("low", "low"),
      trt = c("mi", "mi")
    ))
}

# Data Transformation -----------------------------------------------------

rate_to_prob <- function(x) 1 - exp(-x)

risk_ratio <- function(x, RR) RR * x

convert_time <- function(x, start_interval, end_interval) 
{
  ratio <- end_interval/start_interval
  1 - (1 - x)^(ratio)
}

filter_cohort_data <- function(data, cohort)
{
  data |>
    filter(sex == cohort$sex, ses == cohort$ses)
}


# Variables ---------------------------------------------------------------

# Returns list(model_variables, copula_factory, prob_details) 
set_variables <- function(col_length, variables)
{
  general_tibble <- create_variables_tibble(variables)
  
  model_variables <- general_tibble |>
    filter(!is_column) |>
    select(name, copula_name, random_function)
  
  prob_details <- general_tibble |>
    filter(is_column) |>
    select(name, copula_name, random_function)
  
  copula_factory <- general_tibble |>
    select(copula_name, is_column) |>
    unique() |>
    mutate(length = ifelse(is_column, col_length, 1)) |>
    select(!is_column)
  
  list(model_variables, copula_factory, prob_details)
}

create_variables_tibble <- function(variables)
{
  tibble(
    name = map_chr(variables, ~.x[[1]]),
    copula_name = map_chr(variables, ~.x[[2]]),
    random_function = map(variables, ~.x[[3]]),
    is_column = map_lgl(variables, ~.x[[4]])
  )
}

generate_copula <- function(copula)
{
  copula |>
    mutate(copula = map(length, ~runif(.x))) |>
    select(!length)
}

generate_tornado_copula_unit_variable <- function(copula_factory, name)
{
  res <- copula_factory |>
    mutate(tornado_name = name)
  
  res_low <- res |>
    mutate(copula = map2(copula_name, length, ~ rep(ifelse(.x == name, 0.025, 0.5), .y))) |>
    mutate(bound_type = "low") |>
    select(!length)
  
  res_high <- res |>
    mutate(copula = map2(copula_name, length, ~ rep(ifelse(.x == name, 0.975, 0.5), .y))) |>
    mutate(bound_type = "high") |>
    select(!length)
  
  bind_rows(res_low, res_high)
}

generate_tornado_copula_vector_variable <- function(copula_factory, name, n_sim)
{
  map(1:n_sim, ~copula_factory |>
        mutate(
          tornado_name = name,
          copula = map2(copula_name, length, function(.x, .y) {
            if (name == .x)
              return (runif(.y))
            rep(0.5, .y)
          }),
          sim = .x,
          bound_type = "sim"
        ) |>
        select(!length)) |>
    bind_rows()
}

generate_tornado_copula <- function(copula_factory, n_sim)
{
  map2(copula_factory$copula_name, copula_factory$length, function(name, length) {
    if (length == 1)
      return(generate_tornado_copula_unit_variable(copula_factory, name))
    generate_tornado_copula_vector_variable(copula_factory, name, n_sim)
  }) |>
    bind_rows() |>
    group_split(tornado_name, bound_type)
}

transform_model_variables <- function(variables, current_copula)
{
  variables |>
    left_join(current_copula, by = "copula_name") |>
    mutate(value = map2_dbl(random_function, copula, ~ .x(.y))) |>
    select(name, value)
}

randomise_columns <- function(data, current_copula, prob_details) 
{
  transform_details <- prob_details |>
    left_join(current_copula, by = "copula_name")
  
  for (i in 1:dim(transform_details)[1])
  {
    row <- transform_details[i,]
    col_name <- row$name
    foo <- row$random_function[[1]]
    cop <- row$copula[[1]]
    data[[col_name]] <- foo(data[[col_name]], cop)
  }
  
  data
}

transform_data_cohort <- function(data, cohort, current_model_variables, transformations)
{
  # Load up all 
  setNames(current_model_variables$value, current_model_variables$name) |>
    as.list() |>
    list2env(environment())
  
  for (transformation in transformations)
  {
    environment(transformation) <- environment()
    data <- data |> transformation(cohort)
  }
   
  
  data
}


# Probability Helper Functions --------------------------------------------

beta_confidence <- function(x, p, sd_frac)
{
  mu = x
  s.d = sd_frac * mu
  
  alpha = ((1 - mu)/s.d^2 - 1/mu) * mu^2
  beta = alpha * (1/mu - 1)
  # A bound on beta
  #beta = min(beta, sqrt(mu * (1-mu)))
  
  alpha = ifelse(is.nan(alpha), 0, alpha)
  beta = ifelse(is.nan(beta), 1, beta)
  
  qbeta(p, alpha, beta)
}

beta_low <- \(x, p) beta_confidence(x, p, 0.1)

norm_confidence <- function(x, p, sd_frac)
{
  mu = x
  s.d = mu * sd_frac
  
  qnorm(p, mu, s.d)
}

norm_low <- \(x, p) norm_confidence(x, p, 0.1)


# Markov Model ------------------------------------------------------------

starting_states <- function(data, states)
{
  data <- data |>
    filter(starting_age)
  
  res <- data |>
    select(starts_with(states) & ends_with("prevalence")) |>
    mutate(
      healthy = 1 - rowSums(across(everything())),
      dead = 0
    ) |>
    # Changes names from disease_prevalence to disease
    rename_with(\(x) str_extract(x, "[^_]+"), ends_with("prevalence")) |>
    select(healthy, all_of(states), dead) |>
    unlist() |>
    array(dim = c(dim(data)[1], length(states) + 2))
  
  colnames(res) <- c("healthy", states, "dead")
  res
}

try_pull <- function(data, try_string, default)
{
  if (try_string %in% names(data))
    return(data[[try_string]])
  
  data[[default]]
}

transition_matrix <- function(data, states)
{
  incidence_data <- data |>
    select(starts_with({{states}}) & ends_with(c({{states}}, "incidence")))
  
  mortality_data <- data |>
    select(mortality, starts_with({{states}}) & ends_with("mortality"))
  
  ext_states <- c("healthy", states, "dead")
  res <- array(0, 
               dim = c(length(ext_states), length(ext_states), dim(data)[1]),
               dimnames = list(ext_states, ext_states))
  
  for (i in seq2(1, length(ext_states) - 1))
  {
    src_state <- ext_states[[i]]
    mortality <- mortality_data |>
      try_pull(
        try_string = paste(src_state, "mortality", sep = "_"), 
        default = "mortality"
      ) 
    
    cumsum <- 0
    
    for (j in seq2(i+1, length(ext_states) - 1))
    {
      dest_state <- ext_states[[j]]
      incidence <- incidence_data |>
        try_pull(
          try_string = paste(src_state, dest_state, sep = "_"),
          default = paste(dest_state, "incidence", sep = "_")
        ) * (1 - mortality)
      
      res[src_state, dest_state,] <- incidence
      cumsum <- cumsum + incidence
    }
    
    res[src_state, "dead",] <- mortality
    res[src_state, src_state,] <- 1 - cumsum - mortality
  }
  
  res["dead", "dead",] <- 1
  res
}

markov_function <- function(P, S, s_index, age)
{
  s <- S[s_index,]
  
  res <- array(0, dim = c(dim(P)[3] - s_index + 1, length(s)))
  
  colnames(res) <- colnames(S)
  res[1,] <- s
  
  for (i in seq2(s_index, dim(P)[3] - 1))
  {
    j <- i - s_index
    res[j+2,] <- res[j+1,] %*% P[,,i]
  }
  
  res <- as_tibble(res)
  res$initial_age = age
  res$cycle = 1:(dim(res)[1])
  res
}

apply_markov_model <- function(data, states)
{
  S <- starting_states(data, states)
  P <- transition_matrix(data, states)
  min_age <- data |> pull(age) |> min()
  
  indices <- seq2(1, dim(S)[1])
  indices |>
    lapply(\(x) markov_function(P, S, x, x+min_age-1)) |>
    bind_rows()
}

calculate_disability_daly <- function(data, to_group_by, states)
{
  # Quite inefficient
  data |>
    select(initial_age, cycle, popcount, contains(states)) |>
    rename_with(\(x) paste0(x, "_prop"), states) |>
    pivot_longer(
      cols = !(initial_age:popcount),
      names_to = c("disease", ".value"),
      names_sep = "_"
    ) |>
    mutate(
      daly = prop * disability * popcount,
      .keep = "unused"
    ) |>
    pivot_wider(
      id_cols = initial_age:cycle,
      names_from = disease,
      values_from = daly
    ) |>
    group_by({{to_group_by}}) |>
    summarise(
      across(states, sum)
    ) |>
    rename_with(\(x) paste0(x, "_daly"), states)
}

disability_daly <- function(states)
{
  states <- states
  
  calculate_result <- function(data, state_results)
  {
    to_join <- data |>
      select(age, popcount, starts_with(states) & ends_with("disability"))
    
    state_results |>
      left_join(to_join, by = join_by(initial_age == age)) |>
      calculate_disability_daly(cycle, states)
  }
  
  calculate_result
}

get_prevalence <- function(states)
{
  states <- states
  
  calculate_result <- function(data, state_results)
  {
    data <- data |>
      select(age, popcount)
    
    state_results |>
      left_join(data, by = join_by(initial_age == age)) |>
      mutate(
        alive = 1 - dead,
        across({{states}}, \(x) x / alive)
      ) |>
      group_by(cycle) |>
      summarise(
        totpop = sum(popcount),
        across({{states}}, \(x) sum(x * popcount) / totpop)
      ) |>
      select({{states}}) |>
      rename_with(\(x) paste0(x, "_prevalence"), states)
  }
}

calculate_results <- function(data, states, ...)
{
  state_results <- apply_markov_model(data, states) 
  
  list(...) |>
    map(\(x) x(data, state_results)) |>
    bind_cols() 
}


# Summary Statistics ------------------------------------------------------

get_cumpopcount <- function(data)
{
  data |>
    select(sex, ses, age, popcount) |>
    filter(age <= 84) |>
    group_by(sex, ses) |>
    mutate(cumpopcount = cumsum(popcount)) |>
    slice(1:n(), rep(n(), 110-84)) |>
    select(cumpopcount) |>
    arrange(desc(cumpopcount)) |>
    mutate(cycle = 1:91) |>
    ungroup()
}

group_sex <- function(res, data, across_condition)
{
  res |>
    left_join(get_cumpopcount(data), by = c("sex", "ses", "cycle")) |>
    group_by(cycle, ses, trt) |>
    summarise(
      across({{across_condition}}, \(x) sum(x * cumpopcount) / sum(cumpopcount))
    ) |>
    ungroup()
}

discount_columns <- function(res, columns, discount_rate = 0.03)
{
  discount_rate_tibble <- tibble(cycle = 1:max(res$cycle)) |>
    mutate(dr = (1 - discount_rate)^(cycle - 1))
  
  res |>
    left_join(discount_rate_tibble, by = "cycle") |>
    mutate(across({{columns}}, \(x) x * dr)) |>
    select(!dr)
}

normalise_columns <- function(res, data, columns, groups, per = 1000)
{
  max_pop <- data |>
    select(popcount, {{groups}}) |>
    group_by(pick({{groups}})) |>
    summarise(cohort_popcount = sum(popcount))
    
  res |>
    left_join(max_pop, by = {{groups}}) |>
    mutate(across({{columns}}, \(x) x * per / cohort_popcount))
}

lengthen_data_by_disease <- function(res, cols)
{
  res |>
    pivot_longer(
      cols = ends_with({{cols}}),
      names_to = c("disease", ".value"),
      names_sep = "_"
    )
}

match_cohort <- function(res, cohort_condition, id_keys, columns, cohort_name)
{
  drop_key <- "EcmEnv.y"
  
  res_cohort <- res |>
    filter({{cohort_condition}}) |>
    rename_with(\(x) paste(x, cohort_name, sep = "_"), columns)
  
  res |>
    filter(!{{cohort_condition}}) |>
    left_join(res_cohort, by = id_keys, suffix = c("", "EcmEnv.y")) |>
    select(!ends_with(drop_key))
}

summarise_statistics <- function(res, id_keys, columns, statistics = list(
  mean = \(x) mean(x, na.rm = TRUE),
  sd = \(x) sd(x, na.rm = TRUE),
  lb = \(x) quantile(x, p = 0.025, na.rm = TRUE),
  ub = \(x) quantile(x, p = 0.975, na.rm = TRUE)
))
{
  res |>
    group_by(pick({{id_keys}})) |>
    summarize(across({{columns}}, {{statistics}}))
}