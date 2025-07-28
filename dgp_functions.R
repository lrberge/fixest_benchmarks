# %%
basic_dgp <- function(n = 1000) {
  dum_all = list()
  nb_dum = c(n / 20, floor(sqrt(n)), floor(n**.33))
  N = nb_dum**3
  dum_all[[1]] = sample(nb_dum[1], n, TRUE)
  dum_all[[2]] = sample(nb_dum[2], n, TRUE)
  dum_all[[3]] = sample(nb_dum[3], n, TRUE)

  x1 = rnorm(n)
  x2 = x1**2

  mu = 1 * x1 + 0.05 * x2
  for (m in 1:3) {
    coef_dum = rnorm(nb_dum[m])
    mu = mu + coef_dum[dum_all[[m]]]
  }
  mu = exp(mu)
  y = MASS::rnegbin(mu, theta = 0.5)

  df = data.frame(y, x1, x2, ln_y = log(y + 1))

  for (m in 1:3) {
    df[[paste0("dum_", m)]] = dum_all[[m]]
  }

  return(df)
}

# %%
## This benchmark data set is an adaptation of a benchmark of the authors of the Julia FixedEffectModels.jl software
difficult_dgp <- function(n = 1000) {
  nb_indiv = n / 20
  nb_firm = round(n / 160)
  nb_year = round(n**.3)

  indiv_id = sample(1:nb_indiv, n, TRUE)
  firm_id = pmin(sample(0:20, n, TRUE) + pmax(1, indiv_id %/% 8 - 10), nb_firm)
  year = sample(nb_year, n, TRUE)

  x1 = 5 * cos(indiv_id) + 5 * sin(firm_id) + 5 * sin(year) + runif(n)
  x2 = cos(indiv_id) + sin(firm_id) + sin(year) + rnorm(n)
  y = 3 * x1 + 5 * x2 + cos(indiv_id) + cos(firm_id)^2 + sin(year) + rnorm(n)
  df = data.frame(
    indiv_id = indiv_id,
    firm_id = firm_id,
    year = year,
    x1 = x1,
    x2 = x2,
    y = y,
    abs_y = abs(y),
    binary_y = as.numeric(y > 0),
    ln_y = log(abs(y) + 1)
  )
  return(df)
}

# %%
## Adapated from bipartitepandas python package

#' Simulate Bipartite Labor Market Network
#'
#' Generates panel data for a bipartite network of workers and firms with
#' assortative matching, mobility, and AKM wage structure.
#'
#' @param n_workers Integer. Number of workers (default: 10000)
#' @param n_time Integer. Panel length in time periods (default: 5)
#' @param firm_size Numeric. Average firm size per period (default: 50)
#' @param n_firm_types Integer. Number of firm types (default: 10)
#' @param n_worker_types Integer. Number of worker types (default: 5)
#' @param p_move Numeric. Mobility probability per period (default: 0.5)
#' @param c_sort Numeric. Assortative sorting strength (default: 1)
#' @param c_netw Numeric. Network effect strength (default: 1)
#' @param c_sig Numeric. Sorting/network shock volatility (default: 1)
#' @param alpha_sig Numeric. Worker fixed effect volatility (default: 1)
#' @param psi_sig Numeric. Firm fixed effect volatility (default: 1)
#' @param w_sig Numeric. Wage shock volatility (default: 1)
#' @param x_sig Numeric. Covariate volatility (default: 1)
#' @param y1_beta Numeric. Covariate effect on y1 (default: 0.5)
#' @param y1_sig Numeric. y1 error volatility (default: 1)
#' @param y2_beta Numeric. Covariate effect on y2 (default: 0.3)
#' @param y2_sig Numeric. y2 error volatility (default: 1)
#'
#' @return Data frame with columns:
#'   \item{unit}{Worker identifier}
#'   \item{firm_id}{Firm identifier}
#'   \item{wage}{Original log wage (AKM structure)}
#'   \item{y1}{First outcome variable}
#'   \item{y2}{Second outcome variable}
#'   \item{x}{Covariate}
#'   \item{time}{Time period}
#'   \item{worker_type}{Worker type}
#'   \item{firm_type}{Firm type}
#'   \item{worker_fe}{Worker fixed effect}
#'   \item{firm_fe}{Firm fixed effect}
#'
#' @export
simulate_bipartite <- function(
  n_workers = 10000L,
  n_time = 5L,
  firm_size = 50,
  n_firm_types = 10L,
  n_worker_types = 5L,
  p_move = 0.5,
  c_sort = 1,
  c_netw = 1,
  c_sig = 1,
  alpha_sig = 1,
  psi_sig = 1,
  w_sig = 1,
  x_sig = 1,
  y1_beta = 0.5,
  y1_sig = 1,
  y2_beta = 0.3,
  y2_sig = 1
) {
  # Parameter validation
  stopifnot(
    "n_workers must be positive integer" = n_workers >= 1,
    "n_time must be positive integer" = n_time >= 1,
    "firm_size must be positive" = firm_size > 0,
    "n_firm_types must be positive integer" = n_firm_types >= 1,
    "n_worker_types must be positive integer" = n_worker_types >= 1,
    "alpha_sig must be non-negative" = alpha_sig >= 0,
    "psi_sig must be non-negative" = psi_sig >= 0,
    "w_sig must be non-negative" = w_sig >= 0,
    "c_sig must be non-negative" = c_sig >= 0,
    "p_move must be in [0,1]" = p_move >= 0 && p_move <= 1,
    "x_sig must be non-negative" = x_sig >= 0,
    "y1_sig must be non-negative" = y1_sig >= 0,
    "y2_sig must be non-negative" = y2_sig >= 0
  )

  # ============================================================================
  # GENERATE FIXED EFFECTS
  # ============================================================================

  # Generate fixed effects using inverse normal CDF
  psi <- qnorm(seq_len(n_firm_types) / (n_firm_types + 1)) * psi_sig
  alpha <- qnorm(seq_len(n_worker_types) / (n_worker_types + 1)) * alpha_sig

  # Compute transition matrices
  G <- array(0, dim = c(n_worker_types, n_firm_types, n_firm_types))
  for (l in seq_len(n_worker_types)) {
    for (k_from in seq_len(n_firm_types)) {
      probs <- dnorm((psi - c_netw * psi[k_from] - c_sort * alpha[l]) / c_sig)
      G[l, k_from, ] <- probs / sum(probs)
    }
  }

  # Compute stationary distributions
  H <- matrix(0, nrow = n_worker_types, ncol = n_firm_types)
  for (l in seq_len(n_worker_types)) {
    eig_decomp <- eigen(t(G[l, , ]))
    stationary_idx <- which.min(abs(eig_decomp$values - 1))
    stationary_vec <- Re(eig_decomp$vectors[, stationary_idx])
    stationary_vec <- abs(stationary_vec) / sum(abs(stationary_vec))
    H[l, ] <- stationary_vec
  }

  # ============================================================================
  # SIMULATE MOBILITY
  # ============================================================================

  # Generate worker types
  worker_types <- sample(seq_len(n_worker_types), n_workers, replace = TRUE)

  # Initialize mobility matrices
  firm_types <- matrix(0L, nrow = n_workers, ncol = n_time)
  spell_ids <- matrix(0L, nrow = n_workers, ncol = n_time)
  spell_counter <- 1L

  # Simulate worker mobility over time
  for (i in seq_len(n_workers)) {
    l <- worker_types[i]

    # Initial firm placement
    firm_types[i, 1] <- sample(seq_len(n_firm_types), 1, prob = H[l, ])
    spell_ids[i, 1] <- spell_counter
    spell_counter <- spell_counter + 1L

    # Mobility decisions for subsequent periods
    if (n_time > 1) {
      for (t in 2:n_time) {
        if (runif(1) < p_move) {
          # Worker moves firms
          firm_types[i, t] <- sample(
            seq_len(n_firm_types),
            1,
            prob = G[l, firm_types[i, t - 1], ]
          )
          spell_ids[i, t] <- spell_counter
          spell_counter <- spell_counter + 1L
        } else {
          # Worker stays at same firm
          firm_types[i, t] <- firm_types[i, t - 1]
          spell_ids[i, t] <- spell_ids[i, t - 1]
        }
      }
    }
  }

  # ============================================================================
  # CONSTRUCT PANEL
  # ============================================================================
  # Create long-format panel structure
  indiv_id <- rep(seq_len(n_workers), each = n_time)
  time <- rep(seq_len(n_time), times = n_workers)
  worker_type <- rep(worker_types, each = n_time)
  firm_type <- as.vector(t(firm_types))
  spell <- as.vector(t(spell_ids))

  # Compute spell sizes for firm ID assignment
  spell_summary <- aggregate(
    list(spell_size = spell),
    by = list(spell = spell, firm_type = firm_type),
    FUN = length
  )

  # Assign firm IDs (in_worker_typesined from assign_firm_ids)
  firm_ids <- integer(nrow(spell_summary))

  for (k_type in unique(spell_summary$firm_type)) {
    k_mask <- spell_summary$firm_type == k_type
    k_spells <- spell_summary[k_mask, ]

    # Calculate number of firms needed for this firm type
    total_obs <- sum(k_spells$spell_size)
    n_firms <- max(1L, round(total_obs / (firm_size * n_time)))

    # Randomly assign spells to firms
    assigned_ids <- sample(seq_len(n_firms), nrow(k_spells), replace = TRUE)

    # Create contiguous firm IDs
    firm_mapping <- data.frame(
      temp_id = assigned_ids,
      firm_id = seq_along(unique(assigned_ids))[match(
        assigned_ids,
        unique(assigned_ids)
      )]
    )

    firm_ids[k_mask] <- firm_mapping$firm_id[match(
      assigned_ids,
      firm_mapping$temp_id
    )]
  }

  spell_summary$firm_id <- firm_ids

  # Merge firm IDs back to panel
  panel <- data.frame(
    indiv_id = indiv_id,
    time = time,
    worker_type = worker_type,
    firm_type = firm_type,
    spell = spell
  )
  panel <- merge(panel, spell_summary[c("spell", "firm_id")], by = "spell")

  # Generate fixed effects and outcomes
  panel$worker_fe <- alpha[panel$worker_type]
  panel$firm_fe <- psi[panel$firm_type]
  panel$wage <- panel$worker_fe + panel$firm_fe + rnorm(nrow(panel)) * w_sig

  # Generate covariate and additional outcomes
  panel$x1 <- rnorm(nrow(panel)) * x_sig
  panel$y <- panel$worker_fe +
    panel$firm_fe +
    y1_beta * panel$x +
    rnorm(nrow(panel)) * y1_sig
  panel$y2 <- panel$worker_fe +
    panel$firm_fe +
    y2_beta * panel$x +
    rnorm(nrow(panel)) * y2_sig

  # Sort and return final dataset
  panel <- panel[order(panel$indiv_id, panel$time), ]

  return(panel[c(
    "indiv_id",
    "firm_id",
    "wage",
    "y",
    "y2",
    "x1",
    "time",
    "worker_type",
    "firm_type",
    "worker_fe",
    "firm_fe"
  )])
}
