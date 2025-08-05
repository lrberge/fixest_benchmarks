# %%
library(data.table)
library(fixest)
library(reticulate)
library(JuliaCall)
library(here)
library(nycflights13)
library(tradepolicy)
library(arrow)
source("dgp_functions.R")

# setup R
source("timers/fixest.R")
source("timers/lfe.R")
source("timers/alpaca.R")

# setup python
reticulate::use_virtualenv(here(".venv"), required = TRUE)
pd = reticulate::import("pandas")
pf = reticulate::import("pyfixest")
reticulate::source_python("timers/pyfixest.py")

# setup julia
# This chaos is due to: https://github.com/JuliaInterop/JuliaCall/issues/238
Sys.setenv(JULIA_PROJECT = here())
julia_path <- Sys.which("julia")
if (julia_path != "") {
  julia_bin_cmd <- system("julia -e 'println(Sys.BINDIR)'", intern = TRUE)
  if (length(julia_bin_cmd) > 0 && julia_bin_cmd != "") {
    julia_bin_dir <- julia_bin_cmd[1]
    julia_lib_dir <- file.path(dirname(julia_bin_dir), "lib", "julia")
  } else {
    # Fallback: try to find julia executable
    julia_path <- Sys.which("julia")
    if (julia_path == "") {
      stop("Julia not found")
    }
    julia_bin_dir <- dirname(julia_path)
    julia_lib_dir <- file.path(dirname(julia_bin_dir), "lib", "julia")
  }

  cat("Looking for libunwind in:", julia_lib_dir, "\n")
  # Platform-specific patterns
  if (Sys.info()["sysname"] == "Darwin") {
    # macOS
    libunwind_patterns <- c("libunwind*.dylib", "*unwind*.dylib")
    preload_var <- "DYLD_INSERT_LIBRARIES"
  } else {
    # Linux
    libunwind_patterns <- c("libunwind.so*", "libunwind-*.so*", "*unwind*.so*")
    preload_var <- "LD_PRELOAD"
  }

  # Look for libunwind in Julia's lib directory
  libunwind_path <- NULL
  for (pattern in libunwind_patterns) {
    files <- Sys.glob(file.path(julia_lib_dir, pattern))
    if (length(files) > 0) {
      libunwind_path <- files[1]
      break
    }
  }
  if (!is.null(libunwind_path) && file.exists(libunwind_path)) {
    dyn.load(libunwind_path)
  } else {
    warning("libunwind not found, Julia may have issues")
  }
}

JuliaCall::julia_setup()
JuliaCall::julia_eval('import Pkg; Pkg.activate("."); Pkg.instantiate();')
JuliaCall::julia_source("timers/FixedEffectModels.jl")

# From https://github.com/pachadotdev/capybara/blob/main/dev/benchmarks-no-base.R
ch1_application3 <- tradepolicy::agtpa_applications |>
  as.data.table() |>
  _[year %in% seq(1986, 2006, 4), ] |>
  _[, `:=`(
    exp_year = paste0(exporter, year),
    imp_year = paste0(importer, year),
    year = paste0("intl_border_", year),
    log_trade = log(trade),
    log_dist = log(dist),
    intl_brdr = ifelse(exporter == importer, pair_id, "inter"),
    intl_brdr_2 = ifelse(exporter == importer, 0, 1),
    pair_id_2 = ifelse(exporter == importer, "0-intra", pair_id)
  )] |>
  dcast(... ~ year, value.var = "intl_brdr_2", fill = 0) |>
  _[, sum_trade := sum(trade), by = pair_id]

nyc <- read_parquet(here("data/nyc_taxi.parquet"))


# set seed for (somewhat) reproducibility
set.seed(20250725)

# %%
run_benchmark <- function(
  name = "",
  dgps,
  estimators,
  burn_in = 1L
) {
  res <- NULL
  cat("\n")
  for (dgp_k in seq_len(nrow(dgps))) {
    dgp <- dgps$dgp_function[[dgp_k]]
    n_iters <- dgps$n_iters[dgp_k]
    cat("---- Starting DGP ----\n")
    print(subset(dgps[dgp_k, ], select = -c(dgp_function)))
    cat("\n")

    i = 1L
    while (i <= n_iters + burn_in) {
      cat(".")
      df <- dgp()
      times <- unlist(lapply(seq_len(nrow(estimators)), function(estimator_k) {
        f <- estimators$func[[estimator_k]]
        tryCatch(
          f(df),
          error = function(error) NA_real_
        )
      }))

      if (i > burn_in) {
        res_i <- data.frame(
          iter = i,
          time = times
        )
        res_i <- cbind(res_i, subset(estimators, select = -c(func)))
        res_i <- cbind(res_i, subset(dgps[dgp_k, ], select = -c(dgp_function)))
        res <- rbind(res, res_i)
      }
      i = i + 1
    }
    cat("\n\n")
  }

  return(res)
}

# %%
# fmt: skip
bench_ols <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "simple",    10L, 1e3, list(\() base_dgp(n = 1e3, type = "simple")),
    "difficult", 10L, 1e3, list(\() base_dgp(n = 1e3, type = "difficult")),
    "simple",    10L, 1e4, list(\() base_dgp(n = 1e4, type = "simple")),
    "difficult", 10L, 1e4, list(\() base_dgp(n = 1e4, type = "difficult")),
    "simple",    10L, 1e5, list(\() base_dgp(n = 1e5, type = "simple")),
    "difficult",  5L, 1e5, list(\() base_dgp(n = 1e5, type = "difficult")),
    "simple",    10L, 1e6, list(\() base_dgp(n = 1e6, type = "simple")),
    "difficult",  5L, 1e6, list(\() base_dgp(n = 1e6, type = "difficult"))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | indiv_id"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id)"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | indiv_id
      )
    }),
    "pyfixest.feols", 2L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | indiv_id + year"
      )
    }),
    "FixedEffectModels.reg", 2L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(year)"
      )
    }),
    "lfe::felm", 2L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + year
      )
    }),
    "fixest::feols", 2L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | indiv_id + year
      )
    }),
    "pyfixest.feols", 3L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | indiv_id + year + firm_id"
      )
    }),
    "FixedEffectModels.reg", 3L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(year) + fe(firm_id)"
      )
    }),
    "lfe::felm", 3L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + year + firm_id
      )
    }),
    "fixest::feols", 3L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | indiv_id + year + firm_id
      )
    })
  )
)

# %%
# fmt: skip
bench_akm <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, p_move=, n_iters=, n_obs=, dgp_function=,
    "AKM", 0.01, 5L, 5e5, 
    list(\() simulate_bipartite(n_workers = 2.5e5/5, n_time = 5, p_move = 0.01)),
    "AKM", 0.05, 5L, 5e5, 
    list(\() simulate_bipartite(n_workers = 2.5e5/5, n_time = 5, p_move = 0.05)),
    "AKM", 0.20, 5L, 5e5, 
    list(\() simulate_bipartite(n_workers = 2.5e5/5, n_time = 5, p_move = 0.20)),
    "AKM", 0.40, 5L, 5e5,
    list(\() simulate_bipartite(n_workers = 2.5e5/5, n_time = 5, p_move = 0.40))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 2L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | indiv_id + firm_id"
      )
    }),
    "FixedEffectModels.reg", 2L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id)"
      )
    }),
    "lfe::felm", 2L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id
      )
    }),
    "fixest::feols", 2L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | indiv_id + firm_id
      )
    })
  )
)

# %%
# fmt: skip
bench_poisson <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "simple", 10L, 1e3, list(\() base_dgp(n = 1e3, type = "simple")),
    "simple", 10L, 1e4, list(\() base_dgp(n = 1e4, type = "simple")),
    "simple",  5L, 1e5, list(\() base_dgp(n = 1e5, type = "simple")),
    "simple",  5L, 1e6, list(\() base_dgp(n = 1e6, type = "simple"))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.fepois", 2L, list(\(df) {
      pyfixest_fepois_timer(
        df,
        "ln_y ~ x1 | indiv_id + firm_id"
      )
    }),
    "GLFixedEffectModels Poisson", 2L, list(\(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "ln_y ~ x1 + fe(indiv_id) + fe(firm_id)"
      )
    }),
    "alpaca Poisson", 2L, list(\(df) {
      alpaca_poisson_timer(
        df,
        ln_y ~ x1 | indiv_id + firm_id
      )
    }),
    "fixest::fepois", 2L, list(\(df) {
      fepois_timer(
        df,
        ln_y ~ x1 | indiv_id + firm_id
      )
    }),
    "pyfixest.fepois", 3L, list(\(df) {
      pyfixest_fepois_timer(
        df,
        "ln_y ~ x1 | indiv_id + firm_id + year"
      )
    }),
    "GLFixedEffectModels Poisson", 3L, list(\(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "ln_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)"
      )
    }),
    "alpaca Poisson", 3L, list(\(df) {
      alpaca_poisson_timer(
        df,
        ln_y ~ x1 | indiv_id + firm_id + year      
      )
    }),
    "fixest::fepois", 3L, list(\(df) {
      fepois_timer(
        df,
        ln_y ~ x1 | indiv_id + firm_id + year      
      )
    })
  )
)

# %%
# fmt: skip
bench_logit <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "simple", 10L, 1e3, list(\() base_dgp(n = 1e3, type = "simple")),
    "simple", 10L, 1e4, list(\() base_dgp(n = 1e4, type = "simple")),
    "simple",  5L, 1e5, list(\() base_dgp(n = 1e5, type = "simple")),
    "simple",  5L, 1e6, list(\() base_dgp(n = 1e6, type = "simple"))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    # "pyfixest logit", 2L, list(\(df) {
    #   pyfixest_feglm_logit_timer(
    #     df,
    #     "binary_y ~ x1 | indiv_id + firm_id"
    #   )
    # }),
    "GLFixedEffectModels logit", 2L, list(\(df) {
      julia_call(
        "jl_logit_timer",
        df,
        "binary_y ~ x1 + fe(indiv_id) + fe(firm_id)"
      )
    }),
    "alpaca logit", 2L, list(\(df) {
      alpaca_feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id
      )
    }),
    "fixest logit", 2L, list(\(df) {
      feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id
      )
    }),
    # "pyfixest logit", 3L, list(\(df) {
    #   pyfixest_feglm_logit_timer(
    #     df,
    #     "binary_y ~ x1 | indiv_id + firm_id + year"
    #   )
    # }),
    "GLFixedEffectModels logit", 3L, list(\(df) {
      julia_call(
        "jl_logit_timer",
        df,
        "binary_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)"
      )
    }),
    "alpaca logit", 3L, list(\(df) {
      alpaca_feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id + year      
      )
    }),
    "fixest logit", 3L, list(\(df) {
      feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id + year      
      )
    })
  )
)

# %%
# fmt: skip
bench_ols_multiple_y <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "base_dgp (simple, 2 outcomes)", 10L, 1e3, 
    list(\() base_dgp(n = 1e3, type = "simple")),
    "base_dgp (simple, 2 outcomes)", 10L, 1e4, 
    list(\() base_dgp(n = 1e4, type = "simple")),
    "base_dgp (simple, 2 outcomes)",  5L, 1e5, 
    list(\() base_dgp(n = 1e5, type = "simple")),
    "base_dgp (simple, 2 outcomes)",  5L, 1e6, 
    list(\() base_dgp(n = 1e6, type = "simple"))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y + ln_y ~ x1 | indiv_id + firm_id + year"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "ln_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + year
      ) + 
      lfe_timer(
        df,
        ln_y ~ x1 | indiv_id + firm_id + year
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_timer(
        df,
        c(y, ln_y) ~ x1 | indiv_id + firm_id + year
      )
    })
  )
)

# %%
# fmt: skip
bench_ols_multiple_vcov <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "base_dgp (simple, hc1 + clustered)", 10L, 1e3, 
    list(\() base_dgp(n = 1e3, type = "simple")),
    "base_dgp (simple, hc1 + clustered)", 10L, 1e4, 
    list(\() base_dgp(n = 1e4, type = "simple")),
    "base_dgp (simple, hc1 + clustered)",  5L, 1e5, 
    list(\() base_dgp(n = 1e5, type = "simple")),
    "base_dgp (simple, hc1 + clustered)",  5L, 1e6, 
    list(\() base_dgp(n = 1e6, type = "simple"))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_multiple_vcov_timer(
        df,
        "y ~ x1 | indiv_id + firm_id + year",
        "firm_id"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(year)",
        vcov = "firm_id"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + year
      ) + 
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + year | 0 | firm_id
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_multiple_vcov_timer(
        df,
        y  ~ x1 | indiv_id + firm_id + year,
        cluster = ~firm_id
      )
    })
  )
)

# %%
# fmt: skip
bench_ols_flights <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, n_fe=, dgp_function=,
    "nycflights13", 5L, nrow(nycflights13::flights), 3L, list(\() 
    nycflights13::flights)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, func=,
    "pyfixest.feols", list(\(df) {
      pyfixest_feols_timer(
        df,
        "arr_delay ~ distance | carrier + origin + dest"
      )
    }),
    "FixedEffectModels.reg", list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "arr_delay ~ distance + fe(carrier) + fe(origin) + fe(dest)"
      )
    }),
    "lfe::felm", list(\(df) {
      lfe_timer(
        df,
        arr_delay ~ distance | carrier + origin + dest
      )
    }),
    "fixest::feols", list(\(df) {
      feols_timer(
        df,
        arr_delay ~ distance | carrier + origin + dest
      )
    })
  )
)

# fmt: skip
bench_tradepolicy_ols <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, n_fe=, dgp_function=,
    "tradepolicy", 5L, nrow(ch1_application3), 3L, list(\() ch1_application3)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, func=,
    "pyfixest.feols", list(\(df) {
      pyfixest_feols_timer(
        df,
        "trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr"
      )
    }),
    "FixedEffectModels.reg", list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "trade ~ log_dist + cntg + lang + clny + rta + fe(exp_year) + fe(imp_year) + fe(intl_brdr)"
      )
    }),
    "lfe::felm", list(\(df) {
      lfe_timer(
        df,
        trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr
      )
    }),
    "fixest::feols", list(\(df) {
      feols_timer(
        df,
        trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr
      )
    })
  )
)

# fmt: skip
bench_tradepolicy_ppml <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, n_fe=, dgp_function=,
    "tradepolicy", 5L, nrow(ch1_application3), 3L, list(\() ch1_application3)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, func=,
    "pyfixest.fepois", list(\(df) {
      pyfixest_fepois_timer(
        df,
        "trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr"
      )
    }),
    "GLFixedEffectModels Poisson", list(\(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "trade ~ log_dist + cntg + lang + clny + rta + fe(exp_year) + fe(imp_year) + fe(intl_brdr)"
      )
    }),
    "alpaca Poisson", list(\(df) {
      alpaca_poisson_timer(
        df,
        trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr
      )
    }),
    "fixest::fepois", list(\(df) {
      fepois_timer(
        df,
        trade ~ log_dist + cntg + lang + clny + rta | exp_year + imp_year + intl_brdr
      )
    })
  )
)

# fmt: skip
bench_nyc_taxi_ols <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, n_fe=, dgp_function=,
    "nyc taxi", 2L, nrow(nyc), 3L, list(\() nyc)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, func=,
    "pyfixest.feols", list(\(df) {
      pyfixest_feols_timer(
        df,
        "tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type"
      )
    }),
    "FixedEffectModels.reg", list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "tip_amount ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)"
      )
    }),
    # "lfe::felm", list(\(df) {
    #   lfe_timer(
    #     df,
    #     tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type
    #   )
    # }),
    "fixest::feols", list(\(df) {
      feols_timer(
        df,
        tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type
      )
    })
  )
)

bench_real_data <- rbindlist(
  list(
    bench_ols_flights,
    bench_tradepolicy_ols,
    bench_tradepolicy_ppml,
    bench_nyc_taxi_ols
  ),
  use.names = TRUE,
  fill = TRUE
)


# %%
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
write.csv(
  bench_ols,
  here("results", "bench_ols.csv"),
  row.names = FALSE
)
write.csv(
  bench_poisson,
  here("results", "bench_poisson.csv"),
  row.names = FALSE
)
write.csv(
  bench_logit,
  here("results", "bench_logit.csv"),
  row.names = FALSE
)
write.csv(
  bench_ols_multiple_y,
  here("results", "bench_ols_multiple_y.csv"),
  row.names = FALSE
)
write.csv(
  bench_ols_multiple_vcov,
  here("results", "bench_ols_multiple_vcov.csv"),
  row.names = FALSE
)
write.csv(
  bench_real_data,
  here("results", "bench_ols_real_data.csv"),
  row.names = FALSE
)
