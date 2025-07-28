# %%
library(fixest)
library(reticulate)
library(JuliaCall)
library(here)
library(data.table)
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
    cat("Setting", preload_var, "to:", libunwind_path, "\n")
    do.call(Sys.setenv, setNames(list(libunwind_path), preload_var))
    dyn.load(libunwind_path)
  } else {
    warning("libunwind not found, Julia may have issues")
  }
}

JuliaCall::julia_setup()
JuliaCall::julia_eval('import Pkg; Pkg.activate("."); Pkg.instantiate();')
JuliaCall::julia_source("timers/FixedEffectModels.jl")

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
bench_ols_basic <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "basic_dgp", 10L, 1e3, list(\() basic_dgp(n = 1e3)),
    "basic_dgp", 10L, 1e4, list(\() basic_dgp(n = 1e4)),
    "basic_dgp", 5L, 1e5, list(\() basic_dgp(n = 1e5)),
    "basic_dgp", 5L, 1e6, list(\() basic_dgp(n = 1e6))
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1)"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1
      )
    }),
    "pyfixest.feols", 2L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1 + dum_2"
      )
    }),
    "FixedEffectModels.reg", 2L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1) + fe(dum_2)"
      )
    }),
    "lfe::felm", 2L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1 + dum_2
      )
    }),
    "fixest::feols", 2L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1 + dum_2
      )
    }),
    "pyfixest.feols", 3L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1 + dum_2 + dum_3"
      )
    }),
    "FixedEffectModels.reg", 3L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1) + fe(dum_2) + fe(dum_3)"
      )
    }),
    "lfe::felm", 3L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1 + dum_2 + dum_3
      )
    }),
    "fixest::feols", 3L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1 + dum_2 + dum_3
      )
    })
  )
)


# %%
# fmt: skip
bench_ols_difficult <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "difficult_dgp", 10L, 1e3, list(\() difficult_dgp(n = 1e3)),
    "difficult_dgp", 10L, 1e4, list(\() difficult_dgp(n = 1e4)),
    "difficult_dgp", 5L, 1e5, list(\() difficult_dgp(n = 1e5)),
    "difficult_dgp", 5L, 1e6, list(\() difficult_dgp(n = 1e6))
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
    }),
    "pyfixest.feols", 3L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | indiv_id + firm_id + id_year"
      )
    }),
    "FixedEffectModels.reg", 3L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      )
    }),
    "lfe::felm", 3L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + id_year
      )
    }),
    "fixest::feols", 3L, list(\(df) {
      feols_timer(
        df,
        y ~ x1 | indiv_id + firm_id + id_year
      )
    })
  )
)

# %%
# fmt: skip
bench_akm <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, p_move=, n_iters=, n_obs=, dgp_function=,
    "AKM", 0.01, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.01),
    "AKM", 0.05, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.05),
    "AKM", 0.20, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.20),
    "AKM", 0.40, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.40),
    "AKM", 0.60, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.60),
    "AKM", 0.80, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.80),
    "AKM", 0.95, 5L, 1e5, 
    \() simulate_bipartite(n_workers = 1e5/5, n_time = 5, p_move = 0.95),
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
    "difficult_dgp", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp", 10L, 1e4, \() difficult_dgp(n = 1e4),
    "difficult_dgp", 5L, 1e5, \() difficult_dgp(n = 1e5),
    "difficult_dgp", 5L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.fepois", 2L, list(\(df) {
      pyfixest_fepois_timer(
        df,
        "abs_y ~ x1 | indiv_id + firm_id"
      )
    }),
    "GLFixedEffectModels Poisson", 2L, list(\(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "abs_y ~ x1 + fe(indiv_id) + fe(firm_id)"
      )
    }),
    "alpaca Poisson", 2L, list(\(df) {
      alpaca_poisson_timer(
        df,
        abs_y ~ x1 | indiv_id + firm_id
      )
    }),
    "fixest::fepois", 2L, list(\(df) {
      fepois_timer(
        df,
        abs_y ~ x1 | indiv_id + firm_id
      )
    }),
    "pyfixest.fepois", 3L, list(\(df) {
      pyfixest_fepois_timer(
        df,
        "abs_y ~ x1 | indiv_id + firm_id + id_year"
      )
    }),
    "GLFixedEffectModels Poisson", 3L, list(\(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "abs_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      )
    }),
    "alpaca Poisson", 3L, list(\(df) {
      alpaca_poisson_timer(
        df,
        abs_y ~ x1 | indiv_id + firm_id + id_year      
      )
    }),
    "fixest::fepois", 3L, list(\(df) {
      fepois_timer(
        df,
        abs_y ~ x1 | indiv_id + firm_id + id_year      
      )
    })
  )
)

# %%
# fmt: skip
bench_logit <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "difficult_dgp", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp", 10L, 1e4, \() difficult_dgp(n = 1e4),
    "difficult_dgp", 5L, 1e5, \() difficult_dgp(n = 1e5),
    "difficult_dgp", 5L, 1e6, \() difficult_dgp(n = 1e6)
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
    #     "binary_y ~ x1 | indiv_id + firm_id + id_year"
    #   )
    # }),
    "GLFixedEffectModels logit", 3L, list(\(df) {
      julia_call(
        "jl_logit_timer",
        df,
        "binary_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      )
    }),
    "alpaca logit", 3L, list(\(df) {
      alpaca_feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id + id_year      
      )
    }),
    "fixest logit", 3L, list(\(df) {
      feglm_logit_timer(
        df,
        binary_y ~ x1 | indiv_id + firm_id + id_year      
      )
    })
  )
)

# %%
# fmt: skip
bench_ols_multiple_y <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "difficult_dgp (2 outcomes)", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp (2 outcomes)", 10L, 1e4, \() difficult_dgp(n = 1e4),
    "difficult_dgp (2 outcomes)", 5L, 1e5, \() difficult_dgp(n = 1e5),
    "difficult_dgp (2 outcomes)", 5L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_timer(
        df,
        "y + abs_y ~ x1 | indiv_id + firm_id + id_year"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "abs_y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + id_year
      ) + 
      lfe_timer(
        df,
        abs_y ~ x1 | indiv_id + firm_id + id_year
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_timer(
        df,
        c(y, abs_y) ~ x1 | indiv_id + firm_id + id_year
      )
    })
  )
)

# %%
# fmt: skip
bench_ols_multiple_vcov <- run_benchmark(
  dgps = data.table::rowwiseDT(
    dgp_name=, n_iters=, n_obs=, dgp_function=,
    "difficult_dgp (hc1 + clustered vcov)", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp (hc1 + clustered vcov)", 10L, 1e4, \() difficult_dgp(n = 1e4),
    "difficult_dgp (hc1 + clustered vcov)", 5L, 1e5, \() difficult_dgp(n = 1e5),
    "difficult_dgp (hc1 + clustered vcov)", 5L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = data.table::rowwiseDT(
    est_name=, n_fe=, func=,
    "pyfixest.feols", 1L, list(\(df) {
      pyfixest_feols_multiple_vcov_timer(
        df,
        "y ~ x1 | indiv_id + firm_id + id_year",
        "firm_id"
      )
    }),
    "FixedEffectModels.reg", 1L, list(\(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(indiv_id) + fe(firm_id) + fe(id_year)",
        vcov = "firm_id"
      )
    }),
    "lfe::felm", 1L, list(\(df) {
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + id_year
      ) + 
      lfe_timer(
        df,
        y ~ x1 | indiv_id + firm_id + id_year | 0 | firm_id
      )
    }),
    "fixest::feols", 1L, list(\(df) {
      feols_multiple_vcov_timer(
        df,
        y  ~ x1 | indiv_id + firm_id + id_year,
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
    "nycflights13", 5L, nrow(nycflights13::flights), 3L, \() 
    nycflights13::flights,
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

# %%
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
write.csv(
  bench_ols_basic,
  here("results", "bench_ols_basic.csv"),
  row.names = FALSE
)
write.csv(
  bench_ols_difficult,
  here("results", "bench_ols_difficult.csv"),
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
  bench_ols_real_data,
  here("results", "bench_ols_real_data.csv"),
  row.names = FALSE
)

# %%
library(tinytable)
bench_ols_basic |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_ols_basic), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_ols_difficult |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_ols_difficult), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_akm |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_akm), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_poisson |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_poisson), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_logit |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_logit), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_ols_multiple_vcov |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_ols_multiple_vcov), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")

bench_ols_flights |>
  as.data.table() |>
  _[,
    .(
      mean_time = mean(time, na.rm = TRUE),
      n_failures = sum(is.na(time))
    ),
    by = setdiff(names(bench_ols_real_data), c("iter", "time"))
  ] |>
  tinytable::tt() |>
  print("markdown")
