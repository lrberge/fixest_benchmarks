# %%
library(fixest)
library(reticulate)
library(JuliaCall)
library(here)
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
JuliaCall::julia_setup()
JuliaCall::julia_eval('import Pkg; Pkg.activate(".")') # use project env
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
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "basic_dgp", 10L, 1e3, \() basic_dgp(n = 1e3),
    "basic_dgp", 10L, 1e4, \() basic_dgp(n = 1e4),
    "basic_dgp", 10L, 1e5, \() basic_dgp(n = 1e5),
    "basic_dgp", 10L, 1e6, \() basic_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    "pyfixest.feols", 1L, \(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1"
      )
    },
    "FixedEffectModels.reg", 1L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1)"
      )
    },
    "lfe::felm", 1L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1
      )
    },
    "fixest::feols", 1L, \(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1
      )
    },
    "pyfixest.feols", 2L, \(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1 + dum_2"
      )
    },
    "FixedEffectModels.reg", 2L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1) + fe(dum_2)"
      )
    },
    "lfe::felm", 2L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1 + dum_2
      )
    },
    "fixest::feols", 2L, \(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1 + dum_2
      )
    },
    "pyfixest.feols", 3L, \(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | dum_1 + dum_2 + dum_3"
      )
    },
    "FixedEffectModels.reg", 3L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(dum_1) + fe(dum_2) + fe(dum_3)"
      )
    },
    "lfe::felm", 3L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | dum_1 + dum_2 + dum_3
      )
    },
    "fixest::feols", 3L, \(df) {
      feols_timer(
        df,
        y ~ x1 | dum_1 + dum_2 + dum_3
      )
    }
  )
)


# %%
# fmt: skip
bench_ols_difficult <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "difficult_dgp", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp", 10L, 1e4, \() difficult_dgp(n = 1e4),
    # "difficult_dgp", 3L, 1e5, \() difficult_dgp(n = 1e5),
    # "difficult_dgp", 3L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    "pyfixest.feols", 2L, \(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | id_indiv + id_firm"
      )
    },
    "FixedEffectModels.reg", 2L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(id_indiv) + fe(id_firm)"
      )
    },
    "lfe::felm", 2L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | id_indiv + id_firm
      )
    },
    "fixest::feols", 2L, \(df) {
      feols_timer(
        df,
        y ~ x1 | id_indiv + id_firm
      )
    },
    "pyfixest.feols", 3L, \(df) {
      pyfixest_feols_timer(
        df,
        "y ~ x1 | id_indiv + id_firm + id_year"
      )
    },
    "FixedEffectModels.reg", 3L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      )
    },
    "lfe::felm", 3L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | id_indiv + id_firm + id_year
      )
    },
    "fixest::feols", 3L, \(df) {
      feols_timer(
        df,
        y ~ x1 | id_indiv + id_firm + id_year
      )
    }
  )
)

# %%
# fmt: skip
bench_poisson <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "difficult_dgp", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp", 10L, 1e4, \() difficult_dgp(n = 1e4),
    # "difficult_dgp", 3L, 1e5, \() difficult_dgp(n = 1e5),
    # "difficult_dgp", 3L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    "pyfixest.fepois", 2L, \(df) {
      pyfixest_fepois_timer(
        df,
        "abs_y ~ x1 | id_indiv + id_firm"
      )
    },
    "GLFixedEffectModels Poisson", 2L, \(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "abs_y ~ x1 + fe(id_indiv) + fe(id_firm)"
      )
    },
    "alpaca Poisson", 2L, \(df) {
      alpaca_poisson_timer(
        df,
        abs_y ~ x1 | id_indiv + id_firm
      )
    },
    "fixest::fepois", 2L, \(df) {
      fepois_timer(
        df,
        abs_y ~ x1 | id_indiv + id_firm
      )
    },
    "pyfixest.fepois", 3L, \(df) {
      pyfixest_fepois_timer(
        df,
        "abs_y ~ x1 | id_indiv + id_firm + id_year"
      )
    },
    "GLFixedEffectModels Poisson", 3L, \(df) {
      julia_call(
        "jl_poisson_timer",
        df,
        "abs_y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      )
    },
    "alpaca Poisson", 3L, \(df) {
      alpaca_poisson_timer(
        df,
        abs_y ~ x1 | id_indiv + id_firm + id_year      
      )
    },
    "fixest::fepois", 3L, \(df) {
      fepois_timer(
        df,
        abs_y ~ x1 | id_indiv + id_firm + id_year      
      )
    }
  )
)

# %%
# fmt: skip
bench_logit <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "difficult_dgp", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp", 10L, 1e4, \() difficult_dgp(n = 1e4),
    # "difficult_dgp", 3L, 1e5, \() difficult_dgp(n = 1e5),
    # "difficult_dgp", 3L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    # "pyfixest logit", 2L, \(df) {
    #   pyfixest_feglm_logit_timer(
    #     df,
    #     "binary_y ~ x1 | id_indiv + id_firm"
    #   )
    # },
    "GLFixedEffectModels logit", 2L, \(df) {
      julia_call(
        "jl_logit_timer",
        df,
        "binary_y ~ x1 + fe(id_indiv) + fe(id_firm)"
      )
    },
    "alpaca logit", 2L, \(df) {
      alpaca_feglm_logit_timer(
        df,
        binary_y ~ x1 | id_indiv + id_firm
      )
    },
    "fixest logit", 2L, \(df) {
      feglm_logit_timer(
        df,
        binary_y ~ x1 | id_indiv + id_firm
      )
    },
    # "pyfixest logit", 3L, \(df) {
    #   pyfixest_feglm_logit_timer(
    #     df,
    #     "binary_y ~ x1 | id_indiv + id_firm + id_year"
    #   )
    # },
    "GLFixedEffectModels logit", 3L, \(df) {
      julia_call(
        "jl_logit_timer",
        df,
        "binary_y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      )
    },
    "alpaca logit", 3L, \(df) {
      alpaca_feglm_logit_timer(
        df,
        binary_y ~ x1 | id_indiv + id_firm + id_year      
      )
    },
    "fixest logit", 3L, \(df) {
      feglm_logit_timer(
        df,
        binary_y ~ x1 | id_indiv + id_firm + id_year      
      )
    }
  )
)

# %%
# fmt: skip
bench_ols_multiple_y <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "difficult_dgp (2 outcomes)", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp (2 outcomes)", 10L, 1e4, \() difficult_dgp(n = 1e4),
    # "difficult_dgp (2 outcomes)", 3L, 1e5, \() difficult_dgp(n = 1e5),
    # "difficult_dgp (2 outcomes)", 3L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    "pyfixest.feols", 1L, \(df) {
      pyfixest_feols_timer(
        df,
        "y + abs_y ~ x1 | id_indiv + id_firm + id_year"
      )
    },
    "FixedEffectModels.reg", 1L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "abs_y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      )
    },
    "lfe::felm", 1L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | id_indiv + id_firm + id_year
      ) + 
      lfe_timer(
        df,
        abs_y ~ x1 | id_indiv + id_firm + id_year
      )
    },
    "fixest::feols", 1L, \(df) {
      feols_timer(
        df,
        c(y, abs_y) ~ x1 | id_indiv + id_firm + id_year
      )
    },
  )
)

# %%
# fmt: skip
bench_ols_multiple_vcov <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~dgp_function,
    "difficult_dgp (hc1 + clustered vcov)", 10L, 1e3, \() difficult_dgp(n = 1e3),
    "difficult_dgp (hc1 + clustered vcov)", 10L, 1e4, \() difficult_dgp(n = 1e4),
    # "difficult_dgp (hc1 + clustered vcov)", 3L, 1e5, \() difficult_dgp(n = 1e5),
    # "difficult_dgp (hc1 + clustered vcov)", 3L, 1e6, \() difficult_dgp(n = 1e6)
  ),
  estimators = tibble::tribble(
    ~est_name, ~n_fe, ~func,
    "pyfixest.feols", 1L, \(df) {
      pyfixest_feols_multiple_vcov_timer(
        df,
        "y ~ x1 | id_indiv + id_firm + id_year",
        "id_firm"
      )
    },
    "FixedEffectModels.reg", 1L, \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)"
      ) + 
      julia_call(
        "jl_feols_timer",
        df,
        "y ~ x1 + fe(id_indiv) + fe(id_firm) + fe(id_year)",
        vcov = "id_firm"
      )
    },
    "lfe::felm", 1L, \(df) {
      lfe_timer(
        df,
        y ~ x1 | id_indiv + id_firm + id_year
      ) + 
      lfe_timer(
        df,
        y ~ x1 | id_indiv + id_firm + id_year | 0 | id_firm
      )
    },
    "fixest::feols", 1L, \(df) {
      feols_multiple_vcov_timer(
        df,
        y  ~ x1 | id_indiv + id_firm + id_year,
        cluster = ~id_firm
      )
    },
  )
)

# %%
# fmt: skip
bench_ols_real_data <- run_benchmark(
  dgps = tibble::tribble(
    ~dgp_name, ~n_iters, ~n_obs, ~n_fe, ~dgp_function,
    "nycflights13", 5L, nrow(nycflights13::flights), 3L, \() 
    nycflights13::flights,
  ),
  estimators = tibble::tribble(
    ~est_name, ~func,
    "pyfixest.feols", \(df) {
      pyfixest_feols_timer(
        df,
        "arr_delay ~ distance | carrier + origin + dest"
      )
    },
    "FixedEffectModels.reg", \(df) {
      julia_call(
        "jl_feols_timer",
        df,
        "arr_delay ~ distance + fe(carrier) + fe(origin) + fe(dest)"
      )
    },
    "lfe::felm", \(df) {
      lfe_timer(
        df,
        arr_delay ~ distance | carrier + origin + dest
      )
    },
    "fixest::feols", \(df) {
      feols_timer(
        df,
        arr_delay ~ distance | carrier + origin + dest
      )
    }
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
# library(tidyverse)
# library(tinytable)
#
# bench_ols_basic |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
#
# bench_ols_difficult |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
#
# bench_poisson |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
#
# bench_logit |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
#
# bench_ols_multiple_vcov |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
#
# bench_ols_real_data |>
#   summarize(
#     .by = c(everything(), -iter, -time),
#     mean_time = mean(time, na.rm = TRUE),
#     n_failures = sum(is.na(time))
#   ) |>
#   tinytable::tt() |>
#   print("markdown")
