# %%
source("setup.R")
source("dgp_functions.R")


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
