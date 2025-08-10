# %%
library(nycflights13)
library(tradepolicy)
library(arrow)
source("setup.R")

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
    "nyc taxi", 5L, nrow(nyc), 3L, list(\() nyc)
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
    "lfe::felm", list(\(df) {
      lfe_timer(
        df,
        tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type
      )
    }),
    "fixest::feols", list(\(df) {
      feols_timer(
        df,
        tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type
      )
    })
  )
)

# %%
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

write.csv(
  bench_real_data,
  here("results", "bench_ols_real_data.csv"),
  row.names = FALSE
)

bench_real_data[, .(mean_time = mean(time)), by = .(dgp_name, est_name)]
