# Benchmarks for Fixed-Effect Estimation

This repository contains a set of open-source benchmarks for assessing
the performance of various statistical packages for estimating linear
and generalized linear models. The benchmarks are designed to evaluate
the speed and efficiency of different estimation methods under various
conditions, such as sample size, number of fixed effects, and model
complexity. These benchmarks are scheduled to run weekly on Github
Actions, and the updated results are stored in the `results` directory.

If you would like to contribute to the benchmarks or add alternative
estimators, PRs are welcome. The code `bench.R` is easily adaptable to
include additional estimators. Simply add a timer function to `timers.R`
that returns the time taken to estimate a model using your preferred
package. Then add your estimator to the `bench.R` file and submit a PR
for review.

> [!CAUTION]
>
> As with all benchmarks, these results should be interpreted with
> caution. Performance can vary based on the specific data, model
> specifications, number of cores and RAM available, and other aspects
> of computational environment. Github Actions likely use more limited
> resources than a high-end local machine. But, we use Github Actions to
> ensure that the benchmarks are run in a consistent environment so that
> individuals can use these benchmarks as a potential dev tool.

## Simulation DGP

The code below is used to generate the simulated data for the
benchmarks. The function creates a balanced panel of individuals over a
specified number of years. Additionally workers are assigned to firms.
There are three fixed effects: individual, firm, and year.

How the firms are assigned to individuals can be controlled by the
`type` argument. In the `simple` case, individuals are assigned to firms
randomly. This creates a very “dense” network of individuals and firms
so estimation is relatively fast.

In the `difficult` case, individuals are assigned to firms so as to
create a very “sparse” network. This makes estimation more difficult and
time-consuming. These can be seen as two extreme cases and where a
particular dataset may fall depends on the specific application.


<details closed>
  <summary>DGP Code</summary>
``` r
base_dgp <- function(
  n = 1000,
  nb_year = 10,
  nb_indiv_per_firm = 23,
  type = c("simple", "difficult")
) {
  nb_indiv = round(n / nb_year)
  nb_firm = round(nb_indiv / nb_indiv_per_firm)
  indiv_id = rep(1:nb_indiv, each = nb_year)
  year = rep(1:nb_year, times = nb_indiv)

  if (type == "simple") {
    firm_id = sample(1:nb_firm, n, TRUE)
  } else if (type == "difficult") {
    firm_id = rep(1:nb_firm, length.out = n)
  } else {
    stop("Unknown type of dgp")
  }

  x1 = rnorm(n)
  x2 = x1**2

  firm_fe = rnorm(nb_firm)[firm_id]
  unit_fe = rnorm(nb_indiv)[indiv_id]
  year_fe = rnorm(nb_year)[year]
  mu = 1 * x1 + 0.05 * x2 + firm_fe + unit_fe + year_fe

  df = data.frame(
    indiv_id = indiv_id,
    firm_id = firm_id,
    year = year,
    x1 = x1,
    x2 = x2,
    y = mu,
    negbin_y = MASS::rnegbin(exp(mu), theta = 0.5),
    binary_y = as.numeric(mu > 0),
    ln_y = log(abs(mu) + 1)
  )
  return(df)
}
```
</details>


### OLS Results

![OLS Benchmark Results](results/plot_ols.svg)

### Poisson Results

![Poisson Benchmark Results](results/plot_poisson.svg)

### Logistic Results

![Logistic Benchmark Results](results/plot_logit.svg)

### Real Data

<!-- Real Data -->
| Dataset | Num. obs. | Estimator | Mean Estimation Time |
|----|----|----|----|
| tradepolicy (OLS) | 28566 | pyfixest.feols | 0.262 |
| tradepolicy (OLS) | 28566 | FixedEffectModels.reg | 0.058 |
| tradepolicy (OLS) | 28566 | fixest::feols | 0.037 |
| tradepolicy (Poisson) | 28566 | pyfixest.fepois | 1.894 |
| tradepolicy (Poisson) | 28566 | GLFixedEffectModels Poisson | 4.689 |
| tradepolicy (Poisson) | 28566 | fixest::fepois | 0.243 |
| nycflights13 | 336776 | pyfixest.feols | 0.29 |
| nycflights13 | 336776 | FixedEffectModels.reg | 0.183 |
| nycflights13 | 336776 | fixest::feols | 0.095 |
| Medicare Provider | 9714896 | FixedEffectModels.reg | 13.208 |
| Medicare Provider | 9714896 | fixest::feols | 59.428 |
| nyc taxi | 46099576 | pyfixest.feols | 57.76 |
| nyc taxi | 46099576 | FixedEffectModels.reg | 24.5 |
| nyc taxi | 46099576 | fixest::feols | 43.61 |
<!-- Real Data -->
