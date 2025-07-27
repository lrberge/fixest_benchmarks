library(alpaca)

alpaca_poisson_timer <- function(data, fml, vcov = "iid") {
  start_time <- Sys.time()
  result <- alpaca::feglm(fml, data = data, family = poisson())
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}

alpaca_feglm_logit_timer <- function(data, fml, vcov = "iid") {
  start_time <- Sys.time()
  result <- alpaca::feglm(fml, data = data, family = binomial())
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
