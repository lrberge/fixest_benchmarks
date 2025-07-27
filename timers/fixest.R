library(fixest)

feols_timer <- function(data, fml, vcov = "iid") {
  start_time <- Sys.time()
  result <- feols(fml, data = data, vcov = vcov, notes = FALSE, warn = FALSE)
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
fepois_timer <- function(data, fml, vcov = "iid") {
  start_time <- Sys.time()
  result <- fixest::fepois(
    fml,
    data = data,
    vcov = vcov,
    notes = FALSE,
    warn = FALSE
  )
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
feglm_logit_timer <- function(data, fml, vcov = "iid") {
  start_time <- Sys.time()
  result <- fixest::feglm(
    fml,
    data = data,
    family = "logit",
    vcov = vcov,
    notes = FALSE,
    warn = FALSE
  )
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
feols_multiple_vcov_timer <- function(data, fml, cluster) {
  start_time <- Sys.time()
  result <- feols(fml, data = data, vcov = "hc1", notes = FALSE, warn = FALSE)
  result <- summary(result, cluster = cluster)
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
