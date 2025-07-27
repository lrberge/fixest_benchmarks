library(lfe)

lfe_timer <- function(data, fml) {
  start_time <- Sys.time()
  result <- felm(fml, data = data)
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
