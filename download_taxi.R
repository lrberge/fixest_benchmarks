options(timeout = 300)
bucket = "https://ursa-labs-taxi-data.s3.us-east-2.amazonaws.com"
for (m in sprintf("%02d", 1:3)) {
  fp = file.path(paste0("data/nyc-taxi/year=2012/month=", m))
  dir.create(fp, recursive = TRUE)
  try(
    download.file(
      paste(bucket, "2012", m, "data.parquet", sep = "/"),
      file.path(fp, "data.parquet"),
      mode = "wb"
    ),
    silent = TRUE
  )
}


library(duckdb)
con = dbConnect(duckdb())

nyc = dbGetQuery(
  con,
  "
   FROM 'data/nyc-taxi/**/*.parquet'
   SELECT
      tip_amount, trip_distance, passenger_count,
      vendor_id, payment_type, dropoff_at,
      dayofweek(dropoff_at) AS dofw
   "
)

dbDisconnect(con)
rm(con)

write_parquet(nyc, "data/nyc-taxi.parquet")
unlink("data/nyc-taxi", recursive = TRUE)
