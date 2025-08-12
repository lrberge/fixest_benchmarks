temp_zip <- tempfile(fileext = ".zip")
download.file(
  "https://downloads.cms.gov/files/Medicare-Physician-and-Other-Supplier-PUF.zip",
  temp_zip,
  mode = "wb"
)
unzip(
  temp_zip,
  files = "Medicare_Provider_Util_Payment_PUF_CY2016.txt",
  exdir = "data/"
)
unlink(temp_zip)
