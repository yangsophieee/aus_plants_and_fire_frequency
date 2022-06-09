# https://github.com/traitecoevo/aus_gbif_clean/blob/master/gbif_process.R

library(data.table)
library(curl)
library(zip)
library(tidyverse)
library(CoordinateCleaner)

# GBIF.org (12 August 2021) GBIF Occurrence Download [https://doi.org/10.15468/dl.6rghq7]

curl_download("https://api.gbif.org/v1/occurrence/download/request/0344358-200613084148143.zip",
              "raw_data/gbif.zip",
              quiet = FALSE)
unzip("raw_data/gbif.zip", exdir="raw_data/gbif")

#aus_gbif<-fread("raw_data/gbif/0344358-200613084148143.csv",nrows = 10)

aus_gbif <- fread("raw_data/gbif/0344358-200613084148143.csv", quote = "")

# modified from: https://data-blog.gbif.org/post/gbif-filtering-guide/
data.frame(aus_gbif) %>%
    setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
    filter(occurrencestatus  == "PRESENT")  %>%
    filter(!is.na(decimallongitude)) %>% 
    filter(!is.na(decimallatitude)) %>% 
    filter(establishmentmeans %in% c("INTRODUCED", "INVASIVE", "NATURALISED","")) %>% # Exclude "MANAGED"
    filter(year >= 1900) %>%
    filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
    filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>%
    filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
    filter(!decimallatitude == 0 | !decimallongitude == 0) ->intermediate_check

#adding some additional filters
intermediate_check %>%
filter(!grepl("COUNTRY_COORDINATE_MISMATCH",intermediate_check$issue)) ->i2

i2 %>%
    cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
    cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
    cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
 #   cc_sea() %>% # remove from ocean # doesn't work?
    distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) -> aus_filt

write_csv(aus_filt,"filtered_aus_obs.csv")

select(aus_filt,species,decimalLongitude=decimallongitude,decimalLatitude=decimallatitude,scientificname) %>%
    write_csv("filt_aus_limited_columns_testing.csv")

