
# https://github.com/traitecoevo/aus_gbif_clean/blob/master/make_taxonomy_lookup.R

library(data.table)
library(tidyverse)


#aus_gbif<-fread("raw_data/gbif/0344358-200613084148143.csv",nrows = 1000)
apc_native<-read_csv("processed_data/aus_native_lookup.csv")

#all gbif records for Australia
aus_gbif <-
    fread("raw_data/gbif/0011392-210819072339941.csv", quote = "")

#get entries which 1) started as an APC accepted native species and 2) gbif changed

aus_gbif %>%
    mutate(verbatimSpeciesName=word(verbatimScientificName,1,2)) %>%
    filter(verbatimSpeciesName %in% apc_native$canonicalName) -> apc_starting_name

#only records which changed
apc_starting_name %>%
    filter(species!="") %>%
    dplyr::filter(verbatimSpeciesName!=species) %>%
    select(fullAPCname=verbatimScientificName,binomialAPCname=verbatimSpeciesName,binomialPOWOname=species,POWOinfra=infraspecificEpithet,taxonRank) %>%
    distinct() -> change_log

write_csv(change_log,"processed_data/powo_apc_lookup.csv")

