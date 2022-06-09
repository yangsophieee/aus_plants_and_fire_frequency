

# The Australian Plant Census is available for download here: 
# https://biodiversity.org.au/nsl/services/export/index

apc <- read_csv("raw_data/APC-taxon-2021-06-08-0734.csv")
#apc <- read_csv("raw_data/taxon_list.csv")
apcs <- 
    apc %>% 
    filter(taxonRank == "Species" |
               taxonRank == "Forma" |
               taxonRank == "Varietas" |
               taxonRank == "Subspecies", 
           taxonomicStatus == "accepted")

apcs$nat_count <- str_count(apcs$taxonDistribution, "naturalised")
apcs$comma_count <- str_count(apcs$taxonDistribution, ",")
plot(apcs$nat_count, apcs$comma_count)

apcs$native_somewhere_index <- apcs$nat_count <= apcs$comma_count


hmm <- which(apcs$nat_count <= apcs$comma_count)
apcs$taxonDistribution[hmm][1]
apcs$nat_count[hmm][1]
apcs$comma_count[hmm][1]


sum(apcs$native_somewhere_index, na.rm=T)
dim(apcs)

apcs %>%
    mutate(aus_native = apcs$native_somewhere_index) %>%
    select(taxon_name, aus_native) %>%
    write_csv("processed_data/aus_native_lookup.csv")

