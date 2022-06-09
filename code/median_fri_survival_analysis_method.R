
# Load libraries
library(parallel)
library(raster)
library(tidyverse)
library(lubridate)
library(survival)
library(ciTools)


# Read in fire response and gbif data
fire_response_and_gbif_data <- readRDS(file = "processed_data/fire_response_and_gbif_data.rds")
# From "code/extracting_fire_response_data.Rmd"

# Get list of raster file names from MODIS data
list_of_files <- 
    list.files("raw_data/MODIS_burn_data_australia/Win20", recursive = TRUE, full.names = TRUE) %>% 
    str_subset("burndate") # Filter for the burndate layers not the QA layers 

# A list of geoTIFFs for each year and month
# See pg 15 of the documentation/user guide for MCD64:
#     https://lpdaac.usgs.gov/documents/875/MCD64_User_Guide_V6.pdf 
# Instructions for obtaining the data begins on pg 22.


# Function to find intersecting pixels and extract burn dates
extract_burn_dates <- function(file_path, gbif_coordinates){
    
    # file_path argument must be a string
    raster <- raster(file_path)
    
    # Extract year from file_path
    year <- str_match(file_path, "(?<=.A)[0-9]{4}")
    
    # Extract intersecting pixels
    intersecting_pixels <- raster %>% 
        raster::extract(gbif_coordinates, sp = TRUE)
    
    # Convert to data frame
    df <- 
        intersecting_pixels %>% 
        raster::as.data.frame() %>% 
        rename(burn_date = contains("burndate")) %>% 
        select(decimalLongitude, decimalLatitude, burn_date) %>% 
        # Filter out -2 and 0 values
        filter(burn_date > 0) %>% 
        # Convert Julian date to calendar date
        mutate(burn_date = as.Date(burn_date - 1, origin = str_glue("{year}-01-01"))) 
    
    return(df)
}


# Function for extracting fire return intervals for each coordinate 
calc_fri <- function(burn_dates){
    
    # Add the start and end of the study period as dates
    start_data <- as.Date("2000-10-31")
    end_data <- as.Date("2021-07-01")
    fire_dates <- c(start_data, burn_dates, end_data)
    
    # Calculate fire return intervals
    fire_intervals <- fire_dates - lag(fire_dates)
    
    # Convert to data frame
    fire_intervals <- data.frame(fire_intervals) %>% drop_na()
    
    # Annotate first and last intervals as open with another column
    fire_intervals <-
        fire_intervals %>% 
        mutate(uncensored = case_when(row_number() == 1 | row_number() == n() ~ 0,
                                      TRUE ~ 1))
    
    # Calculate mean fire return interval
    return(fire_intervals)
    
}


# Function for calculating median FRI for a taxon
calculate_median_fri <- function(taxon, data = fire_response_and_gbif_data, list_of_rasters = list_of_files){
    
    # Subset data to taxon
    species_subset <-
        data %>% 
        filter(taxon_name == taxon) %>% 
        select(decimalLongitude, decimalLatitude)
    
    # Round species_subset to resolution of 0.005 decimal degrees so that any coordinates
    # that are approximately within the same MODIS pixel can be removed
    multiple <- 0.005
    species_subset <-
        species_subset %>% 
        mutate(decimalLongitude = multiple*round(decimalLongitude/multiple),
               decimalLatitude = multiple*round(decimalLatitude/multiple)) %>% 
        distinct(.keep_all = TRUE) # Remove duplicates
    
    # Convert to a SpatialPointsDataFrame
    spdf <- SpatialPointsDataFrame(coords = species_subset,
                                   data = species_subset, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
    
    # Extract burn dates for each gbif coordinate
    df <- 
        map_dfr(list_of_rasters, extract_burn_dates, spdf) %>% 
        arrange(decimalLongitude, decimalLatitude, burn_date)
    
    # Find fire return intervals
    fri_df <- 
        df %>% 
        group_by(decimalLongitude, decimalLatitude) %>% 
        summarise(fri = calc_fri(burn_date), .groups = "keep") %>% 
        unpack(cols = fri)
    
    # Join to original gbif coordinates to keep coordinates with no burns
    full_fri_df <-
        species_subset %>% 
        left_join(fri_df, by = c("decimalLongitude", "decimalLatitude"))
    
    # Calculate time interval of the whole MODIS period
    start_data <- as.Date("2000-11-01")
    end_data <- as.Date("2021-07-01")
    modis_period <- end_data - start_data # 7547 days
    
    # Add open intervals for pixels that did not burn in MODIS period
    full_fri_df <-
        full_fri_df %>% 
        replace_na(list(fire_intervals = modis_period, uncensored = 0))
    
    # Make function to capture warning messages 
    quiet_survreg <- quietly(survreg)
    
    # Run survival analysis on interval data including unburnt pixels
    survreg_with_unburnt <- 
        quiet_survreg(Surv(fire_intervals, uncensored) ~ 1, 
                      data = full_fri_df, dist = "weibull")
    
    # Calculate median fire return interval
    med_with_unburnt <-
        (exp(survreg_with_unburnt$result$coefficients)*log(2)^(survreg_with_unburnt$result$scale)) %>% 
        unname()
    
    # Find upper and lower bounds
    med_w_quantiles <- 
        add_quantile(full_fri_df, 
                     survreg_with_unburnt$result, 
                     p = 0.5, confint = TRUE, alpha = 0.05) 
    med_w_confint <- med_w_quantiles[1,c("lcb", "ucb")] %>% unlist() %>% as.vector()
    
    # Run survival analysis on interval data not including unburnt pixels
    survreg_without_unburnt <- 
        quiet_survreg(Surv(fire_intervals, uncensored) ~ 1, 
                      data = fri_df, dist = "weibull")
    
    # Calculate median fire return interval
    med_without_unburnt <-
        (exp(survreg_without_unburnt$result$coefficients)*log(2)^(survreg_without_unburnt$result$scale)) %>% 
        unname()
    
    # Find upper and lower bounds
    med_wo_quantiles <-
        add_quantile(fri_df,
                     survreg_without_unburnt$result,
                     p = 0.5, confint = TRUE, alpha = 0.05)
    med_wo_confint <- med_wo_quantiles[1,c("lcb", "ucb")] %>% unlist() %>% as.vector()
    
    # Return outputs 
    output_df <- data.frame(taxon_name = taxon,
                            med_w = med_with_unburnt,
                            med_w_lwr = med_w_confint[1],
                            med_w_upr = med_w_confint[2],
                            med_w_warnings = ifelse(length(survreg_with_unburnt$warnings) >= 1, 
                                                    paste(survreg_with_unburnt$warnings, 
                                                          collapse = ", "), NA),
                            med_w_messages = ifelse(length(survreg_with_unburnt$messages) >= 1, 
                                                    paste(survreg_with_unburnt$messages, 
                                                          collapse = ", "), NA),
                            med_wo = med_without_unburnt,
                            med_wo_lwr = med_wo_confint[1],
                            med_wo_upr = med_wo_confint[2],
                            med_wo_warnings = ifelse(length(survreg_without_unburnt$warnings) >= 1, 
                                                     paste(survreg_without_unburnt$warnings, 
                                                           collapse = ", "), NA),
                            med_wo_messages = ifelse(length(survreg_without_unburnt$messages) >= 1, 
                                                     paste(survreg_without_unburnt$messages, 
                                                           collapse = ", "), NA),
                            num_unburnt_pixels = nrow(full_fri_df) - nrow(fri_df),
                            num_pixels = nrow(species_subset))
    return(output_df)
    
}


# Extract list of unique taxa 
list_of_taxa <- fire_response_and_gbif_data$taxon_name %>% unique()


# Subset to portion of taxa if needed (process in batches)
#list_of_taxa_1600 <- list_of_taxa[1:1600]


# Find number of cores
num_cores <- detectCores()


# Calculate median FRI for each
median_fri_list <- 
    mclapply(list_of_taxa, 
             calculate_median_fri, 
             mc.cores = num_cores)


# Bind into data frame
median_fri_df <- median_fri_list %>% bind_rows()


# Write to a csv file
#write_csv(median_fri_df, "processed_data/median_fris_from_survival_method.csv")



