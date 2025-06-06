# Load required packages

library(corrplot)
library(lubridate)
library(dplyr)
library(tidyr)
library(rasterVis)
library(viridis)
library(rangeBuilder)
library(sp)
library(raster)
library(blockCV)
library(terra)

# Give more Java memory prior to loading Dismo to run if needed
options(java.parameters = "-Xmx8g" )
library(dismo)

#Set working directory
setwd("C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024")

## STEP 1 - LOAD SPECIES DATA, AOI AND CRS TEMPLATE ######################################################
# Load species records data
species_data <- read.csv("Aggregated_records_2000_2023_ausalbers_yrmonth.csv")
#View(species_data)

# Filter for relevant data columns
species_data <- species_data[,c("ID", "PRESENCE", "OBS.YEAR", "OBS.MONTH", "PRECISION", "X_3577", "Y_3577")]
#View(species_data)

# Set project CRS e.g. EPSG 3577 using Proj4 string code
AA_proj <- CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Load AOI empty raster
template <- raster("./AOI/template_AOI2.tif")

## STEP 2 - SPATIAL RECORDS FILTERING ################################################

# Create a SpatialPoints (SP) object form presence points coords
coordinates <- species_data[, c("X_3577","Y_3577")]
dat_M_s <- SpatialPoints(coordinates)

# Combine SP object with the original data frame
dat_M_s <- SpatialPointsDataFrame(dat_M_s, species_data)

# Ensure AOI template and species data are in matching CRS: e.g. EPSG 3577 allocated to object 'AA_proj'
crs(template)<-AA_proj
crs(dat_M_s)<-AA_proj
plot(template, main = "All data")
points(dat_M_s)

names(dat_M_s@data)[names(dat_M_s@data)=="ID"] <- "ID" 

numbers<-as.numeric(dat_M_s@data$ID)
#print(numbers)

# Paste month and year of observation from species data into new data column month_yr
dat_M_s@data$month_yr<-paste(dat_M_s$OBS.MONTH, dat_M_s$OBS.YEAR)
print(dat_M_s@data$month_yr)

# Create empty output
output <- NULL

# Create object with all unique combinations of months & years
month_year<-unique(dat_M_s@data$month_yr) 

# Spatial filter: Loop through all records to filter out records within 100m of each other for specific month/yr
for (M in unique(dat_M_s@data$month_yr)){ #match unique month and year to M eg "09 2017"
  print(M)
  
  month_subset<-dat_M_s[dat_M_s@data$month_yr==M,] #then take all rows in that specific month
  print(month_subset)
  for (a in 1:nrow(month_subset)){ #then filtering records for that unique month 
    
    d <- spDists(month_subset, month_subset[a,], longlat=F) #get distances from all records in that month to point in row a
    d <- data.frame(d)
    
    print(min(d))
    
    d$ID <- month_subset$ID
    d<-d[(d$d<100) & d$ID != a,] #distance filter to 100m for records in same month/year
    
    
    d <- left_join(d, month_subset@data, by = "ID")[,c(1:2, 6)] #d, id and accuracy location. Choose most precise record if duplicates
    # original:    d <- left_join(d, Pres, by = "ID")[,c(1:2, 8)]
    ptuncert <- month_subset@data$PRECISION[a] #Precision of the point a
    ids<-d$ID
    
    if(nrow(d)>0){ # if there are points within 100m
      if((ptuncert == min(d$PRECISION, na.rm = T)) & # if ptuncert == min of points within 100m
         length(which(ids %in% output$ID))==0) { #And none of the points are already in the output set
        output<-rbind(output,month_subset@data[a,]) #keep point a in final dataset
      } #If two points in 100m proximity of each other collected in the same month + year, only the most accurate or first point will be kept
      
    }else{ # If there are no points within 100m from that time, keep point a
      
      output<-rbind(output,month_subset@data[a,])   
    }
  } 
}

# Check output
output
summary(output)

# Add coordinates to SPDF
output<-dat_M_s[dat_M_s@data$ID %in% output$ID,]

# Set output CRS
crs(output)<-AA_proj

# Plot filtered records
plot(template, main = "filtered spatial clumping")
points(output)

# Reset data object to output
dat_M_fin<-output
View(dat_M_fin@data)

## Export the dataframe as .csv
#write.csv(dat_M_fin@data, file = "filtered_tymp_records.csv", row.names = FALSE)

## STEP 3 - LOAD STATIC EV & EXTRACT VALUES AT POINTS ########################################

# Load static EV rasters
raster_paths <- c("./Enviro predictors/Vegetation/weighted_analysis_native_grasslands.tif",
                  "./Enviro predictors/Vegetation/weighted_analysis_non_native_veg.tif",
                  "./Enviro predictors/Land Use/CLUM_2023_weighted_rasters/Native_grazing.tif",
                  "./Enviro predictors/Land Use/CLUM_2023_weighted_rasters/Native_exot_pasture_mosaic.tif",
                  "./Enviro predictors/Land Use/CLUM_2023_weighted_rasters/Other_minimal_use.tif",
                  "./Enviro predictors/Vegetation/25_veg_height_AOI_100m.tif",
                  "./Enviro predictors/Soil/Nitrogen/Nitrogen_0_30cm.tif",
                  "./Enviro predictors/Static/DEM.tif",
                  "./Enviro predictors/Static/slope.tif",
                  "./Enviro predictors/Static/aspect.tif",
                  "./Enviro predictors/Bioclim_processed/bio_12.tif",
                  "./Enviro predictors/Bioclim_processed/bio_4.tif")

# Stack EV rasters
raster_stack <- raster::stack(raster_paths)

# Set the CRS for raster stack
projection(raster_stack) <- AA_proj

# Check info about the raster stack e.g. number of layers, extent, resolution and CRS
print(raster_stack)

# Rename layers
names(raster_stack)[1] <- c("Nat_grassland")
names(raster_stack)[2] <- c("Non_nat_veg")
names(raster_stack)[3] <- c("Native_grazing")
names(raster_stack)[4] <- c("Mosaic_pasture")
names(raster_stack)[5] <- c("Other_minimal_use")
names(raster_stack)[6] <- c("Veg_height")

# Visualise raster stack layers
plot(raster_stack)  

# Extract static covariates at record points ###################################
static.covs <- raster::extract(raster_stack, dat_M_fin)
# Combine with coords & year data from filtered records dataset
static.covs <- cbind(static.covs, dat_M_fin@data$OBS.YEAR, dat_M_fin@data$X_3577, dat_M_fin@data$Y_3577)
# Rename data columns
colnames(static.covs)[13] <- "Year"
colnames(static.covs)[14] <- "x"
colnames(static.covs)[15] <- "y"
# Check output
View(static.covs)

## STEP 4 - STACK ANNUAL DYNAMIC EVs 2000-2023/2020 & EXTRACT VALUES AT POINTS #################################

# Create new df object with filtered record points e.g. 'filtered Dat_M_fin'
points_df <- dat_M_fin@data

# Set the path to the folder containing dynamic raster files ("/Dynamic NDVI")
folder_path <- "./Enviro predictors/Dynamic NDVI/2000-2023"

# Retrieve all raster files in the folder
tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# Load all raster files into a SpatRaster stack
stacked_raster <- rast(tif_files)

# (If needed) remove the "X" prefix from the layer names in the raster stack
names(stacked_raster) <- gsub("^X", "", names(stacked_raster))
# Print the updated layer names
print("Updated layer names in raster stack:")
print(names(stacked_raster))

# Set unique years in the raster stack
raster_years <- as.numeric(names(stacked_raster))

# Subset points_df to keep only rows with years present in the raster stack
points_df_subset <- points_df[points_df$OBS.YEAR %in% raster_years, ]

# (If needed) Convert OBS.YEAR column to character type in points_df_subset
points_df_subset$OBS.YEAR <- as.character(points_df_subset$OBS.YEAR)

# Check data type matches for names and years in raster stack & points subset
raster_names_type <- class(names(stacked_raster))
print("Data type of names in raster stack:")
print(raster_names_type)
years_type <- class(points_df_subset$OBS.YEAR)
print("Data type of years in points_df_subset:")
print(years_type)

# Extract dynamic values at points #####################################

# Retrieve & display the names of the layers in the raster stack
layer_names <- names(stacked_raster)
print(layer_names)

# Determine the number of dynamic covariates
num_dynamic_covariates <- length(layer_names)
print(num_dynamic_covariates)

# Define the number of dynamic covariates e.g. 24 (one for each yr 2000-2023)
num_dynamic_covariates <- 24

# Initialize an empty matrix to store dynamic covariates
dynamic_cov_matrix <- matrix(NA, nrow = nrow(points_df_subset), ncol = num_dynamic_covariates)

# Extract values for each point based on the corresponding year
for (i in 1:nrow(points_df_subset)) {
  x <- points_df_subset$X_3577[i]
  y <- points_df_subset$Y_3577[i]
  year <- as.character(points_df_subset$OBS.YEAR[i])  # Convert to character to match layer names
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    # Store the extracted value in the matrix
    dynamic_cov_matrix[i, year_index] <- values[[1]]  # Assuming only one value is extracted
    
    # Print the extracted value
    print(paste("Value for point at (", x, ",", y, ") in year", year, ":", values[[1]]))
  }
}

# Convert the output matrix to a data frame
dynamic_cov_df <- as.data.frame(dynamic_cov_matrix)
print(dynamic_cov_df)

# Remove NAs
NDVI <- rowSums(dynamic_cov_df, na.rm = T)
print(NDVI)

covs <- static.covs
#combine static and dynamic covariate values
covs <- cbind(static.covs, NDVI)
View(covs)

# REPEAT FOR ALL DYNAMIC EVs #

### DYNAMIC EV 2 - SOIL MOISTURE ###############################################

# Set the path to the folder containing dynamic raster files ("/Dynamic soil moisture")
folder_path_2 <- "./Enviro predictors/Dynamic soil moisture/2000-2023"

# Define all raster files in the folder
sm_files <- list.files(folder_path_2, pattern = "\\.tif$", full.names = TRUE)

# Load all raster files into a unique SpatRaster stack
stacked_soilm_raster <- rast(sm_files)

# (optional) Visualise raster layers  
plot(stacked_soilm_raster)

# Define unique years in the raster stack
raster_years <- as.numeric(names(stacked_soilm_raster))

# Subset points_df to keep only rows with years present in the raster stack
points_df_subset <- points_df[points_df$OBS.YEAR %in% raster_years, ]

# Convert OBS.YEAR column to character type in points_df_subset
points_df_subset$OBS.YEAR <- as.character(points_df_subset$OBS.YEAR)

# Check data type of names in raster stack and years in subset records matches e.g. 'character'
raster_names_type <- class(names(stacked_soilm_raster))
print("Data type of names in raster stack:")
print(raster_names_type)
years_type <- class(points_df_subset$OBS.YEAR)
print("Data type of years in points_df_subset:")
print(years_type)

# Extract dynamic values at points #####################################

# Retrieve the names of the layers in the raster stack
layer_names <- names(stacked_soilm_raster)

# Display the layer names
print(layer_names)

# Determine the number of dynamic covariates
num_dynamic_covariates <- length(layer_names)
print(num_dynamic_covariates)

# Define the number of dynamic covariates
num_dynamic_covariates <- 24  # Should be same for all dynamic layers e.g. '24'

# Initialize a unique empty matrix to store dynamic covariates
dynamic_cov_matrix_sm <- matrix(NA, nrow = nrow(points_df_subset), ncol = num_dynamic_covariates)

# Extract values for each point based on the corresponding year
for (i in 1:nrow(points_df_subset)) {
  x <- points_df_subset$X_3577[i]
  y <- points_df_subset$Y_3577[i]
  year <- as.character(points_df_subset$OBS.YEAR[i])  # Convert to character to match layer names
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_soilm_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_soilm_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    # Store the extracted value in the matrix
    dynamic_cov_matrix_sm[i, year_index] <- values[[1]]  # Assuming only one value is extracted
    
    # Print the extracted value
    print(paste("Value for point at (", x, ",", y, ") in year", year, ":", values[[1]]))
  }
}

# Convert the matrix to a data frame
dynamic_cov_df_sm <- as.data.frame(dynamic_cov_matrix_sm)
# Remove NAs from other years
Soil_moisture <- rowSums(dynamic_cov_df_sm, na.rm = T)
print(Soil_moisture)

### DYNAMIC EV 3 - DEA LANDCOVER - Herbaceous >15% cover #################################

# Set the path to the folder containing dynamic raster files ("")
folder_path_3 <- "./Enviro predictors/Dynamic_DEA_weighted_av/Herbaceous_over_15"

#Check extents of all rasters match as necessary to create rasterstack:

# List all raster files in the folder
landcover_files <- list.files(folder_path_3, pattern = "\\.tif$", full.names = TRUE)

# IF NEEDED : Loop through each GeoTIFF file and print its extent
for (file in landcover_files) {
  raster_obj <- raster(file)
  file_extent <- extent(raster_obj)
  cat("File:", basename(file), "\n")
  cat("Extent:\n")
  cat("xmin:", xmin(file_extent), "\n")
  cat("xmax:", xmax(file_extent), "\n")
  cat("ymin:", ymin(file_extent), "\n")
  cat("ymax:", ymax(file_extent), "\n\n")
}

# Create a raster stack
stacked_landcover_raster <- terra::rast(landcover_files)

# Extract base filenames without extension
base_names <- tools::file_path_sans_ext(basename(landcover_files))

# Set the names of the raster stack to the base filenames
names(stacked_landcover_raster) <- base_names

# Check layer names
print(names(stacked_landcover_raster))

# Visualise raster layers in stack  
plot(stacked_landcover_raster)

# Get unique years in the raster stack
raster_years <- as.numeric(names(stacked_landcover_raster))

# Subset points_df to keep only rows with years present in the raster stack
points_df_subset <- points_df[points_df$OBS.YEAR %in% raster_years, ]

# Convert OBS.YEAR column to character type in points_df_subset
points_df_subset$OBS.YEAR <- as.character(points_df_subset$OBS.YEAR)

# Check data type of names in raster stack/years in subset records match
raster_names_type <- class(names(stacked_landcover_raster))
print("Data type of names in raster stack:")
print(raster_names_type)
years_type <- class(points_df_subset$OBS.YEAR)
print("Data type of years in points_df_subset:")
print(years_type)

#Extract dynamic values at points#####################################

# Retrieve the names of the layers in the raster stack
layer_names <- names(stacked_landcover_raster)

# Display the layer names
print(layer_names)

# Determine the number of dynamic covariates
num_dynamic_covariates <- length(layer_names)
print(num_dynamic_covariates)

# Define the number of dynamic covariates
num_dynamic_covariates <- 24

# Initialize an empty matrix to store dynamic covariates
dynamic_cov_matrix_landcover <- matrix(NA, nrow = nrow(points_df_subset), ncol = num_dynamic_covariates)

# Check raster extent
raster_extent <- ext(stacked_landcover_raster)
print(raster_extent)

# Loop over each point and extract values
for (i in 1:nrow(points_df_subset)) {
  x <- points_df_subset$X_3577[i]
  y <- points_df_subset$Y_3577[i]
  year <- as.character(points_df_subset$OBS.YEAR[i])
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_landcover_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_landcover_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    # Store a single extracted value in the matrix
    dynamic_cov_matrix_landcover[i, year_index] <- values[[1]]  # Choose the first value
  }
}

# Convert the matrix to a data frame
dynamic_cov_df_landcover <- as.data.frame(dynamic_cov_matrix_landcover)
# Remove NAs for other years
Herbaceous_open <- rowSums(dynamic_cov_df_landcover, na.rm = T)
print(Herbaceous_open)

# Combine static and dynamic covariate values
covs <- cbind(static.covs, NDVI, Soil_moisture, Herbaceous_open)
View(covs)

## (Optional) Export static and dynamic covs
#write.csv(covs, file = "presence_covs.csv")

## STEP 5 - CREATE RANDOM BACKGROUND POINTS AND ASSIGN YEARS ##################################################
# RBG points (no bias accounting)
# # set.seed(133) # If you set this then you will get the same 'random' points when you re-run the code
# rpoints<- randomPoints(template, 10000)
# rpoints<-as.data.frame(rpoints)
# table(is.na(rpoints))
# rpts2<- rpoints[1:10000,] 
# 
# plot(template, main = "random points") # Plot points to check out where they are
# points(rpts2, cex=0.3)
# 
# points_bias_model<-rpts2
# 
# #View(points_bias_model)
# 
# Tymp <- dat_M_fin@data
# #PW.ALL<-subset(PW.ALL, Survey.method!="Acoustic") #removing accoustic monitoring results
# summary(Tymp$OBS.YEAR)
# Year <-sample((Tymp$OBS.YEAR)[!is.na(Tymp$OBS.YEAR)], size=nrow(points_bias_model), replace = TRUE)
# points_bias_model<-as.data.frame(points_bias_model)
# points_bias_model<-cbind(points_bias_model, Year)
# 
# #View(points_bias_model)
# points_bias_model_dynamic_extract<-points_bias_model
# 
# #View(points_bias_model_dynamic_extract)
# 
# ## Print the number of points assigned to each year
# points_per_year <- table(format(points_bias_model_dynamic_extract$Year))
# print(points_per_year)
# 
# #TS- "Random" points are weighted significantly more heavily in 2019-2020 (>1000) against other years, 2005 lightly weighted. REVISIT
# 
# #E.G.
# # 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 
# #  425   86  126  179  121   26  298  226  152  200  465  110  302  366  475  574  798  664  640 1194 1131  902  211  329 
# 
# write.csv(points_bias_model_dynamic_extract, "RBG_points.csv")

### BIAS RBG POINTS #############################################################################
source("C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/Bias_model.R")
###REMEMBER TO RUN BIAS MODEL SCRIPT, OR IMPORT THE .CSV POINTS FROM THE FILE "C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/Bias_RBG_points.csv"
# Using bias SDM background points
# Ensure points_bias_model is properly generated in your script before using this section.
# 
# # # Replace the generation of random background points with bias SDM background points
#  rpts2 <- points_bias_model
#  
#  # Plotting is optional, you may remove this part if not needed
#  plot(template, main = "bias SDM background points") 
#  points(rpts2, cex=0.3)
#  
#  # Assigning years to background points
#  Tymp <- dat_M_fin@data
#  summary(Tymp$OBS.YEAR)
#  Year <- sample((Tymp$OBS.YEAR)[!is.na(Tymp$OBS.YEAR)], size = nrow(rpts2), replace = TRUE)
#  
#  # Combining background points with years
#  rpts2 <- cbind(rpts2, Year)
#  
#  # Optional: View the resulting data frame
#  # View(rpts2)
#  
#  # Print the number of points assigned to each year
#  points_per_year <- table(format(rpts2$Year))
#  print(points_per_year)
#  
#  # Write the bias SDM background points with assigned years to a CSV file
#  write.csv(rpts2, "C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/Bias_RBG_points_15e4.csv")

### EXTRACT STATIC COVARIATES FOR RBG POINTS ######################################  

# Create new object with RBG points (RBG_points.csv)
points_RBG <- read.csv("Bias_RBG_points_15e4.csv")
View(points_RBG)
coordinates(points_RBG)

# Extract static covs ###################################
static_RBG.covs <- raster::extract(raster_stack, coordinates(points_RBG)[,2:3])
static_RBG.covs <- cbind(static_RBG.covs, points_RBG[2:4]) #include years to enable plotting predictions
View(static_RBG.covs)

### EXTRACT DYNAMIC COVARIATES FOR RBG POINTS #####################################

####NDVI##########################

# Initialize an empty matrix to store dynamic covariates
dynamic_cov_bg_matrix <- matrix(NA, nrow = nrow(points_RBG), ncol = num_dynamic_covariates) #num_dynamic_covariates consistent with object already created (n=24)

# Extract values for each point based on the corresponding year
for (i in 1:nrow(points_RBG)) {
  x <- points_RBG$x[i]
  y <- points_RBG$y[i]
  year <- as.character(points_RBG$Year[i])  # Convert to character to match layer names
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    # Store the extracted value in the matrix
    dynamic_cov_bg_matrix[i, year_index] <- values[[1]]  # Assuming only one value is extracted
    
    # Print the extracted value
    print(paste("Value for point at (", x, ",", y, ") in year", year, ":", values[[1]]))
  }
}


# Convert the matrix to a data frame
dynamic_cov_bg_df_NDVI <- as.data.frame(dynamic_cov_bg_matrix)

print(dynamic_cov_bg_df_NDVI)

NDVI_RBG <- rowSums(dynamic_cov_bg_df_NDVI, na.rm = T)
print(NDVI_RBG)

## Save extracted dynamic NDVI RBG values to .csv to save time for future session
#  write.csv(NDVI, file = "C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/dynamic_cov_rbg.csv")

## Read in saved RBG dynamic NDVI values
# NDVI <- read.csv("C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/dynamic_cov_rbg.csv")

#### Soil Moisture #################

#Repeat steps above with Soil moisture raster stack
# Initialize an empty matrix to store dynamic covariates
dynamic_cov_bg_matrix_sm <- matrix(NA, nrow = nrow(points_RBG), ncol = num_dynamic_covariates) #num_dynamic_covariates consistent with object already created (n=24)

# Extract values for each point based on the corresponding year
for (i in 1:nrow(points_RBG)) {
  x <- points_RBG$x[i]
  y <- points_RBG$y[i]
  year <- as.character(points_RBG$Year[i])  # Convert to character to match layer names
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_soilm_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_soilm_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    # Store the extracted value in the matrix
    dynamic_cov_bg_matrix_sm[i, year_index] <- values[[1]]  # Assuming only one value is extracted
    
    # Print the extracted value
    print(paste("Value for point at (", x, ",", y, ") in year", year, ":", values[[1]]))
  }
}

# Convert the matrix to a data frame
dynamic_cov_bg_df_sm <- as.data.frame(dynamic_cov_bg_matrix_sm)

print(dynamic_cov_bg_df_sm)

Soil_moisture_RBG <- rowSums(dynamic_cov_bg_df_sm, na.rm = T)
print(Soil_moisture_RBG) 

#### Landcover - Herbaceous >15% cover ##############

dynamic_cov_bg_matrix_landcover <- matrix(NA, nrow = nrow(points_RBG), ncol = num_dynamic_covariates) #num_dynamic_covariates consistent with object already created (n=24)

# Extract values for each point based on the corresponding year
for (i in 1:nrow(points_RBG)) {
  x <- points_RBG$x[i]
  y <- points_RBG$y[i]
  year <- as.character(points_RBG$Year[i])  # Convert to character to match layer names
  
  # Extract raster values for the given point and year
  values <- terra::extract(stacked_landcover_raster, cbind(x, y))
  
  # Get the index for the specific year
  year_index <- which(names(stacked_landcover_raster) == year)
  
  # Check if a matching year was found
  if (length(year_index) == 0) {
    warning(paste("No matching year found for", year))
  } else {
    dynamic_cov_bg_matrix_landcover[i, year_index] <- values[[1]]
    print(paste("Value for point at (", x, ",", y, ") in year", year, ":", values[[1]]))
  }
}

# Convert the matrix to a data frame
dynamic_cov_bg_df_landcover <- as.data.frame(dynamic_cov_bg_matrix_landcover)

print(dynamic_cov_bg_df_landcover)

Herbaceous_open_RBG <- rowSums(dynamic_cov_bg_df_landcover, na.rm = T)
print(Herbaceous_open_RBG)

## Combine static and dynamic RBG values
RBG_covs <- cbind(static_RBG.covs, NDVI_RBG, Soil_moisture_RBG, Herbaceous_open_RBG)
#RBG_covs <- RBG_covs [,-9] #remove column 9 if auto assigned ID number
colnames(RBG_covs) [16] <- "NDVI" #rename column with RBG dynamic values/yr to "NDVI"/"Soil_moisture"/"Landcover" etc for consistency if necessary
colnames(RBG_covs) [17] <- "Soil_moisture"
colnames(RBG_covs) [18] <- "Herbaceous_open"

View(RBG_covs)

# (optional) Export RBG covs data for future running of model
#write.csv(RBG_covs, file = "RBG_covs.csv")

## STEP 6 - RUN TEMPORAL BLOCK CROSS-VALIDATION TEST ###############

# # (Either) Save time by loading in presaved presence/covs data 
# mod_pres <- read.csv("presence_covs.csv")
# mod_pres <- mod_pres[,-1]

# (Or) Define records covariates as model presence data
mod_pres <- covs
# Create new column (at start of table by default) called "presence" and assign all rows value of 1 (presence record)
mod_pres <- data.frame("presence" = 1, mod_pres)
colnames(mod_pres)
View(mod_pres)

# # (Either) Save time by loading in presaved pseudo-absence/covs data
# mod_abs <- read.csv("RBG_covs.csv")
# mod_abs <- mod_abs[,-1]

# (Or) Define RBG points covariates as pseudo-absence data
mod_abs <- RBG_covs
# Create new column called "presence" and assign all rows value of 0 (non-presence record)
mod_abs <- data.frame("presence" = 0, mod_abs)
View(mod_abs)

# sp_b<-c(rep(1, length(mod_pres[,1])),rep(0,length(mod_bg_b[,1]))) # Code presence as '1' and background as '0'
mod_env_a <-as.data.frame(rbind(mod_pres, mod_abs)) #binding pres and abs together
View(mod_env_a)

sp_b <- mod_env_a[,1]
View(sp_b)

mod_env_b <- mod_env_a[,-1]
View(mod_env_b)

start_yrs<-seq(2000,2023,by=2) #Blocks - leave out these + year after
all_yrs<-seq(2000,2023,1)

k<-length(start_yrs)

blockCV_M <- vector('list', k)

# Training and test data segregation 
###NB -Before executing the training/test data segregation, ensure the working directory is set to the right pathway containing the "Maxent_outputs" folder containing "Species.lambdas" file where it will store the outputs

for (i in 1:k) {
  print(start_yrs[i])
  test_yrs<-c(start_yrs[i],start_yrs[i] + 1)
  
  trainSet_M <- which(!(mod_env_b$Year %in% test_yrs)) # extract the training set indices
  testSet_M <- which(mod_env_b$Year %in% test_yrs) # extract the testing set indices
  table(sp_b[testSet_M], useNA = "always") # added 2025
  
  print(paste0("training data: ",(length(trainSet_M)))) #how many training points
  print(paste0("test data: ",length(testSet_M)) ) #how many test points
  
  colnames(mod_env_b)
  
  #blockcv_mod_maxent_full<-maxent(x=mod_env_b[trainSet_M,-1],#NB edit - exclude the first column, which has the year info. Must edit Also exclude columns with x and y values
  blockcv_mod_maxent_full<-maxent(x=mod_env_b[trainSet_M, c(1:12, 16:18)],#need to modify this line to properly exclude 'x', 'y' and 'Year' columns
                                  p=sp_b[trainSet_M],
                                  args=c("noautofeature", "nothreshold","beta_hinge=1.5","-J","-P"),
                                  removeDuplicates=FALSE,
                                  path = "Maxent_outputs")
  
  #testdat.env_M <- mod_env_b[testSet_M,-1]
  testdat.env_M <- mod_env_b[testSet_M, c(1:12, 16:18)] #need to modify this line to properly exclude 'x', 'y' and 'Year' columns
  testdat.sp_M <- sp_b[testSet_M]
  locust.test.pres_M <- testdat.env_M[testdat.sp_M==1,]
  locust.test.bg_M  <- testdat.env_M[testdat.sp_M==0,]
  
  # print(summary(sp_b[testSet_M]))
  if (nrow(locust.test.pres_M) == 0) {
    warning(paste0("Skipping fold ", i, ": no test presence points."))
    next
  }
  if (nrow(locust.test.bg_M) == 0) {
    warning(paste0("Skipping fold ", i, ": no test background points."))
    next
  }
  if (any(is.na(locust.test.pres_M))) {
    warning(paste0("Skipping fold ", i, ": test presence has NA values."))
    next
  }
  if (any(is.na(locust.test.bg_M))) {
    warning(paste0("Skipping fold ", i, ": test background has NA values."))
    next
  }
  
  # evaluate model with held-out data
  blockCV_M[[i]] <- evaluate(p=locust.test.pres_M,
                             a=locust.test.bg_M,
                             model=blockcv_mod_maxent_full)
}


AUCs_BiasBG <- sapply(blockCV_M, slot, "auc")
round(mean(AUCs_BiasBG),3)

round(sd(AUCs_BiasBG),2)

round(range(AUCs_BiasBG),3)


## CALCULATE VARIABLE IMPORTANCE OF ENVIRONMENTAL PREDICTORS###############################

model <- blockcv_mod_maxent_full
var_imp <- plot(model)
var_imp

## CORRELATION TESTING BETWEEN ENVIRONMENTAL COVARIATES IN TEST MODEL####################

# Remove rows with NA values (analysis won't work with missing data)
cleaned_data <- na.omit(mod_env_b)

# Calculate the correlation matrix
correlation_matrix <- cor(cleaned_data, method = "spearman") #available methods "pearson" (default), "spearman", "kendall" and "partial" in 'cor' function

# Print the correlation matrix
print(correlation_matrix)

## Outwrite correlation as table
# correlation_df <- as.data.frame(correlation_matrix, row.names = rownames(correlation_matrix), col.names = colnames(correlation_matrix))
# write.csv(correlation_df, file = "./Enviro_variables_analysis/spearmans_cor_results_21_05.csv")

## STEP 7 - RUN FULL MAXENT MODEL #########################################

# Train the MaxEnt model using all available data
full_maxent_model <- maxent(x = mod_env_b[, c(1:12, 16:18)],  # Exclude columns 8 to 10 (Year, x, y). Not necessary to exclude the coords but no value added to include
                            p = sp_b,
                            args = c("noautofeature", "nothreshold", "beta_hinge=1.5", "-J", "-P"),
                            removeDuplicates = FALSE,
                            path = "Maxent_outputs")  

# Subset data for the year you want to predict
predict_year <- 2023
# Check data table
View(mod_env_b)
# Define covariates data columns only (exclude year and coords)
predict_data <- mod_env_b[mod_env_b$Year == predict_year, c(1:12, 16:18)]

# Predict habitat suitability using the full MaxEnt model
prediction2023 <- dismo::predict(full_maxent_model, x = predict_data)

#Resample 2023 NDVI (250m) and DEA Landcover rasters (25m) to 100m to stack consistently with static predictors
ndvi250 <- raster(stacked_raster[["2023"]])
ndvi100 <- projectRaster(ndvi250, raster_stack[[1]]) #Resample 2023 layer (change depending on year of prediction) of NDVI raster stack to same res as static environmental variables ie raster_stack (layer 1) 
Herbaceous_open25 <- raster(stacked_landcover_raster[["2023"]])
Herbaceous_open100 <- projectRaster(Herbaceous_open25, raster_stack[[1]])
#create object for 2023 layer in soil moisture & rainfall stack
soil_moisture <- raster(stacked_soilm_raster[["2023"]])

prediction2023 <- stack(raster_stack, ndvi100, soil_moisture, Herbaceous_open100)
plot(prediction2023)
print(names(prediction2023))
names(prediction2023)[13] <- c("NDVI") #change name of dynamic layer prediction from default (eg X2023)
names(prediction2023)[14] <- c("Soil_moisture")
names(prediction2023)[15] <- c("Herbaceous_open")

prediction2023 <- dismo::predict(full_maxent_model, x = prediction2023)

#Visualise prediction
plot(prediction2023)

# Save prediction (outwrite raster)
writeRaster(prediction2023, filename = "C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/SDM/MaxEnt/DS_session_15_01_2024/MaxEnt Predictions/Final_model_preds_2/2023.tif", format = "GTiff", options = "COMPRESS=LZW", overwrite = TRUE)
