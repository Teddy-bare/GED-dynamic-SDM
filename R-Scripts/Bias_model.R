# Bias Maxent model to create RBG points

# Required packages
library(raster)
library(dismo)

setwd("C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/TS HONOURS THESIS/SDM")

# Load AOI if not embedded in main SDM script
template <- raster("./AOI/template_AOI.tif")
plot(template)

# Load filtered records
bias_tymp_filtered<-read.csv("./filtered_tymp_records.csv")
View(bias_tymp_filtered)

Tymp_points_bias <- bias_tymp_filtered[, c("X_3577","Y_3577")]

### Load KDE rasters ##################################################

Tymp_kde_pts_h15e4 <- raster("./Bias_model/kde_h15e4_points.tif")

plot(Tymp_kde_pts_h15e4)

# Bias model setup ######################

# Bias model random background points
B_rpoints<- randomPoints(template, 10000) #Generate 10,000 rb points within the AOI (template)
B_test<-raster::extract(template, B_rpoints) #Extract values at rb points
B_ok<-which(!is.na(B_test)) #Identify rb points with missing values
B_rpoints2<-B_rpoints[B_ok,] #Create subset (B_rpoints2) excluding rb points with missing values
B_rpts2<- B_rpoints2[1:10000,] #Create a subset of 5,000 valid rb points 
plot(template) #Plot the AOI
points(B_rpts2, cex=0.3) #Overlay 5,000 rb points 
B_rpts2<-as.data.frame(B_rpts2) #convert rb points object to data frame for extracting enviro values to points
coordinates(B_rpts2)<-~ x + y #assign coordinates to rb points

# Input KDE as EV
bias_env_inputs<- stack(Tymp_kde_pts_h15e4)
plot(bias_env_inputs)

# Extract KDE to points
values_bias <- raster::extract(bias_env_inputs, Tymp_points_bias)
values_bg_bias<- raster::extract(bias_env_inputs, B_rpts2)

# Full input df
sp_bg<-c(rep(1, length(values_bias[,1])),rep(0,length(values_bg_bias[,1]))) #Assign presence '1' and background/pseudo absence '0'
bias_input_df<-as.data.frame(rbind(values_bias,values_bg_bias)) #Combine presence and pseudo absence data into single frame

#colnames(bias_input_df) #Check columns match layer names

# Run model ######################################

bias_model <- maxent(x=bias_input_df, 
                     p=sp_bg, 
                     path="Maxent_outputs/bias", 
                     args=c('-P', "-J"), 
                     removeDuplicates = FALSE, 
                     progress='text')

response(bias_model) #Plot response curves for each variable in the model

pred_bias_model <- predict(bias_model, bias_env_inputs, filename="./Maxent_outputs/bias/Bias_model_kde15e4.tif", overwrite = TRUE)
plot(pred_bias_model)

# Generate RBG points from bias model ##############################

# Invert the probability values
inverted_prob <- 1 - values(pred_bias_model)

# Normalize the inverted probabilities
normalized_prob <- inverted_prob / sum(inverted_prob, na.rm = TRUE)

# Sample points based on the inverted probabilities
points_bias_model <- as.data.frame(xyFromCell(pred_bias_model, 
                                              sample(which(!is.na(values(pred_bias_model))),
                                                     size = 10000,
                                                     prob = normalized_prob[!is.na(normalized_prob)],
                                                     replace = FALSE)))

# Plot bias model RBG points ######################
plot(template)
points(points_bias_model, cex=0.3)

#write.csv(points_bias_model, "C:/Users/t3dst/OneDrive - The University Of Newcastle/HONOURS/TS HONOURS THESIS/SDM/Bias_model/Bias_model_points_kde15e4_pts.csv")
# # This is the set of background (obtained from the bias model) that you will use for the actual model
