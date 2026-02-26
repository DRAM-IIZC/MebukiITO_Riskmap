
# Sample R code for Spatial distribution of ticks and tick-borne pathogens in central Hokkaido, Japan and associated ecological factors revealed by intensive short-term survey in 2024
# Note: Due to privacy policies regarding sampling coordinates, 
# the original dataset is provided in the manuscript in anonymized format.
# Please replace 'data' with your own dataframe
# The script was developed using R version 4.4.1. 

###################################################################################
# Required packages
library(pROC)
library(MuMIn)
library(dismo)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(mgcv)   
library(caret)
library(performance)

#################################################################################################
#Global and local spatial statistics
#################################################################################################
# simulate spatial weight depends on distance decay parameter (r)
dis <- seq(0, 200, 0.001)
difun <- data.frame( dis = seq(0, 200, 0.001), f1 = exp(-dis/0.1), f2 = exp(-dis), f3 = exp(-dis/10), f4 = exp(-dis/20)) # r = 0.1, 1, 10, 20

ggplot(difun) + geom_line(aes(x = dis, y = f1), col = "#FFC1CA", linewidth = 2)+
                geom_line(aes(x = dis, y = f2), col = "#78A6BD", linewidth = 2)+
                geom_line(aes(x = dis, y = f3), col = "#FDEECC", linewidth = 2)+
                geom_line(aes(x = dis, y = f4), col = "#94654F", linewidth = 2)+  #C9E6F5
                scale_x_continuous(limits = c(0, 100)) +
                scale_y_continuous(limits = c(0, 1)) +
                theme_bw()+
                labs(x="d", y="w", title = "Simulated spatial weight depend on distance decay parameter" )+
                theme(axis.text = element_text(size = 15, family="Times"), axis.title = element_text(size =25, family="Times", face = "italic"), title = element_text(size =25, family="Times"))


it <- 29

moran_p <- data.frame(species = rep(NA, it))
moran_e <- data.frame(species = rep(NA, it))

# sampling site coordinates
coords <- matrix(nrow = nrow(data), ncol = 2)
coords[,1] <- data$Longitude
coords[,2] <- data$Latitude

maxdis <- max(dist(coords))
nb　<-　dnearneigh(coords, d1=0,d2= maxdis) # with all the others
glist　<-　nbdists(nb,coords)

# calculate moran's I and p-value with r = 0.1 to r = 20
#r=0.1~1

for (i in 1:10) {

ik <- i * 100 #m to km

glist_i　<-　lapply(glist,function(x) exp(-x/ik)) # w(i,j)=exp(-d(i,j)/jk)
w　<-　nb2listw(nb,glist=glist_i) 

moran <-moran.test(data$species,listw=w)
moran_p$species[i] <- moran $p
moran_e$species[i] <- moran $estimate[1]

}

#r=2~20

for (j in 2:20) {

jk <- j * 1000 #m to km

glist_j　<-　lapply(glist,function(x) exp(-x/jk)) # w(i,j)=exp(-d(i,j)/jk)
w　<-　nb2listw(nb,glist=glist_j)

moran <-moran.test(data$species,listw=w)

moran_p$species[j+9] <- moran $p
moran_e$species[j+9] <- moran $estimate[1]

}


############################################################################################
#G spatial statistics

maxdis <- max(dist(coords))
nb　<-　dnearneigh(coords, d1=0,d2= maxdis) 
glist　<-　nbdists(nb,coords)
glist_1　<-　lapply(glist,function(x) exp(-x/1000)) # w(i,j)=exp(-d(i,j)/1), km
w　<-　nb2listw(nb,glist= glist_1) 

local_g_species <- localG_perm(data$species, w, zero.policy = TRUE)
data$gi_z_species <- as.vector(local_g_species) 

data <- data %>% 
  mutate(hotspot_category = case_when(
    gi_z_species > 2.58  ~ "Hotspot (99% CI)",
    gi_z_species > 1.96  ~ "Hotspot (95% CI)",
    gi_z_species < -2.58 ~ "Coldspot (99% CI)",
    gi_z_species < -1.96 ~ "Coldspot (95% CI)",
    TRUE              ~ "Not Significant"
  ))


map <- ggplot() +
          geom_sf(data = data, fill = "white", color = "white") + 
          geom_point(data = data, aes(x = Longitude , y = Latitude, color = hotspot_category, size = 5), alpha = 0.8) +   
          scale_color_manual(values = c(
            "Hotspot (99% CI)" = "#d73027",
            "Hotspot (95% CI)" = "#fc8d59",
            "Not Significant" = "grey70",
            "Coldspot (95% CI)" = "#91bfdb",
            "Coldspot (99% CI)" = "#4575b4"
          )) +
          labs(title = "Species name", color = "Category") +
           theme_void() + 
          theme(panel.background = element_rect(fill = "grey90"), legend.position = "none", text = element_text(family="Arial", size = 8, face = "italic")) 

##################################################################################
# Ecological niche modeling
###################################################################################
# Screening potential important variables

model <- glm( species ~ forest+grass+wet+ crop + build+ water + snow+ prec6  + sqprec6 + temp6 + sqtemp6 + suns+ cold + elev+ sqelev + angl+ sqangl + lon + lat + deer+bear+raccoon + tanuki + fox , family="binomial" , data= data) #formula is an example

options(na.action = "na.fail")
selection <- dredge (model, rank="AIC")
comp_models <- get.models(selection,  subset = delta < 2) #dAIC < 2

###################################################################################
#GLM

#data frame for results
data$predicted_value <- -1000

#LOOCV
for(i in 1: nrow(data)){

lo_model <- glm( species ~ forest+ snow + angl+ cold , family="binomial" , data= data[-i,]) #formula is an example

data[i,]$predicted_value <- predict(lo_model, type="response", newdata= data[i,]) 

}

p <- subset(data$predicted_value, data$species == 1) # 1 represents presence
a <- subset(data$predicted_value, data$species == 0) # 0 represents absence 


e <- evaluate(p=p, a=a)
plot(e, 'ROC')     #ROC plot
auc <- e@auc   #AUC

thresh <- threshold(e)$spec_sens  #optimal threshold for TSS

for (j in 1:nrow(data)){
if (data$predicted_value[j] > thresh){
data$predicted_binary[j] <- 1} else {
data$predicted_binary[j] <- 0
}
}

A <-nrow(subset(data, species == 1 & predicted_binary == 1))
B <-nrow(subset(data, species == 1 & predicted_binary == 0))
C <-nrow(subset(data, species == 0 & predicted_binary == 1))
D <-nrow(subset(data, species == 0 & predicted_binary == 0))
sensitivity <- A/(A+B)   
specificity <- D/(C+D)

tss <- sensitivity + specificity-1  #TSS

####################################################################
# check the best model
best_model <- glm( species ~  forest + snow + angl + cold , family="binomial" , data= data) # formula is an example

AIC(best_model)
check_collinearity(best_model)

# check spatial autocorrelation in the residual of model
# grid coordinates
coords<-matrix(nrow = nrow(data), ncol = 2)
coords[,1] <- data$lat
coords[,2] <- data$lon

maxdis <- max(dist(coords))
nb <- dnearneigh(coords, d1=0,d2= maxdis) 
glist <- nbdists(nb,coords)

glist1 <-lapply(glist,function(x) exp(-x/1000)) ## w(i,j)=exp(-d(i,j)), r=1
w1 <-nb2listw(nb,glist=glist1)

residual <-residuals(best_model)
moran <-moran.test( residual ,listw= w1)

####################################################################
#GAM
#LOOCV
for(i in 1: nrow(data)){

lo_model <- gam( species ~  grass + elev  + s(lon) , family="binomial" , data= data[-i,]) #formula is an example

data[i,]$predicted_value <- predict(lo_model, type="response", newdata= data[i,]) 

}

# the same as those for GLM
p <- subset(data$predicted_value, data$species == 1)
a <- subset(data$predicted_value, data$species == 0)


e <- evaluate(p=p, a=a)
plot(e, 'ROC')     #ROC plot
auc <- e@auc   #AUC

thresh <- threshold(e)$spec_sens 

for (j in 1:nrow(data)){
if (data$predicted_value[j] > thresh){
data$predicted_binary[j] <- 1} else {
data$predicted_binary[j] <- 0
}
}

A <-nrow(subset(data, species == 1 & predicted_binary == 1))
B <-nrow(subset(data, species == 1 & predicted_binary == 0))
C <-nrow(subset(data, species == 0 & predicted_binary == 1))
D <-nrow(subset(data, species == 0 & predicted_binary == 0))
sensitivity <- A/(A+B)   
specificity <- D/(C+D)

tss <- sensitivity + specificity-1  #TSS

###################################################################################
#Stratified k-fold Cross-Validation 

k <- 5

# data frame for results
Results <- data.frame(AUC = rep(-100,k), TSS = rep(-100,k))

set.seed(5) 

folds <- createFolds(factor(data$species), k = k, list = TRUE, returnTrain = F)  #Stratified sampling

for (i in 1:k) {

  test_indices <- folds[[i]]
  train_indices <- setdiff(1:nrow(data), test_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
fo_model <- glm( species ~  forest+ snow + angl+ cold , family="binomial" , data= train_data) #formula is an example

test_data$fo_predicted <- predict(fo_model, type="response", newdata= test_data) 

roc_obj <- roc(test_data$species, test_data$fo_predicted)
auc <- auc(roc_obj)                 #AUC


best_coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
sensitivity <- best_coords$sensitivity
specificity <- best_coords$specificity

tss <- sensitivity + specificity - 1 # TSS

Results$AUC[i] <-  auc
Results$TSS[i] <- tss


}



####################################################################
# prediction 

# threashold
for(i in 1: nrow(data)){

lo_model <-  glm( species ~  forest + snow + angl + cold , family="binomial" , data= data[-i,]) #formula is an example of the best model

data[i,]$predicted_value <- predict(lo_model, type="response", newdata= data[i,]) 

}
p <- subset(data$predicted_value, data$species == 1)
a <- subset(data$predicted_value, data$species == 0)
e <- evaluate(p=p, a=a)
thresh <- threshold(e)$spec_sens    #optimal threshold for TSS


data_mapping$probability <- predict(best_model, newdata= data_mapping, type = "response") # best model

# convirt it to binary format
for(i in 1:length(data_mapping $probability)){
if(data_mapping $probability[i] > thresh){
data_mapping $binary[i] <- 1 } else {
data_mapping $binary[i] <- 0}
}

Probability <- ggplot () +
            geom_sf (data= data_mapping, aes(fill = probability, colour = probability)) +　   
            scale_fill_gradient2(low = "green3", mid = "yellow3", high = "red3", limits = c(0,1), midpoint = 0.5) +
            scale_colour_gradient2(low = "green3", mid = "yellow3", high = "red3", limits = c(0,1), midpoint = 0.5) +
            labs(title = "Species name" ) +
            theme_void() +
            theme(legend.position = "none", text = element_text(family="Arial", size = 8, face = "italic"))  

Binary <- ggplot () +
            geom_sf (data= data_mapping, aes(fill = binary, colour = binary)) +　   
            scale_colour_gradient(low = "grey80", high = "chocolate", limits = c(0,1)) +
            scale_fill_gradient(low = "grey80", high = "chocolate", limits = c(0,1)) +
            labs(title = "Species name" ) +
            theme_void() +
            theme(legend.position = "none", text = element_text(family="Arial", size = 8, face = "italic")) 


####################################################################
# Spatial overlap between pathogen and tick
tick <- as.matrix(data_mapping $tick) # predicted probability
pathogen <- as.matrix(data_mapping $pathogen) # predicted probability


# Schoener's D (Schoener, 1968; Warren et al., 2008)
SD <- function(x, y){
  D <- 1 - 0.5 * (sum(abs((x / sum(x))- (y/sum(y))))) 
  return(D)
}


# Warren’s I (Warren et al., 2008)
WI <- function(x, y){
  H <- sqrt(sum((sqrt((x / sum(x)))-sqrt((y/sum(y))))^2)) 
  I <- 1 - 0.5 * H
  return(I)
}


D_orig <- SD(tick, pathogen)
I_orig <- WI(tick, pathogen)


# boot
n_boot <- 1000 
set.seed(123) 

boot_D <- numeric(n_boot)
boot_I <- numeric(n_boot)
    
for(i in 1:n_boot) {
    idx <- sample(1:length(tick), length(tick), replace = TRUE)
    x_boot <- tick[idx]
     y_boot <- pathogen[idx]
      
     boot_D[i] <- SD(x_boot, y_boot)
     boot_I[i] <- WI(x_boot, y_boot)
    }
    
# 95% CI
D_ci <- quantile(boot_D, probs = c(0.025, 0.975), na.rm = TRUE)
I_ci <- quantile(boot_I, probs = c(0.025, 0.975), na.rm = TRUE)
    








