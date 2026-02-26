
# Sample R code for Spatial distribution of ticks and tick-borne pathogens in central Hokkaido, Japan and associated ecological factors revealed by intensive short-term survey in 2024
# Note: Due to privacy policies regarding sampling coordinates, 
# the original dataset is provided in the manuscript in anonymized format.
# Please replace 'data' with your own dataframe
# crs = 2454 in our analysis

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

#################################################################################################
#Global and local spatial statistics
#################################################################################################
# simulate spatial weight depend on distance decay parameter
dis <- seq(0, 200, 0.001)
difun <- data.frame( dis = seq(0, 200, 0.001), f1 = exp(-dis/0.1), f2 = exp(-dis), f3 = exp(-dis/10), f4 = exp(-dis/20))

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

# sampling site coordinate
coords<-matrix(nrow = nrow(data), ncol = 2)
coords[,1] <- data$Longitude
coords[,2] <- data$Latitude

maxdis <- max(dist(coords))
nb　<-　dnearneigh(coords, d1=0,d2= maxdis) # with all the others
glist　<-　nbdists(nb,coords)

# calculate moran's I and p-value with r = 0.1 to r = 20
#r=0.1~1

for (i in 1:10) {

ik <- i * 100

glist_i　<-　lapply(glist,function(x) exp(-x/ik)) # w(i,j)=exp(-d(i,j)/jk)
w　<-　nb2listw(nb,glist=glist_i) 

moiov <-moran.test( species,listw=w)
moran_p$species[i] <- moiov$p
moran_e$species[i] <- moiov$estimate[1]

}

#r=2~20

for (j in 2:20) {

jk <- j * 1000

glist_j　<-　lapply(glist,function(x) exp(-x/jk)) # w(i,j)=exp(-d(i,j)/jk)
w　<-　nb2listw(nb,glist=glist_j)


moiov <-moran.test( species,listw=w)

moran_p$species[j+9] <- moiov$p
moran_e$species[j+9] <- moiov$estimate[1]

}


############################################################################################
#G spatial statistics

maxdis <- max(dist(coords))
nb　<-　dnearneigh(coords, d1=0,d2= maxdis) 
glist　<-　nbdists(nb,coords)
glist_1　<-　lapply(glist,function(x) exp(-x/1000)) # w(i,j)=exp(-d(i,j)/1)
w　<-　nb2listw(nb,glist= glist_1,) 



local_g_species <- localG_perm(data$species, w, zero.policy = TRUE)
data$gi_z_species <- as.vector(local_g_species) 

data <- data %>% 
  mutate(hotspot_category_Iov = case_when(
    gi_z_species > 2.58  ~ "Hotspot (99% CI)",
    gi_z_species > 1.96  ~ "Hotspot (95% CI)",
    gi_z_species < -2.58 ~ "Coldspot (99% CI)",
    gi_z_species < -1.96 ~ "Coldspot (95% CI)",
    TRUE              ~ "Not Significant"
  ))


map <- ggplot() +
          geom_sf(data = data, fill = "white", color = "white") + 
          geom_point(data = data, aes(x = Xmo , y = Ymo, color = hotspot_category_species, size =5), alpha = 0.8) +   
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

#################################################################################################
# Ecological niche modeling
###################################################################################
# Screening potential important variables

model <- glm( iov ~ forest+grass+wet+ crop + build+ water + snow+ prec6  + sqprec6 + temp6 + sqtemp6 + suns+ cold + elev+ sqelev + angl+ sqangl + lon + lat + deer+bear+raccoon + tanuki + fox , family="binomial" , data= data) #formula is an example

options(na.action = "na.fail")
selection <- dredge (model, rank="AIC")
comp_models <- get.models(selection,  subset = delta < 2)

###################################################################################
#GLM

#data frame for results
data$predicted_value <- -1000
aic <- rep(-1000, nrow(data))

#LOOCV
for(i in 1: nrow(data)){

lo_model <- glm( iov ~ forest+ snow + angl+ cold , family="binomial" , data= data[-i,]) #formula is an example

data[i,]$predicted_value <- predict(lo_model, type="response", newdata= data[i,]) 

aic <- AIC (lo_model) #AIC

}

p <- subset(data$predicted_value, data$fact == 1) # 1 represents presence 
a <- subset(data$predicted_value, data$fact == 0) # 0 represents absence 


e <- evaluate(p=p, a=a)
plot(e, 'ROC')     #ROC plot
auc <- e@auc   #AUC


thresh <- threshold(e)$spec_sens  #optimal threshold for TSS

for (j in 1:nrow(data)){
if (data$lo_iov[j] > thresh){
data$PRED[j] <- 1} else {
data$PRED[j] <- 0
}
}

A <-nrow(subset(data, iov == 1 & PRED == 1))
B <-nrow(subset(data, iov == 1 & PRED == 0))
C <-nrow(subset(data, iov == 0 & PRED == 1))
D <-nrow(subset(data, iov == 0 & PRED == 0))
sensitivity <- A/(A+B)   
specificity <- D/(C+D)

tss <- sensitivity + specificity-1  #TSS

####################################################################
# check spatial autocorrelation in the residual of model

#spatial weights matrix
coords<-matrix(nrow = nrow(data), ncol = 2)
coords[,1] <- data$Longitude
coords[,2] <- data$Latitude

maxdis <- max(dist(coords))
nb<-dnearneigh(coords, d1=0,d2= maxdis) 
glist<-nbdists(nb,coords)


glist1 <-lapply(glist,function(x) exp(-x/1000)) ## w(i,j)=exp(-d(i,j)), r=1
w1 <-nb2listw(nb,glist=glist1)

best_model <- glm( iov~  forest+ snow + angl+ cold , family="binomial" , data= data) # formula is an example

residual <-residuals(best_model)
moran <-moran.test( residual ,listw= w1)

####################################################################
#GAM
#LOOCV
for(i in 1: nrow(data)){

lo_model <- gam( ipa ~  grass + elev  + ti (lon, lat, bs = c("tp", "tp")) , family="binomial" , data= data[-i,]) #formula is an example

data[i,]$predicted_value <- predict(lo_model, type="response", newdata= data[i,]) 

}

p <- subset(data$predicted_value, data$fact == 1)
a <- subset(data$predicted_value, data$fact == 0)


e <- evaluate(p=p, a=a)
plot(e, 'ROC')     #ROC plot
auc <- e@auc   #AUC


thresh <- threshold(e)$spec_sens  #spec_sens: best sensitivity+specificity threshold

for (j in 1:nrow(data)){
if (data$lo_iov[j] > thresh){
data$PRED[j] <- 1} else {
data$PRED[j] <- 0
}
}

A <-nrow(subset(data, iov == 1 & PRED == 1))
B <-nrow(subset(data, iov == 1 & PRED == 0))
C <-nrow(subset(data, iov == 0 & PRED == 1))
D <-nrow(subset(data, iov == 0 & PRED == 0))
sensitivity <- A/(A+B)   
specificity <- D/(C+D)

tss <- sensitivity + specificity-1  #TSS

###################################################################################
#Stratified k-fold Cross-Validation 

k <- 5

# data frame for results
Results <- data.frame(AUC = rep(-100,k), TSS = rep(-100,k))

set.seed(5) 

folds <- createFolds(factor(data$species), k = k, list = TRUE, returnTrain = F) 

for (i in 1:k) {

  test_indices <- folds[[i]]
  train_indices <- setdiff(1:nrow(data), test_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
fo_model <- glm( iov~  forest+ snow + angl+ cold , family="binomial" , data= train_data) #formula is an example

test_data$fo_predicted <- predict(fo_model, type="response", newdata= test_data) 

roc_obj <- roc(test_data$iov, test_data$fo_predicted)
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

data$probability <- predict(best_model, newdata= data, type = "response") # best model

# convert it to binary format with
thresh <- threshold(e)$spec_sens  #optimal threshold for TSS

for(i in 1:length(data$probability)){
if(data$probability[i] > thresh){
data$binary[i] <- 1 } else {
data$binary[i] <- 0}
}

Probability <- ggplot () +
            geom_sf (data= data, aes(fill = probability, colour = probability)) +　   
            scale_fill_gradient2(low = "green3", mid = "yellow3", high = "red3", limits = c(0,1), midpoint = 0.5) +
            scale_colour_gradient2(low = "green3", mid = "yellow3", high = "red3", limits = c(0,1), midpoint = 0.5) +
            labs(title = "Species name" ) +
            theme_void() +
            theme(legend.position = "none", text = element_text(family="Arial", size = 8, face = "italic"))  

Binary <- ggplot () +
            geom_sf (data= data, aes(fill = binary, colour = binary)) +　   
            scale_colour_gradient(low = "grey80", high = "chocolate", limits = c(0,1)) +
            scale_fill_gradient(low = "grey80", high = "chocolate", limits = c(0,1)) +
            labs(title = "Species name" ) +
            theme_void() +
            theme(legend.position = "none", text = element_text(family="Arial", size = 8, face = "italic")) 


####################################################################
# Spatial overlap between pathogen and tick
tick <- as.matrix(data$tick)
pathogen <- as.matrix(data$pathogen)


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



#data frame for results
di_results <- data.frame()

D_orig <- SD(tick, pathogen)
I_orig <- WI(tick, pathogen)


# boot
n_boot <- 1000 
set.seed(123) 

boot_D <- numeric(n_boot)
boot_I <- numeric(n_boot)
    
for(i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    x_boot <- tick[idx]
     y_boot <- pathogen[idx]
      
     boot_D[i] <- SD(x_boot, y_boot)
     boot_I[i] <- WI(x_boot, y_boot)
    }
    
# 95% CI
D_ci <- quantile(boot_D, probs = c(0.025, 0.975), na.rm = TRUE)
I_ci <- quantile(boot_I, probs = c(0.025, 0.975), na.rm = TRUE)
    








