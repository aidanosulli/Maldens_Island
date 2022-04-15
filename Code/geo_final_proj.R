library(readr)
library(geoR)
library(gstat)
library(tidyverse)
library(naniar) #for replacing NA's
library(oceanmap) #for the map of buoys south of Malden
library(sp)
library(plotly)



elnino <- read.csv("elnino.csv")


#Find "." entries, assign them to be NA, and then drop them
elnino <- apply(elnino, 2, function(x) {str_replace_all(x, "^\\.", "NA")}) 
elnino <- as.data.frame(elnino)
for (j in 1:dim(elnino)[2]){
    elnino[,j][elnino[,j]=="NA"] <- NA
}
elnino <- elnino %>% drop_na


#Focus Data Around Malden Island
#Malden Left:-3.57329, -155.58559
#Malden Right:-4.5267, -154.32766
elnino <- elnino %>% filter(between(Latitude,-5.2136, -4)) 
elnino <- elnino %>% filter(between(Longitude, -155.3, -154.75))

#Check for Duplicate values and delete them
elnino <- elnino %>% distinct(Longitude, Latitude, .keep_all = TRUE)
elnino <- elnino[-48,] #Drop bouy that strayed from the rest


#Finally, coercing all the columns to their proper type
elnino <- elnino %>% mutate(
  Zonal.Winds = as.numeric(Zonal.Winds), 
  Meridional.Winds = as.numeric(Meridional.Winds), 
  Humidity = as.numeric(Humidity), Air.Temp = as.numeric(Air.Temp), 
  Sea.Surface.Temp = as.numeric(Sea.Surface.Temp),
  Longitude = as.numeric(Longitude),
  Latitude = as.numeric(Latitude)
)

#make final trims to elnino data
elnino <- elnino %>% dplyr::select(Longitude, Latitude, Zonal.Winds, Meridional.Winds, Humidity, Air.Temp, Sea.Surface.Temp)




# sum(between(elnino$Latitude,-5.2136, -3.05505))
# sum(between(elnino$Longitude, -155.3, -154.75))

#create interactive plot to help focus in on long/lat coordinates
elnino2 <- as.data.frame(elnino)
map1 <- ggplot(data = elnino2, aes(x = Longitude, y = Latitude)) + geom_point()
ggplotly(map1)



# library(oceanmap)

#Create Map of Malden Island and Buoys

lon <- c(-155.1, -154.6) 
lat <- c(-5.2, -4)
# figure(width =9, height =5.28)
# plot (elnino$Longitude, elnino$Latitude, xlim =lon , ylim = lat )
# plotmap (main = "Buoys South of Malden Island", add=T)
# points(elnino$Longitude, elnino$Latitude, pch = 1)


#Bubble Plot using ggplot2
Humidity_buoy <- ggplot(data = as.data.frame(elnino), aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(colour = Humidity, size = Humidity/mean(Humidity))) + 
  guides(size = FALSE) +
  scale_colour_gradient(low = "yellow", high = "red")  +
  theme(panel.background = element_rect(fill = "lightblue")) +
  ggtitle("Humidity Across Buoys")


#correlation
subset <- as.data.frame(elnino) %>% select(Zonal.Winds, Meridional.Winds, Humidity, Air.Temp, Sea.Surface.Temp)
upper <- round(cor(subset),2)
upper[upper.tri(upper)]<-"" 
upper_correlation <- as.data.frame(upper) 


#h-scatterplots
# library(gstat)
# library(sp)






####Variograms


###################################### Rose Diagram ##################
#save the data as a geo_data object
library(geoR)
bb <- as.geodata(elnino)

v <- summary(bb)
# v$distances.summary
#max distance in bb is 0.1118

#Compute the variogram for the following directions:
var1 <- variog(bb, dir=pi/2, tol=pi/4, max.dist= 0.05)
var2 <- variog(bb, dir=pi/2.57, tol=pi/4, max.dist=0.05)
var3 <- variog(bb, dir=pi/3.6, tol=pi/4, max.dist=0.05)
var4 <- variog(bb, dir=pi/6, tol=pi/4, max.dist=0.05)
var5 <- variog(bb, dir=pi/18, tol=pi/4, max.dist=0.05)
var6 <- variog(bb, dir=0.944*pi, tol=pi/4, max.dist=0.045)
var7 <- variog(bb, dir=0.833*pi, tol=pi/4, max.dist=0.045)
var8 <- variog(bb, dir=0.722*pi, tol=pi/4, max.dist=0.045)
var9 <- variog(bb, dir=0.611*pi, tol=pi/4, max.dist=0.045)


#plot(var1, main = "Var 1") #range ~ .05
#plot(var2, main = "Var 2") #range ~ .05
# plot(var3, main = "Var 3") #range ~ .05
# plot(var4, main = "Var 4") #range ~ .05
# plot(var5, main = "Var 5") #range ~ .05
# plot(var6, main = "Var 6") #range ~ .045
# plot(var7, main = "Var 7") #range ~ .045
# plot(var8, main = "Var 8") #range ~ .045
# plot(var9, main = "Var 9") #range ~ .045

###################################### End Rose Diagram ####################


#Fit Model Variogram using Normal Estimator
a <- as.data.frame(elnino)
names(a) <- c("x", "y", "Zonal.Winds", "Meridional.Winds", "Humidity", "Air.Temp", "Sea.Surface.Temp")

gstat_ob <- gstat(id = "Humidity", formula = Humidity ~1, locations = ~x+y, data = a)
vario_normal <- variogram(gstat_ob, cutoff = 0.05)
best_fit_normal <- fit.variogram(vario_normal, vgm(15,c("Sph", "Lin", "Exp"), 0.035, nugget = 10), fit.method = 2)


#Fit Model Variogram using Robust Estimator
vario_robust <- variogram(gstat_ob, cutoff = 0.05, cressie = T)
best_fit_robust <- fit.variogram(vario_normal, vgm(15,c("Sph", "Lin", "Exp"), 0.035, nugget = 10), fit.method = 2)


#################################### Kriging ########################
#Create grid of prediction points

#create grid of points to be predicted
grd <- expand.grid(x=seq(from= -155.03, to= -154.925, by=0.001), y=seq(from=-5.03, to=-4.96, by=0.001))
el_ok <- krige(id="Humidity", Humidity~1, locations=~x+y, model=best_fit_normal, data=a, newdata=grd)


#Co Kriging

g1 <- gstat(id="Humidity", formula = Humidity~1, locations = ~x+y, data = a)
g1 <- gstat(g1,id="Air.Temp", formula = Air.Temp~1, locations = ~x+y, data = a) 

vm <- variogram(g1)
vm.fit <- fit.lmc(vm, g1, model=best_fit_normal) 
predicted_ck <- predict(vm.fit, grd)



############################## Cross Validation ####################


#Ordinary vs. Universal Kriging Cross Validation

#Cross validation using Normal Estimator
#Ordinary kriging
cv_normal_ok<- krige.cv(Humidity~1, data=a, locations=~x+y,model=best_fit_normal,nfold=nrow(a)) 
press_normal_OK <- sum(cv_normal_ok$residual^2)/nrow(a)
#Universal kriging:
cv_normal_uk <- krige.cv(Humidity~x+y,data=a, locations=~x+y, model=best_fit_normal,nfold=nrow(a)) 
press_normal_UK <- sum(cv_normal_uk$residual^2)/nrow(a)

#Cross Validation using Robust Estimaor

#Ordinary Kriging
cv_robust_ok <- krige.cv(Humidity~1, data = a, locations = ~x+y, model =best_fit_robust, nfold = nrow(a))
press_robust_OK <- sum(cv_robust_ok$residual^2)/nrow(a)

#Universal
cv_robust_uk <- krige.cv(Humidity ~1, data = a, locations = ~x+y, model = best_fit_robust, nfold = nrow(a))
press_robust_UK <- sum(cv_robust_uk$residual^2)/nrow(a)

answer <- matrix(c(press_normal_OK, press_normal_UK, press_robust_OK, press_robust_UK), nrow = 4, byrow = F) 
rownames(answer) <- c("PRESS Normal, Ordinary", "PRESS Normal, Universal", "PRESS Robust, Oridinary", "PRESS Robust, Universal")


#Co Kriging Cross Validation
cv_ck <- gstat.cv(vm.fit)
press_co_krig <- sum(cv_ck$residual^2)/nrow(a)

answer2 <- matrix(c(press_co_krig, press_normal_OK), nrow =2)
rownames(answer2) <- c("PRESS Co-Kriging", "PRESS Normal, Ordinary")


######################################################



# plot(variogram(g))


#Extra Plots

library(lattice)


# levelplot(predicted_ck$Humidity.pred~x+y, predicted_ck, aspect = "iso")
# contourplot(el_ok$Humidity.pred~x+y, predicted_ck, aspect = "iso", add = T)

# lon <- c(-155.3, -154.6) 
# lat <- c(-5.2, -4)
# figure(width =9, height =5.28)
# plot (elnino$Longitude, elnino$Latitude, xlim =lon , ylim = lat )
# plotmap (main = "Buoys South of Malden Island", add=T)
# points(elnino$Longitude, elnino$Latitude, pch = 1)

