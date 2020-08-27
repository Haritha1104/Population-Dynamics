#Load the necessary packages
library(ape)
library(deSolve) 
library(primer)
library(sf)
library(plyr)


#Import the data
data      <-  read.csv("Allsiteallspeciesoccurencedata2015citizenscience.csv")
head(data)
     
     
#Remove missing values from data
data  <-  data[-which(is.na(data$X)),]
data <- data[-which(is.na(data$No)) , ]
     
     

#Plot the number of species detected in every month
mm = subset(data , data$Species=="Martes martes")
sv = subset(data , data$Species=="Sciurus vulgaris")
sc = subset(data , data$Species=="Sciurus carolinensis")

plot(count(mm$Month),col="green" , xlim=c(1,8),ylim=c(0,400),type="b", pch=19 ,xlab = "Months" , ylab ="Population")
points(count(sc$Month) , col = "grey" , type="b",pch=19)
points(count(sv$Month) , col = "red", type="b",pch=19)
legend("topleft" ,legend =  c("red squirrel" , "grey squirrel" , "pine marten") , col = c("red" , "grey" , "green") , pch = 19 , cex=0.5)


#Subset of the data containing only the species of interest
species_data = data[which(data$Species=="Sciurus carolinensis" | data$Species=="Sciurus vulgaris" | data$Species=="Martes martes"  ),]



#Extract the population size corresponding to each site
T = table(species_data$Site , species_data$Species)

#Subset the species of interest
T = T[,c(8,13,14)]

#Remove the sites which has 0 detection for all the three species
T = T[which(T[,1] !=0 | T[,2] != 0 | T[,3] != 0),]

#Find the unique regions in the dataset
x_values = unique(species_data$X)
y_values = unique(species_data$Y)
regions = unique(species_data$Site)
id_of_uniq_regions = match(regions , species_data$Site)

#Form a dataset with the unique sites, their coordinates and count of species
coords_of_uniq_regions = cbind(species_data$X[id_of_uniq_regions], species_data$Y[id_of_uniq_regions])
rownames(coords_of_uniq_regions) <- regions
my_data =cbind(coords_of_uniq_regions , T)
colnames(my_data)[c(1,2)] <- c("X" , "Y" )



#Spatial correlation tests

#Compute the distance matrix of coordinates 
d = as.matrix(dist(coords_of_uniq_regions, method="euclidean"))

#Moran's test (Auto correlation)

#inverse distance matrix
d_inv <- 1/d
diag(d_inv) <- 0
d_inv[is.infinite(d_inv)] <- 0



#Moran's test for pine marten
Moran.I(my_data[,3] , d_inv)


#Moran's test for grey squirrel
Moran.I(my_data[,4] , d_inv)


#Moran's test for red squirrel
Moran.I(my_data[,5] , d_inv)




#Mantel's test (Crosscorrelation)

#distance matrices of species
mm_dist = as.matrix(dist(my_data[,3]))
sc_dist = as.matrix(dist(my_data[,4]))
sv_dist = as.matrix(dist(my_data[,5]))


#Mantel's test for pine marten and grey squirrel
mantel.test(mm_dist, sc_dist)

#Mantel's test for red squirrel and grey squirrel
mantel.test(sv_dist, sc_dist)

#Mantel's test for pine marten and red squirrel
mantel.test(mm_dist, sv_dist)




#Plot of Northern Ireland with the spatial distribution of species

ire_data <- read_sf("1km grid.shp")
ire_data_sp <- as(ire_data, Class = "Spatial")
coords_ire = coordinates(ire_data_sp)

color = rep(NA, length=length(species_data$Species))
color[which(species_data$Species=="Sciurus vulgaris")] = "red"
color[which(species_data$Species=="Sciurus carolinensis")] = "grey"
color[which(species_data$Species=="Martes martes")] = "green"

par(mar=c(0,0,0,0))
plot(ire_data_sp, col="#f2f2f2", bg="skyblue", lwd=0.25, border=0 )
points(coords ,pch = 19, col=color)
legend("topleft" ,legend =  c("red squirrel" , "grey squirrel" , "pine marten") , col = c("red" , "grey" , "green") , pch = 19)



#Lotka-Volterra models

#Model 1 - Red squirrel


LV1 <-function (t, n, parms) 
{
        with(as.list(parms), {
                dn1dt <- r1 * n[1] * (1 - ( n[1] / K) )
                
                list(c(dn1dt))
        })
}

parms <- c(r1 = 10, K = 0.5) 
initial <- c(2.8) 
out <- ode(y = initial, times = 1:50,func = LV1, parms = parms) 
matplot(out[, 1], out[, -1], type = "l",lwd=2,ylim=c(0,3),xlab="time", ylab="population size",col="red", main="Population growth") 
legend("topright" ,legend =  c("red squirrel") , col = c("red" ) ,lwd=2, lty = 1 , cex=0.6)


#Model 2 - Red squirrel vs Grey squirrel


LV2 <-function(t, n, parms)
{
        with(as.list(parms), {
                dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a11*a12 * n[2] )
                dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a22*a21 * n[1] )
                list(c(dn1dt, dn2dt))
        })
}

parms4 <- c(r1 = 1, r2 = 1.2 , a11 = 1/0.5, a21 =2.55, a22 = 1/1.25, a12 = 1.65 ) 
initial <- c(2.8, 0.07) 
out4 <- ode(y = initial, times = 1:400, func = LV2, parms = parms4) 
matplot(out4[, 1], out4[, -1], col=c("red" , "gray"),lwd=2,lty=c(1,2), type = "l",xlab="time", ylab="population size", main="Population growth")
legend("topright" ,legend =  c("red squirrel","grey squirrel") ,lty=c(1,2), col = c("red","grey" ) , cex=0.5)


#Zero-growth isoclines
plot(x=-1:5 , y = (0.5 - x)/1.65 , type="l", col="red" , ylab="Grey squirrel population" , xlab="Red squirrel population",main="Zero-growth isoclines",xlim=c(0.09,3) , ylim=c(0.09,3))
lines(x=-1:5 , y = 1.25-0.61*x , col="dark grey" )
points(x=0,y=1.25, cex=1.5)
legend("topright" ,legend =  c("red squirrel isocline" , "grey squirrel isocline") , col = c("red" , "dark grey" ) , lty = c(1,1) , cex=0.7)


#Model 3 - Red squirrel vs Grey squirrel vs Pine marten


LV3 <-function (t, n, parms) 
{
        with(as.list(parms), {
                dn1dt <- n[1] * (r1 - a11 * n[1] - a11 *a12 * n[2] - a13 * n[3]) 
                dn2dt <- n[2] * (r2 - a22 * n[2] - a22 *a21 * n[1] - a23 * n[3]) 
                dn3dt <- n[3] * (-r3 + a32 * n[2] + a31 * n[1] )
                list(c(dn1dt, dn2dt, dn3dt))
        })
}


parms <- c(r1 = 1, r2 = 1.2, r3 = 0.384, 
           a11 = 1/0.5, a12 = 1.65, a13 = 0.44, 
           a21 = 0.61, a22 = 1/1.25, a23 = 2.2,
           a32 =1.2, a31=0.5) 
initial <- c(2.8, 0.07, 1) 
out <- ode(y = initial, times = 1:100, func = LV3, parms = parms) 
matplot(out[, 1], out[, -1],type = "l",lwd =2,lty=c(3,1,2),xlab="time", ylab="population size", col=c("red", "grey","dark green")) 
legend("topright" ,legend =  c("red squirrel" , "grey squirrel","pine marten") , lwd =2,lty=c(3,1,2),col = c("red", "grey","dark green") , cex=0.7)



#Model 4 - Red squirrel vs Juvenile Grey squirrel vs Adult Grey squirrel vs Pine martem


LV4 <-function (t, n, parms) 
{
        with(as.list(parms), {
                dn1dt <- n[1] * (r1 - a11 * n[1] - a11 *a13 * n[3] - a14 * n[4])#RS 
                dn2dt <-  -1*((1-s2) +s2 * p) *n[2]  +b* n[3]/2+ n[2] * (r2 - a22 * n[2] - a22 *a21 * n[1] - a24 * n[4])#GS-juv
                dn3dt <-  s2 * p *n[2] - s3 * n[3] + n[3] * (r3 - a33 * n[3] - a33 *a31 * n[1] - a34 * n[4])#GS-adult
                dn4dt <- n[4] * (-s +a41 * n[1]  +a42 * n[2] + a43 * n[3])#PM
                list(c(dn1dt, dn2dt, dn3dt , dn4dt))
        })
}


parms <- c(r1 = 1, r2 = 1.2, r3 = 1.2, s=0.384, 
           a11 = 1/0.5, a12 = 0, a13 = 1.65, a14 = 0.44, 
           a21 = 0, a22 = 1/1.25 ,a24 =3.4,
           a31 =0.61 , a33= 1/1.25 , a34 =2.2 ,
           a41=0.5, a42=1.5, a43=1.2, 
           s2 =0.26 , s3 =0.48, p =0.5, b =  3.4
) 
initial <- c(2.8, 0.14,0.07, 1) 
out <- ode(y = initial, times = 1:150, func = LV4, parms = parms) 
matplot(out[, 1], out[, -1], ylim=c(0,0.4),type = "l",lwd=2,xlab="time", ylab="population size", lty = c(6,1,4,1),col=c("red", "grey" ,"brown","dark green")) 
legend("topright" ,legend =  c("red squirrel" , "grey squirrel juvenile","grey squirrel adult","pine marten") , lwd=2,col = c("red", "grey", "brown","dark green") , lty = c(6,1,4,1) , cex=0.5)


