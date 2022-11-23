
SODA_Kriging=function(Data,X,Y,Expand,Margin)
{

suppressPackageStartupMessages({
require(sf)
require(sp)
require(automap)
require(tidyverse)
require(erer)
})

my_mesh=expand.grid(seq(min(Data[,X]-Margin),max(Data[,X]+Margin),l=Expand),
		seq(min(Data[,Y]-Margin),max(Data[,Y]+Margin),l=Expand))
my_mesh=data.frame(my_mesh)
colnames(my_mesh) <- c('Longitude','Latitude')

pts <- SpatialPoints(coords=my_mesh)
grd <- as(pts, "SpatialPixels")

Grd=st_crs(my_mesh)
sf_data <- st_as_sf(Data, coords = c("Longitude","Latitude"),crs=Grd)

Intrp=matrix(NA,dim(my_mesh)[1],(dim(Data)[2])-2)
colnames(Intrp)=names(Data)[1:(dim(Data)[2]-2)]

P=list()
i=1
while(i<=(dim(Data)[2]-2))
{
Dat=as.data.frame(sf_data[i])[,1]
Nam=names(as.data.frame(sf_data[i]))[1]
# Automatized Kriging
fit_KRIG <- automap::autoKrige(
  formula =Dat ~ 1,
  input_data = as(sf_data, "Spatial"),grd
)
fit=fit_KRIG %>%
  .$krige_output %>%
  as.data.frame() %>%
  dplyr::select(Lon=Longitude,Lat=Latitude, Z= var1.pred)

Intrp[,i]=fit$Z
print(i)
i=i+1

}

Intrp2=data.frame(Longitude=Data[,X],Latitude=Data[,Y],Intrp)
return(Intrp2)

}
