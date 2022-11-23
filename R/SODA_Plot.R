
SODA_Plot=function(Data,Var,X,Y,Type=c("3D","2D","Interactive"),Res=500)
{

suppressPackageStartupMessages({
require(plot3D)
require(RColorBrewer)
require(ggplot2)
require(rgl)
require(magick)
require(plotly)
})

Main=paste0(Var,":","Interpolated Z-Score")

if(Type=="3D")
{
scatter3D(Data[,X],Data[,Y],Data[,Var],surface=FALSE,
		xlab="Longitude",ylab="Latitude",zlab=Var,
		main=Main,
		col = ramp.col(c("blue", "yellow", "red")),
		ticktype = "detailed",type = "h",bty = "b2",
		 pch =20, cex.lab=1.0,phi =30,theta=30,revolutions=Res)
}else if(Type=="2D"){

ggplot(Data, aes(Data[,X],Data[,Y])) +
    		geom_point(aes(color =Data[,Var])) +
    		scale_color_distiller(palette="RdYlBu",limits=c(-2,2.5))+
    		xlab("Longitude")+ ylab("Latitude")+
		ggtitle(Main)+
		labs(color="Z-score")+theme_classic()

}else if(Type=="Interactive"){

fig <- plot_ly(Data, x = ~Data[,X], y = ~Data[,Y], z = ~Data[,Var],
               marker = list(color = ~Data[,Var],colorscale ="Rainbow", 
		   showscale = TRUE))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Longitude'),
                                   yaxis = list(title = 'Latitude'),
                                   zaxis = list(title =Main))
                      )
fig 
}

}

