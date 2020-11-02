#  ------------------------------------------------------------------------
### Modelos no gausiano lognormal
require(INLA)
require(maps)
require(fields)
require(splancs)
library(gridExtra)
library(lattice)
require(fields)
require(mapdata)
require(geoR)

rm(list = ls())
borders = read.csv("borders3.csv")
listSurveys = list.files(path = "data")
mainDir = "C:/Users/moroncog/Documents/GitHub/AnchovySpatialIndicators"

for(i in seq_along(listSurveys)){

setwd(dir = mainDir)  
  
selectSurvey = listSurveys[i]

dataca = read.csv(file.path("data", selectSurvey))

nameSurvey = sub(pattern = ".csv", replacement = "", x = selectSurvey)
nameSurvey = sub(pattern = "Acoustic_Data_", replacement = "", x = nameSurvey)
outNameSurvey = paste0("Cr", nameSurvey)

dir.create(outNameSurvey, showWarnings = FALSE)
setwd(file.path(mainDir, outNameSurvey))

## Cargando los datos
names(dataca) = tolower(names(dataca))

dataca = dataca[, c("lon_m", "lat_m", "ancp")]
dataca = dataca[complete.cases(dataca),]
dataca = subset(dataca, subset = dataca$lat_m <= -4.7 & dataca$lat_m >= -16)
names(dataca) = c("lon_m", "lat_m", "anc")


png(filename = paste0(outNameSurvey, "-ubm.png"), 
    width = 500, height = 600, res = 120)
par(mfrow = c(1,1), mar=c(4,4,0.5,0.5), xaxs = "i", yaxs = "i")
plot(x = dataca$lon_m, y = dataca$lat_m, pch = ".", col = "gray80", 
     xlab = "longitude", ylab = "latitude", axes = FALSE,
     xlim = c(-83, -74), ylim = c(-16, -5))
points(x = dataca$lon_m, y = dataca$lat_m, pch = 1,
       cex = (dataca$anc/max(dataca$anc))*5)
map("worldHires", add=T, fill=T, col=8)
axis(1, at = seq(-82, -74, by = 2),  labels = 
       seq(-82, -74, by = 2))
axis(2, at = seq(-16, -4, 2) ,  labels = 
       seq(-16, -4, 2), las = 2)
box()
dev.off()


## No le tomo logaritmo pero si le sumo 1 
dataca$anc = dataca$anc + 1

#Extrayendo las posiciones lon, lat
coords = cbind(dataca$lon_m, dataca$lat_m)

## Mesh construction
prdomain = inla.nonconvex.hull(as.matrix(coords),
                               convex = -0.075, concave = -0.5,
                               resolution=c(30,30))

#malla triangulada en 2d
mesh = inla.mesh.2d(loc=coords,boundary=prdomain, min.angle = c(20,20),
                    max.edge=c(0.5,5), cutoff=5/60, offset=c(0.6,0.6))

spde =  inla.spde2.matern(mesh, alpha=2)  #modelo spd2 para funcion matern

png(filename = paste0(outNameSurvey, "-mesh.png"), width = 500, height = 600, res = 140)
par(mar=c(1,1,1,1))
plot(mesh, main="")
points(dataca$lon_m, dataca$lat_m, cex=0.2, lwd=2, col="red")
dev.off()

png(filename = paste0(outNameSurvey, "-mesh_in.png"), width = 500, height = 600, res = 140)
par(mar=c(1,1,1,1))
plot(mesh, main="", xlim = c(-80.5, -79), ylim = c(-9,-7.5))
points(dataca$lon_m, dataca$lat_m, cex=0.2, lwd=2, col="red")
dev.off()

##Matriz de prediccion para la malla 
A = inla.spde.make.A(mesh, loc=coords)

## Preparando los datos
stk.dat = inla.stack(data=list(y=dataca$anc),A=list(A,1), tag='dat',
                     effects=list(list(i=1:spde$n.spde),
                                  data.frame(Intercept=rep(1, length(dataca$anc)))))


#---
#construyendo los modelos
  ## primer modelo - sin covariables

f.lat = y ~ 0 + Intercept + f(i, model=spde)

r.lat = inla(f.lat, family='lognormal', control.compute=list(dic=TRUE),
            data=inla.stack.data(stk.dat),
            control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE))

write.table(r.lat$dic$dic, paste0(outNameSurvey, "-DIC.txt"), row.names = F, col.names = "DIC")
write.csv(r.lat$summary.fixed, paste0(outNameSurvey, "-Intercept.csv"))
write.csv(r.lat$summary.hy, paste0(outNameSurvey, "-Hyper.csv"))


##obteniendo la distribucion marginal
    r.f = inla.spde2.result(r.lat,'i', spde, do.transf = 'TRUE')

    #posterior - media de varianza marginal
    MargVarNom = inla.emarginal(function(x) x, r.f$marginals.variance.nominal[[1]])
    #rango practico
    MargRangNom = inla.emarginal(function(x) x, r.f$marginals.range.nominal[[1]])

MarginalDat = data.frame(VarNom = MargVarNom, RangNom = MargRangNom)

write.csv(MarginalDat, paste0(outNameSurvey, "-MarginalData.csv"), row.names = F)

## Observando la distribucion marginal posteriori de alfaz, alfay, phi, beta
png(filename = paste0(outNameSurvey, "-posteriorDist.png"), width = 700, height = 500, res = 140)
par(mfrow=c(2,2), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0)
plot(r.lat$marginals.fix[[1]], type='l', xlab='Intercept',ylab='Density')
# plot(r.lat$summary.random[[1]][,1:2], type='l',xlab='Latitude', ylab='Effect');
plot(r.lat$marginals.hy[[1]], type='l', ylab='Density', xlab=expression(phi))
# plot.default(inla.tmarginal(function(x) 1/exp(x), r.lat$marginals.hy[[4]]),type='l',
#           xlab=expression(kappa), ylab='Density')
plot.default(r.f$marginals.variance.nominal[[1]], type='l',
              xlab=expression(sigma[x]^2), ylab='Density')
plot.default(r.f$marginals.range.nominal[[1]], type='l',
              xlab='Practical range', ylab='Density')
dev.off()

## Prediccion del campo aleatorio
gridsize = 1

  j = 1
  stepsize = j/60 #grilla 3mn
  nxy = round(c(diff(range(borders[,1])), diff(range(borders[,2])))/stepsize)
  
    #Projeccion
    projgrid = inla.mesh.projector(mesh, xlim=range(borders[,1]),
                                  ylim=range(borders[,2]), dims=nxy)
  
    xmean = inla.mesh.project(projgrid, r.lat$summary.random$i$mean) #posterior mean
    xsd   = inla.mesh.project(projgrid, r.lat$summary.random$i$sd) #posterior sd
  
    
    #inside the borders 
  #esto lo puedes cambiar en funcion a la varianza
  
  xy.in = inout(projgrid$lattice$loc,
                      cbind(borders[,1], borders[,2]))
  
  xmean[!xy.in] = xsd[!xy.in] = NA  #doy valores de N.A fuera del borde

  write.csv(xmean, paste0(outNameSurvey,"_", j, "mn", "-xmean.csv"), row.names = F)
  write.csv(xsd, paste0(outNameSurvey,"_", j, "mn", "-xsd.csv"), row.names = F)
  
  
        #plot prediction      
  jet.colors <- colorRampPalette(c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow',
                                         '#FF7F00','red','#7F0000') )
        
  png(filename = paste0(outNameSurvey, "_", j, "mn", "-mapPred_Var.png"), width = 500, height = 750, res = 140)
  do.call('grid.arrange',
        lapply(list(xmean, xsd),
           levelplot, col.regions=jet.colors(200),
           xlab='', ylab='', scales=list(draw=FALSE)))
  dev.off()
        
        Longitud = unique((projgrid$x))
        Latitud = unique((projgrid$y))
        predict = matrix(data=xmean, ncol=length(Latitud), nrow=length(Longitud))
        varianza = matrix(data=xsd, ncol=length(Latitud), nrow=length(Longitud))
  
  write.csv(predict, paste0(outNameSurvey,"_", j, "mn", "-predict_raw.csv"), row.names = F)
  write.csv(varianza, paste0(outNameSurvey,"_", j, "mn", "-variance_raw.csv"), row.names = F)
  
  predict[which(predict < 0)] = 0
  ztemp = predict * 10
  maxi = round(max(ztemp, na.rm = T))
  scale.color = designer.colors(maxi, c('white', 'lightblue','green', 'yellow', 'red','black'))
  
  write.csv(predict, paste0(outNameSurvey, "-predict.csv"), row.names = F)
  write.csv(Longitud, paste0(outNameSurvey, "-Longitud.csv"), row.names = F)
  write.csv(Latitud, paste0(outNameSurvey, "-Latitud.csv"), row.names = F)
  
  z_temp = as.vector(predict)
  coordenadas = merge(x = Longitud, y = Latitud, all = T)
  outMatrix = data.frame(lon = coordenadas$x, lat = coordenadas$y,
                         var = z_temp)
  write.csv(outMatrix, paste0(outNameSurvey,"_", j, "mn", "-predData.csv"), row.names = F)

  
  png(filename = paste0(outNameSurvey,"_", j, "mn", "-mapPred_show.png"), 
      width = 500, height = 600, res = 120)
  par(mfrow = c(1,1), mar=c(4,4,0.5,0.5))
  image(Longitud, Latitud, predict, col=scale.color, main="",
             zlim = c(0,ceiling(max(predict, na.rm = T))), 
             axes=F, xlim = c(-83, -74), ylim = c(-16, -5), 
        xlab = "longitude", ylab = "latitude")
  map("worldHires", add=T, fill=T, col=8)
  axis(1, at = seq(-82, -74, by = 2),  labels = 
         seq(-82, -74, by = 2))
  axis(2, at = seq(-16, -4, 2) ,  labels = 
         seq(-16, -4, 2), las = 2)
  legend.krige(c(-82.5,-82.1),c(-15.5,-9.5), 
               predict, vertical=T, col=scale.color)
  box()
  dev.off()

  print(i)

}