#  ------------------------------------------------------------------------
mainDir = "C:/Users/moroncog/Documents/GitHub/AnchovySpatialIndicators"
setwd(dir = mainDir)  

require(fields)
require(mapdata)
require(maps)
require(geoR)
library(sp)
require(rgeos)
library('raster')
source("Spatial_indicators_functions.R")
source("distCoast.R")

infTh = 0.05

listSurveys = list.dirs(path = "surveys")
listSurveys = listSurveys[-1]
nameSurvey = sub(pattern = "surveys/", replacement = "", x = listSurveys)

std_level = 4.86 # calculado con thr level de 0.05

levelhot = numeric(length(listSurveys))
areahot = numeric(length(listSurveys))
areahot_std = numeric(length(listSurveys))
arealow = numeric(length(listSurveys))
saplus = numeric(length(listSurveys))
satot = numeric(length(listSurveys))
npatchhot = numeric(length(listSurveys))
npatchhot_std = numeric(length(listSurveys))
npatchlow = numeric(length(listSurveys))
isotropy = numeric(length(listSurveys))
hotpatch_dc = numeric(length(listSurveys))
hotpatchstd_dc = numeric(length(listSurveys))
for(i in seq_along(listSurveys)){

  selectSurvey = listSurveys[i]
  
  lonlat = read.csv(file.path(selectSurvey, paste0(nameSurvey[i], "_1mn-predData.csv")))
  predict = read.csv(file.path(selectSurvey, paste0(nameSurvey[i], "_1mn-predict_raw.csv")))
  predict = as.matrix(predict)

  predict2 = t(predict)
  predict3 = predict2[(nrow(predict2):1),]
  
  rr = raster(predict3, xmn = -83.88621, xmx = -73.7257, ymn = -16.20002, ymx = -4.564909) 
  
  maxVal = ceiling(max(predict, na.rm = T))
  
  Level = quantile(lonlat$var[which(lonlat$var > infTh)], probs = 0.95, na.rm = T)
  levelhot[i] = Level
  
  # limites siempre los mismos a no ser que cambie el area de prediccion
  rc = cut(rr, breaks = c(Level, maxVal))
  pols = rasterToPolygons(rc, dissolve=T)
  areahot[i] = pols@polygons[[1]]@area*3600
  npatchhot[i] = length(pols@polygons[[1]]@Polygons)

  rcstd = cut(rr, breaks = c(std_level, maxVal))
  polsstd = rasterToPolygons(rcstd, dissolve=T)
  areahot_std[i] = polsstd@polygons[[1]]@area*3600
  npatchhot_std[i] = length(polsstd@polygons[[1]]@Polygons)
  
    
  rc2 = cut(rr, breaks = c(infTh, maxVal))
  pols2 = rasterToPolygons(rc2, dissolve=T)
  arealow[i] = pols2@polygons[[1]]@area*3600
  npatchlow[i] = length(pols2@polygons[[1]]@Polygons)
    
  pos = which(lonlat$var >= infTh)
  saplus[i] = mean(lonlat$var[pos])
  satot[i] = sum(lonlat$var[pos])
    
  plot.new()
  isotropy[i] = cgi(x = lonlat$lon[pos], y = lonlat$lat[pos], z = lonlat$var[pos])$Iso
  
  tmp_dc = numeric(length(pols@polygons[[1]]@Polygons))
  tmp_areas = numeric(length(pols@polygons[[1]]@Polygons))
  for(k in seq_along(pols@polygons[[1]]@Polygons)){
    tmp_areas[k] = pols@polygons[[1]]@Polygons[[k]]@area
    tmp_dc[k] = mean(distCoast(lon = pols@polygons[[1]]@Polygons[[k]]@coords[,1], 
                       lat = pols@polygons[[1]]@Polygons[[k]]@coords[,2]))
  }
  
  hotpatch_dc[i] = sum(tmp_dc*tmp_areas)/sum(tmp_areas)

  
  tmpstd_dc = numeric(length(polsstd@polygons[[1]]@Polygons))
  tmpstd_areas = numeric(length(polsstd@polygons[[1]]@Polygons))
  for(k in seq_along(polsstd@polygons[[1]]@Polygons)){
    tmpstd_areas[k] = polsstd@polygons[[1]]@Polygons[[k]]@area
    tmpstd_dc[k] = mean(distCoast(lon = polsstd@polygons[[1]]@Polygons[[k]]@coords[,1], 
                               lat = polsstd@polygons[[1]]@Polygons[[k]]@coords[,2]))
  }
  
  hotpatchstd_dc[i] = sum(tmpstd_dc*tmpstd_areas)/sum(tmpstd_areas)
  
    
  print(i)
  
}

# f1 = 1
# l1 = 35
# f2 = 36
# l2 = 43

outDat = data.frame(survey = nameSurvey,
                    Level_agreg = levelhot,
                    Area_agreg = areahot,
                    Area_agreg_std = areahot_std,
                    Npatch_agreg = npatchhot,
                    Npatch_agreg_std = npatchhot_std,
                    DistCoast_agreg = hotpatch_dc,
                    DistCoast_agreg_std = hotpatchstd_dc,
                    Area_tot = arealow,
                    sAplus = saplus,
                    sAtot = satot,
                    Isotropia = isotropy,
                    Npatch_low = npatchlow)

write.csv(outDat, "SpatialIndicators.csv", row.names = F)

# AREA ESTANDAR : 84578.35

