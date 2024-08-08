

################################################################################################
### Title: Past Human Settlements in the Tropical Forest of Borneo Reveal Distinct Climate, Soil, and Settlement Patterns                   
### Journal: Plants People Planet August, 2024
### Author: Nathalia Perez Cardenas et al.
### Institution: Gepgraphy Department, University of Zurich
### Summary of the code and functions used in the development of the manuscript specified above    
### the original data is not provided and may be accessed following the Data Availability statement
####################################################################################################


library(terra)
library(elevatr)
library(geodata)
library(dplyr)
library(raster)
library(sf)
library(MultiscaleDTM)    ####### Pending
library(sp)

rm(list=ls());gc(verbose=T);gcinfo(FALSE)

setwd("/Users/your environment")

# Establish the path of the documents

pth_bor <- "Borneo_mask.shp"
pth_utm <- "UTM_Borneo_Merge_WGS84.shp"
pth_riv <- "River_Borneo_clip.shp"
pth_sites <- "sites_6ka_CALIBRATED.shp"

# read the documents using the sf package 

bor <- st_read(pth_bor)
bor <- st_cast(bor, "MULTILINESTRING")
utm <- st_read(pth_utm)
riv <- st_read(pth_riv)
arc_sites <- st_read(pth_sites)
names6ka <- arc_sites$name

#######################################################
#### 1. Get Environmental data for archaeological sites ##
#######################################################

sites6ka <- arc_sites

# Elevation with Geodata package

elev <- elevation_global(res = 2.5, path="Rdata/")
values.elev <- extract(elev, sites6ka)
sites.elev <- cbind.data.frame(st_coordinates(sites6ka), values.elev)
sites.elev <- cbind.data.frame(sites.elev, names6ka)

rm(values.elev, elev)

##########################################
# Get Soil texture: Sand, Silt and Clay  #####
##########################################

sand <- soil_world(var="sand", depth=15, stat="mean", path="Rdata/")
values_snd <- terra::extract(sand, sites6ka)
values_snd <- dplyr::select(values_snd, "sand_5-15cm")
rm(sand)

silt <- soil_world(var="silt", depth=15, stat="mean", path="Rdata/")
values_silt <- terra::extract(silt, sites6ka)
values_silt <- dplyr::select(values_silt, "silt_5-15cm")
rm(silt)

clay <- soil_world(var="clay", depth=15, stat="mean", path="Rdata/")
values_clay <- terra::extract(clay, sites6ka)
values_clay <-  dplyr::select(values_clay, "clay_5-15cm")
rm(clay)

#####################################
## Get Climatic data ###############
####################################

## Annual Mean Temperature (bio1) ##########
############################################

bio = geodata::worldclim_global("bio", res = 2.5,  path=tempdir(), version = "2.1")

tmp_av <- bio[[1]]
valuestmp <- terra::extract(tmp_av, sites6ka)  
valuestmp <- valuestmp %>% rename(bio1 = wc2.1_2.5m_bio_1) %>%
  dplyr::select("bio1") 
rm(tmp_av)

## Annual Precipitation (bio12)   ##########
############################################

pcp_av <- bio[[12]]
valuespcp_av <- terra::extract(pcp_av, sites6ka)  
valuespcp_av <- valuespcp_av %>% rename(bio12 = wc2.1_2.5m_bio_12) %>%
  dplyr::select("bio12")
rm(pcp_av)

## Precipitation Seasonality (Coefficient of Variation) (bio15)  ######
#######################################################################

pcp_ses <- bio[[15]]
valuespcp_ses <- terra::extract(pcp_ses, sites6ka)  
valuespcp_ses <- valuespcp_ses %>% rename(bio15 = wc2.1_2.5m_bio_15) %>%
  dplyr::select("bio15")
rm(pcp_ses)

evm_data <- cbind(sites.elev, values_snd, values_silt, values_clay, valuestmp, valuespcp_av, valuespcp_ses)


# Distance to rivers and distance to coast ######################################################
############################################################################

# read the documents using the sf package 

for(i in 1:nrow(utm))
{
  tmp_riv <- st_intersection(riv, utm[i,])
  tmp_bor <- st_intersection(bor, utm[i,])
  tmp_sites <- st_intersection(arc_sites, utm[i,])
  
  tmp_id <- as.data.frame(utm[i,11])
  tmp_id <- as.character(tmp_id[1,1])
  
  if(substr(tmp_id, 3, 3)=="N")
  {
    epsg_code <- 32600+as.numeric(substr(tmp_id, 1, 2))
    rm(tmp_id)
  }else{
    epsg_code <- 32700+as.numeric(substr(tmp_id, 1, 2))
    rm(tmp_id)
  }
  
  tmp_riv <- st_transform(tmp_riv, crs=epsg_code)
  tmp_bor <- st_transform(tmp_bor, crs=epsg_code)
  tmp_sites <- st_transform(tmp_sites, crs=epsg_code)
  
  dist_tmp <- st_distance(tmp_sites, tmp_riv)
  dist_min <- apply(dist_tmp, 1, min)
  rm(dist_tmp)
  
  dist_bor <- st_distance(tmp_sites, tmp_bor)
  min_bor <- apply(dist_bor, 1, min)
  rm(dist_bor)
  
  tmp_sites$min_dist <- dist_min/1000
  tmp_sites$min_dist_bor <- min_bor/1000
  rm(dist_min, min_bor)
  tmp_sites <- as.data.frame(tmp_sites[,c("ID", "name", "min_dist", "min_dist_bor")]) #,c(2,13,15,16)
  
  if(i==1)
  {
    tmp_df <- tmp_sites
    rm(tmp_sites)
  }else{
    tmp_df <- rbind(tmp_df, tmp_sites)
    rm(tmp_sites)
  }
}

dist_dt <- rename(tmp_df, names6ka = name) 
env_arc_sites <- merge(evm_data, dist_dt, by=c("names6ka"))       ######################################## all_sites6ka <- 

all_arc_sites <- st_as_sf(env_arc_sites, coords = c("X", "Y"), crs = 4326) #   

##################################################################
### Get palaeoclimatic data for archaeological sites #############
##################################################################


pth_bio1_0 <- "CHELSA_TraCE21k_bio01_20_V1.0.tif"
pth_bio1_2 <- "CHELSA_TraCE21k_bio01_0_V1.0.tif"                  
pth_bio1_4 <- "CHELSA_TraCE21k_bio01_-20_V1.0.tif"
pth_bio1_6 <- "CHELSA_TraCE21k_bio01_-40_V1.0.tif"
pth_bio12_0 <- "CHELSA_TraCE21k_bio12_20_V1.0.tif"
pth_bio12_2 <- "CHELSA_TraCE21k_bio12_0_V1.0.tif"
pth_bio12_4 <- "CHELSA_TraCE21k_bio12_-20_V1.0.tif"
pth_bio12_6 <- "CHELSA_TraCE21k_bio12_-40_V1.0.tif"
pth_bio15_0 <- "CHELSA_TraCE21k_bio15_20_V1.0.tif"
pth_bio15_2 <- "CHELSA_TraCE21k_bio15_0_V1.0.tif"
pth_bio15_4 <- "CHELSA_TraCE21k_bio15_-20_V1.0.tif"
pth_bio15_6 <- "CHELSA_TraCE21k_bio15_-40_V1.0.tif"


bio1_0 <- rast(pth_bio1_0)
bio1_2 <- rast(pth_bio1_2)
bio1_4 <- rast(pth_bio1_4)
bio1_6 <- rast(pth_bio1_6)
bio12_0 <- rast(pth_bio12_0)
bio12_2 <- rast(pth_bio12_2)
bio12_4 <- rast(pth_bio12_4)
bio12_6 <- rast(pth_bio12_6)
bio15_0 <- rast(pth_bio15_0)
bio15_2 <- rast(pth_bio15_2)
bio15_4 <- rast(pth_bio15_4)
bio15_6 <- rast(pth_bio15_6)


############# Annual Mean Temperature (bio1)  #######
#####################################################

tem_0 <- terra::extract(bio1_0, sites6ka)  
tem_2 <- terra::extract(bio1_2, sites6ka)
tem_4 <- terra::extract(bio1_4, sites6ka)
tem_6 <- terra::extract(bio1_6, sites6ka)

paleotem <- cbind(tem_0, tem_2, tem_4, tem_6) %>%
  dplyr::select(c(1,2,4,6,8)) %>%
  rename(tem_0 = "CHELSA_TraCE21k_bio01_20_V1.0", tem_2 = "CHELSA_TraCE21k_bio01_0_V1.0", tem_4 = "CHELSA_TraCE21k_bio01_-20_V1.0", tem_6 = "CHELSA_TraCE21k_bio01_-40_V1.0")

############# Annual Mean precipitation (bio12)  ######
#####################################################

pre_0 <- terra::extract(bio12_0, sites6ka)  
pre_2 <- terra::extract(bio12_2, sites6ka)
pre_4 <- terra::extract(bio12_4, sites6ka)
pre_6 <- terra::extract(bio12_6, sites6ka)

paleopre <- cbind(pre_0, pre_2, pre_4, pre_6) %>%
  dplyr::select(c(1,2,4,6,8)) %>%
  rename(pre_0 = "CHELSA_TraCE21k_bio12_20_V1.0", pre_2 = "CHELSA_TraCE21k_bio12_0_V1.0", pre_4 = "CHELSA_TraCE21k_bio12_-20_V1.0", pre_6 = "CHELSA_TraCE21k_bio12_-40_V1.0")

############# Seasonality (bio15)  ######################
#####################################################

sas_0 <- terra::extract(bio15_0, sites6ka)  
sas_2 <- terra::extract(bio15_2, sites6ka)
sas_4 <- terra::extract(bio15_4, sites6ka)
sas_6 <- terra::extract(bio15_6, sites6ka)

paleosas <- cbind(sas_0, sas_2, sas_4, sas_6) %>%
  dplyr::select(c(1,2,4,6,8)) %>%
  rename(sas_0 = "CHELSA_TraCE21k_bio15_20_V1.0", sas_2 = "CHELSA_TraCE21k_bio15_0_V1.0", sas_4 = "CHELSA_TraCE21k_bio15_-20_V1.0", sas_6 = "CHELSA_TraCE21k_bio15_-40_V1.0")


paleo_arc_clm <- cbind(paleotem, paleopre, paleosas) %>%
  dplyr::select("ID" , "tem_0", "tem_2", "tem_4", "tem_6", "pre_0", "pre_2", "pre_4", "pre_6",  "sas_0", "sas_2", "sas_4", "sas_6")


############################################################################################################################
#  2. Create Random locations and extract values of environmental variables ###################################################### 
##############################################################################################################################

library(sf)
library(geodata)
library(terra)
library(dplyr)
library(raster)

#  Establish the path of the documents

pth_polygon <- "Bounding_polygon.shp"
          
pol <- st_read(pth_polygon)

num_iteraciones <- c(1:1000)


bb <- st_bbox(pol)

sand <- soil_world(var="sand", depth=15, stat="mean", path="Rdata/")
silt <- soil_world(var="silt", depth=15, stat="mean", path="Rdata/")
clay <- soil_world(var="clay", depth=15, stat="mean", path="Rdata/")

bio = geodata::worldclim_global("bio", res = 2.5,  path=tempdir(), version = "2.1")

tmp_av <- bio[[1]]
pcp_av <- bio[[12]]
pcp_ses <- bio[[15]]


for (i in num_iteraciones) 
{
  # Create random points   ###############################################################
  
  vec_x <- runif(500, 110.1173, 118.6195)
  vec_y <- runif(500, -2.808650, 5.814714)
  
  xny <- data.frame(cbind(vec_x, vec_y))
  
  names(xny) <- c("Lon", "Lat")
  
  tmp_sf <- st_as_sf(xny, coords = c("Lon", "Lat"), crs = 4326)
  pnt_in_pol <- st_intersection(tmp_sf, pol)
  rm(xny, vec_x, vec_y)
  
  t1=Sys.time()
  pnt_in_pol_rand <- pnt_in_pol[sample(1:nrow(pnt_in_pol), 47),]
  pnt_in_pol_rand$ID <- 1:nrow(pnt_in_pol_rand)
  
  # Then extract the values of the points: Elevation  #####################################################################
  
  sites <- pnt_in_pol_rand  
  rm(pnt_in_pol_rand, tmp_sf)
  
  # Elevation with Geodata package
  
  elev <- elevation_global(res = 2.5, path="Rdata/")
  values.elev <- terra::extract(elev, sites)
  sites.elev <- cbind.data.frame(st_coordinates(sites), values.elev)
  
  # Geomorphology: Sand, Silt and Clay  ##########################################################
  
  values_snd <- terra::extract(sand, sites)
  values_snd <- dplyr::select(values_snd, "sand_5-15cm")
  
  values_silt <- terra::extract(silt, sites)
  values_silt <- dplyr::select(values_silt, "silt_5-15cm")
  
  values_clay <- terra::extract(clay, sites)
  values_clay <-  dplyr::select(values_clay, "clay_5-15cm")
  
  sites.elev <- cbind(sites.elev, values_snd, values_silt, values_clay)
  names(sites.elev) <- c("X", "Y", "ID", "Elev", "Sand", "Silt", "Clay")
  rm(values_snd, values_silt, values_clay)
  
  ## Annual Mean Temperature (bio1)  ##########################################################
  valuestmp <- terra::extract(tmp_av, sites)  
  #valuestmp = valuestmp/10
  valuestmp <- valuestmp %>% rename(bio1 = wc2.1_2.5m_bio_1) %>%
    dplyr::select("bio1")
  
  ## Annual Precipitation (bio12)   ##########################################################
  
  valuespcp_av <- terra::extract(pcp_av, sites)  
  valuespcp_av <- valuespcp_av %>% rename(bio12 = wc2.1_2.5m_bio_12) %>%
    dplyr::select("bio12")
  
  ## Precipitation Seasonality (Coefficient of Variation) (bio15)  ##########################################################
  
  valuespcp_ses <- terra::extract(pcp_ses, sites)  
  valuespcp_ses <- valuespcp_ses %>% rename(bio15 = wc2.1_2.5m_bio_15) %>%
    dplyr::select("bio15")
  
  site.elev <- cbind(sites.elev, valuestmp,valuespcp_av, valuespcp_ses)
  rm(valuestmp,valuespcp_av, valuespcp_ses)
  
  site.elev$Grp_ID <- paste0("GRP_", i)
  
  if(i==1)
  {
    rndm_data_env <- site.elev
  }else{
    rndm_data_env <- rbind(rndm_data_env, site.elev)
  }
  t2=Sys.time()
  print(i)
  print(t2-t1)
}

all_sites <- st_as_sf(rndm_data_env, coords = c("X", "Y"), crs = 4326)

rm(bio,clay, elev, pcp_av, pcp_ses, pnt_in_pol,sand, silt, sit, site.elev,sites,sites.elev, tmp_av, values.elev, t1, t2)

# Distance to rivers and distance to coast  #########################################################################################


for(j in 1:nrow(utm))
{
  t1=Sys.time()
  tmp_bor <- st_intersection(bor, utm[j,])
  tmp_riv <- st_intersection(riv, utm[j,])
  tmp_sites <- st_intersection(all_sites, utm[j,])
  
  tmp_id <- as.data.frame(utm[j,11])
  tmp_id <- as.character(tmp_id[1,1])
  
  if(substr(tmp_id, 3, 3)=="N")
  {
    epsg_code <- 32600+as.numeric(substr(tmp_id, 1, 2))
    rm(tmp_id)
  }else{
    epsg_code <- 32700+as.numeric(substr(tmp_id, 1, 2))
    rm(tmp_id)
  }
  
  tmp_bor <- st_transform(tmp_bor, crs=epsg_code)
  tmp_riv <- st_transform(tmp_riv, crs=epsg_code)
  tmp_sites <- st_transform(tmp_sites, crs=epsg_code)
  
  dist_riv <- st_distance(tmp_sites, tmp_riv)
  min_riv <- apply(dist_riv, 1, min)
  rm(dist_riv)
  
  dist_bor <- st_distance(tmp_sites, tmp_bor)
  min_bor <- apply(dist_bor, 1, min)
  rm(dist_bor)
  
  tmp_sites$min_dist_river <- min_riv
  tmp_sites$min_dist_bor <- min_bor
  rm(min_riv, min_bor)
  tmp_sites <- as.data.frame(tmp_sites[,c("ID","Grp_ID" ,"layer", "min_dist_river", "min_dist_bor")])
  tmp_sites <- tmp_sites[,-6]
  
  if(j==1)
  {
    tmp_min <- tmp_sites
    rm(tmp_sites)
  }else{
    tmp_min <- rbind(tmp_min, tmp_sites)
    rm(tmp_sites)
  }
  rm(tmp_riv, tmp_bor)
  t2=Sys.time()
  print(j)
  print(t2-t1)
}

rndm_data_env <- merge(rndm_data_env, tmp_min, all.x=TRUE, by=c("ID", "Grp_ID"))
rm(tmp_min, all_sites)
all_rndm_sites <- st_as_sf(rndm_data_env, coords = c("X", "Y"), crs = 4326)


############################################################################################################################
#  Extract palaeoclimate data for Random locations  ###################################################### 
##############################################################################################################################


#pht_rnd_sites <- "M:/group/ess_data/PerezNathalia/Archaeological-settlements/Worshop_project/Outcome/Random_Sites_Env_Variables_final6ka47.shp"   
rnd_sites <- all_rndm_sites  #vect(pht_rnd_sites)

rnd_sites_dt <- read.csv("Random_Sites_Env_Variables_final6ka47.csv")   #### or env_arc_sites 
rnd_sites_dt_sbst <- dplyr::select(rnd_sites_dt, 1,2)

############# Annual Mean Temperature (bio1)  ##########################################################
#################################################################################################


tem_0 <- terra::extract(bio1_0, rnd_sites)  
tem_2 <- terra::extract(bio1_2, rnd_sites)
tem_4 <- terra::extract(bio1_4, rnd_sites)
tem_6 <- terra::extract(bio1_6, rnd_sites)

paleotem <- cbind(tem_0, tem_2, tem_4, tem_6, rnd_sites_dt_sbst) %>%
  dplyr::select(c(2,4,6,8:10)) %>%
  rename(tem_0 = "CHELSA_TraCE21k_bio01_20_V1.0", tem_2 = "CHELSA_TraCE21k_bio01_0_V1.0", tem_4 = "CHELSA_TraCE21k_bio01_-20_V1.0", tem_6 = "CHELSA_TraCE21k_bio01_-40_V1.0")

############# Annual Mean precipitation (bio12)  ##########################################################
#################################################################################################

pre_0 <- terra::extract(bio12_0, rnd_sites)  
pre_2 <- terra::extract(bio12_2, rnd_sites)
pre_4 <- terra::extract(bio12_4, rnd_sites)
pre_6 <- terra::extract(bio12_6, rnd_sites)

paleopre <- cbind(pre_0, pre_2, pre_4, pre_6) %>%
  dplyr::select(c(1,2,4,6,8)) %>%
  rename(pre_0 = "CHELSA_TraCE21k_bio12_20_V1.0", pre_2 = "CHELSA_TraCE21k_bio12_0_V1.0", pre_4 = "CHELSA_TraCE21k_bio12_-20_V1.0", pre_6 = "CHELSA_TraCE21k_bio12_-40_V1.0")

############# Seasonality (bio15)  ##########################################################
#################################################################################################

sas_0 <- terra::extract(bio15_0, rnd_sites)  
sas_2 <- terra::extract(bio15_2, rnd_sites)
sas_4 <- terra::extract(bio15_4, rnd_sites)
sas_6 <- terra::extract(bio15_6, rnd_sites)

paleosas <- cbind(sas_0, sas_2, sas_4, sas_6) %>%
  dplyr::select(c(1,2,4,6,8)) %>%
  rename(sas_0 = "CHELSA_TraCE21k_bio15_20_V1.0", sas_2 = "CHELSA_TraCE21k_bio15_0_V1.0", sas_4 = "CHELSA_TraCE21k_bio15_-20_V1.0", sas_6 = "CHELSA_TraCE21k_bio15_-40_V1.0")


paleo_rnd_clm <- cbind(paleotem, paleopre, paleosas) %>%
  dplyr::select("ID" , "Grp_ID", "tem_0", "tem_2", "tem_4", "tem_6", "pre_0", "pre_2", "pre_4", "pre_6",  "sas_0", "sas_2", "sas_4", "sas_6")


##################################################################
#### 3. Randomization test   ######################################
###############################################################

## Random mean and mean of archaeological sites

pth_rnd_pnt <- "Arc_Sites_Env_Variables_6ka_dist_cec_47000.CSV" #Random_Sites_Env_Variables_final_2.csv"
pth_arc_pnt <- "Arc_sites_env_variables_6ka_dist_cec_calib.csv" #Environmental_dt_for_points.csv"

rnd_pnt <- read.csv(pth_rnd_pnt)
arc_pnt <- read.csv(pth_arc_pnt)

 #### rndm_data_env  ### This table has been created above and can we replace in the line 477. It does'nt count with the cec column. 
 #### env_arc_sites  ### This table has been created above and can we replace in the line 478. It does'nt count with the cec column.

## Random mean and all distribution archaeological values plot

# Harmonize the data tables

rnd_pnt <- rnd_pnt %>% 
  mutate("Group" = "rdm") %>%
  mutate(mn_dst_r = mn_dst_r/1000,
         mn_dst_b = mn_dst_b/1000) %>%
  dplyr::select("Grp_ID", "Elev", "Sand",  "Silt", "Clay", "bio1", "bio12",  "bio15",  "mn_dst_r", "mn_dst_b", "merge_all_cec_borneo") %>% 
  rename(min_dist_river = mn_dst_r,
         min_dist_bor = mn_dst_b,
         cec_value  =  merge_all_cec_borneo) %>%
  na.omit()


mn_tbl <- rnd_pnt %>% group_by(Grp_ID) %>% 
  summarize(Elev = mean(Elev), Sand = mean(Sand), Silt = mean(Silt), Clay = mean(Clay), bio1 = mean(bio1), bio12 = mean(bio12), bio15 = mean(bio15),
            min_dist_river = mean(min_dist_river), min_dist_bor = mean(min_dist_bor), cec_value = mean(cec_value)) %>%
  dplyr::select("Elev", "Sand",  "Silt", "Clay", "bio1", "bio12",  "bio15",  "min_dist_river", "min_dist_bor", "cec_value") 

#########

arc_pnt <- arc_pnt %>% 
  mutate(Grp_ID = "Obs") %>%
  rename(cec_value = merge_all_cec_borneo,
         min_dist_river = min_dst,
         min_dist_bor = mn_dst_) %>%
  dplyr::select("Grp_ID", "Elev", "Sand",  "Silt", "Clay", "bio1", "bio12",  "bio15",  "min_dist_river", "min_dist_bor", "cec_value")


arch_mn_tbl <- arc_pnt %>% group_by(Grp_ID) %>% 
  summarize(Elev = mean(Elev), Sand = mean(Sand), Silt = mean(Silt), Clay = mean(Clay), bio1 = mean(bio1), bio12 = mean(bio12), bio15 = mean(bio15),
            min_dist_river = mean(min_dist_river), min_dist_bor = mean(min_dist_bor), cec_value = mean(cec_value)) %>%
  select(2:11)


mn_tbl <- rbind(mn_tbl, arch_mn_tbl )

#### -test

variables <- c("Elev"  , "Sand"  , "Silt"  , "Clay" , "bio1"  , "bio12" , "bio15" , "min_dist_river" ,"min_dist_bor" , "cec_value")
prop_exc <- data.frame("Variable" =  variables, "pexc_value" = numeric(length(variables)))


prop_exc$pexc_value[prop_exc$Variable == "Elev"] <- sum(mn_tbl$Elev[1:1000]<= mn_tbl$Elev[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "Sand"] <- 1-(sum(mn_tbl$Sand[1:1000]<= mn_tbl$Sand[1001])/1000)  
prop_exc$pexc_value[prop_exc$Variable == "Silt"] <- 1-sum(mn_tbl$Silt[1:1000]<= mn_tbl$Silt[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "Clay"] <- sum(mn_tbl$Clay[1:1000]<= mn_tbl$Clay[1001])/1000  
prop_exc$pexc_value[prop_exc$Variable == "bio1"] <- 1-sum(mn_tbl$bio1[1:1000]<= mn_tbl$bio1[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "bio12"] <- sum(mn_tbl$bio12[1:1000]<= mn_tbl$bio12[1001])/1000  
prop_exc$pexc_value[prop_exc$Variable == "bio15"] <- 1-sum(mn_tbl$bio15[1:1000]<= mn_tbl$bio15[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "min_dist_river"] <- sum(mn_tbl$min_dist_river[1:1000]<= mn_tbl$min_dist_river[1001])/1000 
prop_exc$pexc_value[prop_exc$Variable == "min_dist_bor"] <- sum(mn_tbl$min_dist_bor[1:1000]<= mn_tbl$min_dist_bor[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "cec_value"] <- sum(mn_tbl$cec_value[1:1000]<= mn_tbl$cec_value[1001])/1000 


##################################################################################
###### Paleoclimate Random mean and mean of Paleoclimate archaeological site 
#####################################################################################

#pth_paleo_rnd_pnt <- "paleo_rnd_clm47.CSV" #Random_Sites_Env_Variables_final_2.csv"
#pth_paleo_arc_pnt <- "paleo_arc_clmCALIB.csv" #Environmental_dt_for_points.csv"

#rnd_pnt <- read.csv(pth_paleo_rnd_pnt)  # paleo_rnd_clm
#arc_pnt <- read.csv(pth_paleo_arc_pnt)  # paleo_arc_clm

rnd_pnt <- paleo_rnd_clm  
arc_pnt <- paleo_arc_clm

# Harmonize the data tables

rnd_pnt <- rnd_pnt %>% 
  dplyr::select( "Grp_ID", "tem_0" , "tem_2" , "tem_4" , "tem_6" , "pre_0" , "pre_2" , "pre_4" , "pre_6" , "sas_0" , "sas_2" , "sas_4"  ,"sas_6") %>% 
  na.omit()


mn_tbl <- rnd_pnt %>% group_by(Grp_ID) %>% 
  summarize(tem_0 = mean(tem_0), tem_2 = mean(tem_2), tem_4 = mean(tem_4), tem_6 = mean(tem_6), pre_0 = mean(pre_0), pre_2 = mean(pre_2), pre_4 = mean(pre_4),
            pre_6 = mean(pre_6), sas_0 = mean(sas_0), sas_2 = mean(sas_2), sas_4 = mean(sas_4), sas_6 = mean(sas_6)) %>%
  dplyr::select("tem_0" , "tem_2" , "tem_4" , "tem_6" , "pre_0" , "pre_2" , "pre_4" , "pre_6" , "sas_0" , "sas_2" , "sas_4"  ,"sas_6") 

#########

arc_pnt <- arc_pnt %>% 
  mutate(Grp_ID = "Obs") %>%
  dplyr::select("Grp_ID", "tem_0" , "tem_2" , "tem_4" , "tem_6" , "pre_0" , "pre_2" , "pre_4" , "pre_6" , "sas_0" , "sas_2" , "sas_4"  ,"sas_6")


arch_mn_tbl <- arc_pnt %>% group_by(Grp_ID) %>% 
  summarize(tem_0 = mean(tem_0), tem_2 = mean(tem_2), tem_4 = mean(tem_4), tem_6 = mean(tem_6), pre_0 = mean(pre_0), pre_2 = mean(pre_2), pre_4 = mean(pre_4),
            pre_6 = mean(pre_6), sas_0 = mean(sas_0), sas_2 = mean(sas_2), sas_4 = mean(sas_4), sas_6 = mean(sas_6))  %>%
  select(2:13)


mn_tbl <- rbind(mn_tbl, arch_mn_tbl )

#### -test

variables <- c("tem_0" , "tem_2" , "tem_4" , "tem_6" , "pre_0" , "pre_2" , "pre_4" , "pre_6" , "sas_0" , "sas_2" , "sas_4"  ,"sas_6")
prop_exc <- data.frame("Variable" =  variables, "pexc_value" = numeric(length(variables)))


prop_exc$pexc_value[prop_exc$Variable == "tem_0"] <- 1-sum(mn_tbl$tem_0[1:1000]<= mn_tbl$tem_0[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "tem_2"] <- 1-sum(mn_tbl$tem_2[1:1000]<= mn_tbl$tem_2[1001])/1000  
prop_exc$pexc_value[prop_exc$Variable == "tem_4"] <- 1-sum(mn_tbl$tem_4[1:1000]<= mn_tbl$tem_4[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "tem_6"] <- 1-sum(mn_tbl$tem_6[1:1000]<= mn_tbl$tem_6[1001])/1000  
prop_exc$pexc_value[prop_exc$Variable == "pre_0"] <- sum(mn_tbl$pre_0[1:1000]<= mn_tbl$pre_0[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "pre_2"] <- sum(mn_tbl$pre_2[1:1000]<= mn_tbl$pre_2[1001])/1000  
prop_exc$pexc_value[prop_exc$Variable == "pre_4"] <- sum(mn_tbl$pre_4[1:1000]<= mn_tbl$pre_4[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "pre_6"] <- sum(mn_tbl$pre_6[1:1000]<= mn_tbl$pre_6[1001])/1000 
prop_exc$pexc_value[prop_exc$Variable == "sas_0"] <- 1-sum(mn_tbl$sas_0[1:1000]<= mn_tbl$sas_0[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "sas_2"] <- 1-sum(mn_tbl$sas_2[1:1000]<= mn_tbl$sas_2[1001])/1000 
prop_exc$pexc_value[prop_exc$Variable == "sas_4"] <- 1-sum(mn_tbl$sas_4[1:1000]<= mn_tbl$sas_4[1001])/1000
prop_exc$pexc_value[prop_exc$Variable == "sas_6"] <- 1-sum(mn_tbl$sas_6[1:1000]<= mn_tbl$sas_6[1001])/1000 


#############################################################################################
########### Distribution values plot
#############################################################################################


#################  Current conditions plot
#######################################################

labels <- c(bio1_mn = "MAT Current", bio12_mn= "AP Current", bio15_mn = "PS Current", 
            clay_mn = "% Clay",   sand_mn = "% Sand", cec_mn = "CEC", 
            elev_mn = "Elevation",  riv_mn_km = "River Distance (km)", 
            coast_mn_km = "Coast Distance (km)", silt_mn = "% Silt")

density <- ggplot(data_long_rnd, aes(x= value, fill = type_variable)) +  #colour = type_variable,
  geom_density(alpha=0.3) + 
  scale_fill_manual(values = c("#52854C", "#4E84C4", "#440154FF")) +
  facet_wrap(~ factor(name, c("bio1_mn", "clay_mn", "coast_mn_km", "silt_mn", "bio12_mn",  "sand_mn",  "riv_mn_km", "elev_mn", "bio15_mn", "cec_mn")), scales = "free", labeller = as_labeller(labels), ncol = 4) + 
  geom_vline(data = data_long_arc,
             aes(xintercept = value), color ="brown1", linetype="dashed", linewidth=0.9) +
  labs(y = "Density") +
  theme(strip.text = element_text(size=13))


density + geom_vline(data= mean_values, aes(xintercept = mean_value, color = type_variable), linetype = "dashed", size = 0.9) +
  scale_color_manual(values = c("#52854C", "#4E84C4", "#440154FF")) #+
#theme(legend.position = c(1.00, 0.13)) +
#theme(legend.key.size = unit(1.3, 'cm')) +
#theme(legend.text = element_text(size=12)) +
#theme(legend.title = element_text(size=12)) 


#ggsave("evn_pnts_line6katest_HOR_CALIB.png",   width=18, height=12, dpi = 300)


#############################   Palaeoclimate conditios Plot
##############################################################

labels <- c(tem_0_mn = "MAT 0BP", tem_2_mn = "MAT 2ka BP", tem_4_mn = "MAT 4ka BP", tem_6_mn = "MAT 6ka BP", 
            pre_0_mn= "AP 0BP", pre_2_mn= "AP 2ka BP", pre_4_mn= "AP 4ka BP", pre_6_mn= "AP 6ka BP",
            sas_0_mn = "PS 0BP", sas_2_mn = "PS 2ka BP", sas_4_mn = "PS 4ka BP", sas_6_mn = "PS 6ka BP")

density <- ggplot(data_long_rnd, aes(x= value, fill = type_variable)) +  #colour = type_variable,
  geom_density(alpha=0.3) + 
  scale_fill_manual(values = c("#52854C", "#4E84C4", "#440154FF")) +
  facet_wrap(~ factor(name, c("tem_6_mn", "tem_4_mn",  "tem_2_mn" , "tem_0_mn", "pre_6_mn", "pre_4_mn" ,  "pre_2_mn" , "pre_0_mn" ,"sas_6_mn" , "sas_4_mn" , "sas_2_mn" , "sas_0_mn")), scales = "free", labeller = as_labeller(labels), ncol = 4) + 
  geom_vline(data = data_long_arc,
             aes(xintercept = value), color ="brown1", linetype="dashed", linewidth=0.9) +
  labs(y = "Density") +
  theme(strip.text = element_text(size=13))


density + geom_vline(data= mean_values, aes(xintercept = mean_value, color = type_variable), linetype = "dashed", size = 0.9) +
  scale_color_manual(values = c("#52854C", "#4E84C4", "#440154FF")) #+
#theme(legend.position = c(0.78, 0.13)) +
#theme(legend.key.size = unit(1.3, 'cm')) +
#theme(legend.text = element_text(size=12)) +
#theme(legend.title = element_text(size=12)) 


#ggsave("paleo_evn_pnts_line6ka_HORCALIB.png",  width=18, height=12, dpi = 300)


## plot supplementary material histogramas #################################

arc_pnt <- arc_pnt  %>% 
  select(1:8, 11:13)


labels <- c(bio1 = "MAT Current", bio12= "AP Current", bio15 = "PS Current", 
            Clay = "% Clay",  Silt = "% Silt", Sand = "% Sand", cec_arc_po = "Cec", 
            Elev = "Elevation",  min_dist_river  = "River Distance (km)", 
            min_dist_bor = "Coast Distance (km)")


p <- arc_pnt %>% 
  pivot_longer(-Grp_ID) %>% 
  ggplot(aes(x = value, fill=name))+
  geom_histogram(bins = 47) +
  facet_wrap(~ factor(name, c("bio1", "bio12", "bio15", "Clay", "Silt", "Sand", "cec_arc_po", "min_dist_river", "min_dist_bor", "Elev")), scales = "free", labeller = as_labeller(labels), ncol = 3)


p + scale_fill_brewer(palette = "Spectral", labels = labels) +
  theme(strip.text = element_text(size = 11))


#ggsave("hist_evn_arc_pnts6kaCALIB.png",  width=12, height=18, dpi = 300)


################################################################################################################################################
############# 4. Comparion with History database of the Global Environment HYDE ###################################################
###################################################################################################################################


# Libraries

library(dplyr)
library(ggplot2)
library(sf)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(terra)


###########  Data sites for join with the index 

pth_sites <- "Archaeological_sites.shp"
arc_sites <- st_read(pth_sites)
subset_arc_sites <- arc_sites[, c("Name", "geometry")]
df_arc_sites <- as.data.frame(subset_arc_sites)

df_arc_sites$long <- st_coordinates(df_arc_sites$geometry)[, 1]
df_arc_sites$lat <- st_coordinates(df_arc_sites$geometry)[, 2]
df_arc_sites <- dplyr::select(df_arc_sites, c(1,3,4))

df_arc_sites <- df_arc_sites %>%
  mutate(has_dot = grepl("\\.", Name)) %>%
  transmute(
    number = ifelse(has_dot, as.integer(substring(Name, 1, regexpr("\\.", Name) - 1)), NA),
    name = ifelse(has_dot, trimws(substring(Name, regexpr("\\.", Name) + 1, nchar(Name))), trimws(Name)),
    long,
    lat
  ) %>%
  dplyr::select(number, name, long, lat)

df_arc_sites$ID <- 1:73

df_arc_sites <- dplyr::select(df_arc_sites, c("ID", "name", "number", "long", "lat")) %>%
  rename(old_ID = number)


#################################### Archaeological sites 6 ka BP

pth_sites6ka <- "sites_6ka_CALIBRATED.csv"
sites6ka <- read.csv(pth_sites6ka)
sites6ka <- rename(sites6ka, name = site2, old_ID = ID)

sites6kacoor <- merge(sites6ka, df_arc_sites, by="old_ID")

sites6kacoor <- dplyr::select(sites6kacoor, c("old_ID", "name.x", "presencia", "cantidad", "long", "lat", "ID")) %>%
  rename(name = name.x)

############################################## Data antiquity index for mapping Hyde and observations


conv_pth <- "converted_rangeland_archaeological_points.xlsx"
graz_pth <- "grazing_archaeological_points.xlsx"
past_pth <- "pasture_archaeological_points.xlsx"
popd_pth <- "popdensity_archaeological_points.xlsx"

conv_arc <- read_xlsx(conv_pth)
graz_arc <- read_xlsx(graz_pth)
past_arc <- read_xlsx(past_pth)
popd_arc <- read_xlsx(popd_pth)


new_names <- c("ID" , "Name" , "10000BC", "9000BC",  "8000BC",  "7000BC",  "6000BC" , "5000BC" , "4000BC" , "3000BC" , "2017AD" , "2016AD" , "2015AD" , "2014AD", "2013AD" ,
               "2012AD" , "2011AD" , "2010AD" , "2009AD",  "2008AD" , "2007AD" , "2006AD" , "2005AD" , "2004AD" , "2003AD" , "2002AD" ,"2001AD" , "2000BC",  "2000AD",  "1990AD", 
               "1980AD" , "1970AD" , "1960AD" , "1950AD" , "1940AD" , "1930AD",  "1920AD",  "1910AD" , "1900AD" , "1890AD",  "1880AD" , "1870AD"  ,"1860AD" , "1850AD" , "1840AD" ,
               "1830AD",  "1820AD" , "1810AD" , "1800AD" , "1790AD" , "1780AD" , "1770AD"  ,"1760AD" , "1750AD" , "1740AD" , "1730AD" , "1720AD" , "1710AD" , "1700AD" , "1600AD" ,
               "1500AD",  "1400AD" , "1300AD" , "1200AD" , "1100AD" , "1000BC" , "1000AD",  "900AD"  , "800AD"  , "700AD" ,  "600AD" ,  "500AD" ,  "400AD" ,  "300AD" ,  "200AD"  ,
               "100AD" , "0AD", "landcover")

conv_arc <- conv_arc %>%
  mutate(landcover = "converted_rangeland") 
colnames(conv_arc) <- new_names

conv_arc <- conv_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

conv_arc <- conv_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

conv_arc$year_n <- as.numeric(conv_arc$year_n)
conv_arc$ID <- as.factor(conv_arc$ID)

conv_arc <- conv_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

conv_arc_ant <- conv_arc %>%
  filter(value > 0) %>%
  filter(before_present < 6000) %>%
  group_by(ID) %>%
  filter(before_present == max(before_present)) %>%
  ungroup() %>%
  arrange(desc(before_present)) 

#antq_index <- merge(df_arc_sites, conv_arc_ant, by= "ID")

antq_index <- merge(x=sites6kacoor, y=conv_arc_ant, by = "ID", all.x=T)

#comvertir a st y visualizar en qgis

antq_index <- st_as_sf(antq_index, coords = c("long", "lat"), crs = 4326)
#st_write(antq_index, "conv_antq_index6kav2.shp")

################################

past_arc <- past_arc %>%
  mutate(landcover = "pasture") 
colnames(past_arc) <- new_names

past_arc <- past_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

past_arc <- past_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

past_arc$year_n <- as.numeric(past_arc$year_n)
past_arc$ID <- as.factor(past_arc$ID)

past_arc <- past_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

past_arc_ant <- past_arc %>%
  filter(value > 0) %>%
  filter(before_present < 6000) %>%
  group_by(ID) %>%
  filter(before_present == max(before_present)) %>%
  ungroup() %>%
  arrange(desc(before_present)) 

#past_antq_index <- merge(df_arc_sites, past_arc_ant, by= "ID")

past_antq_index <- merge(x=sites6kacoor, y=past_arc_ant, by = "ID", all.x=T)
#comvertir a st y visualizar en qgis

past_antq_index <- st_as_sf(past_antq_index, coords = c("long", "lat"), crs = 4326)
#st_write(past_antq_index, "past_antq_index6kav2.shp")


###############################################################################################

graz_arc <- graz_arc %>%
  mutate(landcover = "grazing") 
colnames(graz_arc) <- new_names

graz_arc <- graz_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

graz_arc <- graz_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

graz_arc$year_n <- as.numeric(graz_arc$year_n)
graz_arc$ID <- as.factor(graz_arc$ID)

graz_arc <- graz_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

graz_arc_ant <- graz_arc %>%
  filter(value > 0) %>%
  filter(before_present < 6000) %>%
  group_by(ID) %>%
  filter(before_present == max(before_present)) %>%
  ungroup() %>%
  arrange(desc(before_present)) 

#graz_antq_index <- merge(df_arc_sites, graz_arc_ant, by= "ID")

graz_antq_index <- merge(x=sites6kacoor, y=graz_arc_ant, by = "ID", all.x=T)

graz_antq_index <- st_as_sf(graz_antq_index, coords = c("long", "lat"), crs = 4326)
#st_write(graz_antq_index, "graz_antq_index6kav2.shp")


######################################################################################################

popd_arc <- popd_arc %>%
  mutate(landcover = "popdensity") 
colnames(popd_arc) <- new_names

popd_arc <- popd_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

popd_arc <- popd_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

popd_arc$year_n <- as.numeric(popd_arc$year_n)
popd_arc$ID <- as.factor(popd_arc$ID)

popd_arc <- popd_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

popd_arc_ant <- popd_arc %>%
  filter(value > 0) %>%
  filter(before_present < 6000) %>%
  group_by(ID) %>%
  filter(before_present == max(before_present)) %>%
  ungroup() %>%
  arrange(desc(before_present)) 

#popd_antq_index <- merge(df_arc_sites, popd_arc_ant, by= "ID")

popd_antq_index <- merge(x=sites6kacoor, y=popd_arc_ant, by = "ID", all.x=T)

popd_antq_index <- st_as_sf(popd_antq_index, coords = c("long", "lat"), crs = 4326)
#st_write(popd_antq_index, "popd_antq_index6kav2.shp")


##############################################

obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)

# Iterar sobre los nombres de las columnas desde la cuarta hasta la 82
for (i in 4:82) {
  # Reemplazar "X" con una cadena vacía en cada nombre de columna
  column_names[i] <- gsub("X", "", column_names[i])
}

# Asignar los nuevos nombres de columna al dataframe
names(obs_arc) <- column_names

obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  dplyr::select(c(2,3, 21:80,82, 83)) %>%
  pivot_longer(cols = c(3:63),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)

obs_arc_ant <- obs_arc %>%
  #filter(year <= 10000) %>%
  filter(value != 0) %>%
  group_by(nw_id) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  arrange(desc(year)) %>%
  rename(old_id = ID) %>%
  rename(ID = nw_id)

obs_antq_index <- merge(sites6kacoor, obs_arc_ant, by= "ID") %>%
  dplyr::select("ID","old_ID", "name" , "presencia" ,"cantidad" , "long" , "lat" ,"year"  ,"value") %>%
  rename(bfr_prs = year)

obs_antq_index <- st_as_sf(obs_antq_index, coords = c("long", "lat"), crs = 4326)
#st_write(obs_antq_index, "obs_antq_index6kav2Calib.shp")


##############################################################################################################
############################################## Data Ocupation index for mapping Hyde and observations
##############################################################################################################


conv_pth <- "converted_rangeland_archaeological_points.xlsx"
graz_pth <- "grazing_archaeological_points.xlsx"
past_pth <- "pasture_archaeological_points.xlsx"
popd_pth <- "popdensity_archaeological_points.xlsx"

conv_arc <- read_xlsx(conv_pth)
graz_arc <- read_xlsx(graz_pth)
past_arc <- read_xlsx(past_pth)
popd_arc <- read_xlsx(popd_pth)

conv_arc <- conv_arc %>%
  mutate(landcover = "converted_rangeland") 
colnames(conv_arc) <- new_names

conv_arc <- conv_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

conv_arc <- conv_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

conv_arc$ID <- as.factor(conv_arc$ID)
conv_arc$year_n <- as.numeric(conv_arc$year_n)

conv_arc <- conv_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

# conv_arc_occ <- conv_arc %>%
#   filter(year_n < 5000) %>%
#   mutate(occ = ifelse(value > 0, 1, 0)) %>%
#   group_by(ID) %>%
#   summarise(occ_tot = sum(occ))

# conv_occ_index <- merge(sites6kacoor, conv_arc_occ, by= "ID")
# 
# conv_occ_index <- st_as_sf(conv_occ_index, coords = c("long", "lat"), crs = 4326)
# st_write(conv_occ_index, "conv_occ_index6kav2.shp")

TEST <- conv_arc %>%
  filter(value != 0) %>%
  filter(before_present < 6000) %>%
  arrange(ID, desc(before_present)) %>%
  group_by(ID)

data <- TEST %>%
  arrange(ID, before_present)

# Difference between adjacent rows per site (ID)
data <- data %>%
  group_by(ID) %>%
  mutate(time_diff = lead(before_present, default = last(before_present)) - before_present) %>%
  ungroup()

timetot <- data %>%
  group_by(ID) %>%
  summarise(total_time_habitated = sum(time_diff, na.rm = TRUE))

# Print the result
print(timetot)

conv_occ_index <- merge(sites6kacoor, timetot, by= "ID")

conv_occ_index <- st_as_sf(conv_occ_index, coords = c("long", "lat"), crs = 4326)
st_write(conv_occ_index, "conv_occ_index6kav2.shp")


################################################################################################


graz_arc <- graz_arc %>%
  mutate(landcover = "converted_rangeland") 
colnames(graz_arc) <- new_names

graz_arc <- graz_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

graz_arc <- graz_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

graz_arc$year_n <- as.numeric(graz_arc$year_n)
graz_arc$ID <- as.factor(graz_arc$ID)

graz_arc <- graz_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

# graz_arc_occ <- graz_arc %>%
#   filter(year_n < 5000) %>%
#   mutate(occ = ifelse(value > 0, 1, 0)) %>%
#   group_by(ID) %>%
#   summarise(occ_tot = sum(occ))

graz_arc_occ <- graz_arc %>%
  filter(value != 0) %>%
  filter(before_present < 6000) %>%
  arrange(ID, desc(before_present)) %>%
  group_by(ID)

graz_arc_occ <- graz_arc_occ %>%
  arrange(ID, before_present)

# Difference between adjacent rows per site (ID)
graz_arc_occ <- graz_arc_occ %>%
  group_by(ID) %>%
  mutate(time_diff = lead(before_present, default = last(before_present)) - before_present) %>%
  ungroup()

timetot <- graz_arc_occ %>%
  group_by(ID) %>%
  summarise(total_time_habitated = sum(time_diff, na.rm = TRUE))

# Print result
print(timetot)

graz_occ_index <- merge(sites6kacoor, timetot, by= "ID")

graz_occ_index <- st_as_sf(graz_occ_index, coords = c("long", "lat"), crs = 4326)
st_write(graz_occ_index, "graz_occ_index6kav2.shp")


#########################################

 past_arc <- past_arc %>%
   mutate(landcover = "converted_rangeland") 
 colnames(past_arc) <- new_names
 
 past_arc <- past_arc %>% 
   select(c(1,3:78)) %>%
   pivot_longer(cols = c(2:76), 
                names_to = "year", 
                values_to = "value") 
 
 past_arc <- past_arc %>% 
   separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 
 
 past_arc$ID <- as.factor(past_arc$ID)
 
# past_arc_occ <- past_arc %>%
#   filter(year_n < 5000) %>%
#   mutate(occ = ifelse(value > 0, 1, 0)) %>%
#   group_by(ID) %>%
#   summarise(occ_tot = sum(occ))
# 
 past_occ_index <- merge(sites6kacoor, past_arc_occ, by= "ID")
 
 past_occ_index <- st_as_sf(past_occ_index, coords = c("long", "lat"), crs = 4326)
 st_write(past_occ_index, "past_occ_index6kav2.shp")


#########################################

popd_arc <- popd_arc %>%
  mutate(landcover = "converted_rangeland") 
colnames(popd_arc) <- new_names

popd_arc <- popd_arc %>% 
  dplyr::select(c(1,3:78)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

popd_arc <- popd_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 


popd_arc$year_n <- as.numeric(popd_arc$year_n)
popd_arc$ID <- as.factor(popd_arc$ID)

popd_arc <- popd_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

# 
# popd_arc_occ <- popd_arc %>%
#   filter(year_n < 5000) %>%
#   mutate(occ = ifelse(value > 0, 1, 0)) %>%
#   group_by(ID) %>%
#   summarise(occ_tot = sum(occ))

popd_arc_occ <- popd_arc %>%
  filter(value != 0) %>%
  filter(before_present < 6000) %>%
  arrange(ID, desc(before_present)) %>%
  group_by(ID)

popd_arc_occ <- popd_arc_occ %>%
  arrange(ID, before_present)

# Difference between adjacent rows per site (ID)
popd_arc_occ <- popd_arc_occ %>%
  group_by(ID) %>%
  mutate(time_diff = lead(before_present, default = last(before_present)) - before_present) %>%
  ungroup()

timetot <- popd_arc_occ %>%
  group_by(ID) %>%
  summarise(total_time_habitated = sum(time_diff, na.rm = TRUE))

# Print the result
print(timetot)


popd_occ_index <- merge(sites6kacoor, timetot, by= "ID")

popd_occ_index <- st_as_sf(popd_occ_index, coords = c("long", "lat"), crs = 4326)
st_write(popd_occ_index, "popd_occ_index6kav2.shp")


##############################################


obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)

# Iterate over the names of the columns from the 4-82 
for (i in 4:82) {
  # Replace "X" with an empty string in each column name 
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign new names of columns in the dataframe 
names(obs_arc) <- column_names


obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  dplyr::select(c(2,3, 21:80,82, 83)) %>%
  pivot_longer(cols = c(3:63),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)
obs_arc$nw_id <- as.factor(obs_arc$nw_id)

# obs_arc_occ <- obs_arc %>%
#   group_by(nw_id) %>%
#   summarise(occ_tot = sum(value)) %>%
#   rename(ID = nw_id)


obs_arc_occ <- obs_arc %>%
  filter(value != 0) %>%
  arrange(nw_id, desc(year)) %>%
  group_by(nw_id)

obs_arc_occ <- obs_arc_occ %>%
  arrange(nw_id, year)

# Difference between adjacent rows per site (ID)
obs_arc_occ <- obs_arc_occ %>%
  group_by(nw_id) %>%
  mutate(time_diff = lead(year, default = last(year)) - year) %>%
  ungroup() %>%
  filter(time_diff <= 200)

timetot <- obs_arc_occ %>%
  group_by(nw_id) %>%
  summarise(total_time_habitated = sum(time_diff, na.rm = TRUE)) %>%
  rename(ID = nw_id)

# Print result
print(timetot)

obs_occ_index <- merge(sites6kacoor, timetot, by= "ID")

obs_occ_index <- st_as_sf(obs_occ_index, coords = c("long", "lat"), crs = 4326)
#st_write(obs_occ_index, "obs_occ_index6kav2.shp")



#######################################################################################################
############################################## Data total by time for mapping Hyde and observations
##########################################################################################################


conv_arc <- "converted_rangeland_archaeological_points.xlsx"
graz_arc <- "grazing_archaeological_points.xlsx"
past_arc <- "pasture_archaeological_points.xlsx"
popd_arc <- "popdensity_archaeological_points.xlsx"

conv_arc <- read_xlsx(conv_arc)
graz_arc <- read_xlsx(graz_arc)
past_arc <- read_xlsx(past_arc)
popd_arc <- read_xlsx(popd_arc)

new_names <- c("ID" , "Name" , "10000BC", "9000BC",  "8000BC",  "7000BC",  "6000BC" , "5000BC" , "4000BC" , "3000BC" , "2017AD" , "2016AD" , "2015AD" , "2014AD", "2013AD" ,
               "2012AD" , "2011AD" , "2010AD" , "2009AD",  "2008AD" , "2007AD" , "2006AD" , "2005AD" , "2004AD" , "2003AD" , "2002AD" ,"2001AD" , "2000BC",  "2000AD",  "1990AD", 
               "1980AD" , "1970AD" , "1960AD" , "1950AD" , "1940AD" , "1930AD",  "1920AD",  "1910AD" , "1900AD" , "1890AD",  "1880AD" , "1870AD"  ,"1860AD" , "1850AD" , "1840AD" ,
               "1830AD",  "1820AD" , "1810AD" , "1800AD" , "1790AD" , "1780AD" , "1770AD"  ,"1760AD" , "1750AD" , "1740AD" , "1730AD" , "1720AD" , "1710AD" , "1700AD" , "1600AD" ,
               "1500AD",  "1400AD" , "1300AD" , "1200AD" , "1100AD" , "1000BC" , "1000AD",  "900AD"  , "800AD"  , "700AD" ,  "600AD" ,  "500AD" ,  "400AD" ,  "300AD" ,  "200AD"  ,
               "100AD" , "0AD")


colnames(conv_arc) <- new_names

conv_arc <- conv_arc %>% 
  select(c(1,3:77)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

conv_arc <- conv_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

conv_arc$year_n <- as.numeric(conv_arc$year_n)

conv_arc <- conv_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

conv_arc <- merge(sites6kacoor, conv_arc, by= "ID")

conv_arc_tot_time2 <- conv_arc %>%
  filter(before_present <= 6000) %>%
  mutate(occ = ifelse(value > 0, 1, 0)) %>%
  group_by(before_present) %>%  # ID
  summarise(conv_arc = sum(occ))


##########################

colnames(graz_arc) <- new_names

graz_arc <- graz_arc %>% 
  select(c(1,3:77)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

graz_arc <- graz_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

graz_arc$year_n <- as.numeric(graz_arc$year_n)

graz_arc <- graz_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

graz_arc <- merge(sites6kacoor, graz_arc, by= "ID")

graz_arc_tot_time <- graz_arc %>%
  filter(before_present <= 6000) %>%
  mutate(occ = ifelse(value > 0, 1, 0)) %>%
  group_by(before_present) %>%
  summarise(graz_arc = sum(occ))

psc_tot <- merge(conv_arc_tot_time2, graz_arc_tot_time, by = "before_present")

##################################


colnames(past_arc) <- new_names

past_arc <- past_arc %>% 
  select(c(1,3:77)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

past_arc <- past_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

past_arc$year_n <- as.numeric(past_arc$year_n)

past_arc <- past_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

past_arc <- merge(sites6kacoor, past_arc, by= "ID")

past_arc_tot_time <- past_arc %>%
  filter(before_present <= 6000) %>%
  mutate(occ = ifelse(value > 0, 1, 0)) %>%
  group_by(before_present) %>%
  summarise(past_arc = sum(occ))

psc_tot <- merge(psc_tot, past_arc_tot_time, by = "before_present")

#####################################

colnames(popd_arc) <- new_names

popd_arc <- popd_arc %>% 
  select(c(1,3:77)) %>%
  pivot_longer(cols = c(2:76), 
               names_to = "year", 
               values_to = "value") 

popd_arc <- popd_arc %>% 
  separate(year, into = c("year_n", "timescale"), sep = "(?<=\\d)(?=[A-Z])", convert = TRUE) 

popd_arc$year_n <- as.numeric(popd_arc$year_n)

popd_arc <- popd_arc %>% 
  mutate(before_present = ifelse(timescale == "BC", year_n + 1950, 1950 - year_n)) %>%
  filter(timescale %in% c("BC", "AD"))

popd_arc <- merge(sites6kacoor, popd_arc, by= "ID")

popd_arc_tot_time <- popd_arc %>%
  mutate(occ = ifelse(value > 0, 1, 0)) %>%
  group_by(before_present) %>%
  summarise(popd_arc = sum(occ))

psc_tot <- merge(psc_tot, popd_arc_tot_time, by = "before_present")


####################################
####################################

obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)
# # Iterate over the names of the columns from the 4-82
for (i in 4:82) {
  # Replace "X" with an empty string in each column_name
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign the new column names to the dataframe
names(obs_arc) <- column_names


obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  select(c(2,3, 21:80,82, 83)) %>%
  pivot_longer(cols = c(3:63),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)

obs_arc_tot <- obs_arc %>%
  group_by(year) %>%
  summarise(obs = sum(value)) %>%
  rename(before_present = year) ``

psc_tot_tst <- full_join(psc_tot, obs_arc_tot, by = "before_present")
#write.csv(psc_tot_tst, "count_hyde_obs6kav2.csv")

###############################################

psc_tot_tst_long <- psc_tot_tst %>%
  pivot_longer(cols = -c(before_present, past_arc) , names_to = "name", values_to = "value")

psc_tot_tst_long <- psc_tot_tst_long %>%
  filter(!is.na(value))  

psc_tot_tst_long$name <- factor(psc_tot_tst_long$name,
                                levels = c("obs", "conv_arc", "graz_arc", "popd_arc"),
                                labels = c("Observations", "Converted Rangeland", "Grazing", "Population Density"))

# Plot
ggplot(psc_tot_tst_long, aes(x = before_present, y = value, color = name)) +
  geom_line(size = 1.2, linetype = "solid") +
  geom_point(size = 2) +
  labs(x = "Years BP", y = "Number of archaeological sites inhabited", color = "Legend") +
  scale_x_continuous(trans = "reverse") +  # Reverse the x-axis
  scale_color_manual(values = c("Observations" = "#DF536B", "Converted Rangeland" = "#61D04F", "Grazing" = "#2297E6",  "Population Density" = "gray62")) +
  facet_wrap(~ name, ncol = 1, scales = "free_y")



#ggsave("simultaneous_inhab_1.png", width=12, height=18,  dpi = 300) 



###########################################################

# Exclusion data from 6.0 ka BP. 


obs_arc <- "arch_obs_liter_new.xlsx"

obs_arc <- read_xlsx(obs_arc)

column_names <- names(obs_arc)
# Iterate over the names of the columns from the 4-82
for (i in 4:82) {
  # Replace "X" with an empty string in each column_name
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign new names of columns of the dataframe
names(obs_arc) <- column_names

obs_arc1 <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  select(c(2,3, 21:80, 82, 83)) %>%
  pivot_longer(cols = c(3:63),
               names_to = "year",
               values_to = "value") %>%
  group_by(nw_id, site2) %>%
  summarize(cantidad = sum(value ==1, na.rm= TRUE)) %>%
  rename(ID = nw_id) 

st_afte_6 <- filter(obs_arc1, cantidad >0)
write.csv(st_afte_6, "sites_after6.csv")

sitios_con_cero <- sum(obs_arc1$cantidad == 0)


obs_dsd_6 <- merge(df_arc_sites, obs_arc1, by= "ID")

obs_dsd_6 <- st_as_sf(obs_dsd_6, coords = c("long", "lat"), crs = 4326)
#st_write(obs_dsd_6, "obs_dsd_6v2.shp")



########################################################
# Sites by time interval
########################################################


####### 6000 to 4000
##################################

obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)

# Iterate over the names of the columns from the 4-82
for (i in 4:82) {
  # Replace "X" with an empty string in each column_name
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign new names of columns of the dataframe
names(obs_arc) <- column_names

obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  dplyr::select(c(2,3, 21:40, 83)) %>%
  pivot_longer(cols = c(3:22),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)

obs_arc_ant <- obs_arc %>%
  #filter(year <= 10000) %>%
  filter(value != 0) %>%
  group_by(nw_id) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  arrange(desc(year)) %>%
  rename(old_id = ID) %>%
  rename(ID = nw_id)

obs_antq_index6_4 <- merge(sites6kacoor, obs_arc_ant, by= "ID") %>%
  dplyr::select("ID","old_ID", "name" , "presencia" ,"cantidad" , "long" , "lat" ,"year"  ,"value") %>%
  rename(bfr_prs = year)

obs_antq_index6_4 <- st_as_sf(obs_antq_index6_4, coords = c("long", "lat"), crs = 4326)
#st_write(obs_antq_index6_4, "obs_antq_index6_4Calib.shp")


####### 4000 to 2000
##################################

obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)

# Iterate over the names of the columns from the 4-82
for (i in 4:82) {
  # Replace "X" with an empty string in each column_name
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign new names of columns of the dataframe
names(obs_arc) <- column_names

obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  dplyr::select(c(2,3, 41:60, 83)) %>%
  pivot_longer(cols = c(3:22),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)

obs_arc_ant <- obs_arc %>%
  #filter(year <= 10000) %>%
  filter(value != 0) %>%
  group_by(nw_id) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  arrange(desc(year)) %>%
  rename(old_id = ID) %>%
  rename(ID = nw_id)

obs_antq_index4_2 <- merge(sites6kacoor, obs_arc_ant, by= "ID") %>%
  dplyr::select("ID","old_ID", "name" , "presencia" ,"cantidad" , "long" , "lat" ,"year"  ,"value") %>%
  rename(bfr_prs = year)

obs_antq_index4_2 <- st_as_sf(obs_antq_index4_2, coords = c("long", "lat"), crs = 4326)
#st_write(obs_antq_index4_2, "obs_antq_index4_2Calib.shp")


####### 2000 to 0
##################################

obs_arc <- "matrix_arch_obs_by_year_CALIBRATED.csv"

obs_arc <- read.csv(obs_arc)

column_names <- names(obs_arc)

# Iterate over the names of the columns from the 4-82
for (i in 4:82) {
  # Replace "X" with an empty string in each column_name
  column_names[i] <- gsub("X", "", column_names[i])
}

# Assign new names of columns of the dataframe
names(obs_arc) <- column_names

obs_arc <- obs_arc %>%
  filter(Mapped == "yes") %>%
  mutate(nw_id = 1:73) %>% 
  dplyr::select(c(2,3, 61:82, 83)) %>%
  pivot_longer(cols = c(3:22),
               names_to = "year",
               values_to = "value")

obs_arc$year <- as.numeric(obs_arc$year)

obs_arc_ant <- obs_arc %>%
  #filter(year <= 10000) %>%
  filter(value != 0) %>%
  group_by(nw_id) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  arrange(desc(year)) %>%
  rename(old_id = ID) %>%
  rename(ID = nw_id)

obs_antq_index2_0 <- merge(sites6kacoor, obs_arc_ant, by= "ID") %>%
  dplyr::select("ID","old_ID", "name" , "presencia" ,"cantidad" , "long" , "lat" ,"year"  ,"value") %>%
  rename(bfr_prs = year)

obs_antq_index2_0 <- st_as_sf(obs_antq_index2_0, coords = c("long", "lat"), crs = 4326)
#st_write(obs_antq_index2_0, "obs_antq_index2_0Calib.shp")

















