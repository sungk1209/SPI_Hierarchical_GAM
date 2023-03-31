# *------------------------------------------------------------------
# | PROGRAM NAME: Analyze a multi-centennial meteorological drought trend in North America
# | FILE NAME: 2.download_inst.R
# | DATE: 
# | CREATED BY:  Kay Sung       
# *----------------------------------------------------------------
# | PURPOSE:  download Gridmet and CRU 
# |          3-months precipitation
# |        
# |
# *------------------------------------------------------------------
require(tidyverse)
require(here)
require(lubridate)
require(zoo)
require(ggplot2)
require(ncdf4)
require(dplyr)
require(fitdistrplus)
require(RANN)
require(mgcv)
require(furrr)
require(parallel)
select <- dplyr::select

data_path <- "../data/"
output_path <- "../output/"

output.p <- file.path(output_path, "test")
dir.create(output.p, recursive = FALSE)

write_data_path <- file.path(output.p, "data")
dir.create(write_data_path, recursive=FALSE, showWarnings = FALSE)

naspa_grid<- readRDS(paste0(write_data_path,"/naspa_grid.rds"))

download_inst <- function(j){
  
  loc <- data.frame(index = j,lat = naspa_grid[j,]$lat, lon = naspa_grid[j,]$lon)
 print(loc)
 
#######################################################################
##		 2.1 download GRIDMET
###################################################################

lat_col <- which.min(abs(lat_list_gridmet - loc$lat))
lon_col <- which.min(abs(lon_list_gridmet - loc$lon))

#count <- lat_col[length(lat_col)]-lat_col[1] 

var_data_j <- ncvar_get(nc_file_gridmet, as.character(var_prcp_gridmet), 
                        start=c(lon_col[1], lat_col[1], 1), 
                        count=c(1,1,-1)) 

### Remove missing values
var_data_j[var_data_j ==32767] <- NA

####Gridmet has fine resolution, need to average over all grid cells
#Unit of gridcell: daily accumulated 
#a <- colMeans(var_data_j, dim = 2)

yup <- tibble(site =  loc$index, date = date_list_gridmet, 
              variable = as.character(var_prcp_gridmet), value = var_data_j)

accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  group_by(year,month)%>%
  summarise(pr_mm_month = sum(value)) %>%
  ungroup()

n_roll <-3
n_roll_min <-2
gridmet_df <- accum_df %>%
  mutate(pr_3M_ave = rollmeanr(x=pr_mm_month, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(units = "mm/month") %>%
  mutate(date = as.Date(paste0(year,'-',month,'-01'))) %>%
  mutate(site = loc$index) %>%
  mutate(model = "Gridmet") %>%
  mutate(variable = as.character(ncdf_df_gridmet$short_name)) %>%
  select(date, site, year, variable, model, pr_3M_ave,units)

### Close the nc file
#nc_close(nc_file)

#######################################################################
## 2.2 Download CRU 
#####################################################################

#Example names as the row names, column 1 as longitude, and column 2 as latitude
lat_col <- which.min(abs(lat_list_cru - loc$lat))
lon_col <- which.min(abs(lon_list_cru - loc$lon))

var_data_j <- ncvar_get(nc_cru, varid= var_prcp_cru, start = c(lon_col,lat_col,1), 
                        count = c(1,1,-1))

var_data_j[var_data_j < -9.96920996838687e+36] <- NA

yup <- tibble(site = loc$index, date = date_list_cru, 
              variable = var_prcp_cru, value = var_data_j)

cru_df <- yup %>%
  arrange(date) 

n_roll <- 3
n_roll_min <- 2
cru_df <- cru_df %>%
  mutate( pr_3M_ave = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE))  %>%
  #  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  #  mutate(precip = case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
  #                              TRUE ~ NA_real_)) %>%
  mutate(year = year(date)) %>%
  mutate(units = "mm/month") %>% 
  mutate(model = "CRU") %>%
  select(date, site, year, variable,model, pr_3M_ave,units)
#nc_close(nc_cru)


################################################################
##########Combine instrumental data 
###############################################################

instrument_df <- rbind(gridmet_df,cru_df)

instrument_df <- instrument_df %>%
  arrange(date) %>%
  drop_na(pr_3M_ave) %>%
  mutate(scenario = "observed") %>%
  mutate(emission = "historical")

saveRDS(instrument_df,file = paste0(write_data_path,"/",loc$index,"instrument_df.rds"))

}

##########################################################3
###Variables for Gridmet 
###########################################################
ncdf_df_gridmet <- data.frame(short_name = c("pr"))

ncdf_df_gridmet <- ncdf_df_gridmet %>%
  mutate(url = paste0('http://thredds-prod.nkn.uidaho.edu:8080/thredds/dodsC/agg_met_',short_name, 
                      '_1979_CurrentYear_CONUS.nc#fillmismatch'))


### Open the NCDF file
nc_file_gridmet <- nc_open(ncdf_df_gridmet$url, verbose = FALSE)
var_prcp_gridmet <- attributes(nc_file_gridmet$var)$names[4]

### Extract lat, lon, and time info
lat_list_gridmet <- ncvar_get(nc_file_gridmet, "lat")
lon_list_gridmet <- ncvar_get(nc_file_gridmet, "lon")
date_list_gridmet <- ncvar_get(nc_file_gridmet, "day") + as.Date("1900-01-01")

#######################################################################
## 2.2 Variables for CRU 
#####################################################################

nc_cru <- nc_open(filename = "../data/cru_ts4.05.1901.2020.pre.dat.nc")
attributes(nc_cru$var)$names
var_prcp_cru <- attributes(nc_cru$var)$names[[1]]
var_unit <- ncatt_get(nc_cru, attributes(nc_cru$var)$names[[1]], "units")

lat_list_cru <- ncvar_get(nc_cru,"lat")
lon_list_cru <- ncvar_get(nc_cru,"lon")

time_scale <- ncatt_get(nc_cru, "time", "units")$value %>%
  str_split_fixed("since",n =2)
date_list_cru <- ncvar_get(nc_cru, "time") + as.Date(time_scale[2])

results<- map(c(1:100), download_inst)
                   