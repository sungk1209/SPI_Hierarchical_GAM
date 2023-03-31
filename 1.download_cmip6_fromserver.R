# *------------------------------------------------------------------
# | PROGRAM NAME: Analyze a multi centennial meteorological drought trend in North America
# | FILE NAME: 1.download_cmip6_fromserver.R
# | DATE: 
# | CREATED BY:  Kay Sung       
# *----------------------------------------------------------------
# | PURPOSE: Download MIROC, MRI of Pmip, historical and CMIP6 scenario.
# |           (I also added gfdl just in case we need in the future)
# |          and construct 3 months average precipitation              
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

select <- dplyr::select
data_path <- "./data/"
output_path <- "../output/"

output.p <- file.path(output_path, "test")
dir.create(output.p, recursive = FALSE)

write_data_path <- file.path(output.p, "data")
dir.create(write_data_path, recursive=FALSE, showWarnings = FALSE)

###########################################################################
###  Prepare a list of PMIP files
###########################################################################
download_cmip <- function(j){
  
  loc <- data.frame(index = j,lat = naspa_grid[j,]$lat, lon = naspa_grid[j,]$lon)
  print(loc)
  
  ###########################################################################
  ###  Function Download CMIP6 data
  #########################################################################
  
  for (i in seq(1,dim(ncdf_df)[1])){
    
    ### Choose a single NCDF file
    ncdf_info <- ncdf_df[i,]
    
    ### Open the NCDF file
    nc_file <- nc_open(ncdf_info$url,verbose = FALSE)
    #45,94
    ### Extract lat, lon, and time info
    lat_list <- ncvar_get(nc_file, "lat")
    lon_list <- ncvar_get(nc_file, "lon")
    #lon_list <- ifelse(lon_list > 180, -(360 - lon_list), lon_list)
    
    begin_date <- ncatt_get(nc_file, "time", attname = "units")$value
    begin_date <- substring(begin_date,11,nchar(begin_date))
    
    date_list <- ncvar_get(nc_file, "time") + as.Date(begin_date)
    
    lat_col <- which.min(abs(lat_list - loc$lat))
    if(loc$lon < 0) {
      lon_col <- which.min(abs(lon_list - (360 + loc$lon)))
    }else{
      lon_col <- which.min(abs(lon_list - loc$lon))}
    
    ### Extract variable data
    var_data_j <-  ncvar_get(nc_file, short_name, 
                             start=c(lon_col,lat_col,1), 
                             count=c(1,1,-1)) 
    ### Remove missing values
    var_data_j[var_data_j > 1.00e+20] <- NA
    
    ### Extract variable
    yup <- tibble(site = loc$index, date = date_list, 
                  variable = as.character(short_name), 
                  value = var_data_j, model = ncdf_info$model, 
                  scenario = ncdf_info$scenario,
                  emission = ncdf_info$emission)
    
    ### Close the nc file
    nc_close(nc_file)
    
    
    if (i == 1) {
      clim_df <- yup
    } else {
      clim_df <- rbind(clim_df,yup)
    }
    
    print(i)
  }
  
  ###########################################################################
  ###  Do some processing
  ###########################################################################
  ### Convert from kg/m2s to mm/month
  
  clim_df <- clim_df %>%
    mutate(month = month(date)) %>%
    #mutate(precip_mm_day = value * 24*60*60)
    mutate(precip_mm_month = ifelse(month %in% c(1,3,5,7,8,10,12), value*24*60*60*31,
                                    ifelse(month == 2, value*24*60*60*28,
                                           value*24*60*60*30)))
  
  clim_df <- clim_df %>%
    arrange(date)
  
  # p <- ggplot(clim_df, aes(x=date, y= precip_mm_month, colour= interaction(scenario,model))) + geom_line(size= 1,alpha = 0.5) +
  #   theme_bw(base_size = 15) +
  #   scale_y_continuous(name="month_ave_Prcp (mm/day)") +
  #   scale_x_date(breaks = "100 years",date_labels = "%Y", name = "year")
  # p
  #############################################################
  ######calculate 3 months average precipitation ##############
  #############################################################
  #make monthly sum
  n_roll <-3
  n_roll_min <-2
  cmip_df <- clim_df %>%
    group_by(model,scenario) %>%
    mutate(pr_3M_ave = rollmeanr(x=precip_mm_month, k=n_roll, fill=NA, na.rm=TRUE)) %>%
    #mutate(roll_mean_3_notna = rollsumr(x=!is.na(precip_mm_month), k=n_roll, fill=NA, na.rm=TRUE)) %>%
    # mutate(pr_3m_month = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
    #                          TRUE ~ NA_real_)
    #) %>%
    mutate(units = "mm/month") %>%
    mutate(site = loc$index) %>%
    mutate(year = year(date), month = month(date)) %>%
    ungroup() %>%
    select(date, year,site, variable, model, pr_3M_ave,units, scenario,month, emission)
  
  # p <- ggplot(cmip_df %>% filter(month(date) == 7), 
  #             aes(x=date, y=pr_3M_ave, color = interaction(scenario,model))) +
  #   geom_line(alpha = 0.5) + 
  #   theme_bw() + scale_y_continuous(name="3M-ave.Precip (mm/day)") +
  #   scale_x_date(breaks = "100 years",date_labels = "%Y", name = "year") 
  # p
  
  saveRDS(cmip_df, file = file.path(paste0(write_data_path,"/",j,"_pr_cmip6.rds") ))
  
  print(paste0("downloaded_",j))
  
  
}


#################################################################
###Download List
#################################################

#setup global variables
short_name <- "pr"

### Create an object to hold the data links
mri_past <- tibble(begin_date = c(NA,NA),
                   end_date = c(NA,NA), 
                   scenario = c("PMIP","historical"),
                   emission = "historical",
                   model="MRI")
mri_past <- mri_past%>%
  mutate(url = c(
    paste0(data_path,"/cmip/","pr_Amon_MRI-ESM2-0_past1000_r1i1p1f1_gn_085001-184912.nc"),
    paste0(data_path,"/cmip/","pr_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc")
  ))

mri_ssp <- tibble(begin_date = c("201501","201501"),
                  end_date = c("210012","210012"), 
                  scenario = c("ssp126","ssp585"), 
                  emission = c("ssp126","ssp585"), 
                  #num_sce = c(126,585),
                  model="MRI")
mri_ssp <- mri_ssp %>% 
  mutate(url = 
           paste0(data_path,"/cmip/","pr_Amon_MRI-ESM2-0_",scenario,"_r1i1p1f1_gn_",begin_date,"-",end_date,".nc")
  ) #%>%
# select(-num_sce)

ncdf_mri <- mri_past %>%
  bind_rows(mri_ssp)
#link for miroc
mrc_past <- tibble(begin_date = paste0(sprintf('%0.4d',seq(0850, 1650, 200)),"01"), 
                   end_date = paste0(seq(1049,1849, 200),"12"), 
                   scenario = "PMIP", 
                   emission = "historical",
                   model="Miroc")
### Add URL
mrc_past <- mrc_past %>% 
  mutate(url = 
           paste0(data_path,"/cmip/","pr_Amon_MIROC-ES2L_past1000_r1i1p1f2_gn_",begin_date,"-",end_date,".nc")
  )

#download historical data for MIROC

mrc_hist <- tibble(begin_date = NA, 
                   end_date =   NA, 
                   scenario = "historical", 
                   emission = "historical",
                   model="Miroc")
### Add URL
mrc_hist <- mrc_hist  %>% 
  mutate(url = paste0(data_path,"/cmip/","pr_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc"))
mrc_hist

mrc_ssp <- tibble(begin_date = "201501", 
                  end_date = "210012", 
                  scenario = c("ssp126", "ssp585"), 
                  emission = c("ssp126", "ssp585"), 
                  model="Miroc")

mrc_ssp <- mrc_ssp %>%
  mutate(url = c(paste0(data_path,"/cmip/","pr_Amon_MIROC6_ssp126_r1i1p1f1_gn_201501-210012.nc"), 
                 paste0(data_path,"/cmip/","pr_Amon_MIROC6_ssp585_r1i1p1f1_gn_201501-210012.nc")))

ncdf_miroc <- mrc_past %>%
  bind_rows(mrc_hist) %>%
  bind_rows(mrc_ssp)

# gfdl_historical <- tibble(begin_date = c("185001", "195001"), 
#                           end_date = c("194912", "201412"), 
#                           scenario = "historical", 
#                           emission = "historical",
#                           model="GFDL-ESM4")
# ### Add URL
# gfdl_historical <- gfdl_historical %>% 
#   mutate(url = paste0('http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Amon/pr/gr1/v20190726/pr_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_',begin_date, '-', end_date,'.nc'))
# 
# gfdl_historical
# 
# ### Create an object to hold the data links
# gfdl_585 <- tibble(   begin_date = NA, 
#                       end_date = NA, 
#                       scenario = "ssp585", 
#                       emission = "ssp585",
#                       model="GFDL-ESM4")
# ### Add URL
# gfdl_585 <- gfdl_585 %>% 
#   mutate(url = paste0("http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp585_r1i1p1f1_gr1_201501-210012.nc"))
# 
# ncdf_585

# gfdl_126 <- tibble(  begin_date = NA, 
#                      end_date = NA, 
#                      scenario = "ssp126", 
#                      emission = "ssp126", 
#                      model="GFDL-ESM4")
### Add URL
# gfdl_126 <- gfdl_126 %>% 
#   mutate(url = paste0('http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp126/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp126_r1i1p1f1_gr1_201501-210012.nc'))
# 
# gfdl_126
# 
# ### Merge
# gfdl_df <- gfdl_historical %>%
#   bind_rows(gfdl_585)  %>%
#   bind_rows(gfdl_126)
### Check
ncdf_df <- ncdf_miroc %>%
  bind_rows(ncdf_mri) 

### Check
ncdf_df

###################################################################
####   Determine the spatial coverage based on Naspa
#####################################################################

fn_wm <- "NASPA_WARM_SPI.nc"

nc_wm <- nc_open(paste0(data_path,fn_wm))
#print(nc_wm)

var1 <- attributes(nc_wm$var)
lat_list_naspaw <- ncvar_get(nc_wm, "lat")
lon_list_naspaw <- ncvar_get(nc_wm, "lon")

#naspa_grid <- expand_grid(lat = lat_list, lon = lon_list)
var_data_k <- ncvar_get(nc_wm, varid= var1$names[1], start = c(1400,1,1), 
                        count = c(1,-1,-1))

#####make grid covers of NASPA
# I used this only one time for determine the cells where have NASPA.
#Made file naspa_grid.rds and ran GAM model for those cells.

grids <- data.frame(which(!is.na(var_data_k), arr.ind = TRUE))
naspa_grid <- data.frame(lat = lat_list_naspaw[grids$row], 
                         lon = lon_list_naspaw[grids$col])

#I saved this in "naspa_grid.rds"

#############################################################
#### Run model
##################################################
#numCores<- parallel::detectCores() - 1

# warning: it is hard to run in parallel, as you keep accessing the webpage at the same time.
# Downloading the whole file instead of accessing server extracting the data is a lot faster.


results<- purrr::map(c(4509:6638), download_cmip)
# 
# temp <- cmip_df %>% filter(scenario == "historical")
# 
# p <- ggplot(temp %>% filter(month(date) == 7), 
#             aes(x=date, y=pr_3M_ave, color = interaction(scenario,model))) +
#      geom_line(alpha = 0.5) + 
#     theme_bw() + scale_y_continuous(name="3M-ave.Precip (mm/day)") +
#       scale_x_date(breaks = "100 years",date_labels = "%Y", name = "year") 
#   p
#                    