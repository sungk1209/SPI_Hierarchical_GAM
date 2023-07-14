# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: 3.gpcc_convert_spi3.R
# | DATE: Oct.02.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: -
# | 
# | 
# *--------------------------------------------------------------------------
require(lubridate)
require(tidyverse)
require(dplyr)
require(ncdf4)
require(fitdistrplus)
require(zoo)
library(RANN)
select <- dplyr::select

data_path <- "../data/"
output_path <- "../output/"

naspa_download <- function(loc){
  
fn_wm <- "NASPA_WARM_SPI.nc"
short_name <- "pr"

ncdf_df <- data.frame(short_name = c("pr"), 
                      var_name = c("precipitation_amount"))
nc_wm <- nc_open(paste0(data_path,fn_wm))
#print(nc_wm)

var1 <- attributes(nc_wm$var)

lat_list <- ncvar_get(nc_wm, "lat")
lon_list <- ncvar_get(nc_wm, "lon")
date_list <- ncvar_get(nc_wm, "time")
tm_orig <- as.Date("0000-08-01")
date_list <- tm_orig %m+% years(date_list)

lat_col <- which.min(abs(lat_list - loc$lat))
lon_col <- which.min(abs(lon_list - loc$lon))

var_data_k <- ncvar_get(nc_wm, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                        count = c(-1,1,1))

yuc <- tibble(date = date_list,spi3 = var_data_k) %>%
  drop_na()

begin.y <- min(year(yuc$date))
end.y <- max(year(yuc$date)) 
yrs <- list(begin.y, end.y)

### Add a month column
naspa_spi3 <- yuc %>% 
  mutate(year = year(date)) %>%
  mutate(month = 7) %>%
  complete(year = seq(begin.y,end.y), month = seq(1,12)) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-01"))) %>%
  mutate(date = as.Date(ceiling_date(date, "month")-1)) %>%
  arrange(date)

fn_co <- "NASPA_COOL_SPI.nc"

nc_co <- nc_open(paste0(data_path,fn_co))
print(nc_co)

var1 <- attributes(nc_co$var)
lat_list <- ncvar_get(nc_co, "lat")
lon_list <- ncvar_get(nc_co, "lon")
date_list <- ncvar_get(nc_co, "time")
tm_orig <- as.Date("0000-04-30")
date_list <- tm_orig %m+% years(date_list)

lat_col <- which.min(abs(lat_list - loc$lat))
lon_col <- which.min(abs(lon_list - loc$lon))

var_data_k <- ncvar_get(nc_co, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                        count = c(-1,1,1))

yuc <- tibble(date = date_list, 
              spi5 = var_data_k) %>%
  drop_na()

naspa_spi5 <- yuc %>% 	
  mutate(month = 4) %>%
  mutate(year = year(date)) %>%
  complete(year = seq(begin.y,end.y), month = seq(1,12)) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-01"))) %>%
  mutate(date = as.Date(ceiling_date(date, "month")-1)) %>%
  arrange(date)

naspa_spi3 <- naspa_spi3 %>%
  mutate(spi5 = naspa_spi5$spi5)

saveRDS(naspa_spi3,file = paste0(write_data_path,"/naspa_spi.rds"))

#############################################################
### GPCC
##################################################################

n_days <- 12
n_roll <- 3

gpcc <- "../data/precip.mon.total.v7.nc"

nc_gpcc <- nc_open(gpcc,verbose = FALSE)

var1 <- attributes(nc_gpcc$var)
lat_list <- ncvar_get(nc_gpcc, "lat")
lon_list <- ncvar_get(nc_gpcc, "lon")
date_list <- ncvar_get(nc_gpcc, "time")
tm_orig <- as.Date("1800-01-01")
date_list <- date_list + tm_orig

lat_col <- which.min(abs(lat_list - loc$lat))
if(loc$lon < 0) {
  lon_col <- which.min(abs(lon_list - (360 + loc$lon)))
}else{
  lon_col <- which.min(abs(lon_list - loc$lon))}


var_data_j <- ncvar_get(nc_gpcc, varid= var1$names[2], start = c(lon_col,lat_col,1), 
                        count = c(1,1,-1))
yup <- tibble(site = loc$index, date = date_list, 
              value = as.numeric(var_data_j))

gpcc_df <- yup %>%
  mutate(model = "GPCC") %>%
  mutate(units = "mm/month") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  select(date,site, year,month, model, value,units)

n_roll <- 3
n_roll_min <- 2

gpcc_df <- gpcc_df %>%
  mutate(roll_mean_3 = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(pr_3M_ave = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                            TRUE ~ NA_real_)
  ) %>%
  mutate(date = ceiling_date(date, "month") - 1) %>%
  select(date, site, model,pr_3M_ave,units, value) 

n_roll <- 5
n_roll_min <- 4

gpcc_df <- gpcc_df %>%
  mutate(roll_mean_5 = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(pr_5M_ave = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_5,
                              TRUE ~ NA_real_)
  ) %>%
  select(-roll_mean_3_notna, -roll_mean_5,-value) %>%
  filter(pr_3M_ave >0 & pr_5M_ave >0)

gpcc_df <- gpcc_df %>%
  group_by(month(date)) %>%
  mutate(shape3 = fitdist(pr_3M_ave, "gamma")$estimate[[1]],
         rate3 = fitdist(pr_3M_ave, "gamma")$estimate[[2]]) %>%
  mutate(prob3 = pgamma(pr_3M_ave, shape = shape3, rate = rate3)) %>%
  mutate(spi3 = qnorm(prob3,mean = 0, sd = 1)) %>%
  mutate(shape5 = fitdist(pr_5M_ave, "gamma")$estimate[[1]],
         rate5 = fitdist(pr_5M_ave, "gamma")$estimate[[2]]) %>%
  mutate(prob5 = pgamma(pr_5M_ave, shape = shape5, rate = rate5)) %>%
  mutate(spi5 = qnorm(prob5,mean = 0, sd = 1)) %>%
  ungroup()

saveRDS(gpcc_df,file = paste0(write_data_path,"/gpcc_df.rds"))
return(yrs)
}
