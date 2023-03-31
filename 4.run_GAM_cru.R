# *------------------------------------------------------------------
# | PROGRAM NAME: Analyze a multicentennial meteorological drought trend in AMerica
# | FILE NAME: 4.run_GAM_cru.R
# | DATE: 
# | CREATED BY:  Kay Sung       
# *----------------------------------------------------------------
# | PURPOSE: Generalize Additive Model to analyze trends of SPI-3  
# |          in NASPA, MIROC, MRI, CRU, Gridmet(only available in US)
# |          This should work in pararell. 
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

data.p <- file.path(output_path, "test1")
dir.create(data.p, recursive = FALSE)

write_data_path <- file.path(data.p, "data")
dir.create(write_data_path, recursive=FALSE, showWarnings = FALSE)

naspa_grid<- readRDS(paste0(data.p,"/naspa_grid.rds"))
output.p <- file.path(output_path, "Gam_result")
###########################################################################
### 1. Download CMIP6 data
#########################################################################
GI_model_grid <- function(j){
  
  loc <- data.frame(index = j,lat = naspa_grid[j,]$lat, lon = naspa_grid[j,]$lon)
 print(loc)
 
##############################################################
####          Call precipitation data                       #####
##############################################################
 instrument_df <- readRDS(paste0(write_data_path,"/",j,"instrument_df.rds"))
 ds_naspa_df <- readRDS(paste0(write_data_path,"/",j,"ds_naspa.rds"))
 pred_df<- readRDS(paste0(write_data_path,j,"_pr_cmip6.rds") )
 
 instrument_df <- instrument_df %>%
   mutate(month = month(date))
 
 ds_naspa_df <- ds_naspa_df %>%
   mutate(month = month(date))
 
 prcp_df <- instrument_df %>%
   bind_rows(ds_naspa_df) %>%
   bind_rows(pred_df) 
 #I downloaded GFDL, but did not used for Gam modeling because it didn't have PMIP
 prcp_df <- prcp_df %>%
   filter(model !="GFDL-ESM4") %>%
   drop_na() %>%
   arrange(date) 
 
 prcp_pos <- prcp_df %>%
   filter(pr_3M_ave >0)
prcp_pos$model <- as.factor(prcp_pos$model) 
#p <- ggplot(prcp_pos %>% filter(month == 12)) +
#   geom_line(aes(x = date, y = pr_3M_ave, color = model, linetype = scenario)) + theme_classic()
#p
# temp <- prcp_pos %>% filter(model == "CRU" | model == "Naspa")
# p <- ggplot(data=temp %>%filter(month == 1 & year %in% seq(1900,2020))) +
#   geom_line(aes(x = date, y = pr_3M_ave, group= model,colour = model)) +
#   geom_vline(xintercept = as.Date("2100-01-31")) +
#   labs(title = "",
#        x="year",y="precip(mm/Month)", size = 20) +
#   scale_x_date() +
#   theme_classic(base_size = 20)

# p

 scenarios <- c("ssp126","ssp585")
for (i in c(1,2)) {
  prcp_tmp <- prcp_pos %>%
    filter(emission == scenarios[i] | emission == "historical") %>%
    filter(year > 800 & year < 2100) %>%
    arrange(date)
  
   print(paste0("start GAM model_",i) )
  #Original Model with weights
  start_time <- Sys.time()
  gam_fit <- gam(list(
    pr_3M_ave ~ model + s(month, by = model, bs = "cc", k = 6) + 
      te(year, month, bs = c("cr","cc"), k = c(18,6)) , 
    ~model + s(month, by = model, bs = "cc", k = 6) + 
      te(year, month, bs = c("cr","cc"), k = c(18,6))),
    family=gammals, link=list("identity","log"), 
    data = prcp_tmp)
  end_time <- Sys.time()
  dur <- end_time - start_time
  print(dur)
  
  saveRDS(gam_fit, file = paste0(output.p,"/",loc$index,"gam_fit_",i,".rds"))
  
  ### make a new combined_df with Gridmet as model
  min_cru <- min((prcp_tmp %>% filter(model == "CRU"))$year) 
  max_cru <- max((prcp_tmp %>% filter(model == "CRU"))$year)
  min_naspa <- min((prcp_tmp %>% filter(model == "Naspa"))$year) 
  max_naspa <- max((prcp_tmp %>% filter(model == "Naspa"))$year)
  min_miroc <- min((prcp_tmp %>% filter(model == "Miroc"))$year) 
  max_miroc <- max((prcp_tmp %>% filter(model == "Miroc"))$year)
  min_mri <- min((prcp_tmp %>% filter(model == "MRI"))$year) 
  max_mri <- max((prcp_tmp %>% filter(model == "MRI"))$year)
  
  if(sum(is.na(prcp_tmp$model == "Gridmet")) ==0){
 
    new_ts_df <-data.frame(expand.grid(year = seq(min_naspa,max_naspa),
                                     month = seq(1,12)), model = "CRU", plot_model = "Naspa") %>% 
    bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),
                                     month = seq(1,12)), model = "CRU", plot_model = "CRU")) %>% 
    bind_rows(data.frame(expand.grid(year = seq(min_miroc,max_miroc),
                                     month = seq(1,12)), model = "CRU", plot_model = "Miroc")) %>% 
    bind_rows(data.frame(expand.grid(year = seq(min_mri,max_mri),
                                     month = seq(1,12)), model = "CRU", plot_model = "MRI")) 
  }else {
    min_gridm <- min((prcp_tmp %>% filter(model == "Gridmet"))$year) 
    max_gridm <- max((prcp_tmp %>% filter(model == "Gridmet"))$year)
    
    new_ts_df <-data.frame(expand.grid(year = seq(min_naspa,max_naspa),
                                       month = seq(1,12)), model = "CRU", plot_model = "Naspa") %>% 
      bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),
                                       month = seq(1,12)),model = "CRU", plot_model = "Gridmet")) %>% 
      bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),
                                       month = seq(1,12)), model = "CRU", plot_model = "CRU")) %>% 
      bind_rows(data.frame(expand.grid(year = seq(min_miroc,max_miroc),
                                       month = seq(1,12)), model = "CRU", plot_model = "Miroc")) %>% 
      bind_rows(data.frame(expand.grid(year = seq(min_mri,2100),
                                       month = seq(1,12)), model = "CRU", plot_model = "MRI")) 
    
  }
  ### Make predictions based on this
  GAM_predict <- predict(gam_fit, newdata = new_ts_df, se.fit = TRUE, type = "response")
  
  GAM_predict  <- GAM_predict %>%
    data.frame() %>%
    as_tibble() %>%
    rename(est_mean = 1) %>%
    rename(est_shape = 2) 
  
  GAM_predict <- GAM_predict %>%
    mutate(est_shape = 1/exp(est_shape)) %>%
    mutate(est_scale = est_mean/est_shape) %>%
    mutate(est_sd = (est_shape * est_scale^2)^0.5)
  
  temp_model_ts <- transform(new_ts_df,
                             modGI = GAM_predict$est_mean,
                             shapeGI = GAM_predict$est_shape,
                             scaleGI = GAM_predict$est_scale,
                             lower2 = qgamma(0.025, shape = GAM_predict$est_shape,
                                            scale = GAM_predict$est_scale),
                             upper2 = qgamma(0.975, shape = GAM_predict$est_shape,
                                            scale = GAM_predict$est_scale),
                             #lowcrit = GAM_predict$est_mean + crit * GAM_predict$se.fit.1,
                             #upcrit = GAM_predict$est_mean - crit * GAM_predict$se.fit.1,
                             #lowp = GAM_predict$est_mean -  2* GAM_predict$se.fit.1,
                             #upp = GAM_predict$est_mean +  2* GAM_predict$se.fit.1,
                             est_sd = GAM_predict$est_sd,
                             scenario = scenarios[i])
  
  temp_model_ts <- temp_model_ts %>%
    mutate(date = as.Date(paste0(year,"-",month,"-01")))
  
  if (i == 1) {modeled_df <- temp_model_ts
  }else{modeled_df <- modeled_df %>% bind_rows(temp_model_ts)}
  
}
modeled_df <- modeled_df %>%
   arrange(date)
saveRDS(modeled_df, file = paste0(output.p,"/",loc$index,"_modeled_ts.rds"))

return(temp_model_ts$lower15)
}

numCores<- parallel::detectCores() - 1  
  
results<- mclapply(naspa_grid, GI_model_grid,numCores)


                   
