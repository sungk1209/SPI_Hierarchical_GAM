# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: 4. naspa_downscale_knn.R
# | DATE: Oct.02.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
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

### Read in
knn <- function(loc){
  
gpcc_df <- readRDS(file = paste0(write_data_path,"/gpcc_df.rds"))
naspa_spi3<- readRDS(file = paste0(write_data_path,"/naspa_spi.rds"))

paras <- gpcc_df %>% 
  select(date,shape3, rate3) %>%
  group_by(month(date)) %>%
  summarise(shape = mean(shape3), rate = mean(rate3)) %>%
  rename(month =`month(date)`)

n_library <- length(gpcc_df$pr_3M_ave)
### Convert library into a dataframe
library_df <- data.frame(i = seq(1,n_library), spi_thisjuly = as.numeric(gpcc_df$spi3)) %>%
  mutate(spi_nextmay = dplyr::lead(gpcc_df$spi5, n=9,  default = NA))%>% 
  mutate(spi_nextjuly = dplyr::lead(gpcc_df$spi3, n=12,  default = NA))%>% 
  drop_na()

#library_df[1:50,]
library_df <- library_df %>% select(-i)

n_neighbors <- 10
begin.y <- yrs[[1]][1]
end.y <- yrs[[2]][1] -1
predicted_df <- data.frame(year = NA)

for (year_i in c(begin.y:end.y)){ 

#### Eventually put this in a loop through years
date_subset <- seq(as.Date(paste0(year_i,"-08-01")),
                   as.Date(paste0(year_i+1,"-08-01")), 
                     by = "month")-1
  
naspa_subset <- naspa_spi3 %>%
	filter(date %in% date_subset)

naspa_points <- data.frame(spi_thisjuly = naspa_subset$spi3[[1]],
                           spi_nextapr = naspa_subset$spi5[[10]],
                           spi_nextjuly = naspa_subset$spi3[[13]])

if (sum(is.na(naspa_points))>0) {
  next
}
### Find the k closest points

closest <- nn2(data= library_df,
               query = naspa_points, k=n_neighbors)

#want to improve: pick only July 
for(k in seq(1,n_neighbors)){
	
	closest_k <- closest$nn.idx[[k]]
	fragment_k <- library_df[seq(closest_k, closest_k+12),]
	
	naspa_subset_k <- naspa_subset %>%
		mutate(iter = paste0("iter_", k)) %>%
		mutate(spi3 = fragment_k$spi_thisjuly)

	if(k == 1){
		naspa_subset_iter <- naspa_subset_k
		
	} else {
		naspa_subset_iter <- naspa_subset_iter %>% bind_rows(naspa_subset_k)
	}
}
  #smoothed <- loess(spi3 ~ as.numeric(date),data = naspa_subset_iter)
  #predicted_k <- data.frame(date= date_subset,
  #                         spi3 = predict(smoothed, newdata = naspa_subset$date))

if(is.na(predicted_df$year[1]) == TRUE){
  predicted_df <- naspa_subset_iter
    } else {
  predicted_df <- predicted_df %>% bind_rows(naspa_subset_iter)
  }

}

predicted_df <- predicted_df %>%
  right_join(paras, by= "month") %>%
  mutate(prob = pnorm(spi3)) %>%
  mutate(precip = qgamma(prob, shape = shape, rate = rate)) 

ds_naspa_df <- predicted_df %>%
  group_by(date) %>%
  summarise(pr_3M_ave = mean(precip, na.rm = TRUE))

ds_naspa_df <- ds_naspa_df %>%
  mutate(year = year(date)) %>%
  mutate(variable = "pr") %>%
  mutate(model = "Naspa")%>%
  mutate(units = "mm/month") %>%
  mutate(site = loc$index) %>%
  select(date,site, year, variable, model , pr_3M_ave, units)

saveRDS(ds_naspa_df, file = paste0(write_data_path,"/ds_naspa.rds"))

return()

}
###########################################################
###                     plotting                      #####
###########################################################
predicted_df$iter <- as.factor(predicted_df$iter)
temp <- gpcc_df %>% 
  filter(date %in% date_subset) %>%
  mutate(iter = "GPCC") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
 select(year,month,spi3,date, iter,spi5)

temp2 <- rbind(naspa_subset_iter,temp)

p <-ggplot(temp2, aes(x=date, y=spi3)) +
  geom_point(aes(colour= iter)) + 
    geom_smooth(alpha = 0.5) +
  geom_line(data = temp2 %>% filter(iter == "GPCC"), size = 2) +
  labs(title = paste0(year_i)) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  theme_classic(base_size = 18)

p
ggsave(p,filename = paste0(output_path,year_i,loc$site[1],"1930_5OKC.png"))

p <-ggplot(predicted_df %>% filter(year(date) > 1900 & year(date) < 1930), 
           aes(x=date, y=spi3)) +
  geom_line(aes(colour = "black")) + 
  geom_line(data = gpcc_df %>% filter(year(date) > 1900 & year(date) < 1930),
            aes(y = spi3), color = "blue" ) +
 
  labs(title = paste0("1900 - 1930",loc$site),
       color = "data") 
p

ggsave(p,filename = paste0(write_figures_path,"/ds_naspa.png"),width = 6, height = 4)

###########################################################
###                  plotting all NNs and neighbors   #####
###########################################################
p <-ggplot(predicted_df %>% 
             filter(year(date) > 1900 & year(date) < 1910),
           aes(x=date, y= precip)) +
  geom_line(aes(group = iter), color = "skyblue3") +
 # scale_color_brewer(type = "seq") +
  geom_line(data = gpcc_df %>% filter(year(date) > 1900 & year(date) < 1910),
            aes(x = date, y = pr_3M_ave),size = 1.0) +
  #geom_point(naspa_spi3 %>% filter(year > 1990 & year  < 2000 & month == 7), 
  #           mapping = aes(y = spi3, color = "naspa_spi3")) +
  #geom_point(naspa_spi5 %>% filter(year > 1990 & year  < 2000 & month == 4), 
  #          mapping = aes(y = spi5, color = "naspa_spi5")) +
  labs(title = paste0("3months ave.mean and all NNs"),
       y = "3-months Precip(mm/m)",
       color = "data") +
  theme_classic(base_size = 18)

p <- p + stat_summary(fun = "mean", colour = "red", size = 1, geom = "line")
p

ggsave(p,filename = paste0(write_figures_path,loc$site[1],"ds_naspa.png"),width = 6, height = 4 )

#readRDS(file = paste0(output_path,"predictedpreip.rds"))
#Need naspa_df from the file"building_naspa_precip"
p <- ggplot(predicted_df %>% 
              filter(year(date) > 1990 & year(date) < 2020 & month == 7),
            aes(x = date, y= precip)) + 
  geom_line(aes(group = iter), color = "skyblue3") +
  #geom_point(naspa_df %>%
  #            filter(year(date) > 900 & year(date) < 1000), 
  #           mapping = aes(y = precip, color = "naspa")) +
  geom_line(data = gpcc_df %>%filter(year(date) > 1990, month(date) == 7), 
            aes(y = pr_3M_ave,color = "gpcc"),size = 0.5) +
  labs(title ="precip, naspa and predicted",
       y = "3months ave.prcp(mm/m)") +
  scale_color_brewer(palette = "Set1") +
  theme_classic(base_size = 18)
p <- p + stat_summary(fun = "mean", colour = "black", size = 1, geom = "line")
p
ggsave(p,filename = paste0(output_path,loc$site[1],"10.png"),width = 6, height = 4 )


