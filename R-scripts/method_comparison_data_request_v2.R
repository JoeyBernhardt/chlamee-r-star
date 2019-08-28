
########

## Scratch code for extracting high quality time series for conducting sensitivity analysis

########

# NOTE: this will require installing/using the development version of growthTools:

# remove any old versions of the package before installing
# remove.packages('growthTools')
# devtools::install_github("ctkremer/growthTools@develop",auth_token = '6af60cb683ef4cbab1d6b25e0ec2ed6b925c831b',
# 						 build_vignettes = F,upgrade_dependencies=F,force=TRUE)

library(growthTools)
library(tidyverse)
library(cowplot)

# load data set of time series
# fdat <- read.csv("./data/derived_data/201902_all_data_processed_fluor.csv",stringsAsFactors = F)
fdat <- read_csv("data-raw/chlamee-tpc/chlamee-acclimated-rfu-time.csv") %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = days) %>% 
	filter(round == "repeat") %>% 
	filter(population != "COMBO")

fdat <- read_csv("data-raw/chlamee-tpc/chlamee-acute-rfu-time.csv") %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = days) %>% 
	filter(round == "repeat") %>% 
	filter(population != "COMBO")

library(here)
fdat <- read_csv(here("data-processed", "light-rstar-rfus-time.csv")) %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = days) %>% 
	filter(population != "COMBO")
	
fdat <- read_csv('data-raw/chlamee-tpc/globe-chlamy-acclimated-RFU-time.csv') %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = days) %>% 
	filter(population != "COMBO") %>% 
	filter(round == "repeat") %>% 
	filter(RFU != 0)

fdat <- read_csv('data-raw/chlamee-tpc/globe-chlamy-acute-rfu-time.csv') %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = days) %>% 
	filter(population != "COMBO") %>% 
	filter(round == "repeat") %>% 
	filter(RFU != 0)


fdat <- read_csv('data-raw/chlamee-salt-old/chlamee-salt-rfus-time.csv') %>% 
	mutate(log.fluor = log(RFU)) %>% 
	rename(dtime = time) %>% 
	filter(population != "COMBO") %>% 
	filter(RFU != 0) %>% 
	rename(well_plate = unique_id)



# 
# # create single id column for each unique time series
# fdat <- fdat %>%
# 	mutate(id=paste(experiment.ID,isolate.id,media,temperature,dilution,replicate,sep='-')) 

# First fit via AICc, the default
growth.rates.aicc <- fdat %>%
	# filter(temperature %in% c("22", "28", "34", "16", "8")) %>% 
	# filter(well_plate == "B02_32") %>% 
	group_by(well_plate) %>%
	do(grs=get.growth.rate(x=.$dtime,y=.$log.fluor,id=.$well_plate,plot.best.Q=F,model.selection = 'AICc'))


# Process the results:
growth.rates.clean.aicc <- growth.rates.aicc %>%
	summarise(well_plate, mu=grs$best.slope,r2=grs$best.model.rsqr,
                                                           se=grs$best.se,best.model=grs$best.model,
                                                           best.model.n=grs$best.model.slope.n,
                                                           best.model.r2=grs$best.model.slope.r2,
                                                           pre.n=grs$best.model.pre.ns,post.n=grs$best.model.post.ns,
                                                           pre.r2=grs$best.model.pre.rs,post.r2=grs$best.model.post.rs,
                                                           B1=coef(grs$best.model.contents[[1]])['B1'][[1]],
                                                           B2=coef(grs$best.model.contents[[1]])['B2'][[1]])
growth.rates.clean.aicc<-as.data.frame(growth.rates.clean.aicc)
head(growth.rates.clean.aicc)

# save intermediate data (optional)
# write.csv(growth.rates.clean.aicc,"./data/derived_data/All_data_growth_rate_estimates_AICc_2019-02-28.csv",row.names=F)


# Re-fit via AIC
growth.rates.aic <- fdat %>%
	group_by(well_plate) %>%
	# filter(temperature %in% c("22", "28", "34", "16", "8")) %>% 
	do(grs=get.growth.rate(x=.$dtime,y=.$log.fluor,id=.$well_plate,plot.best.Q=F,fpath=fpath2,model.selection = 'AIC'))

# Process the results:
growth.rates.clean.aic <- growth.rates.aic %>%
	summarise(well_plate,
                                                         mu=grs$best.slope,r2=grs$best.model.rsqr,
                                                         se=grs$best.se,best.model=grs$best.model,
                                                         best.model.n=grs$best.model.slope.n,
                                                         best.model.r2=grs$best.model.slope.r2,
                                                         pre.n=grs$best.model.pre.ns,post.n=grs$best.model.post.ns,
                                                         pre.r2=grs$best.model.pre.rs,post.r2=grs$best.model.post.rs,
                                                         B1=coef(grs$best.model.contents[[1]])['B1'][[1]],
                                                         B2=coef(grs$best.model.contents[[1]])['B2'][[1]])
growth.rates.clean.aic<-as.data.frame(growth.rates.clean.aic)
head(growth.rates.clean.aic)

# save intermediate data (optional)
# write.csv(growth.rates.clean.aic,"./data/derived_data/All_data_growth_rate_estimates_AIC_2019-02-28.csv",row.names=F)
# 
# 
# # either load the intermediate data files...
# gdat1<-read.csv("./data/derived_data/All_data_growth_rate_estimates_AIC_2019-02-28.csv")
# gdat2<-read.csv("./data/derived_data/All_data_growth_rate_estimates_AICc_2019-02-28.csv")

# ... or just redefine data frame names and move on
gdat1<-growth.rates.clean.aic
gdat2<-growth.rates.clean.aicc

# find time series that are best fit as gr.lagsat via both IC methods:
id1<-gdat1$well_plate[gdat1$best.model=='gr.lagsat']
id2<-gdat2$well_plate[gdat2$best.model=='gr.lagsat']

# targets:
targets <- gdat1 %>% 
	# filter(well_plate %in% id1[id1 %in% id2] & pre.n>=3 & post.n>=3) %>%
	filter(well_plate %in% id1[id1 %in% id2]) %>% 
	filter(pre.n >=3 , post.n >=3) %>% 
	select(well_plate)
targets<-unlist(targets)

# pull out time series on this list of targets
# that have at least 3 observations in both the lag (pre.n) and saturated (post.n) portions
subset.fdat<-fdat %>% filter(well_plate %in% targets)
subset_fdat_temperature <- subset.fdat
# save output
write_csv(subset.fdat,"data-processed/high_quality_acclimated_chlamee_time_series.csv")
write_csv(subset.fdat,"data-processed/high_quality_chlamee_light_time_series.csv")


subset.fdat %>% 
	ggplot(aes(x = dtime, y = RFU)) + geom_point()


fdat %>% 
	filter(light_level == "10") %>% 
	ggplot(aes(x = dtime, y = RFU)) + geom_point() +
	facet_wrap( ~ well_plate) 
ggsave("figures/light_subsample_time_series.png", width = 30, height = 30 )


subset.fdat %>% 
	ggplot(aes(x = dtime, y = log.fluor)) + geom_point() +
	facet_wrap(~well_plate)
ggsave("figures/lag_sat_time_series.pdf", width = 12, height = 8)
