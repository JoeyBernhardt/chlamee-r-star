library(minpack.lm)
library(broom)
library(tidyverse)
library(rootSolve)
### light growth direct

light_data <- read_csv("data-processed/light-rstar-rfus-time.csv")
treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))

ldata <- light_data %>% 
	mutate(percentage = ifelse(percentage == "0.5-0.7", "0.6", percentage)) %>% 
	mutate(percentage = as.numeric(percentage)) %>% 
	mutate(percentage = percentage/100) %>% 
	mutate(light = percentage*250) %>% 
	group_by(population, light, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(population != "COMBO") %>% 
	filter(!is.na(RFU)) %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0)) 
	

ldata %>% 
	mutate(exponential = case_when(light > 20 & days < 2 ~ "yes",
								   light < 20 & days < 3 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") %>% 
	ggplot(aes(x = days, y = RFU, group = well_plate, color = factor(light))) + geom_point() +
	scale_color_viridis_d() + facet_wrap( ~ population) + geom_line()
ggsave("figures/light-exponential.png", width = 20, height = 15)


light_exp <- ldata %>% 
	mutate(exponential = case_when(light > 20 & days < 2 ~ "yes",
								   light < 20 & days < 3 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0)) %>% 
	left_join(., treatments)


fitting_window <- function(x) {
	
	growth_rates <- ldata %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(light, population, well_plate) %>%
		do(tidy(nls(RFU ~ N0 * exp(r*days),
					data= .,  start=list(r=0.01),
					control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

fitting_window_log_linear <- function(x) {
	growth_rates <- ldata %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(light, population, well_plate) %>% 
		do(tidy(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}



windows <- seq(3,7, by = 1)
multi_fits_log <- windows %>% 
	map_df(fitting_window_log_linear)


multi_fits <- windows %>% 
	map(fitting_window)

multi_fits %>% View

#+ fig.width = 8, fig.height = 6
multi_fits_log %>% 
	filter(term == "days") %>%
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate, color = light)) + geom_point() + geom_line() +
	facet_wrap( ~ population)

exp_fits_top <- multi_fits_log %>%
	filter(term == "days") %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate) 

write_csv(exp_fits_top, "data-processed/light-sequential-fits.csv")

exp_fits_top <- read_csv("data-processed/light-sequential-fits.csv")

light_split <- ldata %>% 
	arrange(days) %>% 
	split(.$well_plate) 


l2 <- light_split %>% 
	map_df(rownames_to_column, var = "time_point") 

light_exponential <- left_join(ldata, treatments, by = "population") %>% 
	left_join(., l2) %>% 
	left_join(., exp_fits_top) %>% 
	filter(time_point <= number_of_points)


results_light <- light_exponential %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(tidy(nlsLM(RFU ~ N0_mean * exp(((light/((1/(alpha * 
												   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light +
											   	(1/alpha))))*(days)),
				  data= .,  
				  start= c(alpha = 0.2, eopt = 100, ps = 1),
				  lower = c(alpha = 0, eopt = 50, ps = 0.1),
				  upper = c(alpha = 4, eopt = 400, ps = 5),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	))


write_csv(results_light, "data-processed/light-EP-params-direct.csv")
results_light <- read_csv("data-processed/light-EP-params-direct.csv")


### try fitting the light data the Monod way

light_exponential %>% 
	filter(population == 8) %>% 
	ggplot(aes(x = days, y = RFU, group = well_plate, color = factor(light))) + geom_line() +
	scale_color_viridis_d() + facet_wrap( ~ light)

results_light_monod <- light_exponential %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(tidy(nlsLM(RFU ~ N0_mean * exp((umax*(light/ (ks+ light)))*(days)),
				  data= .,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 1),
				  upper = c(umax = 6, ks = 200),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))))

write_csv(results_light_monod, "data-processed/light-monod-params-direct.csv")

umax*(light/ (ks+ light))
umax <- 1.468292
ks <- 0.25
light <- 10

### draw the light monod curves
source(here("R-scripts", "predict_monod.R"))

results_light_monod <- read_csv("data-processed/light-monod-params-direct.csv")

bs_split <- results_light_monod %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_l <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_predsl_2 <- left_join(all_preds_l, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE) %>% 
	rename(light = nitrate_concentration.x)

all_predsl_2 %>% 
	ggplot(aes(x = light, y = growth_rate, group = population)) + geom_line() +
	facet_wrap( ~ population)

growth_16 <- all_predsl_2 %>% 
	select(population, light, growth_rate) %>% 
	spread(key = light, value = growth_rate, 2:3) %>% 
	ungroup() %>% 
	select(-population)

pca16 <- rda(growth_16, scale=TRUE)

loadings16 <- scores(pca16,choices=c(1,2))
summary(eigenvals(pca16))

lights <- unique(all_predsl_2$light)
pcs16 <- as_data_frame((loadings16[[1]]))
pc1_16 <- pcs16 %>% 
	mutate(light = lights) 
pc1_16 %>% 
	filter(light != 0) %>% 
	ggplot(aes(x = light, y = PC1)) + geom_line() +
	geom_hline(yintercept = 0) +
	# xlab(bquote('Light ('*mu * mol~ m^-2~s^-1*')')) +
	xlab(expression(Light ~ (mu * mol ~ m^{-2} * s^{-1}))) +
	geom_line() + xlim(0, 250)
ggsave("figures/light-pca.png", width = 6, height = 4)



### bootstrap the Monod light curves

nsplit <- light_exponential %>% 
	split(.$population)

results3<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp((umax*(light/ (ks+ light)))*(days)),
				  data= nsplit[[i]],  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results3 <- bind_rows(results3, hold)
}
m <- 0.56
write_csv(results3, "data-processed/light_monod_params_bootstrapped.csv")
results3 <- read.csv("data-processed/light_monod_params_bootstrapped.csv") %>% 
	mutate(population_index = as.character(population_index))

unique(results3$population_index)
ancestors_i <-  results3 %>%
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	filter(treatment == "none") %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(anc_rstar_q2.5=quantile(rstar, probs=0.025),
			  anc_rstar_q97.5=quantile(rstar, probs=0.975),
			  anc_rstar_mean = mean(rstar)) %>% 
	ungroup() %>% 
	select(-treatment, - population)

rstar_CI_i <- results3 %>%
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	left_join(., ancestors_i) %>% 
	mutate(change_pstar = rstar - anc_rstar_mean) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(change_rstar_lower=quantile(change_pstar, probs=0.025),
			  change_rstar_upper=quantile(change_pstar, probs=0.975),
			  change_rstar_mean = mean(change_pstar)) %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) 


cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")

light_monod <- read_csv("data-processed/light-monod-params-direct.csv") %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(rstar = ks*m/(umax-m))

light_monod_ancestors <- light_monod %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, rstar) %>% 
	rename(ancestral_rstar = rstar) 

all_light_monod <- left_join(light_monod, light_monod_ancestors) %>% 
	mutate(change_rstar = rstar - ancestral_rstar) %>% 
	mutate(treatment = str_replace(treatment, "none", "Ancestors")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS")))


rstar_CI_i %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
	# geom_point(position=position_dodge(width=0.5)) +
	geom_pointrange(size = .2, aes(x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper), position=position_jitter(width=0.1)) +
	geom_hline(yintercept = 0) +
	ylab("Trait change relative to ancestor") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position = "none") 
# +
# 	geom_point(aes(x = treatment, y = change_rstar), data = all_light_monod, color = "red")
ggsave("figures/i-star-trait-change-95CI.png", width = 6, height = 4)

write_csv(rstar_CI_i, "data-processed/change-istar-monod-boot-0.56.csv")

#### bootstrap the EP curves

nsplit <- light_exponential %>% 
	# filter(population %in% c(11)) %>% 
	filter(!is.na(N0_mean), !is.na(light), !is.na(days), !is.na(RFU), !is.na(population)) %>% 
	split(.$population)


results4<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp(((light/((1/(alpha * 
												   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))))*(days)),
				  data= nsplit[[i]],  
				  start= c(alpha = 0.2, eopt = 100, ps = 1),
				  lower = c(alpha = 0, eopt = 50, ps = 0.1),
				  upper = c(alpha = 4, eopt = 400, ps = 5),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results4 <- bind_rows(results4, hold)
}

write_csv(results4, "data-processed/light-ep-bootstrap-params.csv")

results_light_boot_split <- results4 %>% 
	mutate(unique_curve = paste(replicate, population_index, sep = "_")) %>% 
	split(.$unique_curve)

# Now find I* from root of EP ---------------------------------------------

epcurve_m <- function(light, alpha, eopt, ps){
	res <-  (light/((1/(alpha * 
							eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha)))-0.1
	res
}

df <- ep_fits_split[[1]]
df <- results_light_boot_split[[1]]
get_rstar_i <- function(df){
	roots <- uniroot.all(function(x) epcurve_m(x, alpha = df$estimate[[1]], eopt = df$estimate[[2]],
											   ps = df$estimate[[3]]), interval = c(0, 100))
	rstar <- roots[[1]]
	rstar2 <- enframe(x = rstar, value = "rstar")
	return(rstar2)
}

ep_fits_split <- results_light %>%
	# select(1:5) %>% 
	# spread(key = term, value = estimate) %>% 
	split(.$population)


rstars_ep <- ep_fits_split %>% 
	map_df(get_rstar_i, .id = "population") %>% 
	left_join(., treatments, by = "population")

write_csv(rstars_ep, "data-processed/rstars-light-direct.csv")

#### bootstrap I*

library(rootSolve)

get_rstar_i_boot <- function(df){
	roots <- uniroot.all(function(x) epcurve_m(x, alpha = df$alpha[[1]], eopt = df$eopt[[1]],
											   ps = df$ps[[1]]), interval = c(0, 100))
	rstar <- roots[[1]]
	rstar2 <- enframe(x = rstar, value = "rstar")
	return(rstar2)
}

rstars_ep_bootstrap <- results_light_boot_split %>% 
	map_df(get_rstar_i_boot, .id = "unique_curve") 

rstars_ep_bootstrap2 <- rstars_ep_bootstrap %>% 
	separate(unique_curve, into = c("replicate", "population")) %>% 
	left_join(., treatments, by = "population") %>% 
	rename(i_star = rstar)

write_csv(rstars_ep_bootstrap2, "data-processed/light-ep-params-bootstrap.csv")

library(plotrix)
rstars_ep %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	ylab("I* (umols/m2/s)") + xlab("Selection treatment")
ggsave("figures/light-istar-direct.png", width = 8, height = 6)


# draw the EP curves ------------------------------------------------------

results_light <- read_csv("data-processed/light-EP-params-direct.csv")


