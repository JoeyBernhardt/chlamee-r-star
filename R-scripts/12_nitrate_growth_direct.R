library(readxl)
library(tidyverse)
library(minpack.lm)
library(here)
library(janitor)
library(cowplot)
library(broom)

treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")

nitrate_abundances <- read_csv("data-processed/nitrate-RFUs-time-exponential-phase.csv") %>% 
	filter(population != "COMBO") 

nitrate_split <- nitrate_abundances %>% 
	arrange(days) %>% 
	split(.$well_plate) 


n2 <- nitrate_split %>% 
	map_df(rownames_to_column, var = "time_point") 

# t() %>% 
# 	as.data.frame() %>% 
# 	rownames_to_column(var = "well_plate") %>%
# 	filter(well_plate != "well_plate") %>% 
# 	rename(number_time_points = V1)
	


exponential_fits <- read_csv("data-processed/nitrate-sequential-fits.csv")

nitrate <- left_join(nitrate_abundances, treatments, by = "population") %>% 
	left_join(., n2) %>% 
	left_join(., exponential_fits) %>% 
	filter(time_point <= number_of_points)
	

nitrate %>% 
	ggplot(aes(x = days, y = RFU, group = well_plate, color = factor(nitrate_concentration))) + geom_point() +
	scale_color_viridis_d() + facet_wrap( ~ population) + geom_line()
ggsave("figures/nitrate-exponential.png", width = 20, height = 15)

unique(nitrate$population)


nitrate_n0 <- nitrate %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0))



results <- nitrate_n0 %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(tidy(nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ (ks+ nitrate_concentration)))*(days)),
				  data= .,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))))

results_alpha <- nitrate_n0 %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(tidy(nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ ((umax/alpha)+ nitrate_concentration)))*(days)),
				  data= .,  
				  start= c(umax = 1, alpha = 0.01),
				  lower = c(umax = 0, alpha= 0),
				  upper = c(umax = 3, alpha = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))))

results_alpha %>% 
	ggplot(aes(x = treatment, y = estimate)) + geom_point() +
	facet_wrap( ~ term, scales = "free")


write_csv(results, "data-processed/nitrate-monod-params-direct.csv")
results <- read_csv("data-processed/nitrate-monod-params-direct.csv")

results2 <- nitrate_n0 %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(fit = nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ (ks+ nitrate_concentration)))*(days)),
				  data= .,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))) %>% 
	as.data.frame(nlsBoot(.$fit))





m <- 0.1
res_wide <- results %>% 
	select(ancestor_id, treatment, population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(rstar = ks*m/(umax-m))

### plot R* vs 1/umax

write_csv(res_wide, "data-processed/nitrate-rstars-direct.csv")

n_stars <- read_csv("data-processed/nitrate-rstars-direct.csv")

n_stars %>% 
	ggplot(aes(x = rstar, y = 1/(umax- m))) + geom_point() +
	geom_smooth(method = "lm")

library(plotrix)
res_wide %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar, umax) %>% 
	ggplot(aes(x = reorder(treatment, rstar_mean), y = rstar_mean)) + geom_point() +
	geom_errorbar(aes(ymin = rstar_mean - rstar_std.error, ymax = rstar_mean + rstar_std.error), width = 0.1) +
	ylab("R* (uM N)") + xlab("Selection treatment")
ggsave("figures/nitrate-rstar-direct.png", width = 8, height = 6)


res_wide %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar, umax) %>% 
	ggplot(aes(x = reorder(treatment, umax_mean), y = umax_mean)) + geom_point() +
	geom_errorbar(aes(ymin = umax_mean - umax_std.error, ymax = umax_mean + umax_std.error), width = 0.1) +
	ylab("umax (per day)") + xlab("Selection treatment")
ggsave("figures/nitrate-umax-direct.png", width = 8, height = 6)


mod1 <- lm(umax ~ treatment, data = res_wide)
summary(mod1)


res_wide %>% 
	group_by(treatment) %>% 
	ggplot(aes(x = ancestor_id, y = rstar, color = treatment)) + geom_point() +
	# geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	ylab("R* (uM N)") + xlab("Selection treatment")



# plot the fits -----------------------------------------------------------

predict_abundance <- function(df){

	growth_fun <- function(x){
	abundance <-  df$N0_mean[[1]] * exp((df$umax[[1]]*(df$nitrate_concentration[[1]]/ (df$ks[[1]]+ df$nitrate_concentration[[1]])))*(x))
	abundance
}

pred <- function(x) {
	y <- growth_fun(x)
}

x <- seq(0, 2, by = 0.1)

preds <- sapply(x, pred)
preds <- data.frame(x, preds) %>% 
	rename(days.x = x, 
		   abundance = preds)
}


nitrate_setup <- nitrate_n0 %>% 
	group_by(ancestor_id, treatment, population, well_plate) %>% 
	distinct(nitrate_concentration, N0_mean) 

nitrate_setup2 <- left_join(nitrate_setup, res_wide)

nitrate_split <- nitrate_setup2 %>% 
	split(.$well_plate)

abundances <- nitrate_split %>% 
	map_df(predict_abundance, .id = "well_plate")

abundances2 <- left_join(abundances, nitrate_setup)

abundances2 %>% 
	ggplot(aes(x = days.x, y = abundance, group = well_plate, color = factor(nitrate_concentration))) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = factor(nitrate_concentration)), data = nitrate_n0) +
	facet_grid(ancestor_id ~ treatment, scales = "free_y") + scale_color_viridis_d() +
	xlab("Days") + ylab("RFU")
ggsave("figures/nitrate-direct-fits-population.png", width = 15, height = 15)

abundances2 %>% 
	ggplot(aes(x = days.x, y = abundance, group = well_plate, color = factor(nitrate_concentration))) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = factor(nitrate_concentration)), data = nitrate_n0) +
	facet_wrap( ~ well_plate, scales = "free_y") + scale_color_viridis_d() +
	xlab("Days") + ylab("RFU")
ggsave("figures/nitrate-direct-fits.png", width = 45, height = 35)


#### bootstrap



fitc <- nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ (ks+ nitrate_concentration)))*(days)),
				  data= nsplit[[30]],  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
nls_boot_c <- nlsBoot(fitc, niter = 999)
boots30 <- data.frame(nls_boot_c$coefboot)


fitcb <- nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ ((umax/alpha)+ nitrate_concentration)))*(days)),
			  data= nsplit[[3]],  
			  start= c(umax = 1, alpha = 0.1),
			  lower = c(umax = 0, alpha = 0.00001),
			  upper = c(umax = 3, alpha = 10),
			  control = nls.control(maxiter=1024, minFactor=1/204800000))
nls_boot_c <- nlsBoot(fitcb, niter = 999)
boots30 <- data.frame(nls_boot_c$coefboot)

summary(fitc)
best_fit_c <- coef(fitc)
nlsResiduals(fitc)

nls_boot_coefs_c <- as_data_frame(nls_boot_c$coefboot)
nls_boot_coefs_c %>% 
	summarise_each(funs(mean, std.error), umax, ks) %>% View
nls_boot_coefs_c %>% 
	summarise(ks_q2.5=quantile(ks, probs=0.025),
		  ks_q97.5=quantile(ks, probs=0.975),
		  ks_mean = mean(ks)) %>% View


df <- nsplit[[12]]

get_bootstrap <- function(df){
	nls_boot_c <- nlsBoot(nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ (ks+ nitrate_concentration)))*(days)),
				  data= df,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000)), niter = 999)
	# nls_boot_coefs_c <- as_tibble(nls_boot_c$coefboot)
	# nls_boot_coefs_c <- nls_boot_coefs_c %>% 
	# 	mutate(population = unique(df$population))
	
}


df <- nsplit[[1]]

library(nlstools)
boot_fits <- nsplit %>% 
	map(get_bootstrap)

boot_fits

nsplit <- nitrate_n0 %>% 
	# filter(population %in% c(11)) %>% 
	filter(!is.na(N0_mean), !is.na(nitrate_concentration), !is.na(days), !is.na(RFU), !is.na(population)) %>% 
	split(.$population)


unique(nitrate_n0$population)
View(nsplit[[5]])
i <-5
library(nlstools)
results3<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ (ks+ nitrate_concentration)))*(days)),
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

results3b<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp((umax*(nitrate_concentration/ ((umax/alpha)+ nitrate_concentration)))*(days)),
				  data= nsplit[[i]],  
				  start= c(umax = 1, alpha = 0.2),
				  lower = c(umax = 0, alpha= 0.00001),
				  upper = c(umax = 3, alpha = 10),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results3b <- bind_rows(results3b, hold)
}

results3b
unique(results3$population_index)
write_csv(results3, "data-processed/nitrate-monod-params-bootstrap.csv")
write_csv(results3b, "data-processed/nitrate-monod-alpha-params-bootstrap.csv")

rstar_CI_n_alpha <- results3b %>%
	mutate(rstar = (umax/alpha)*m/(umax-m)) %>% 
	group_by(population_index) %>% 
	summarise(rstar_lower=quantile(rstar, probs=0.025),
			  rstar_upper=quantile(rstar, probs=0.975),
			  rstar_mean = mean(rstar),
			  alpha_lower = quantile(alpha, probs=0.025),
			  alpha_upper=quantile(alpha, probs=0.975),
			  alpha_mean = mean(alpha)) %>% 
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	filter(treatment %in% c("Ancestors", "N"))


ggplot() +
	geom_pointrange(aes(color = ancestor_id, x = treatment, y = alpha_mean, ymin = alpha_lower, ymax = alpha_upper),
					data = rstar_CI_n_alpha, position=position_jitter(width=0.1)) +
	ylim(0, 0.2)

ggplot() +
	geom_pointrange(aes(x = treatment, y = alpha_mean, ymin = alpha_lower, ymax = alpha_upper, color = ancestor_id),
					data = rstar_CI_n_alpha) +
	geom_line(aes(x = treatment, y = alpha_mean, color = ancestor_id, group = ancestor_id), data = rstar_CI_n_alpha)


results3 <- read.csv("data-processed/nitrate-monod-params-bootstrap.csv") %>% 
	mutate(population_index = as.character(population_index))
unique(results3c$population_index)
m <- 0.56
rstar_CI_n <- results3 %>%
	mutate(unique_id = paste(population_index, replicate, sep = "_")) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(population_index) %>% 
	summarise(rstar_q2.5=quantile(rstar, probs=0.025),
			  rstar_q97.5=quantile(rstar, probs=0.975),
			  rstar_mean = mean(rstar)) 

unique(results3$population_index)
ancestors_n <-  results3 %>% 
	dplyr::rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	filter(treatment == "none") %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	dplyr::group_by(treatment, ancestor_id, population) %>% 
	summarise(anc_rstar_q2.5=quantile(rstar, probs=0.025),
			  anc_rstar_q97.5=quantile(rstar, probs=0.975),
			  anc_rstar_mean = mean(rstar)) %>%
	ungroup() %>% 
	select(-treatment, - population)

rstar_CI_n <- results3 %>%
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	left_join(., ancestors_n) %>% 
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
rstar_CI_n %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
	# geom_point(position=position_dodge(width=0.5)) +
	geom_pointrange(size = .2, aes(x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper), position=position_jitter(width=0.1)) +
	geom_hline(yintercept = 0) +
	ylab("Trait change relative to ancestor") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position = "none")
ggsave("figures/n-star-trait-change-95CI.png", width = 6, height = 4)

write_csv(rstar_CI_n, "data-processed/change-nstar-monod-boot-0.56.csv")
predict_monod <- function(df) {
	
	monodcurve <- function(x){
		growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 1000, by = 1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(nitrate_concentration.x = x, 
			   growth_rate = preds)
}


monod_split <- res_wide %>% 
	split(.$population)

monod_split <- results3 %>% 
	mutate(unique_id = paste(population_index, replicate, sep = "_")) %>% 
	split(.$unique_id)


predictions_nitrate <- monod_split %>% 
	map_df(predict_monod, .id = "unique_id") 
## ended here

tail(predictions_nitrate)

predictions_nitrate <- predictions_nitrate %>% 
	separate(unique_id, into = c("population", "replicate"))

predictions_nitrate <- predictions_nitrate %>%  
	left_join(., treatments)



predictions_nitrate %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate)) + geom_line(alpha = 0.1) +
	facet_grid(ancestor_id ~ treatment) + ylab("Growth rate (per day)") + xlab("Nitrate (uM N)")
ggsave("figures/nitrate-monod-bootstrap.png", width = 12, height = 8)

names(predictions_nitrate)
predictions_summary <- predictions_nitrate %>% 
	group_by(ancestor_id, treatment, population, nitrate_concentration.x) %>% 
	summarise(growth_rate_q2.5=quantile(growth_rate, probs=0.025),
			  growth_rate_q97.5=quantile(growth_rate, probs=0.975),
			  growth_rate_mean = mean(growth_rate)) 

write_csv(predictions_summary, "data-processed/nitrate-bootstrap-monod-summary.csv")

predictions_summary <- read_csv("data-processed/nitrate-bootstrap-monod-summary.csv")

cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")

predictions_summary %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate_mean, group = population)) + geom_line() +
	geom_ribbon(aes(ymin = growth_rate_q2.5, ymax = growth_rate_q97.5), fill = "grey70", alpha = 0.5) +
	facet_grid(ancestor_id ~ treatment) +ylab("Growth rate (per day)") + xlab("Nitrate (uM N)")
ggsave("figures/nitrate-monod-bootstrap-CI.png", width = 13, height = 7)



predictions_summary %>% 
	filter(treatment == "none") %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate_mean, group = ancestor_id, color = ancestor_id)) + geom_line(size = 1) +
	# facet_grid(ancestor_id ~ treatment) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (per day)") + xlab("Nitrate (uM N)") + xlim(0, 15) + ylim(0, 0.7)
ggsave("figures/ancestors-nitrate-monod-low-N.png", width = 6, height = 4)


predictions_summary %>% 
	filter(treatment == "none") %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate_mean, group = ancestor_id, color = ancestor_id)) + geom_line(size = 1) +
	# facet_grid(ancestor_id ~ treatment) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (per day)") + xlab("Nitrate (uM N)") + xlim(900, 1000) + ylim(1.7, 1.9)

predictions_summary %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>%
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	# filter(treatment == "none") %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate_mean, group = ancestor_id, color = treatment)) + geom_line() +
	geom_ribbon(aes(ymin = growth_rate_q2.5, ymax = growth_rate_q97.5, fill = treatment), alpha = 0.2, linetype = "blank") +
	facet_wrap(~ treatment, nrow = 2) +ylab("Growth rate (per day)") + xlab("Nitrate (uM N)") + 
	scale_color_manual(values = cols) +
	scale_fill_manual(values = cols)
ggsave("figures/nitrate-monod-bootstrap.png", width = 12, height = 6)



predictions_summary %>%
	group_by(population, nitrate_concentration.x, ancestor_id) %>% 
	summarise(growth_rate_mean = mean(growth_rate_mean)) %>% 
	# filter(treatment == "none", ancestor_id %in% c("anc2")) %>% 
	ggplot(aes(x = nitrate_concentration.x, y = growth_rate_mean, group = population, color = ancestor_id)) + geom_line() +
	# facet_wrap(ancestor_id ~ treatment) +
	# geom_ribbon(aes(ymin = growth_rate_q2.5, ymax = growth_rate_q97.5), fill = "grey70", alpha = 0.5) +
ylab("Growth rate (per day)") + xlab("Nitrate (uM N)")
ggsave("figures/monod-nitrate-ancestors.png", width = 8, height = 6)


### do PCA on these curves

growth_16 <- predictions_summary %>% 
	select(population, nitrate_concentration.x, growth_rate_mean) %>% 
	group_by(population, nitrate_concentration.x) %>% 
	summarise(growth_rate_mean = mean(growth_rate_mean)) %>% 
	# filter(temperature %in% c(10, 16, 22, 28, 34, 40)) %>% 
	spread(key = nitrate_concentration.x, value = growth_rate_mean, 2:3) %>% 
	ungroup() %>% 
	select(-population)

unique(predictions_summary$population)

cov_mat <- cov(growth_16)
pca_res <- prcomp(cov_mat, center = TRUE,scale. = TRUE)
pca_res <- prcomp(cov_mat, center = TRUE,scale. = TRUE)
summary(pca_res)
ggbiplot(pca_res, labels= tpc_temps)
ggsave("figures/PCA-non-acclimated-biplot.pdf", width = 6, height = 6)


library(vegan)
pca16 <- rda(growth_16)
summary(pca16)


pca16 <- rda(growth_16, scale=TRUE)

loadings16 <- scores(pca16,choices=c(1,2))
summary(eigenvals(pca16))

nitrate_conc <- unique(predictions_summary$nitrate_concentration.x)
pcs16 <- as_data_frame((loadings16[[1]]))
pc1_16 <- pcs16 %>% 
	mutate(nitrate = nitrate_conc) 
pc1_16 %>% 
	filter(nitrate != 0) %>% 
	ggplot(aes(x = nitrate, y = PC1)) + geom_line() +
	xlim(0, 1000) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Nitrate concentration (uM N)") + geom_line()
ggsave("figures/nitrate-pca.png", width = 8, height = 6)


## is there a ks umax trade-off?

res_wide %>% 
	ggplot(aes(x = ks, y = umax)) + geom_point() +
	geom_smooth(method = "lm", color = "grey") + ylab("Max growth rate (umax)") + xlab("Half saturation constant (ks)")
ggsave("figures/nitrate-umax-ks-trade-off.png", width = 8, height = 6)

