

listviewer::jsonedit(subset_f10)

?get.growth.rate


all_preds2 %>% 
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
	geom_point(aes(x = days, y = y, shape = exponential)) +
	facet_wrap(~ well_plate) + ylab("Ln(RFU)") +xlab("Days")
# ggsave(here("figures", "all_growth_tools_plots_well_plate.pdf"), width = 45, height = 45)
ggsave(here("figures", "all_growth_tools_plots_well_plate_AIC.png"), width = 45, height = 45)


growth_sum_params2 <- left_join(growth_sum_params, key, by = "well_plate")

write_csv(growth_sum_params2, here("data-processed", "nitrate_exp_growth_w_growthtools_AIC.csv"))



all_growth_n <- read_csv(here("data-processed", "nitrate_exp_growth_w_growthtools_AIC.csv")) %>% 
	filter(population != "COMBO")

all_growth_n_c <- read_csv(here("data-processed", "nitrate_exp_growth_w_growthtools.csv")) %>% 
	filter(population != "COMBO")

monod_fits_c <- all_growth_n_c %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	rename(estimate = mu) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))




treatments <- read_csv(here("data-processed", "nitrate-treatments-processed.csv"))

prediction_function <- function(df) {
	
	
	monodcurve<-function(x){
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

bs_split <- monod_fits %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)

bs_splitc <- monod_fits_c %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_n <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(prediction_function, .id = "population")

all_preds_nc <- bs_splitc %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(prediction_function, .id = "population")
all_predsn_2 <- left_join(all_preds_n, treatments, by = c("population"))
all_predsn_2c <- left_join(all_preds_nc, treatments, by = c("population"))



Plot the Monod fits

```{r, eval = FALSE}


filter(all_predsn_2, nitrate_concentration < 500, ancestor_id == "anc5", treatment == "N") %>% View

all_growth_n2 <- left_join(all_growth_n, treatments)

all_growth_n2c <- left_join(all_growth_n_c, treatments)

all_growth_n2c %>% 
	mutate(estimate = mu) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	filter(nitrate_concentration < 1500, ancestor_id == "anc4", treatment == "C") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data=filter(all_predsn_2, nitrate_concentration < 1500, ancestor_id == "anc4", treatment == "C"), aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") + xlim(0, 110)
ggsave(here("figures", "nitrate_monod_growthTools.png"), width = 12, height = 8)


monod_wide <-  monod_fits %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
monod_curve_mortality <- function(nitrate_concentration, umax, ks){
	res <- (umax* (nitrate_concentration / (ks+ nitrate_concentration))) - 0.1
	res
}	


## find the R*
rstars <- monod_wide %>% 
	mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629") 

rstars2 <- left_join(rstars, treatments, by = "population")
rstars2 %>%
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id), size = 2, data = rstars2, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-means-AIC.png"), width = 6, height = 4)

```