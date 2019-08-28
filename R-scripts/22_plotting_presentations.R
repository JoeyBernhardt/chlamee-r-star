

library(tidyverse)
library(cowplot)
library(broom)
library(readxl)
library(plotrix)


all_stars3 <- read_csv("data-processed/all-rstars.csv")


#### plots for Anita
cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")


all_stars4 <- all_stars3 %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS")))

all_stars4_sum <- all_stars4 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), p_star, n_star, i_star)

### phosphorus plots

all_stars4_sum %>% 
	ggplot() +
	# geom_point(aes(x = treatment, y = p_star, color = treatment), size = 2, alpha = 0.5, data = all_stars4) +
	geom_point(aes(x = treatment, y = p_star_mean, color = treatment), size = 2, data = all_stars4_sum) +
	geom_errorbar(aes(x = treatment, ymin = p_star_mean - p_star_std.error, ymax = p_star_mean + p_star_std.error, color = treatment),
				  data = all_stars4_sum, width = 0.1) +
	theme_classic() + 
	theme(legend.position="none") +
	scale_color_manual(values=cols) +
	ylab("P* (uM P)") + xlab("Selection treatment")
ggsave("figures/p_star_dots.png", width = 6, height = 4)
ggsave("figures/p_star_dots_means_only.png", width = 6, height = 4)


### nitrogen plots

all_stars4_sum %>% 
	ggplot() +
	# geom_point(aes(x = treatment, y = n_star, color = treatment), size = 2, alpha = 0.5, data = all_stars4) +
	geom_point(aes(x = treatment, y = n_star_mean, color = treatment), size = 2, data = all_stars4_sum) +
	geom_errorbar(aes(x = treatment, ymin = n_star_mean - n_star_std.error, ymax = n_star_mean + n_star_std.error, color = treatment),
				  data = all_stars4_sum, width = 0.1) +
	theme_classic() + 
	theme(legend.position="none") +
	scale_color_manual(values=cols) +
	ylab("N* (uM N)") + xlab("Selection treatment")
ggsave("figures/n_star_dots.png", width = 6, height = 4)
ggsave("figures/n_star_dots_means_only.png", width = 6, height = 4)



### light plots

all_stars4_sum %>% 
	ggplot() +
	# geom_point(aes(x = treatment, y = i_star, color = treatment), size = 2, alpha = 0.5, data = all_stars4) +
	geom_point(aes(x = treatment, y = i_star_mean, color = treatment), size = 2, data = all_stars4_sum) +
	geom_errorbar(aes(x = treatment, ymin = i_star_mean - i_star_std.error, ymax = i_star_mean + i_star_std.error, color = treatment),
				  data = all_stars4_sum, width = 0.1) +
	theme_classic() + 
	theme(legend.position="none") +
	scale_color_manual(values=cols) +
	ylab(bquote('I* ('*mu~'mol'~ m^-2~s^-1*')')) + xlab("Selection treatment")
ggsave("figures/i_star_dots.png", width = 6, height = 4)
ggsave("figures/i_star_dots_means_only.png", width = 6, height = 4)



# Trade-off plots ---------------------------------------------------------

changes_all <- read_csv("data-processed/changes_all.csv") %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS")))
changes_sum <- changes_all %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_n_star, change_i_star, change_p_star) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS")))

changes_all %>% 
	ggplot(aes(x = change_n_star, y = change_p_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.09, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.7, xmax=0, ymin=0, ymax=0.09),
			  color="transparent", alpha=0.5, fill = "grey") + 
	geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	theme_classic() + 
	theme(legend.position="none") +
	geom_point(aes(x = change_n_star_mean, y = change_p_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	geom_point(aes(x = change_n_star_mean, y = change_p_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	geom_point(size = 3) +
	geom_point(size = 3, shape = 1, color = "black") +
	geom_pointrange(aes(x = change_n_star_mean, y = change_p_star_mean, ymin = change_p_star_mean - change_p_star_std.error, ymax = change_p_star_mean + change_p_star_std.error, color = treatment), data = changes_sum) +
	geom_segment(aes(x = change_n_star_mean - change_n_star_std.error, xend = change_n_star_mean + change_n_star_std.error,  y = change_p_star_mean, yend = change_p_star_mean, color = treatment), data = changes_sum) +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	scale_color_manual(values=cols, name = "Selection treatment") +
	ylim(-0.09, 0.09) + xlim(-1.7, 1.7) +
	ylab("Change in P* (um P)") + xlab("Change in N* (um N)")
ggsave("figures/n-star-p-star-trade-off-T-color-treatment.png", width = 6.5, height = 5)



changes_all %>% 
	ggplot(aes(x = change_i_star, y = change_p_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.09, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.7, xmax=0, ymin=0, ymax=0.09),
			  color="transparent", alpha=0.5, fill = "grey") +
	# geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	# tableau_color_pal(palette = "Tableau 10") +
	# geom_point(aes(x = change_i_star_mean, y = change_p_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	# geom_point(aes(x = change_i_star_mean, y = change_p_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	# geom_point(size = 3) +
	theme_classic() + 
	theme(legend.position="none") +
	# geom_point(size = 3, shape = 1, color = "black") +
	# scale_color_manual(values=cols, name = "Selection treatment") +
	# geom_pointrange(aes(x = change_i_star_mean, y = change_p_star_mean,  ymin = change_p_star_mean - change_p_star_std.error, ymax = change_p_star_mean + change_p_star_std.error, color = treatment), data = changes_sum) +
	# geom_segment(aes(x = change_i_star_mean - change_i_star_std.error, xend = change_i_star_mean + change_i_star_std.error,  y = change_p_star_mean, yend = change_p_star_mean, color = treatment), data = changes_sum) +
	# # geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	ylab("Change in P* (um P)") + 
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')')) +
	ylim(-0.09, 0.09) + xlim(-1.7, 1.7) 
ggsave("figures/p-star-i-star-trade-off-T-color-treatment-no-points.png", width = 6.5, height = 5)



changes_all %>% 
	ggplot(aes(x = change_i_star, y = change_n_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.9, ymin=-0.9, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.9, xmax=0, ymin=0, ymax=0.9),
			  color="transparent", alpha=0.5, fill = "grey") +
	# geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	theme_classic() + 
	theme(legend.position="none") +
	# tableau_color_pal(palette = "Tableau 10") +
	# geom_point(aes(x = change_i_star_mean, y = change_n_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	# geom_point(aes(x = change_i_star_mean, y = change_n_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	# geom_point(size = 3) +
	# geom_point(size = 3, shape = 1, color = "black") +
	# geom_pointrange(aes(x = change_i_star_mean, y = change_n_star_mean,  ymin = change_n_star_mean - change_n_star_std.error, ymax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	# geom_segment(aes(x = change_i_star_mean - change_i_star_std.error, xend = change_i_star_mean + change_i_star_std.error,  y = change_n_star_mean, yend = change_n_star_mean, color = treatment), data = changes_sum) +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	scale_color_manual(values=cols, name = "Selection treatment") +
	ylab("Change in N* (um N)") + xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')'))
ggsave("figures/n-star-i-star-trade-off-T-color-treatment-no-points.png", width = 6.5, height = 5)



### salt
library(here)
library(janitor)
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")


salt <- read.csv("data-processed/salt_growth_rates.csv") %>% 
	separate(salt_level, into = c("S", "concentration"), sep = 1) %>% 
	mutate(concentration = as.numeric(concentration)) %>%
	mutate(population = as.character(population)) %>% 
	left_join(., treatments, by = "population") %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS")))


salt %>% 
	filter(treatment == "A") %>% 
	ggplot(aes(x = concentration, y = mu, color = treatment)) + geom_point(alpha = 0.2, size = 2) +
	geom_point(shape = 1, aes(color = treatment), size = 2) +
	# facet_wrap( ~ treatment) +
	geom_smooth(aes(color = treatment), se = FALSE) +
	scale_color_manual(values=cols, name = "Selection treatment") +
	ylab("Growth rate (/day)") + xlab("Salt concentration (g/L)")
ggsave("figures/salt-growth-rates.png", width = 8, height = 5)



### fit some sort of decreasing function to the data

unique(salt$population)

salt %>% 
	filter(is.na(population)) %>% View


# salt_tolerances <- salt %>% 
# 	group_by(treatment, ancestor_id, population) %>% 
# 	do(tidy(nls(mu ~ a/(b + concentration),
# 				data= .,  start=list(a = 1, b = 1), algorithm="port", lower=list(c=0.01, d=0),
# 				control = nls.control(maxiter=500, minFactor=1/204800000))))


salt_tolerances <- salt %>% 
	filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	do(tidy(nls(mu ~ a/(1 + exp(-b * (concentration-c))),
				data= .,  start=list(a = 1, b = 0.5, c = 3), algorithm="port",
				control = nls.control(maxiter=500, minFactor=1/204800000))))

salt %>% 
	filter(treatment == "A") %>% 
	ggplot(aes(x = concentration, y = mu, color = treatment)) + geom_point()

library(nls.multstart)

## a simplified form of the common logistic function
## f(x) = a/(1 + exp(-b*x + c))

fits <- salt %>% 
	# filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(mu ~ a/(1 + exp(-b * (concentration-c))),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 1, b = -3, c = 3),
												  start_upper = c(a = 2, b = 1, c = 5),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.1, b = -10, c = 1),
												  upper = c(a = 10, b = 20, c = 100),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

fits2 <- salt %>% 
	# filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(mu ~ a/(1 + exp(-b*concentration + c)),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 1, b = -3, c = 3),
												  start_upper = c(a = 2, b = 1, c = 5),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.001, b = -100, c = -100),
												  upper = c(a = 10, b = 200, c = 100),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

params2 <- fits2 %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 
params <- fits %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 

p1 <- params %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(salt_tolerance = c)
	

p2  <- params2 %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(salt_tolerance = c/b) 

params %>% 
	filter(term == "c") %>% 
	ggplot(aes(x = treatment, y = estimate)) + geom_point()


ssplit <- params %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$population)

ssplit2 <- params2 %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$population)

## now we need to plot the fits
predict_growth <- function(df){
	
	predict_salt <- function(x){
		growth.rate <- df$a[[1]]/(1 + exp(-df$b[[1]] * (x-df$c[[1]])))
		return(growth.rate)
	}
	
	pred <- function(x) {
		y <- predict_salt(x)
	}
	
	x <- seq(0, 60, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(salt_concentration = x, 
			   growth_rate = preds)
	
}

predict_growth2 <- function(df){
	
	predict_salt <- function(x){
		growth.rate <- df$a[[1]]/(1 + exp(-df$b[[1]]*x + df$c[[1]]))
		return(growth.rate)
	}
	
	pred <- function(x) {
		y <- predict_salt(x)
	}
	
	x <- seq(0, 60, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(salt_concentration = x, 
			   growth_rate = preds)
	
}

df <- ssplit[[1]]
predictions_salt <- ssplit %>% 
	map_df(predict_growth, .id = "population")

predictions_salt2 <- ssplit2 %>% 
	map_df(predict_growth2, .id = "population")

salt %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt) +
	geom_vline(aes(xintercept = salt_tolerance), data = p1, color = "lightblue") + 
	facet_wrap( ~ population, scales = "free") + xlim(0, 10) + ylab("Population growth rate") + xlab("Salt concentration (g/L)")
ggsave("figures/predictions-salt-sigmoid.png", width = 15, height = 10)

salt %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt2, color = "lightblue") +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt, color = "orange") +
	geom_vline(aes(xintercept = salt_tolerance), data = p2, color = "lightblue") + 
	geom_vline(aes(xintercept = salt_tolerance), data = p1, color = "orange") + 
	facet_wrap( ~ population, scales = "free") + xlim(0, 10) + ylab("Population growth rate") + xlab("Salt concentration (g/L)")
ggsave("figures/predictions-salt-sigmoid-logistic-comparison.png", width = 15, height = 10)

sl2 <- salt_tolerances %>% 
	filter(term == "b") %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	rename(salt_tolerance = estimate)

write_csv(sl2, "data-processed/salt_tolerances.csv")



ssplit <- salt_tolerances %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$population)


predict_growth <- function(df){

predict_salt <- function(x){
	growth.rate <- df$a[[1]]/(df$b[[1]] + x)
	return(growth.rate)
}
	
	pred <- function(x) {
		y <- predict_salt(x)
	}
	
	x <- seq(0, 60, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(salt_concentration = x, 
			   growth_rate = preds)
	
}

df <- ssplit[[1]]
predictions_salt <- ssplit %>% 
	map_df(predict_growth, .id = "population")

salt %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt) +
	facet_wrap( ~ population, scales = "free") + xlim(0, 10)
ggsave("figures/salt-predictions.png", width = 12, height = 8)



salt_summary %>% 
	ggplot(aes(x = treatment, y = mean, color = treatment)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# geom_point(aes(x = treatment, y = estimate, color = treatment), data = sl2, alpha = 0.5, size = 2) +
	scale_color_manual(values=cols, name = "Selection treatment") +ylab("Salt tolerance (g/L)") + 
	xlab("Selection treatment") +
	theme_classic() + 
	theme(legend.position="none") 
ggsave("figures/salt_tol_dots.png", width = 6, height = 4)
ggsave("figures/salt_tol_dots-means.png", width = 6, height = 4)

salt_summary <- salt_tolerances %>% 
	filter(term == "b") %>%
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), estimate) 


library(beyonce)
all_changes <- read_csv("data-processed/changes_allb.csv")

cols2 <- c("#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")



levels(allc$trait)

all_changes %>% 
	ggplot(aes(x = n_star, y = n_ks)) + geom_point() +
	geom_smooth(method = "lm")


all_changes %>% 
	filter(treatment != "none") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	select(ancestor_id, treatment, population, contains("change")) %>%
	mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
	mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
	mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
	mutate(`ks_n` = scale(change_n_ks, center = FALSE)) %>% 
	mutate(`ks_p` = scale(change_p_ks, center = FALSE)) %>% 
	mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
	gather(key = trait, value = change, 13:18) %>% 
	mutate(trait = as.factor(trait)) %>% 
	mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*", "ks_n", "ks_p"))) %>% 
	# group_by(trait) %>% 
	# arrange(change) %>% 
	ggplot(aes(x = treatment, y = change, color = treatment)) + geom_jitter(width = 0.2, size =2, alpha = 0.5) +
	xlab("") + ylab("Standardized trait change") + ylim(-4.3, 4.3) +
	geom_hline(yintercept = 0) +
	scale_color_manual(values=cols2, name = "Selection treatment") +
	facet_wrap( ~ trait) +
	# theme(strip.background = element_blank(),
	# strip.text.y = element_blank()) +
	# geom_point(aes(x = treatment, y = mean), data = all_changes_sum, color = "grey", shape = 5, size = 5) +
	geom_pointrange(aes(x = treatment, y = mean, ymin = mean + std.error, ymax = mean - std.error, color = treatment), data = all_changes_sum, size = 0.75)
	ggsave("figures/standardized-trait-change-facet-means-center-grey.png", width = 10, height = 6)

	
	
	all_changes_sum2 <- all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>% 
		select(ancestor_id, treatment, population, contains("change")) %>% 
		mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
		mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
		mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
		mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		group_by(treatment, trait) %>% 
		summarise_each(funs(mean, std.error), change) %>% 
		mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*"))) 
	
	all_changes_sum3 <- all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>% 
		select(ancestor_id, treatment, population, contains("change")) %>% 
		mutate(`N*` = change_n_star) %>% 
		mutate(`P*` = change_p_star) %>% 
		mutate(`I*` = change_i_star) %>% 
		mutate(`Salt tolerance` = change_salt_tol) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		group_by(treatment, trait) %>% 
		summarise_each(funs(mean, std.error), change) %>% 
		mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*"))) 
	
	
	cols_no_anc  <- c(tableau_color_pal(palette = "Summer", direction = -1)(7))	
	all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>%
		select(ancestor_id, treatment, population, contains("change")) %>%
		mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
		mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
		mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
		# mutate(`ks_n` = scale(change_n_ks, center = FALSE)) %>% 
		# mutate(`ks_p` = scale(change_p_ks, center = FALSE)) %>% 
		mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		mutate(trait = as.factor(trait)) %>% 
		mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*"))) %>% 
		ggplot(aes(x = treatment, y = change, color = treatment)) + geom_jitter(width = 0.2, size =2, alpha = 0.5) +
		xlab("") + ylab("Standardized trait change") + ylim(-4.3, 4.3) +
		geom_hline(yintercept = 0) +
		scale_color_manual(values=cols_no_anc, name = "Selection treatment") +
		facet_wrap( ~ trait) +
		theme(strip.background = element_blank(),
		strip.text.y = element_blank(),
		strip.text.x = element_text(size=20, family = "Arial Rounded MT Bold", color = "black")) +
		geom_point(aes(x = treatment, y = mean), data = all_changes_sum2, color = "black", shape = 1, size = 3.9) +
		geom_pointrange(aes(x = treatment, y = mean, ymin = mean + std.error, ymax = mean - std.error, color = treatment), data = all_changes_sum2, size = 0.75) +
		theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 1.5)  
	ggsave("figures/standardized-trait-changes-tableau.pdf", width = 8, height = 6)
	
## paper plot	
	
	tableau_cols <- tableau_color_pal(palette = "Summer", direction = -1)(7)
	
	
	cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")	
	all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>%
		select(ancestor_id, treatment, population, contains("change")) %>%
		mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
		mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
		mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
		mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		mutate(trait = as.factor(trait)) %>% 
		mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*"))) %>% 
		ggplot(aes(x = treatment, y = change, color = treatment)) + geom_jitter(width = 0.2, size =2, alpha = 0.5) +
		xlab("") + ylab("Standardized trait change relative to ancestor") + ylim(-4.3, 4.3) +
		geom_hline(yintercept = 0) +
		scale_color_manual(values=cols_no_anc, name = "Selection treatment") +
		facet_wrap( ~ trait) +
		theme(strip.background = element_blank(),
			  strip.text.y = element_blank(),
			  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
		geom_point(aes(x = treatment, y = mean), data = all_changes_sum2, color = "black", shape = 1, size = 3.9) +
		geom_pointrange(aes(x = treatment, y = mean, ymin = mean + std.error*1.96, ymax = mean - std.error*1.96, color = treatment), data = all_changes_sum2, size = 0.75) +
		theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=14, family = "Arial", color = "black"),
			axis.title=element_text(size=14, family = "Arial", color = "black"),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 1) 
	ggsave("figures/standardized-trait-changes-tableau.png", width = 8, height = 6)
	
	all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>%
		select(ancestor_id, treatment, population, contains("change")) %>%
		mutate(`N*` = change_n_star) %>% 
		mutate(`P*` = change_p_star) %>% 
		mutate(`I*` = change_i_star) %>% 
		mutate(`Salt tolerance` = change_salt_tol) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		mutate(trait = as.factor(trait)) %>% 
		mutate(trait = factor(trait,
							  levels=c("Salt tolerance", "P*", "N*",
							  		 "I*"))) %>% 
		ggplot(aes(x = treatment, y = change, color = treatment)) + geom_jitter(width = 0.2, size =2, alpha = 0.5) +
		xlab("") + ylab("Trait change relative to ancestor") +
		geom_hline(yintercept = 0) +
		scale_color_manual(values=cols_no_anc, name = "Selection treatment") +
		facet_wrap( ~ trait, scales = "free") +
		theme(strip.background = element_blank(),
			  strip.text.y = element_blank(),
			  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
		geom_point(aes(x = treatment, y = mean), data = all_changes_sum3, color = "black", shape = 1, size = 3.9) +
		geom_pointrange(aes(x = treatment, y = mean, ymin = mean + std.error, ymax = mean - std.error, color = treatment),
						data = all_changes_sum3, size = 0.75) +
		theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=14, family = "Arial", color = "black"),
			axis.title=element_text(size=14, family = "Arial", color = "black"),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 1) 
	ggsave("figures/unstandardized-trait-changes-tableau-se-error.png", width = 8, height = 6)
	
	
	all_changes %>% 
		filter(treatment != "none") %>% 
		mutate(treatment = factor(treatment,
								  levels=c("C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>%
		select(ancestor_id, treatment, population, contains("change")) %>%
		mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
		mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
		mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
		mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
		gather(key = trait, value = change, 13:16) %>% 
		mutate(change_direction = ifelse(change > 0, "positive", "negative")) %>% View
		group_by(trait, treatment) %>% 
		count(change_direction) %>% View
	
all_changes %>% 
	filter(treatment != "none") %>% 
	select(ancestor_id, treatment, population, contains("change")) %>% 
	mutate(change_n_star_scale = scale(change_n_star)) %>% 
	mutate(change_p_star_scale = scale(change_p_star)) %>% 
	mutate(change_i_star_scale = scale(change_i_star)) %>% 
	mutate(change_salt_tol_scale = scale(change_salt_tol)) %>% 
	gather(key = trait, value = change, 11:14) %>% 
	group_by(treatment, trait) %>% 
	summarise_each(funs(mean, std.error), change) %>% 
	ggplot(aes(x = trait, treatment, y = mean, color = trait)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.3, alpha = 0.5) +
	xlab("Selection treatment") + ylab("Standardized trait change") +
	geom_hline(yintercept =  0)




### is there a trade-off between the ks for each nutrient?

all_changes %>% 
	ggplot(aes(x = n_ks, y = p_ks)) + geom_point() +
	geom_smooth(method = 'lm') +
	ylab("Half saturation constant (um P)") +
	xlab("Half saturation constant (um N)")
ggsave("figures/ks-trade-off.png", width = 8, height = 6)




all_changes_sum <- all_changes %>% 
	filter(treatment != "none") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	select(ancestor_id, treatment, population, contains("change")) %>% 
	mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
	mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
	mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
	mutate(`ks_n` = scale(change_n_ks, center = FALSE)) %>% 
	mutate(`ks_p` = scale(change_p_ks, center = FALSE)) %>% 
	mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
	gather(key = trait, value = change, 13:18) %>% 
	group_by(treatment, trait) %>% 
	summarise_each(funs(mean, std.error), change) %>% 
	mutate(trait = factor(trait,
						  levels=c("Salt tolerance", "P*", "N*",
						  		 "I*", "ks_n", "ks_p"))) 


all_changes %>% 
	filter(treatment != "none") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	select(ancestor_id, treatment, population, contains("change")) %>% 
	mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
	mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
	mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
	mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
	ggplot(aes(x = `I*`, y = `Salt tolerance`, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=2.3, ymax=4.1, ymin=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymax=0, ymin=-4.1),
			  color="transparent", alpha=0.5, fill = "grey") +
	# geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	theme_classic() + 
	# theme(legend.position="none") +
	# tableau_color_pal(palette = "Tableau 10") +
	geom_point(aes(x = `I*_mean`, y = `Salt tolerance_mean`, color = treatment), data = salt_changes_sum, shape = 18, size = 8) +
	geom_point(aes(x = `I*_mean`, y = `Salt tolerance_mean`), data = salt_changes_sum, shape = 5, size = 5, color = 'black') +
	
	# geom_point(aes(x = change_i_star_mean, y = change_n_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	# geom_point(aes(x = change_i_star_mean, y = change_n_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	geom_point(size = 3) +
	geom_point(size = 3, shape = 1, color = "black") +
	geom_pointrange(aes(x = `I*_mean`, y = `Salt tolerance_mean`,  ymin = `Salt tolerance_mean` - `Salt tolerance_std.error`, ymax = `Salt tolerance_mean` + `Salt tolerance_std.error`, color = treatment), data = salt_changes_sum) +
	geom_segment(aes(x = `I*_mean` - `I*_std.error`, xend = `I*_mean` + `I*_std.error`,  y = `Salt tolerance_mean`, yend = `Salt tolerance_mean`, color = treatment), data = salt_changes_sum) +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = salt_changes_sum) +
	scale_color_manual(values=cols2, name = "Selection treatment") +
	ylab("Change in salt tolerance (g/L)") + xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')')) +
	ylim(-4.1, 4.1) + xlim(-2.3, 2.3) 
ggsave("figures/salt-tol-star-i-star-trade-off-T-color-treatment.png", width = 6.5, height = 5)



salt_changes_sum <- all_changes %>% 
	filter(treatment != "none") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	select(ancestor_id, treatment, population, contains("change")) %>% 
	mutate(`N*` = scale(change_n_star, center = FALSE)) %>% 
	mutate(`P*` = scale(change_p_star, center = FALSE)) %>% 
	mutate(`I*` = scale(change_i_star, center = FALSE)) %>% 
	mutate(`Salt tolerance` = scale(change_salt_tol, center = FALSE)) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), `N*`, `P*`, `I*`, `Salt tolerance`) 
