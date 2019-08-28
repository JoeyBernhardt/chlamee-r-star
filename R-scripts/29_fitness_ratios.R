


library(tidyverse)
library(cowplot)


all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv")
combos <- read_csv("data-processed/chlamee-unique-combos2.csv") %>% 
	gather(key = competitor, value = identity, comp1, comp2)
all_combos_treatment <- read_csv("data-processed/all_combos_treatment.csv")


treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")

nppops <- treatments %>% 
	select(-note) %>% 
	filter(treatment %in% c("N", "P")) %>% 
	rename(comp1 = population, 
		   comp2 = ancestor_id) %>% 
	mutate(combination = rownames(.)) %>% 
	gather(key = competitor, value = population, comp1, comp2) 
	
write_csv(nppops, "data-processed/nppops.csv")





npops_star <- nppops %>% 
	left_join(., all_rstars, by = "population")


npops_star %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color = competitor)) + geom_point() +
	geom_segment(aes(color = competitor, x = max_n_star, y = max_p_star, xend = max_n_star + nc*9, yend = max_p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_wrap( ~ combination) + 
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")

coexist_or_not <- function(combos){
	snippet <- all_rstars %>% 
		filter(population %in% c(combos$comp1[[1]], combos$comp2[[1]]))## pull out on set of competitors
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	n_winner <- snippet$population[snippet$n_star == min(snippet$n_star)]
	p_winner <- snippet$population[snippet$p_star == min(snippet$p_star)]
	steeper_zngi <- snippet$cons_vector[snippet$population == n_winner] > snippet$cons_vector[snippet$population == p_winner]
	zngis_cross <- n_winner != p_winner
	coexist <- zngis_cross == TRUE & steeper_zngi == TRUE
	output <- data.frame(pop1 = pop1, pop2 = pop2, coexist = coexist, zngis_cross = zngis_cross, steeper_zngi = steeper_zngi)
	
}

splits <- nppops %>% 
	spread(key = competitor, value = population) %>% 
	split(.$combination) %>% 
	map_df(coexist_or_not, .id = "combination")


outcomes_all_together <- nppops %>% 
	spread(key = competitor, value = population) %>% 
	split(.$combination) %>% 
	map_df(coexist_or_not, .id = "combination") %>% 
	left_join(., nppops) %>% 
	mutate(competition_outcome = case_when(zngis_cross == TRUE & steeper_zngi == TRUE ~ "stable coexistence",
										   zngis_cross == TRUE & steeper_zngi == FALSE ~ "unstable coexistence",
										   zngis_cross == FALSE ~ "exclusion"))

nested_combos <- outcomes_all_together %>% 
	filter(zngis_cross == FALSE)

crossing_combos <- outcomes_all_together %>% 
	filter(zngis_cross == TRUE)


find_alphas_nested <- function(combos, SN1, SP1){
	snippet <- all_rstars %>% 
		filter(population %in% c(combos$comp1[[1]], combos$comp2[[1]])) 
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	c1P <- snippet$pc[snippet$population == pop1]
	c2P <- snippet$pc[snippet$population == pop2]
	c1N <- snippet$nc[snippet$population == pop1]
	c2N <- snippet$nc[snippet$population == pop2]
	R1P <- snippet$p_star[snippet$population == pop1]
	R2N <- snippet$n_star[snippet$population == pop2]
	R2P <- snippet$p_star[snippet$population == pop2]
	R1N <- snippet$n_star[snippet$population == pop1]
	D <- 0.5
	
	## how can we define the consumption vector lines?
	
	SN <- SN1
	SP <- SP1
	
	
	r1 <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
	r2 <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])
	
	cons_vec1_intercept <- R1P -(c1P/c1N)*R1N
	cons_vec2_intercept <- R2P -(c2P/c2N)*R2N
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1], nitrogen_supply = SN1,
						phosphorus_supply = SP1, zone = zone)
	alphas_calc2 <- bind_cols(alphas_calc, r_s, comps)
	return(alphas_calc2)
}

find_alphas_combos <- function(combos, SN1, SP1){
	snippet <- all_rstars %>% 
		filter(population %in% c(combos$comp1[[1]], combos$comp2[[1]])) 
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	c1P <- snippet$pc[snippet$population == pop1]
	c2P <- snippet$pc[snippet$population == pop2]
	c1N <- snippet$nc[snippet$population == pop1]
	c2N <- snippet$nc[snippet$population == pop2]
	R1P <- snippet$p_star[snippet$population == pop1]
	R2N <- snippet$n_star[snippet$population == pop2]
	R2P <- snippet$p_star[snippet$population == pop2]
	R1N <- snippet$n_star[snippet$population == pop1]
	D <- 0.5
	
	## how can we define the consumption vector lines?
	
	SN <- SN1
	SP <- SP1
	
	
	r1 <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
	r2 <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])
	
	cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
	cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1], nitrogen_supply = SN1,
						phosphorus_supply = SP1, zone = zone, combination = combos$combination[1])
	alphas_calc2 <- bind_cols(alphas_calc, r_s, comps)
	return(alphas_calc2)
}

nested_sympatry <- nppops %>% 
	spread(key = competitor, value = population) %>% 
	filter(combination %in% c(nested_combos$combination)) %>% 
	split(.$combination) 

crossing_sympatry <- nppops %>% 
	spread(key = competitor, value = population) %>% 
	filter(combination %in% c(crossing_combos$combination)) %>% 
	split(.$combination) 


### now ask whether fitness increased in the low resource environments, i.e. 
### fd / fa > 1

ns <- c(1000, 10)
ps <- c(50, 0.5)

results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(nested_sympatry, find_alphas_nested, SN1 = i, SP1 = j, .id = "combination")
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results_symp_nested <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(supply_ratio = n_supply/p_supply) %>% 
	select(-treatment)

results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(crossing_sympatry, find_alphas_combos, SN1 = i, SP1 = j, .id = "combination")
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results_symp_crossed <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(supply_ratio = n_supply/p_supply) %>% 
	select(-treatment)


all_results_symp <- bind_rows(results_symp_nested, results_symp_crossed) %>% 
	mutate(fit_ratio = ifelse(pop1 %in% c("anc2", "anc3", "anc4", "anc5", "cc1690"), 1/fit_ratio, fit_ratio)) %>% 
	mutate(nsupply = ifelse(nitrogen_supply == 1000, "high", "low")) %>% 
	mutate(psupply = ifelse(phosphorus_supply == 50, "high", "low")) %>% 
	mutate(both_supplies = paste(nsupply, psupply, sep = "")) %>% 
	filter(both_supplies %in% c("lowhigh", "highlow")) %>% 
	left_join(., nppops, by = "combination") %>% 
	filter(treatment == "N" & n_supply == 10 | treatment == "P" & p_supply == 0.5)
write_csv(all_results_symp, "data-processed/fitness_ratios_symp.csv")

all_results_symp %>% 
	ggplot(aes(x = fit_ratio, color = treatment, fill = treatment)) + geom_density() +
	facet_wrap( ~ treatment) + geom_vline(xintercept = 1) +
	xlab(expression(paste("Fitness ratio ( ", frac(italic(f[descendant]), italic(f[ancestor])), " )", sep=""))) 
	ggsave("figures/fitness-ratios-desc-anc.png", width = 8, height = 4)
