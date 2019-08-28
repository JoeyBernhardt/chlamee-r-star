
library(tidyverse)
library(cowplot)

# all_rstars <- read_csv("data-processed/all-rstars.csv") ## this has the old stoich data
all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv")
combos <- read_csv("data-processed/chlamee-unique-combos2.csv") %>% 
	gather(key = competitor, value = identity, comp1, comp2)

mylims.x.Chesson <-  1
mylims.y.Chesson  <- 4





### Tilman to Chesson mechanics

## write a function to loop through all the different combinations of pairs and resource conditions
combos <- combos_split[[18]]

find_alphas <- function(combos){
		snippet <- all_rstars %>% 
	filter(ancestor_id %in% c(combos$identity)) %>% 
			filter(treatment == c(combos$treatment[[1]]))

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

SN <- 1000/300
SP <- 140/350

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
comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1])
alphas_calc2 <- bind_cols(alphas_calc, r_s, comps)
return(alphas_calc2)
}




combos_split <- combos %>% 
	split(.$combination)


all_alphas_ancestors <- combos_split %>% 
	map_df(find_alphas, .id = "combination")

all_alphas <- all_alphas_ancestors %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho)


# side bar to get plots ---------------------------------------------------

x.1 = seq(1, 0, by=-0.001)
x.2 = seq(0, -1, by=-0.001)
y1.1 = 1/((1-x.1))
y2.1 = 1-x.1
y1.2 = 1/((1-x.2))
y2.2 = 1-x.2

cf.df.1 = data.frame(my.x = c(x.1, x.1), 
					 my.y = c(y1.1, y2.1), 
					 whichy = c(rep("y1", times = length(seq(1, 0, by=-0.001))), 
					 		   rep("y2", times = length(seq(1, 0, by=-0.001)))))
rib.dims.1 = data.frame(min.dim = y1.1, max.dim = y2.1, x.dim = x.1)

cf.df.2 = data.frame(my.x = c(x.2, x.2), 
					 my.y = c(y1.2, y2.2), 
					 whichy = c(rep("y1", times = length(seq(0, -1, by=-0.001))), 
					 		   rep("y2", times = length(seq(0, -1, by=-0.001)))))
rib.dims.2 = data.frame(min.dim = y1.2, max.dim = y2.2, x.dim = x.2)



ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) +
	geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = treatment), size = 3, data = all_alphas) +
	geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = all_alphas, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 6.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 6.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 6.0, 
			  fontface = "bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 6.0, 
			  fontface = "bold") +

	scale_linetype_manual(values = c("y2" = "solid", 
									 "y1" = "dotted")) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		  axis.title.y = element_text(angle = 90), 
		  axis.title = element_text(size = 14)) 
		
ggsave("figures/Chesson-plot-all.png", width = 7, height = 5)

### ok so now we want to generate the combination lists for all the pairwise combinations in different
### treatments from each ancestor


library(readxl)
library(here)
library(janitor)
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")

anc2 <- treatments %>% 
	filter(ancestor_id == "anc2")

anc3 <- treatments %>% 
	filter(ancestor_id == "anc3")
anc4 <- treatments %>% 
	filter(ancestor_id == "anc4")
anc5 <- treatments %>% 
	filter(ancestor_id == "anc5")
cc1690 <- treatments %>% 
	filter(ancestor_id == "cc1690")


anc2_combos <- combn(anc2$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_id = "anc2")

anc3_combos <- combn(anc3$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_id = "anc3")

anc4_combos <- combn(anc4$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_id = "anc4")

anc5_combos <- combn(anc5$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_id = "anc5")

cc1690_combos <- combn(cc1690$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_id = "cc1690")

all_combos <- bind_rows(anc2_combos, anc3_combos, anc4_combos, anc5_combos, cc1690_combos) %>% 
	mutate(combination = rownames(.)) %>% 
	split (.$combination)


all_combos_all <- bind_rows(anc2_combos, anc3_combos, anc4_combos, anc5_combos, cc1690_combos) %>% 
	mutate(combination = rownames(.)) 

write_csv(all_combos_all, "data-processed/all_combos_allopatry.csv")

### now the 'after evolution combos'

anc2_3 <- treatments %>% 
	filter(ancestor_id %in% c("anc2", "anc3"))

anc2_3_combos <- combn(anc2_3$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc2_3")

anc2_4 <- treatments %>% 
	filter(ancestor_id %in% c("anc2", "anc4"))

anc2_4_combos <- combn(anc2_4$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc2_4")

anc2_5 <- treatments %>% 
	filter(ancestor_id %in% c("anc2", "anc5"))

anc2_5_combos <- combn(anc2_5$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc2_5")

anc2_cc1690 <- treatments %>% 
	filter(ancestor_id %in% c("anc2", "cc1690"))

anc2_cc1690_combos <- combn(anc2_cc1690$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc2_cc1690")

anc3_4 <- treatments %>% 
	filter(ancestor_id %in% c("anc3", "anc4"))

anc3_4_combos <- combn(anc3_4$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc3_4")

anc3_5 <- treatments %>% 
	filter(ancestor_id %in% c("anc3", "anc5"))

anc3_5_combos <- combn(anc3_5$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc3_5")

anc3_cc1690 <- treatments %>% 
	filter(ancestor_id %in% c("anc3", "cc1690"))

anc3_cc1690_combos <- combn(anc3_cc1690$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc3_cc1690")

anc4_5 <- treatments %>% 
	filter(ancestor_id %in% c("anc4", "anc5"))

anc4_5_combos <- combn(anc4_5$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc4_5")

anc4_cc1690 <- treatments %>% 
	filter(ancestor_id %in% c("anc4", "cc1690"))

anc4_cc1690_combos <- combn(anc4_cc1690$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc4_cc1690")

anc5_cc1690 <- treatments %>% 
	filter(ancestor_id %in% c("anc5", "cc1690"))

anc5_cc1690_combos <- combn(anc5_cc1690$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(ancestor_combo = "anc5_cc1690")


all_before_after_combos <- bind_rows(anc5_cc1690_combos, anc4_cc1690_combos, anc4_5_combos,
									 anc3_cc1690_combos, anc3_5_combos, anc3_4_combos, anc2_cc1690_combos,
									 anc2_5_combos, anc2_4_combos, anc2_3_combos)

write_csv(all_before_after_combos, "data-processed/all-before-after-combos.csv")


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

ns <- c(seq(1*(min(all_rstars$n_star) - 0.01), 1*(max(all_rstars$n_star) + 0.01), length.out = 4), 
		seq(2*(min(all_rstars$n_star) - 0.01), 2*(max(all_rstars$n_star) + 0.01), length.out = 4),
		seq(3*(min(all_rstars$n_star) - 0.01), 3*(max(all_rstars$n_star) + 0.01), length.out = 4))
ps <- c(seq(1*(min(all_rstars$p_star) - 0.001), 1*(max(all_rstars$p_star) + 0.001), length.out = 4),
		seq(2*(min(all_rstars$p_star) - 0.001), 2*(max(all_rstars$p_star) + 0.001), length.out = 4),
		seq(3*(min(all_rstars$p_star) - 0.001), 3*(max(all_rstars$p_star) + 0.001), length.out = 4))

write_csv(data.frame(n_supply = ns), "data-processed/n-supplies.csv")
write_csv(data.frame(p_supply = ps), "data-processed/p-supplies.csv")

summary(ns)

### now try for many supply points

all_combos_allop <- read_csv("data-processed/all_combos_allopatry.csv") %>% 
	split(.$combination)

results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(all_combos_allop, find_alphas_combos, SN1 = i, SP1 = j)
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results3 <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(location = "evolved in allopatry") %>% 
	mutate(supply_ratio = n_supply/p_supply)

write_csv(results3, "data-processed/chesson-params-allopatry.csv")


results3 %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	# top_n(n = 1, wt = n_supply) %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 3)) %>% 
	filter(combination == "10_1") %>% 
	filter(supply_ratio_round == 10.681) %>% 
ggplot(aes(x = supply_ratio, y = stabil_potential)) + geom_point() +
	scale_color_viridis_c()

### ok now plot this on the ZNGI plot to see if it looks right

results3 %>% 
	group_by(coexist) %>% tally()

snip <- results3 %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 3)) %>% 
	filter(combination == "10_1") %>% 
	filter(supply_ratio_round == 10.681) 

all_rstars %>% 
	filter(population %in% c(10, 1)) %>% 
	ggplot(aes(x = n_star, y = p_star)) + geom_point(color = "grey") +
	geom_segment(aes(x = max(n_star), y = max(p_star), xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, color = "grey") +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = nitrogen_supply, y = phosphorus_supply, color = zone), data = snip) +
	# geom_point(aes(x = 1000/500, y = 50/500, color = zone), color = "pink", size = 3) +
	# geom_point(aes(x = 10/500, y = 50/500, color = zone), color = "purple", size = 3) +
	# geom_point(aes(x = 1000/500, y = 0.5/500, color = zone), color = "orange", size = 3) +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/allopatry-example-22-6.png", width = 8, height = 6)


	ggplot()+
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = 1000/500, y = 50/500, color = zone), color = "pink", size = 3) +
	geom_point(aes(x = 10/500, y = 50/500, color = zone), color = "purple", size = 3) +
	geom_point(aes(x = 1000/500, y = 0.5/500, color = zone), color = "orange", size = 3) +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 



limitation_slopes <- all_rstars %>% 
	mutate(cons_slope = pc/nc) %>% 
	mutate(combo_slope = 50/1000) %>% 
	mutate(P_limitation_slope = 0.5/1000) %>% 
	mutate(N_limitation_slope = 50/10) %>% 
	mutate(limiting_resource = ifelse(cons_slope > combo_slope, "P", "N")) %>% 
	mutate(P_limiting_resource = ifelse(cons_slope > P_limitation_slope, "P", "N")) %>% 
	mutate(N_limiting_resource = ifelse(cons_slope > N_limitation_slope, "P", "N")) 
## ok so it likes in full COMBO, everybody is P limited


results3 <- read_csv("data-processed/chesson-params-allopatry.csv")
results3b <- results2b %>% 
	filter(supply_ratio < 50)

results3b <- results3 %>% 
	filter(supply_ratio < 50)

results2 %>% 
	ggplot(aes(x = supply_ratio)) + geom_histogram()


ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) +
	geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = supply_ratio), size = 2, data = results3b, alpha = 0.5) +
	# geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = all_alphas_results, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	scale_colour_viridis_c(name = "N:P supply ratio") +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 4.0, 
			  fontface = "bold") +
	# xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
	xlab(expression(paste("Niche difference ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		axis.title.y = element_text(angle = 90), 
		axis.title = element_text(size = 14))
ggsave("figures/Chesson-plot-all-allopatry-supplies-less50-reduced.png", width = 8, height = 5)

ggplot() +
geom_point(aes(x = supply_ratio, y = stabil_potential, color = supply_ratio), size = 2,
		   data = results3b, alpha = 0.5) + scale_color_viridis_c()

results3bb <- results3b %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 1)) %>% 
	filter(p_supply < 0.095) %>% 
	filter(n_supply < 2.1) %>% 
	group_by(supply_ratio_round) %>% 
	mutate(minp = min(p_supply)) %>% 
	mutate(minn = min(n_supply)) %>% 
	filter(p_supply == minp & n_supply == minn) %>% 
	top_n(n = 1, wt = n_supply) %>% 
	top_n(n = 1, wt = p_supply)
ggplot() +
	geom_point(aes(x = p_supply, y = n_supply, color = supply_ratio), data = results3b) +
	geom_point() +
	scale_colour_viridis_c()
ggsave("figures/supply_points.png", width = 8, height = 6)

results3b %>% 
	ggplot(aes(x = supply_ratio)) + geom_density()

results3bb %>% 
	group_by(supply_ratio_round) %>% 
	tally() %>% View

all_alphas_results %>%
	ggplot(aes(x = stabil_potential)) + geom_histogram(bins = 50) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) 
ggsave("figures/niche-differences-allopatry.png", width = 5, height = 3)



ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) +
	geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = treatment), size = 3, data = all_alphas_results) +
	geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = all_alphas_results, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 4.0, 
			  fontface = "bold") +
	
	scale_linetype_manual(values = c("y2" = "solid", 
									 "y1" = "dotted")) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		axis.title.y = element_text(angle = 90), 
		axis.title = element_text(size = 14)) 
ggsave("figures/Chesson-plot-all-allopatry-0.01.png", width = 7, height = 5)

all_alphas_results %>%
	ggplot(aes(x = stabil_potential)) + geom_histogram(bins = 50) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) 
ggsave("figures/niche-differences-allopatry.png", width = 5, height = 3)

all_alphas_results %>%
	ggplot(aes(x = 1, y = stabil_potential)) + geom_jitter(width = 0.1, alpha = 0.5) 


#### now let's get the differences within treatments -- sympatry

n_treat <- treatments %>% 
	filter(treatment == "N")

p_treat <- treatments %>% 
	filter(treatment == "P")

none_treat <- treatments %>% 
	filter(treatment == "none")

l_treat <- treatments %>% 
	filter(treatment == "L")

b_treat <- treatments %>% 
	filter(treatment == "B")

bs_treat <- treatments %>% 
	filter(treatment == "BS")

s_treat <- treatments %>% 
	filter(treatment == "S")

c_treat <- treatments %>% 
	filter(treatment == "C")

n_combos <- combn(n_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "N")
s_combos <- combn(s_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "S")

p_combos <- combn(p_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "P")
l_combos <- combn(l_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "L")

b_combos <- combn(b_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "B")

bs_combos <- combn(bs_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "BS")

c_combos <- combn(c_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "C")

none_combos <- combn(none_treat$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(treatment = "none")

all_combos_treatment <- bind_rows(n_combos, p_combos, l_combos, s_combos, bs_combos, b_combos, none_combos) %>% 
	mutate(combination = rownames(.)) %>% 
	split (.$combination)

write_csv(all_combos_treatment, "data-processed/all_combos_treatment.csv")

comp_results_treatment <- all_combos_treatment %>% 
	map_df(find_alphas_combos)


all_alphas_results_treaments <- comp_results_treatment %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(location = "evolved in sympatry")


ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) +
	geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = treatment), size = 3, data = all_alphas_results_treaments) +
	geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = all_alphas_results_treaments, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 4.0, 
			  fontface = "bold") +
	
	scale_linetype_manual(values = c("y2" = "solid", 
									 "y1" = "dotted")) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		axis.title.y = element_text(angle = 90), 
		axis.title = element_text(size = 14)) 
ggsave("figures/Chesson-plot-all-sympatry.png", width = 7, height = 5)

all_alphas_results_treaments %>%
	ggplot(aes(x = stabil_potential)) + geom_histogram(bins = 50) +
	xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) 
ggsave("figures/niche-differences-sympatry.png", width = 5, height = 3)


all_alphas_ancestors <- all_alphas_results_treaments %>%
	filter(treatment == "none") %>% 
	mutate(location = "ancestors")

all_comp_locations <- bind_rows(all_alphas_results, all_alphas_results_treaments, all_alphas_ancestors, results3) %>% 
	mutate(location = factor(location, levels = c("ancestors", "evolved in sympatry", "evolved in allopatry")))

all_comp_locations_summary <- all_comp_locations %>% 
	group_by(location) %>% 
	summarise(mean_stabil = mean(stabil_potential),
			  stabil_se = std.error(stabil_potential))

all_comp_locations %>% 
	ggplot(aes(x = location, y = stabil_potential, color = log(supply_ratio))) + geom_jitter(width = 0.1, alpha = 0.5) +
	geom_point(aes(x = location, y = mean_stabil), data= all_comp_locations_summary, color = "purple") +
	geom_errorbar(aes(ymin = mean_stabil - stabil_se, ymax = mean_stabil + stabil_se), width = 0.1, color = "purple", data = all_comp_locations_summary) +
	ylab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) +
	xlab("Location") +
	geom_hline(yintercept = 0) + scale_color_viridis_c(begin = 0.9, end = 0.1)
ggsave("figures/stabil-potential-location.png", width = 8, height = 6)

#395199

results3 %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	# top_n(n = 1, wt = n_supply) %>% 
	filter(combination == "10_1", supply_ratio > 12.5, supply_ratio < 13) %>% 
	ggplot(aes(x = supply_ratio, y = stabil_potential)) + geom_line() +
	scale_color_viridis_c()
ggsave("figures/stabil_potential_supply_ratio.png", width = 8, height = 6)

results3b %>%
	mutate(combination = paste(pop1, pop2, sep = "_")) %>% 
	group_by(supply_ratio) %>% 
	# top_n(n = 1, wt = n_supply) %>% 
	mutate(supply_ratio_round = round(supply_ratio, digits = 3)) %>% 
	filter(combination == "10_1") %>% 
	filter(supply_ratio_round == 10.681) %>% View
ggplot(aes(x = supply_ratio, y = stabil_potential)) + geom_line() +
	scale_color_viridis_c()
ggsave("figures/stabil_potential_supply_ratio.png", width = 8, height = 6)


### ok so the new plan is to find the species which cannot coexist because they don't have crossing ZNGI's and then apply the different 
## nutrient limitation lines. 

outcomes_allopatry <- read_csv("data-processed/tilman-competition-outcomes-allopatry.csv")
nested <- outcomes_allopatry %>% 
	filter(zngis_cross == FALSE)
steeper <- outcomes_allopatry %>% 
	filter(steeper_zngi == TRUE)

all_combos_allo_nested <- read_csv("data-processed/all_combos_allopatry.csv") %>% 
	filter(combination %in% c(nested$combination)) %>% 
	split(.$combination)





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





results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(all_combos_allo_nested, find_alphas_nested, SN1 = i, SP1 = j)
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results_allo_nested <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(location = "evolved in allopatry") %>% 
	mutate(supply_ratio = n_supply/p_supply)

## ok this is weird because we are getting some cases where coexistence is predicted to be possible, but it cannot be

allo_long <- outcomes_allopatry %>% 
	gather(key = competitor, value = population, pop1, pop2) %>% 
	left_join(., all_rstars)

## plot it! 
all_rstars %>% 
	filter(population %in% c(29, 27)) %>% 
	ggplot(aes(x = n_star, y = p_star)) + geom_point() +
	geom_segment(aes(color = population, x = n_star, y = p_star, xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = population, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = population, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = 0.3439753, y = 0.008908737), color ="orange") +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")  
ggsave("figures/population27-29-coexist.png", width = 8, height = 6)

## ok now let's join the two results chunks -- for the ones that are nested and the ones that are crossed.

allop_chesson <- read_csv("data-processed/chesson-params-allopatry.csv")

allop_crossed <- allop_chesson %>% 
	mutate(nested_crossed = ifelse(combination %in% c(nested$combination), "nested", "crossed")) %>% 
	filter(nested_crossed == "crossed") %>% 
	mutate(coexistence_outcome = ifelse(combination %in% c(steeper$combination), "stable coexistence", "unstable coexistence"))

allop_nested <- results_allo_nested %>% 
	mutate(coexistence_outcome = "exclusion") %>% 
	mutate(nested_crossed = "nested")

all_outcomes <- bind_rows(allop_crossed, allop_nested) 
write_csv(all_outcomes, "data-processed/all_chesson_outcomes_allopatry.csv")

all_outcomes <- read_csv("data-processed/all_chesson_outcomes_allopatry.csv")

all_outcomes %>% 
	ggplot(aes(x = nitrogen_supply, y = phosphorus_supply, color  = supply_ratio)) + geom_point() + 
	scale_color_viridis_c()

low_n_supply_ratio <- 50/10
low_p_supply_ratio <- 0.5/1000


all_outcomes2 <- all_outcomes %>% 
	filter(supply_ratio < 25, supply_ratio > 15)

ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) + geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = supply_ratio), size = 3, data = all_outcomes2) +
	# geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = all_outcomes, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	scale_color_viridis_c() +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 4.0, 
			  fontface = "bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 4.0, 
			  fontface = "bold") +
	
	# scale_linetype_manual(values = c("y2" = "solid", 
	# 								 "y1" = "dotted")) +
	xlab(expression(paste("Niche differences ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		axis.title.y = element_text(angle = 90), 
		axis.title = element_text(size = 14)) 
ggsave("figures/coexistence-outcomes-chesson-supply-point.png", width = 10, height = 6)
	

