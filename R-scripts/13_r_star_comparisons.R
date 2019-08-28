### Compare Rstars and consumption vectors

library(tidyverse)
library(cowplot)
library(plotrix)
library(here)
library(readxl)
library(janitor)

treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")



p_star <- read_csv("data-processed/phosphate-rstars-direct.csv") %>% 
	rename(p_star = rstar) %>% 
	rename(p_ks = ks) %>% 
	rename(p_umax = umax)
n_star <- read_csv("data-processed/nitrate-rstars-direct.csv") %>% 
	rename(n_star = rstar) %>% 
	rename(n_ks = ks) %>% 
	rename(n_umax = umax)
i_star <- read_csv("data-processed/rstars-light-direct.csv") %>% 
	rename(i_star = rstar) 

results_light <- read_csv("data-processed/light-EP-params-direct.csv") %>% 
	filter(term == "ps") %>% 
	rename(light_umax = estimate) 

all_stars <- left_join(p_star, n_star) %>% 
	left_join(., i_star) %>% 
	left_join(., results_light)

p_star %>% 
	filter(treatment %in% c("none", "P", "N")) %>% 
	ggplot(aes(x = treatment, y = rstar)) + geom_point()

i_star %>% 
	filter(treatment %in% c("none", "L")) %>% 
	ggplot(aes(x = treatment, y = i_star, color = ancestor_id, group = ancestor_id)) + geom_point() +
	geom_line()

i_star %>% 
	# filter(treatment %in% c("none", "L")) %>% 
	ggplot(aes(x = treatment, y = i_star, color = ancestor_id, group = ancestor_id)) + geom_point() +
	geom_line()

n_star %>% 
	filter(treatment %in% c("none", "P", "N")) %>% 
	ggplot(aes(x = treatment, y = rstar)) + geom_point()

i_star %>% 
	# filter(treatment %in% c("none", "P", "N", "L")) %>% 
	ggplot(aes(x = reorder(treatment, rstar), y = rstar)) + geom_point() +
	xlab("Selection treatment") + ylab("I star")

n_star %>% 
	# filter(treatment %in% c("none", "P", "N", "L")) %>% 
	ggplot(aes(x = reorder(treatment, rstar), y = rstar)) + geom_point() +
	xlab("Selection treatment") + ylab("N star")

p_star %>% 
	# filter(treatment %in% c("none", "P", "N", "L")) %>% 
	ggplot(aes(x = reorder(treatment, rstar), y = rstar)) + geom_point() +
	xlab("Selection treatment") + ylab("P star")

p_star %>% 
	filter(treatment %in% c("none", "P")) %>% 
	ggplot(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id)) + geom_point() +
	geom_line() +
	xlab("Selection treatment") + ylab("P star")

p_star %>% 
	filter(treatment %in% c("none", "P", "C")) %>% 
	ggplot(aes(x = reorder(ancestor_id, rstar), y = rstar, color = treatment)) + geom_point() +
	geom_line() +
	xlab("Selection treatment") + ylab("P star")

n_star %>% 
	filter(treatment %in% c("none", "N", "C")) %>% 
	ggplot(aes(x = reorder(ancestor_id, rstar), y = rstar, color = treatment)) + geom_point() +
	geom_line() +
	xlab("Selection treatment") + ylab("N star")


p_star %>% 
	filter(treatment %in% c("none", "P", "C")) %>% 
	lm(rstar ~ treatment, data = .) %>% summary()




# zngi plots --------------------------------------------------------------

### now for essential resources


## define parameters

m_1 <- 0.2 
m_2 <- 0.2
k_1P <- 0.5
k_1N <- 0.7
k_2N <- 0.8
k_2P <- 0.4
r_1 <- 1.5
r_2 <- 1.7

### resource 1 is N, resource 2 is P

## R star for species 1, resource 2
R_1P <- (m_1*k_1P)/(r_1 - m_1)

## R star for species 1, resource 1
R_1N <- (m_1*k_1N)/(r_1 - m_1)

## R star for species 2, resource 1
R_2N <- (m_2*k_2N)/(r_2 - m_2)

## R star for species 2, resource 2
R_2P <- (m_2*k_2P)/(r_2 - m_2)

c1N = 2; c1P = 4; c2N = 4; c2P = 2 ### cij = per capita consumption of consumer i on resource j

mylims.x <-  3
mylims.y <-  3
ZNGI.df <-  data.frame(orange = c(5), purple = c(6))

D <- 0.25

# N_1_star <- (D / c_12) * (S_2 - R_12) - (c_22 / c_12) * N_2_star


SN <- 1000 #0.248
SP <- 50  #0.245


a11 <- c1P / (D * (SP - R_1P))
a12 <- c2P / (D * (SP - R_1P))
a21 <- c1N / (D * (SN - R_2N))
a22 <- c2N / (D * (SN - R_2N))


rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
stabil_potential <- 1 - rho #stabilizing potential
fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r_2-D)/(r_1-D) #fitness ratio 
coexist <- rho < fit_ratio &  fit_ratio < 1/rho



#' plot it!
ggplot(ZNGI.df, aes(x=orange, y=purple)) +
	# geom_point(x = S1, y = S2, size = 3) +
	coord_cartesian(expand = 0, ylim=c(0, mylims.y), xlim = c(0, mylims.x)) +
	geom_abline(slope = 1, intercept = 0, color = "grey") +
	
	## species 1 ZNGIs
	geom_segment(aes(x = R_1P, y = R_1N, xend = 1, yend = R_1N), size = 0.5, col='purple') +
	geom_segment(aes(x = R_1P, y = R_1N, xend = R_1P, yend = 1), size = 0.5, col='purple') +
	
	## species 2 ZNGIs
	geom_segment(aes(x = R_2N, y = R_2P, xend = 1, yend = R_2P), size = 0.5, col='orange') +
	geom_segment(aes(x = R_2N, y = R_2P, xend = R_2N, yend = 1), size = 0.5, col='orange') +
	
	geom_segment(aes(x = R_2N, y = R_1N, xend = R_2N + c1N, yend = R_1N + c1P),
				 size = 0.5, col='purple', linetype = 2) +
	
	geom_segment(aes(x = R_2N, y = R_1N, xend = R_2N-c1N/100, yend = R_1N-c1P/100),
				 size = 0.5, col='purple', arrow = arrow(type = "closed", length = unit(0.05, "inches"))) +
	
	geom_segment(aes(x = R_2N, y = R_1N, xend = R_2N + c2N, yend = R_1N + c2P), size = 0.5, col='orange', 
				 linetype = 2) + 
	geom_segment(aes(x = R_2N, y = R_1N, xend = R_2N-c2N/100, yend = R_1N-c2P/100), 
				 size = 0.5, col='orange', 
				 arrow = arrow(type = "closed", length = unit(0.05, "inches"))) + 
	xlab("P") + ylab("N") + 
	geom_point(aes(x = 0.1, y = 2)) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 



library(readxl)
library(janitor)
library(tidyverse)

stoich <- read_excel("data-raw/stoichiometry/ChlamEE_140218_stoich_growthrates.xlsx") %>% 
	clean_names() %>% 
	mutate(cN = b_n/b_c) %>% 
	mutate(cP = b_p/b_c) %>% 
	mutate(cNb = b_n/biomass) %>% 
	mutate(cPb = b_p/biomass) %>% 
	mutate(nc = 1/cto_n_molar) %>% 
	mutate(pc = 1/cto_p_molar) %>% 
	mutate(c_diff = abs(cN - cP)) %>% 
	mutate(ancestor = ifelse(ancestor == "CC 1690", "cc1690", ancestor)) %>% 
	select(cN, cP, c_diff, everything()) %>% 
	mutate(selection_treatment = ifelse(is.na(selection_treatment), "none", selection_treatment)) %>% 
	mutate(ancestor = ifelse(is.na(ancestor), strain_id, ancestor))


stoich_sub <- stoich %>% 
	filter(strain_id %in% c(1, 27)) %>% 
	select(strain_id, cN, cP)

c1N <- stoich_sub[1,2][[1]] 
c1P <- stoich_sub[1,3][[1]]  
c2N <- stoich_sub[2,2][[1]] 
c2P <- stoich_sub[2,3][[1]]



### looks for tradeoffs

all_stars %>% 
	ggplot(aes(x = p_star, y = n_star)) + geom_point() +
	geom_smooth(method = "lm")

all_stars %>% 
	ggplot(aes(x = p_star, y = i_star, color = ancestor_id)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = p_star, y = i_star), color = "black") +
	ylab("I* (umol/m2/s)") + xlab("P* (uM P)")
ggsave("figures/istar-pstar-tradeoff.png", width = 8, height = 6)

all_stars %>% 
	ggplot(aes(x = p_star, y = n_star, color = ancestor_id)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = p_star, y = i_star), color = "black") +
	ylab("N* (uM N)") + xlab("P* (uM P)")
ggsave("figures/nstar-pstar-tradeoff.png", width = 8, height = 6)

stoich %>% 
	ggplot(aes(x = cP, y = cN, color = ancestor)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = cP, y = cN), color = "black") +
	ylab("Biomass nitrogen per unit carbon") + xlab("Biomass phosphorus per unit carbon")
ggsave("figures/cP-cN-tradeoff.png", width = 8, height = 6)

stoich %>% 
	ggplot(aes(x = b_p, y = b_n, color = ancestor)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = b_p, y = b_n), color = "black") +
	ylab("Biomass nitrogen") + xlab("Biomass phosphorus")
ggsave("figures/bP-bN-tradeoff.png", width = 8, height = 6)



stoich %>% 
	ggplot(aes(x = cNb, y = cPb, color = ancestor)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = cNb, y = cPb), color = "black") +
	ylab("Biomass nitrogen per unit biomass") + xlab("Biomass phosphorus per unit biomass")
ggsave("figures/bP-bN-tradeoff-biomass.png", width = 8, height = 6)
	

stoich %>% 
	ggplot(aes(x = b_c, y = biomass, color = ancestor)) +
	geom_point() +
	geom_smooth(method = "lm", se = FALSE, color = "black") +
	geom_abline(slope = 1, intercept = 0, color = "red") +
	ylab("Total biomass (mg /L)")  +
	xlab("Biomass C (mg C/ L)") + xlim(0, 12) + ylim(0, 60)

stoich %>% 
	ggplot(aes(x = nc, y = pc, color = ancestor)) + geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	geom_smooth(method = "lm", se = FALSE, aes(x = nc, y = pc), color = "black") +
	ylab("P:C") + xlab("N:C")
ggsave("figures/nutrient-c-tradeoff-biomass.png", width = 8, height = 6)

cs <- stoich %>% 
	select(strain_id, ancestor, selection_treatment, nc, pc) %>% 
	rename(population = strain_id,
		   ancestor_id = ancestor,
		   treatment = selection_treatment) %>% 
	mutate(ancestor_id = case_when(ancestor_id == "Anc 2" ~ "anc2",
								   ancestor_id == "Anc 3" ~ "anc3",
								   ancestor_id == "Anc 4" ~ "anc4",
								   ancestor_id == "Anc 5" ~ "anc5",
								   ancestor_id == "cc1690" ~ "cc1690")) %>% 
	mutate(population = case_when(population == "Anc 2" ~ "anc2",
								  population == "Anc 3" ~ "anc3",
								  population == "Anc 4" ~ "anc4",
								  population == "Anc 5" ~ "anc5",
								  population == "cc1690" ~ "cc1690",
								  TRUE ~ population))


all_stars2 <- left_join(all_stars, cs) 


all_stars2 %>% 
	ggplot(aes(x = n_star, y = p_star, color = ancestor_id)) + geom_point() +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*6, yend = p_star + pc*6, color = ancestor_id),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = ancestor_id),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = 3, yend = p_star, color = ancestor_id),
				 size = 0.5, linetype = 1) +
	xlim(0, 3) + ylim(0, 1) + 
	facet_wrap( ~ ancestor_id)

all_stars2 %>% 
	# filter(treatment != "none") %>% 
	filter(treatment == "C") %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment)) + geom_point() +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = treatment),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = treatment),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = treatment),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_point(aes(x = n_star, y = p_star), color = "black", data= filter(all_stars2, treatment == "none")) +
	# xlim(0, 3) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_grid(~ ancestor_id) +
	coord_cartesian() +
	# geom_point(aes(x = 10/6, y = 50/6), color = "grey", size = 3) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/zngis-ancestor-n-supply.png", width = 12, height = 3)
ggsave("figures/zngis-ancestor-all.png", width = 10, height = 10)


# distance from supply points ---------------------------------------------

all_stars2 %>% 
	mutate(r_ratio = n_star/p_star) %>%
	mutate(lowN_s_ratio = 10/50) %>% 
	mutate(lowP_s_ratio = 1000/0.5) %>% 
	mutate(combo_s_ratio = 1000/50) %>% 
	# filter(treatment != "none") %>% 
	filter(treatment == "N") %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment)) + geom_point() +
	geom_point(aes(x = n_star, y = p_star), color = "black", data= filter(all_stars2, treatment == "none")) +
	# xlim(0, 3) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_grid(~ ancestor_id) +
	coord_cartesian() +
	geom_point(aes(x = 10, y = 50), color = "grey", size = 3) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 


all_stars2 %>% 
	filter(treatment == "none") %>% 
	mutate(r_ratio = p_star/n_star) %>% 
	mutate(lowN_s_ratio = 50/10) %>% 
	mutate(lowP_s_ratio = 0.5/1000) %>% 
	mutate(combo_s_ratio = 50/1000) %>% 
	select(contains("ratio"), ancestor_id) %>% 
	gather(key = ratio_type, value = ratio, 1:4) %>%
	ggplot(aes(x = ratio_type, y = ratio, color = ancestor_id)) +
	geom_point() +
	geom_hline(yintercept = 1) + scale_y_log10() +ylab("P:N") +
	theme(axis.text.x = element_text(angle = 90))
ggsave("figures/p-n-ratios-supplies.png", width = 6, height  = 3)

fold_changes <- all_stars2 %>% 
	filter(treatment == "none") %>% 
	mutate(r_ratio = p_star/n_star) %>% 
	mutate(lowN_s_ratio = 50/10) %>% 
	mutate(lowP_s_ratio = 0.5/1000) %>% 
	mutate(combo_s_ratio = 50/1000) %>% 
	mutate(lowN_fold_change = lowN_s_ratio/r_ratio) %>% 
	mutate(lowP_fold_change = lowP_s_ratio/r_ratio) %>% 
	mutate(COMBO_fold_change = combo_s_ratio/r_ratio) %>% 
	select(ancestor_id, contains("change"))

all_stars3 <- left_join(all_stars2, fold_changes)

write_csv(all_stars3, "data-processed/all-rstars.csv")
all_stars3 <- read_csv("data-processed/all-rstars.csv")

all_stars3 %>% 
	# filter(treatment != "none") %>% 
	filter(treatment == "N") %>% 
	ggplot(aes(x = n_star, y = p_star, color = lowN_fold_change)) + geom_point() +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = lowN_fold_change),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = lowN_fold_change),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = lowN_fold_change),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_point(aes(x = n_star, y = p_star), color = "black", data= filter(all_stars2, treatment == "none")) +
	# xlim(0, 3) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_grid(~ ancestor_id) +
	scale_color_viridis_c(begin = 0.5) +
	coord_cartesian() +
	# geom_point(aes(x = 10/6, y = 50/6), color = "grey", size = 3) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/zngis-ancestor-n-fold-change.png", width = 12, height = 3)

all_stars3 %>% 
	# filter(treatment == "P") %>% 
	ggplot(aes(x = n_star, y = p_star, color = lowP_fold_change)) + geom_point(color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, color = "grey") +
	
	
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = lowP_fold_change),
				 size = 1, linetype = 2, data = filter(all_stars3, treatment == "P")) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = lowP_fold_change),
				 size = 1, linetype = 1, data = filter(all_stars3, treatment == "P")) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = lowP_fold_change),
				 size = 1, linetype = 1, data = filter(all_stars3, treatment == "P")) +
	# geom_point(aes(x = n_star, y = p_star, color = low_P_fold_change), data = filter(all_stars3, treatment == "P")) +
	
	
	
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_point(aes(x = n_star, y = p_star), color = "black", data= filter(all_stars2, treatment == "none")) +
	# xlim(0, 3) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_grid(~ ancestor_id) +
	scale_color_viridis_c(option = "viridis", name = "Fold change in P") +
	coord_cartesian() +
	# geom_point(aes(x = 10/6, y = 50/6), color = "grey", size = 3) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/zngis-ancestor-p-fold-change.png", width = 12, height = 3)



# did they become more similar to one another? ----------------------------

all_stars3 %>% 
	# filter(treatment == "P") %>% 
	ggplot(aes(x = n_star, y = p_star)) + geom_point(color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, color = "grey") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, color = "grey") +
	
	# 
	# geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = lowP_fold_change),
	# 			 size = 1, linetype = 2, data = filter(all_stars3, treatment == "P")) +
	# geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = lowP_fold_change),
	# 			 size = 1, linetype = 1, data = filter(all_stars3, treatment == "P")) +
	# geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = lowP_fold_change),
	# 			 size = 1, linetype = 1, data = filter(all_stars3, treatment == "P")) +
	# geom_point(aes(x = n_star, y = p_star, color = low_P_fold_change), data = filter(all_stars3, treatment == "P")) +
	
	
	
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
				 size = 0.5, linetype = 2, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1, data = filter(all_stars2, treatment == "none"), color = "black") +
	geom_point(aes(x = n_star, y = p_star), color = "black", data= filter(all_stars2, treatment == "none")) +
	# xlim(0, 3) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_grid(~ treatment) +
	scale_color_viridis_c(option = "viridis", name = "Fold change in P") +
	coord_cartesian() +
	# geom_point(aes(x = 10/6, y = 50/6), color = "grey", size = 3) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/n-p-stars-by-treatment.png", width = 12, height = 4)

all_stars3 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_p = p_star_sd/p_star_mean) %>% 
	mutate(cv_n = n_star_sd/n_star_mean) %>% 
	gather(key = nutrient, value = cv, cv_n, cv_p) %>%
	ggplot(aes(x = reorder(treatment, cv), y = cv, fill = nutrient, group = treatment)) + geom_bar(stat = "identity") +
	facet_wrap( ~ nutrient, scales = "free")
ggsave("figures/cv-by-treatment-r-stars.png", width = 8, height = 4)


### let's try to test the Vasseur and Fox hypothesis -- that there should be convergent trait evolution,
## when grown in sympatry. To do this, we need to estimate the trait variation within each of the 
## treatments, and compare that to the variation across the different environments. 


## start with P

all_stars3 %>% 
	filter(treatment == "P") %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_p = p_star_sd/p_star_mean) %>% View

all_stars3 %>% 
	filter(treatment != "P") %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_p = p_star_sd/p_star_mean) %>% View

n_se <- all_stars3 %>% 
	filter(treatment == "N") %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_n = n_star_sd/n_star_mean) %>% 
	mutate(type = "n_once")

n_other_se <- all_stars3 %>% 
	filter(treatment != "N") %>% 
	filter(treatment != "none") %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_n = n_star_sd/n_star_mean) %>% 
	mutate(type = "not_n")

subn <- all_stars3 %>% 
	filter(treatment == "N") 

subn2 <- bind_rows(subn, subn, subn, subn, subn, subn) %>% 
	summarise_each(funs(mean, std.error, sd), p_star, n_star) %>% 
	mutate(cv_n = n_star_sd/n_star_mean) %>% 
	mutate(type = "n_only")

all_ns <- bind_rows(subn2, n_other_se, n_se)

## ok how to visualize this?

n_se <- all_stars3 %>% 
	filter(treatment == "N") %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star) %>% 
	mutate(cv = n_star_sd/n_star_mean) %>% 
	mutate(type = "n_once") %>% 
	mutate(se = n_star_std.error) %>% 
	mutate(sd = n_star_sd) %>% 
	mutate(var = n_star_var) 

n_other_se <- all_stars3 %>% 
	filter(treatment != "N") %>% 
	filter(treatment != "none") %>% 
	# sample_n(size = 5, replace = FALSE) %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star) %>% 
	mutate(cv = n_star_sd/n_star_mean) %>% 
	mutate(type = "not_n") %>% 
	mutate(se = n_star_std.error) %>% 
	mutate(sd = n_star_sd) %>% 
	mutate(var = n_star_var) 

n_comparisons <- bind_rows(n_se, n_other_se) %>% 
	mutate(resource = "N*") 

p_se <- all_stars3 %>% 
	filter(treatment == "P") %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star) %>% 
	mutate(cv = p_star_sd/p_star_mean) %>% 
	mutate(type = "p_once") %>% 
	mutate(se = p_star_std.error) %>% 
	mutate(sd = p_star_sd) %>% 
	mutate(var = p_star_var) 

p_other_se <- all_stars3 %>% 
	filter(treatment != "P") %>% 
	filter(treatment != "none") %>% 
	# sample_n(size = 5, replace = FALSE) %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star) %>% 
	mutate(cv = p_star_sd/p_star_mean) %>% 
	mutate(type = "not_p") %>% 
	mutate(se = p_star_std.error) %>% 
	mutate(sd = p_star_sd) %>% 
	mutate(var = p_star_var) 

p_comparisons <- bind_rows(p_se, p_other_se) %>% 
	mutate(resource = "P*")

i_se <- all_stars3 %>% 
	filter(treatment == "L") %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star, i_star) %>% 
	mutate(cv = i_star_sd/i_star_mean) %>% 
	mutate(type = "i_once") %>% 
	mutate(se = i_star_std.error) %>% 
	mutate(sd = i_star_sd) %>% 
	mutate(var = i_star_var) 

i_other_se <- all_stars3 %>% 
	filter(treatment != "L") %>% 
	filter(treatment != "none") %>% 
	# sample_n(size = 5, replace = FALSE) %>% 
	summarise_each(funs(mean, std.error, sd, var), p_star, n_star, i_star) %>% 
	mutate(cv = i_star_sd/i_star_mean) %>% 
	mutate(se = i_star_std.error) %>% 
	mutate(sd = i_star_sd) %>% 
	mutate(type = "not_i") %>% 
	mutate(var = i_star_var) 

i_comparisons <- bind_rows(i_se, i_other_se) %>% 
	mutate(resource = "I*")

## now let's add in the salt

salt_tolerances <- read.csv("data-processed/salt-tolerances-bootstrapped.csv") %>% 
	rename(population = population_index) %>% 
	group_by(population) %>% 
	summarise(salt_tolerance = mean(c)) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = ifelse(treatment == "A", "none", treatment)) 



all_tol_traits <- left_join(all_stars3, salt_tolerances)

s_se <- salt_tolerances %>% 
	filter(treatment == "S") %>% 
	summarise_each(funs(mean, std.error, sd, var), salt_tolerance) %>% 
	mutate(cv = sd/mean) %>% 
	mutate(type = "s_once") %>% 
	mutate(se = std.error) %>% 
	mutate(sd = sd) %>% 
	mutate(var = var) 

s_other_se <- salt_tolerances %>% 
	filter(treatment != "S") %>% 
	filter(treatment != "none") %>% 
	# sample_n(size = 5, replace = FALSE) %>% 
	summarise_each(funs(mean, std.error, sd, var), salt_tolerance) %>% 
	mutate(cv = sd/mean) %>% 
	mutate(type = "not_s") %>% 
	mutate(se = std.error) %>% 
	mutate(sd = sd) 

s_comparisons <- bind_rows(s_se, s_other_se) %>% 
	mutate(resource = "Salt tolerance")


all_comps <- bind_rows(n_comparisons, p_comparisons, i_comparisons, s_comparisons) %>% 
	mutate(location = ifelse(grepl("once", type), "Resource limited or high salt", "All other environments"))

all_comps %>% 
	ggplot(aes(x = resource, y = cv, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("CV of minimum resource requirement")
ggsave("figures/cv-r-stars-with-salt.png", width = 8, height = 6)

all_comps %>% 
	ggplot(aes(x = resource, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("Sd of minimum resource requirement")
ggsave("figures/cv-r-stars-with-salt.png", width = 8, height = 6)

all_comps %>% 
	ggplot(aes(x = location, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("SD of trait value") +
	facet_wrap( ~ resource, scales = "free", ncol = 4) +
	theme(axis.title.x=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank())
ggsave("figures/sd-r-stars-with-salt.png", width = 10, height = 4)

### new version of the graph with axis labels etc


plotP <- all_comps %>% 
	filter(resource == "P*") %>% 
	mutate(location = ifelse(location == "Resource limited or high salt", "P limited", "All other environments")) %>% 
	mutate(location = factor(location, levels = c("P limited", "All other environments"))) %>% 
	ggplot(aes(x = location, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("SD of trait value") +
	# facet_wrap( ~ resource, scales = "free", ncol = 4) + 
	theme(legend.position = "none") + xlab("") + ylab("SD of P*")

plotN <- all_comps %>% 
	filter(resource == "N*") %>% 
	mutate(location = ifelse(location == "Resource limited or high salt", "N limited", "All other environments")) %>% 
	mutate(location = factor(location, levels = c("N limited", "All other environments"))) %>% 
	ggplot(aes(x = location, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("SD of trait value") +
	# facet_wrap( ~ resource, scales = "free", ncol = 4) + 
	theme(legend.position = "none") + xlab("") + ylab("SD of N*")

plotI <- all_comps %>% 
	filter(resource == "I*") %>% 
	mutate(location = ifelse(location == "Resource limited or high salt", "Light limited", "All other environments")) %>% 
	mutate(location = factor(location, levels = c("Light limited", "All other environments"))) %>% 
	ggplot(aes(x = location, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("SD of trait value") +
	# facet_wrap( ~ resource, scales = "free", ncol = 4) + 
	theme(legend.position = "none") + xlab("") + ylab("SD of I*")

plotS <- all_comps %>% 
	filter(resource == "Salt tolerance") %>% 
	mutate(location = ifelse(location == "Resource limited or high salt", "High salt", "All other environments")) %>% 
	mutate(location = factor(location, levels = c("High salt", "All other environments"))) %>% 
	ggplot(aes(x = location, y = sd, fill = location)) + geom_bar(stat = "identity", position = "dodge") + 
	scale_fill_viridis_d(begin = 0.5, end = 0.9) + ylab("SD of trait value") +
	# facet_wrap( ~ resource, scales = "free", ncol = 4) + 
	theme(legend.position = "none") + xlab("") + ylab("SD of salt tolerance")

multi_plotS <- plot_grid(plotP, plotN, plotI, plotS, labels = c("A", "B", "C", "D"), align = "h", nrow = 1)
save_plot("figures/sd-traits.png", multi_plotS,
		  ncol = 4, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1
)



# Plot stars vs each other ------------------------------------------------

## no apparent trade-offs between the Rstars
library(beyonce)

plota <- all_stars3 %>% 
	ggplot(aes(x = n_star, y = i_star, color = ancestor_id)) + geom_point(size = 3) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	geom_smooth(method = "lm", color = "black") + ylab("I* (umol/m2/s2)") + xlab("N* (uM N)")
ggsave("figures/istar-vs-nstar.png", width = 6, height = 4)

plotb <- all_stars3 %>% 
	ggplot(aes(x = p_star, y = i_star, color = ancestor_id)) + geom_point(size = 3) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	geom_smooth(method = "lm", color = "black") + ylab("I* (umol/m2/s2)") + xlab("P* (uM P)")
ggsave("figures/istar-vs-pstar.png", width = 6, height = 4)

plotc <- all_stars3 %>% 
	ggplot(aes(x = n_star, y = p_star, color = ancestor_id)) + geom_point(size = 3) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	geom_smooth(method = "lm", color = "black") + ylab("P* (uM P)") + xlab("N* (uM N)")
ggsave("figures/pstar-vs-nstar.png", width = 6, height = 4)

legend <- get_legend(plota)

multi_plot_stars <- plot_grid(plota + theme(legend.position="none"), plotb + theme(legend.position="none"), plotc + theme(legend.position="none"), legend, labels = c("A", "B", "C"), align = "h", nrow = 1)
save_plot("figures/rstars-relative-to-eachother.png", multi_plot_stars,
		  ncol = 3, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.5
)


## anc3, 4
## anc 4, cc1690
	
combos <- read_csv("data-processed/chlamee-unique-combos2.csv") %>% 
	gather(key = competitor, value = identity, comp1, comp2)
## combo5, 9



cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
supply_vec <- SP/SN

cons_vec1_fun <- function(x){
	y <- (c1P/c1N)*x + cons_vec1_intercept
}

cons_vec2_fun <- function(x){
	y <- (c2P/c2N)*x + cons_vec2_intercept
}

supply_vec_fun <- function(x){
	y <- (SP/SN)*x
}

SN <- 1000/300
SP <- 140/350

all_stars3 %>% 
	filter(treatment %in% c("none", "P"), ancestor_id %in% combos$identity[combos$combination==5]) %>% 
	mutate(treatment = ifelse(treatment == "none", "ancestors", "low P selected")) %>% 
	ggplot(aes(x = n_star, y = p_star, color = ancestor_id)) + geom_point() +
	geom_segment(aes(linetype = treatment, x = max(n_star), y = max(p_star), xend = max(n_star) + nc*7, yend = max(p_star) + pc*7, color = ancestor_id),
				 size = 0.5) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = ancestor_id, linetype = treatment),
				 size = 0.5) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = ancestor_id, linetype = treatment),
				 size = 0.5) +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = max(n_star), y = max(p_star)), color = "black") +
	# geom_abline(slope = cons_vec1, intercept = cons_vec1_intercept, color = "purple") +
	# geom_abline(slope = cons_vec2, intercept = cons_vec2_intercept, color = "green") +
	# geom_abline(slope = supply_vec, intercept = 0) +
	# facet_grid(~ ancestor_id) +
	scale_color_viridis_d(option = "viridis", end = 0.5) +
	coord_cartesian(xlim = c(0, 4), ylim = c(0,1)) +
	# geom_abline(slope = 140/1000, intercept = 0, color = "grey") +
	geom_point(aes(x = SN, y = SP), color = "grey", size = 3, shape = 8) +
	panel_border(colour = "black") 
ggsave("figures/competition-combin5-corrSP-lowP.png", width = 8, height =6)

## now take this combination and translate to niche and fitness differences


# translate to mct terms --------------------------------------------------
## species i, resource j (resource 1 = N, resource 2 = P)
## species 1 is limited by R2, species 2 is limited by R1
## species cc1690 is limited more by P, species anc4 is more limited by N

snippet <- all_stars3 %>% 
	filter(treatment %in% c("none", "P"), ancestor_id %in% combos$identity[combos$combination==9]) %>% 
	filter(treatment == "none")
pop1 <- snippet$population[[2]]
pop2 <- snippet$population[[1]]

c1P <- snippet$pc[snippet$population == pop1]
c2P <- snippet$pc[snippet$population == pop2]
c1N <- snippet$nc[snippet$population == pop1]
c2N <- snippet$nc[snippet$population == pop2]
R1P <- snippet$p_star[snippet$population == pop1]
R2N <- snippet$n_star[snippet$population == pop2]
R2P <- snippet$p_star[snippet$population == pop2]
R1N <- snippet$n_star[snippet$population == pop1]
D <- 0.5
SN <- 1000/300
SP <- 140/300

a11 <- c1P / (D * (SP - R1P))
a12 <- c2P / (D * (SP - R1P))
a21 <- c1N / (D * (SN - R2N))
a22 <- c2N / (D * (SN - R2N))
r1 <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
r2 <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])

rho_anc <- rho
fit_ratio_anc <- fit_ratio

rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
stabil_potential <- 1 - rho #stabilizing potential
fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D) #fitness ratio
coexist <- rho < fit_ratio &  fit_ratio < 1/rho

ratio_evolved <- data.frame(rho = rho, fit_ratio = fit_ratio, treatment = "low P selected")
ratio_anc <- data.frame(rho = rho, fit_ratio = fit_ratio, treatment = "ancestors")
all_ratios <- bind_rows(ratio_evolved, ratio_anc)

mylims.x.Chesson = 1
mylims.y.Chesson = 4
ggplot() + geom_point(aes(x = rho, y = fit_ratio, shape = treatment), color = "black", size = 2, data = all_ratios) +
	# geom_point(aes(x = rho_anc, y = fit_ratio_and), color = "black", size = 2) +
	geom_line(aes(x = seq(0, 1, by = 0.01), y = seq(0, 1, by = 0.01))) +
	geom_line(aes(x = seq(0, 1, by = 0.01), y = 1/seq(0, 1, by = 0.01))) +
	geom_ribbon(aes(x = seq(0, 1, by = 0.01),
					ymin = seq(0, 1, by = 0.01),
					ymax = 1/seq(0, 1, by = 0.01)),
				fill = "grey",
				alpha = 0.3) +
	scale_shape_manual(values = c(16, 1)) +
	geom_hline(yintercept = 1, linetype = 2, color = "grey") +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(0, mylims.x.Chesson), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	panel_border(colour = "black") +
	xlab(expression(paste("Niche overlap (", rho, ")", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep="")))
ggsave("figures/chesson-plot-combin9.png", width = 7, height = 5)
	

# figure out what zone we’re in -------------------------------------------


zone_middle <- cons_vec2_fun(3) < supply_vec_fun(3) & supply_vec_fun(3) < cons_vec1_fun(3)
zone_bottom <- cons_vec2_fun(3) > supply_vec_fun(3) & supply_vec_fun(3) < cons_vec1_fun(3)
zone_top <- cons_vec2_fun(3) < supply_vec_fun(3) & supply_vec_fun(3) > cons_vec1_fun(3)

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

alphas(zone)

## zone 1 alphas (zone_top)
## ## species i, resource j (resource 1 = N, resource 2 = P)

a11 <- c1N / (D * (SN - R1N))
a12 <- c2N / (D * (SN - R1N))
a21 <- c1N / (D * (SN - R2N))
a22 <- c2N / (D * (SN - R2N))

### zone 2 alphas (zone_middle)
a11 <- c1P / (D * (SP - R1P))
a12 <- c2P / (D * (SP - R1P))
a21 <- c1N / (D * (SN - R2N))
a22 <- c2N / (D * (SN - R2N))

### zone 3 alphas (zone_bottom)

a11 <- c1P / (D * (SP - R1P))
a12 <- c2P / (D * (SP - R1P))
a21 <- c1P / (D * (SP - R2P))
a22 <- c2P / (D * (SP - R2P))

### now we need to figure out what zone we are in
SN <- 1000/300
SP <- 140/300

## how can we define the consumption vector lines?

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

cons_vec1_fun(3)
cons_vec2_fun(3)
supply_vec_fun(3)


ggplot() +
	geom_abline(slope = cons_vec1, intercept = 0, color = "purple") +
	geom_abline(slope = cons_vec2, intercept = 0, color = "green") +
	geom_abline(slope = supply_vec, intercept = 0) 
	


zone_middle <- cons_vec2_fun(3) < supply_vec_fun(3) & supply_vec_fun(3) < cons_vec1_fun(3)
zone_bottom <- cons_vec2_fun(3) > supply_vec_fun(3) & supply_vec_fun(3) < cons_vec1_fun(3)
zone_top <- cons_vec2_fun(3) < supply_vec_fun(3) & supply_vec_fun(3) > cons_vec1_fun(3)

### now ask whether the drop in P* is related to the difference in the P:N ratio of the ancestor relative to the supply point

all_stars3 %>% 
	filter(treatment %in% c("P", "none")) %>% 
	select(ancestor_id, treatment, p_star, lowP_fold_change) %>% 
	spread(key = treatment, value = p_star, 3:4) %>% 
	mutate(change_p_star = P - none) %>% 
	ggplot(aes(x = lowP_fold_change, y = change_p_star)) + geom_point() +
	geom_smooth(method = "lm")
ggsave("figures/p-star-change.png", width = 6, height = 4)

all_stars3 %>% 
	filter(treatment %in% c("N", "none")) %>% 
	select(ancestor_id, treatment, n_star, lowN_fold_change) %>% 
	spread(key = treatment, value = n_star, 3:4) %>% 
	mutate(change_n_star = N - none) %>% 
	ggplot(aes(x = lowN_fold_change, y = change_n_star)) + geom_point() +
	geom_smooth(method = "lm")
ggsave("figures/n-star-change.png", width = 6, height = 4)



# plot the change in N* vs change in P* -----------------------------------

p_change <- all_stars3 %>% 
	filter(treatment %in% c("P", "none")) %>% 
	select(ancestor_id, treatment, p_star, lowP_fold_change) %>% 
	spread(key = treatment, value = p_star, 3:4) %>% 
	mutate(change_p_star = P - none)  %>% 
	select(ancestor_id, change_p_star, lowP_fold_change)

n_change <- all_stars3 %>% 
	filter(treatment %in% c("N", "none")) %>% 
	select(ancestor_id, treatment, n_star, lowN_fold_change) %>% 
	spread(key = treatment, value = n_star, 3:4) %>% 
	mutate(change_n_star = N - none) %>% 
	select(ancestor_id, change_n_star, lowN_fold_change)

i_change <- all_stars3 %>% 
	filter(treatment %in% c("L", "none")) %>% 
	select(ancestor_id, treatment, i_star) %>% 
	spread(key = treatment, value = i_star, 2:3) %>% 
	mutate(change_i_star = L - none) %>% 
	select(ancestor_id, change_i_star)

all_changes <- left_join(p_change, n_change, by = "ancestor_id") %>% 
	left_join(., i_change)

all_changes %>% 
	ggplot(aes(x = change_n_star, y = change_i_star)) + geom_point() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

all_changes %>% 
	ggplot(aes(x = change_p_star, y = change_n_star)) + geom_point() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_smooth(method = "lm")


salt_tolerances <- read_csv("data-processed/salt_tolerances.csv") %>% 
	mutate(treatment = ifelse(treatment == "A", "none", treatment))

all_tol_traits <- left_join(all_stars3, salt_tolerances)

ancestors <- all_tol_traits %>% 
	filter(treatment == "none") %>% 
	rename(ancestral_n_star = n_star,
		   ancestral_p_star = p_star,
		   ancestral_i_star = i_star,
		   ancestral_n_ks = n_ks,
		   ancestral_p_ks = p_ks,
		   ancestral_salt_tolerance = salt_tolerance) %>% 
	select(ancestor_id, ancestral_n_star, ancestral_p_star, ancestral_i_star, ancestral_salt_tolerance, ancestral_n_ks, ancestral_p_ks)

write_csv(ancestors, "data-processed/ancestral-traits.csv")

changes_all <- all_tol_traits %>% 
	left_join(., ancestors) %>% 
	group_by(ancestor_id) %>% 
	# mutate(ancestral_n_star = filter(n_star, treatment == "none")) %>% View
	# filter(treatment %in% c("N", "none")) %>% 
	# select(ancestor_id, treatment, n_star, lowN_fold_change) %>% 
	# spread(key = treatment, value = n_star, 3:4) %>% 
	mutate(change_n_star =  n_star - ancestral_n_star) %>% 
	mutate(change_i_star =  i_star - ancestral_i_star) %>% 
	mutate(change_p_star =  p_star - ancestral_p_star) %>% 
	mutate(change_n_ks =  n_ks - ancestral_n_ks) %>% 
	mutate(change_p_ks =  p_ks - ancestral_p_ks) %>% 
	mutate(change_salt_tol =  salt_tolerance - ancestral_salt_tolerance) ## change to add salt

write_csv(changes_all, "data-processed/changes_allb.csv")

changes_all <- read_csv("data-processed/changes_allb.csv")
changes_all %>% 
	ggplot(aes(x = change_n_star, y = change_p_star, color = ancestor_id)) + geom_point() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	# geom_smooth(method = "lm", color = "black") +
	ylab("Change in P* (um P)") + xlab("Change in N* (um N)")
ggsave("figures/change-n-star-pstar.png", width = 8, height = 6)


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

changes_all %>% 
	ggplot(aes(x = change_i_star, y = change_p_star, color = ancestor_id)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.10, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.5, xmax=0, ymin=0, ymax=0.06),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	ylim(-0.1, 0.06) +
	coord_cartesian() +
	ylab("Change in P* (um P)") + xlab("Change in I* (umol/m2/s)")
ggsave("figures/p-star-i-star-trade-off.png", width = 8, height = 6)


library(vegan)
data(dune)

changes <- changes_all %>%
	ungroup() %>% 
	select(change_n_star, change_i_star, change_p_star)

dune.pca <- rda(changes)
uscores <- data.frame(dune.pca$CA$u)
uscores1 <- inner_join(rownames_to_column(changes_all), rownames_to_column(data.frame(uscores)), type = "right", by = "rowname")
vscores <- data.frame(dune.pca$CA$v)

p1 <- ggplot(uscores1, aes(x = PC1, y = PC2, col = treatment)) + 
	geom_point(size = 3) +
	# scale_color_manual(values=cbPalette) +
	# scale_fill_manual(values=cbPalette) +
	# scale_shape_manual(values = c(21:25)) +
	theme_bw() +
	theme(strip.text.y = element_text(angle = 0))
p1



biplot(dune.pca)
changes_all %>% 
	ggplot(aes(x = change_n_star, y = change_i_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-1.8, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.5, xmax=0, ymin=0, ymax=1),
			  color="transparent", alpha=0.5, fill = "grey") + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +

	# geom_smooth(method = "lm", se = FALSE, color = "black") +
	ylab("Change in I* (umol/m2/s)") + xlab("Change in N* (um N)")
ggsave("figures/n-star-i-star-trade-off.png", width = 8, height = 6)

cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")

library(plotrix)
library(ggthemes)
changes_sum <- changes_all %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_n_star, change_i_star, change_p_star)
	

# plot1 <- 
# Qualitative color schemes by Paul Tol
tol1qualitative=c("#4477AA")
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative_black=c("black", "#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

pal(tol7qualitative)	


cols  <- c("black", tableau_color_pal(palette = "Summer", direction = -1)(7))

	changes_all %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
								  levels=c("A", "C", "L", "N",
								  		 "P", "B", "S", "BS"))) %>%
	ggplot(aes(x = change_n_star, y = change_p_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.09, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.5, xmax=0, ymin=0, ymax=0.07),
			  color="transparent", alpha=0.5, fill = "grey") + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = cols) +
		geom_point(size = 3) +
	geom_pointrange(aes(x = change_n_star_mean, y = change_p_star_mean,  ymin = change_p_star_mean - change_p_star_std.error, ymax = change_p_star_mean + change_p_star_std.error, color = treatment),
					data = changes_sum, size = 1) +
	geom_segment(aes(x = change_n_star_mean - change_n_star_std.error, xend = change_n_star_mean + change_n_star_std.error,  y = change_p_star_mean, yend = change_p_star_mean, color = treatment),
				 data = changes_sum, size = 1) +
		geom_point(size = 3, shape = 1, color = "black") +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	ylab("Change in P* (um P)") + xlab("Change in N* (um N)") + geom_point(size = 3, shape = 1, color = "black") +
		geom_point(aes(x = change_n_star_mean, y = change_p_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
		geom_point(aes(x = change_n_star_mean, y = change_p_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
		theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 1.5)  
ggsave("figures/n-star-p-star-trade-off-color-treatment-tableau.pdf", width = 6, height = 5)
ggsave("figures/n-star-p-star-trade-off-color-treatment.png", width = 8, height = 6)


changes_all %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	ggplot(aes(x = change_i_star, y = change_p_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.9, ymin=-0.09, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.9, xmax=0, ymin=0, ymax=0.07),
			  color="transparent", alpha=0.5, fill = "grey") + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = cols) +
	geom_point(size = 3) +
	geom_pointrange(aes(x = change_i_star_mean, y = change_p_star_mean,  ymin = change_p_star_mean - change_p_star_std.error, ymax = change_p_star_mean + change_p_star_std.error, color = treatment),
					data = changes_sum, size = 1) +
	geom_segment(aes(x = change_i_star_mean - change_i_star_std.error, xend = change_i_star_mean + change_i_star_std.error,  y = change_p_star_mean, yend = change_p_star_mean, color = treatment),
				 data = changes_sum, size = 1) +
	geom_point(size = 3, shape = 1, color = "black") +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	ylab("Change in P* (um P)") + xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')')) + geom_point(size = 3, shape = 1, color = "black") +
	geom_point(aes(x = change_i_star_mean, y = change_p_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	geom_point(aes(x = change_i_star_mean, y = change_p_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1.5)  
 ggsave("figures/i-star-p-star-trade-off-color-treatment-tableau.pdf", width = 6, height = 5)

 multi_plot <- plot_grid(plot1, plot2, align = "h", nrow = 1)
 save_plot("figures/all_changes_rel_ancestor_poster.pdf", multi_plot,
 		  ncol = 2, # we're saving a grid plot of 2 columns
 		  nrow = 1, # and 2 rows
 		  # each individual subplot should have an aspect ratio of 1.3
 		  base_aspect_ratio = 1.2
 )


paleta_col <- ggthemes::tableau_color_pal("Tableau 10")(10)


plot2 <- changes_all %>% 
	# filter(treatment != "none") %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	ggplot(aes(x = change_i_star, y = change_p_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.09, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.5, xmax=0, ymin=0, ymax=0.07),
			  color="transparent", alpha=0.5, fill = "grey") + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_brewer(type = "qual") +
	# tableau_color_pal(palette = "Tableau 10") +
	geom_point(aes(x = change_i_star_mean, y = change_p_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	geom_point(aes(x = change_i_star_mean, y = change_p_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	geom_point(size = 3) +
	geom_point(size = 3, shape = 1, color = "black") +
	geom_pointrange(aes(x = change_i_star_mean, y = change_p_star_mean,  ymin = change_p_star_mean - change_p_star_std.error, ymax = change_p_star_mean + change_p_star_std.error, color = treatment), data = changes_sum) +
	geom_segment(aes(x = change_i_star_mean - change_i_star_std.error, xend = change_i_star_mean + change_i_star_std.error,  y = change_p_star_mean, yend = change_p_star_mean, color = treatment), data = changes_sum) +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	ylab("Change in P* (um P)") + xlab("Change in I* (umol/m2/s)")
ggsave("figures/i-star-p-star-trade-off-color-treatment.png", width = 8, height = 6)


plot3 <- changes_all %>% 
	# filter(treatment != "none") %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	ggplot(aes(x = change_i_star, y = change_n_star, color = treatment)) +
	geom_rect(aes(xmin=0, xmax=1.7, ymin=-0.9, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-1.9, xmax=0, ymin=0, ymax=0.9),
			  color="transparent", alpha=0.5, fill = "grey") + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_brewer(type = "qual") +
	# tableau_color_pal(palette = "Tableau 10") +
	geom_point(aes(x = change_i_star_mean, y = change_n_star_mean, color = treatment), data = changes_sum, shape = 18, size = 8) +
	geom_point(aes(x = change_i_star_mean, y = change_n_star_mean), data = changes_sum, shape = 5, size = 5, color = "black") +
	geom_point(size = 3) +
	geom_point(size = 3, shape = 1, color = "black") +
	geom_pointrange(aes(x = change_i_star_mean, y = change_n_star_mean,  ymin = change_n_star_mean - change_n_star_std.error, ymax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	geom_segment(aes(x = change_i_star_mean - change_i_star_std.error, xend = change_i_star_mean + change_i_star_std.error,  y = change_n_star_mean, yend = change_n_star_mean, color = treatment), data = changes_sum) +
	# geom_errorbarh(aes(xmin = change_n_star_mean - change_n_star_std.error, xmax = change_n_star_mean + change_n_star_std.error, color = treatment), data = changes_sum) +
	ylab("Change in N* (um N)") + xlab("Change in I* (umol/m2/s)")
ggsave("figures/i-star-n-star-trade-off-color-treatment.png", width = 8, height = 6)


multi_plot <- plot_grid(plot1, plot2, plot3, labels = c("A", "B", "C"), align = "h", nrow = 1)
save_plot("figures/all_changes_rel_ancestor.png", multi_plot,
		  ncol = 3, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.5
)

all_stars3 %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment)) + geom_point() +
	ylab("P* (uM P)") + xlab("N* (uM N)")
ggsave("figures/n-star-pstar.png", width = 8, height = 6)



# calculate selection differentials ---------------------------------------


## beta is the slope of the regression of fitness against the trait value. So in this case, 
## it'd be the slope of the value of population growth rate at the resource level.
## for phosphorus, the growth rate at the P supply point (0.5 uM P)

### which is the limiting resource?

all_stars3 %>% 
	mutate(limiting_resource = ifelse(n_star > p_star, "nitrogen", "phosphorus")) %>% View


predict_r <- function(nutrient_concentration = 0.5){
	r <- umax*(nutrient_concentration/ (ks+ nutrient_concentration))
}

all_stars3 %>% 
	mutate(growth_at_supply_p = p_umax*(0.5/(p_ks + 0.5))) %>% 
	mutate(growth_at_supply_n = n_umax*(10/(n_ks + 10))) %>% 
	filter(treatment == "none") %>% 
	ggplot(aes(x = p_star, y = growth_at_supply_p, color = ancestor_id)) + geom_point() +
	# facet_wrap( ~ ancestor_id) +
	geom_smooth(method = "lm")
ggsave("figures/selection-gradient-p.png", width = 6, height = 4)

p_monod_boot <- read.csv("data-processed/phosphate-monod-params-bootstrap.csv") %>% 
	left_join(., treatments, by = c("population_index" = "population"))

all_rstars4_anc <- all_stars3 %>% 
	mutate(growth_at_supply_p50 = p_umax*(50/(p_ks + 50))) %>%  
	mutate(growth_at_supply_p5 = p_umax*(5/(p_ks + 5))) %>% 
	mutate(growth_at_supply_p0.5 = p_umax*(0.5/(p_ks + 0.5))) %>%
	mutate(growth_at_supply_n = n_umax*(10/(n_ks + 10))) %>% 
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 18:20)



p_monod_boot %>% 
	mutate(growth_at_supply_p50 = umax*(50/(ks + 50))) %>%  
	mutate(growth_at_supply_p5 = umax*(5/(ks + 5))) %>% 
	mutate(growth_at_supply_p0.5 = umax*(0.5/(ks + 0.5))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 8:10) %>%
	ggplot(aes(x = p_star, y = growth_rate, color = ancestor_id)) + geom_point(alpha =0.5) +
	facet_wrap( ~ supply_level, scales = "free") + geom_smooth(method = "lm") +
	geom_point(aes(x = p_star, y = growth_rate), data = all_rstars4_anc, color = "black") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (/day)") + xlab("P* (um P)")
ggsave("figures/selection-gradient-p-bootstrap-supply-all-P.png", width = 8, height = 4)

p_monod_boot %>% 
	mutate(growth_at_supply_p50 = umax*(50/(ks + 50))) %>%  
	mutate(growth_at_supply_p5 = umax*(5/(ks + 5))) %>% 
	mutate(growth_at_supply_p0.5 = umax*(0.5/(ks + 0.5))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 8:10) %>%
	ggplot(aes(x = p_star, y = growth_rate, color = ancestor_id)) + geom_point(alpha =0.5) +
	facet_wrap( ~ supply_level, scales = "free") + geom_smooth(method = "lm") +
	geom_point(aes(x = p_star, y = growth_rate), data = all_rstars4_anc, color = "black") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (/day)") + xlab("P* (um P)")
ggsave("figures/selection-gradient-p-bootstrap-supply-all-P-scales.png", width = 10, height = 4)

library(broom)
p_monod_boot %>% 
	mutate(growth_at_supply_p = umax*(5/(ks + 5))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_p ~ p_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "p_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/p-selection-gradient-slopes.png", width = 6, height = 4)

p_monod_boot %>% 
	mutate(growth_at_supply_p = umax*(50/(ks + 50))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_p ~ p_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "p_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/p-selection-gradient-slopes-control.png", width = 6, height = 4)



all_stars3 %>% 
	mutate(growth_at_supply_p = p_umax*(0.5/(p_ks + 0.5))) %>% 
	mutate(growth_at_supply_n = n_umax*(10/(n_ks + 10))) %>% 
	# filter(treatment == "none") %>% 
	ggplot(aes(x = n_star, y = growth_at_supply_n, color = ancestor_id)) + geom_point()
ggsave("figures/selection-gradient-n.png", width = 6, height = 4)



#### nitrate selection gradient

n_monod_boot <- read.csv("data-processed/nitrate-monod-params-bootstrap.csv") %>% 
	left_join(., treatments, by = c("population_index" = "population"))


all_rstars5_anc <- all_stars3 %>% 
	mutate(growth_at_supply_n1000 = n_umax*(1000/(n_ks + 1000))) %>%  
	mutate(growth_at_supply_n100 = n_umax*(100/(n_ks + 100))) %>% 
	mutate(growth_at_supply_n10 = n_umax*(10/(n_ks + 10))) %>%  
	# mutate(growth_at_supply_n = n_umax*(10/(n_ks + 10))) %>% 
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 23:25)


n_monod_boot %>% 
	mutate(growth_at_supply_n1000 = umax*(1000/(ks + 1000))) %>%  
	mutate(growth_at_supply_n100 = umax*(100/(ks + 100))) %>% 
	mutate(growth_at_supply_n10 = umax*(10/(ks + 10))) %>%  
	mutate(n_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 8:10) %>% 
	ggplot(aes(x = n_star, y = growth_rate, color = ancestor_id)) + geom_point(alpha =0.5) +
	facet_wrap( ~ supply_level, scales = "free") + geom_smooth(method = "lm") +
	geom_point(aes(x = n_star, y = growth_rate), data = all_rstars5_anc, color = "black") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (/day)") + xlab("N* (um N)")
ggsave("figures/selection-gradient-n-bootstrap-supply-all-N.png", width = 10, height = 4)



n_monod_boot %>% 
	mutate(growth_at_supply_n = umax*(10/(ks + 10))) %>%  
	mutate(n_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_n ~ n_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "n_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/n-selection-gradient-slopes.png", width = 6, height = 4)

n_monod_boot %>% 
	mutate(growth_at_supply_n = umax*(1000/(ks + 1000))) %>%  
	mutate(n_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_n ~ n_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "n_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/n-selection-gradient-slopes-control.png", width = 6, height = 4)


#### light selection gradient

results_light <- read_csv("data-processed/light-EP-params-direct.csv") %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate)

light_boot <- read.csv("data-processed/light-ep-bootstrap-params.csv") %>% 
	left_join(., treatments, by = c("population_index" = "population")) %>% 
	rename(population = population_index)

istars_boot <- read_csv("data-processed/light-ep-params-bootstrap.csv") %>% 
	left_join(., light_boot)


(light/((1/(alpha * 
				eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha)))

all_rstars6_anc <- results_light %>% 
	mutate(growth_at_supply_l100 = 100/((1/(alpha * 
											  	eopt^2)) * 100^2 + (1/ps - 2/(alpha * eopt)) * 100 + (1/alpha))) %>%  
	mutate(growth_at_supply_l20 = 20/((1/(alpha * 
												eopt^2)) * 20^2 + (1/ps - 2/(alpha * eopt)) * 20 + (1/alpha))) %>%  
	mutate(growth_at_supply_l5 = 5/((1/(alpha * 
										  	eopt^2)) * 5^2 + (1/ps - 2/(alpha * eopt)) * 5 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 7:9)


istars_boot %>% 
	mutate(growth_at_supply_l100 = 100/((1/(alpha * 
												eopt^2)) * 100^2 + (1/ps - 2/(alpha * eopt)) * 100 + (1/alpha))) %>%  
	mutate(growth_at_supply_l20 = 20/((1/(alpha * 
										  	eopt^2)) * 20^2 + (1/ps - 2/(alpha * eopt)) * 20 + (1/alpha))) %>%  
	mutate(growth_at_supply_l5 = 5/((1/(alpha * 
											eopt^2)) * 5^2 + (1/ps - 2/(alpha * eopt)) * 5 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	gather(key = supply_level, value = growth_rate, 11:13) %>% 
	ggplot(aes(x = i_star, y = growth_rate, color = ancestor_id)) + geom_point(alpha =0.5) +
	facet_wrap( ~ supply_level, scales = "free") + geom_smooth(method = "lm") +
	# geom_point(aes(x = n_star, y = growth_rate), data = all_rstars6_anc, color = "black") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (/day)") + xlab("I* (umol m2/s)")
ggsave("figures/selection-gradient-light-bootstrap-supply-all-light.png", width = 10, height = 4)



istars_boot %>% 
	mutate(growth_at_supply_light = 5/((1/(alpha * 
										   	eopt^2)) * 5^2 + (1/ps - 2/(alpha * eopt)) * 5 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_light ~ i_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "i_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/light-selection-gradient-slopes.png", width = 6, height = 4)

istars_boot %>% 
	mutate(growth_at_supply_light = 100/((1/(alpha * 
										   	eopt^2)) * 100^2 + (1/ps - 2/(alpha * eopt)) * 100 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_light ~ i_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "i_star") %>% 
	ggplot(aes(x = reorder(ancestor_id, estimate), y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")
ggsave("figures/light-selection-gradient-slopes-control.png", width = 6, height = 4)


#### figure with all the selection gradients

n_gradients_lowN <- n_monod_boot %>% 
	mutate(growth_at_supply_n = umax*(10/(ks + 10))) %>%  
	mutate(n_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_n ~ n_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "n_star") %>% 
	mutate(resource = "Nitrogen") %>% 
	mutate(resource_level = "Low supply")

n_gradients_highN <- n_monod_boot %>% 
	mutate(growth_at_supply_n = umax*(1000/(ks + 1000))) %>%  
	mutate(n_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_n ~ n_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "n_star") %>% 
	mutate(resource = "Nitrogen") %>% 
	mutate(resource_level = "High supply")


p_gradients_highP <- p_monod_boot %>% 
	mutate(growth_at_supply_p = umax*(50/(ks + 50))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_p ~ p_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "p_star") %>% 
	mutate(resource = "Phosphorus") %>% 
	mutate(resource_level = "High supply")

p_gradients_lowP <- p_monod_boot %>% 
	mutate(growth_at_supply_p = umax*(0.5/(ks + 0.5))) %>%  
	mutate(p_star = ks*0.1/(umax-0.1)) %>% 
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_p ~ p_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "p_star")%>% 
	mutate(resource = "Phosphorus") %>% 
	mutate(resource_level = "Low supply")


i_gradients_highI  <- istars_boot %>% 
	mutate(growth_at_supply_light = 100/((1/(alpha * 
											 	eopt^2)) * 100^2 + (1/ps - 2/(alpha * eopt)) * 100 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_light ~ i_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "i_star") %>% 
	mutate(resource = "Light") %>% 
	mutate(resource_level = "High supply")

i_gradients_lowI  <- istars_boot %>% 
	mutate(growth_at_supply_light = 5/((1/(alpha * 
											 	eopt^2)) * 5^2 + (1/ps - 2/(alpha * eopt)) * 5 + (1/alpha))) %>%  
	filter(treatment == "none") %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_light ~ i_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "i_star") %>% 
	mutate(resource = "Light") %>% 
	mutate(resource_level = "Low supply")


all_gradients <- bind_rows(i_gradients_highI, i_gradients_lowI, n_gradients_highN, n_gradients_lowN, p_gradients_highP, p_gradients_lowP)


all_gradients %>% 
	ggplot(aes(x = ancestor_id, y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + ylab("Selection gradient") +
	xlab("Ancestor")+
	facet_wrap(resource ~ resource_level, scales = "free", nrow = 1) +
	theme(axis.text.x = element_text(angle = 90)) +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank())
ggsave("figures/selection-gradients-all.png", width = 10, height = 3)



changes_all %>% 
	filter(treatment == "C") %>% 
	select(change_p_star, everything()) %>% View


library(broom)
all_stars3 %>% 
	mutate(growth_at_supply_p = p_umax*(0.5/(p_ks + 0.5))) %>% 
	mutate(growth_at_supply_n = n_umax*(10/(n_ks + 10))) %>% 
	group_by(ancestor_id) %>% 
	do(tidy(lm(growth_at_supply_p ~ p_star, data = .), conf.int = TRUE)) %>% 
	filter(term == "p_star") %>% 
	ggplot(aes(x = ancestor_id, y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1)



all_stars2 %>% 
	filter(treatment == "none", ancestor_id %in% c("anc4", "cc1690")) %>% 
	ggplot(aes(x = n_star, y = p_star, color = ancestor_id)) + geom_point() +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = ancestor_id),
				 size = 0.5, linetype = 3) +
	geom_segment(aes(x = max(n_star), y = max(p_star), xend = max(n_star) + nc*7, yend = max(p_star) + pc*7, color = ancestor_id),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = ancestor_id),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = ancestor_id),
				 size = 0.5, linetype = 1) +
	geom_point(aes(x = max(n_star), y = max(p_star)), color = "black") + 
	# xlim(0, 3) + ylim(0, 0.5) +
	ylab("P (uM)") + xlab("N (uM)") +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")
ggsave("figures/zngis-ancestors-only.png", width = 6, height = 4)


all_stars2 %>% 
	filter(treatment %in% c("L"), ancestor_id %in% c("anc4", "cc1690")) %>% 
	ggplot(aes(x = n_star, y = p_star, color =treatment)) + geom_point() +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7, color = treatment),
				 size = 0.5, linetype = 3) +
	geom_segment(aes(x = max(n_star), y = max(p_star), xend = max(n_star) + nc*7, yend = max(p_star) + pc*7, color = treatment),
				 size = 0.5, linetype = 1, alpha = 0.5) +
	geom_segment(aes(x = n_star, y = p_star, xend = n_star, yend = 1, color = treatment),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(x = n_star, y = p_star, xend = 4, yend = p_star, color = treatment),
				 size = 0.5, linetype = 1) +
	# geom_abline(slope = 1, intercept = 0) + 
	# geom_point(aes(x = max(n_star), y = max(p_star)), color = "black") + 
	facet_wrap(~ ancestor_id) +
	# scale_color_manual(values = c("blue", "red")) +
	xlim(1, 4) + ylim(0, 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/zngis-ancestors-only-l-selected.png", width = 6, height = 4)


	


# from Tamminen -----------------------------------------------------------

xls <- read_excel("data-raw/stoichiometry/ChlamEE_phenotypic_measurements.xlsx") %>%
	mutate(Ratio = (2*Cass) / Respiration,
		   NO3toBiomass = NO3 / Biomass,
		   PO4toBiomass = PO4 / Biomass,
		   Selection_Treatment = factor(Selection_Treatment,
		   							 levels=c("A", "C", "L", "N",
		   							 		 "P", "B", "S", "BS")))
xls_long <- gather(xls, key=variable, value=value, -Ancestor, -Selection_Treatment)
cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")

## Define a plotting function to prepare boxplots for each variable
draw_bp <- function(Treatment, output_file)
{
	filter(xls_long, variable == Treatment) %>%
		ggplot(aes(x=Selection_Treatment, y=value, color=Selection_Treatment)) +
		geom_boxplot() +
		geom_point(aes(size=2)) +
		theme_classic() + 
		theme(legend.position="none") +
		scale_color_manual(values=cols)
	ggsave(output_file, plot=last_plot(), useDingbats=FALSE)
}

## Draw the boxplots of the phenotypic measurements
draw_bp("CtoP_Molar", "figures/ctop.pdf")
draw_bp("CtoN_Molar", "figures/cton.pdf")
draw_bp("Cass_Chl", "figures/cass_chl.pdf")
draw_bp("Ratio", "figures/ratio.pdf")
draw_bp("Respiration_Chl", "figures/respiration_chl.pdf")
draw_bp("Biomass", "figures/biomass.pdf")
draw_bp("NO3toBiomass", "figures/no3tobiomass.pdf")
draw_bp("PO4toBiomass", "figures/po4tobiomass.pdf")

## Define a function to test the signficance of the responses
## using two-sided Wilcoxon tests
wilcoxon_test <- function(treatment, stoich, alt)
{
	control <- filter(xls, Selection_Treatment == "A")[[stoich]]
	tmt <- filter(xls, Selection_Treatment == treatment)[[stoich]]
	tidy(wilcox.test(x=control, y=tmt, alternative=alt)) %>%
		mutate(Treatment = treatment, Stoich = stoich) %>%
		dplyr::select(-method, -alternative, -statistic)
}

## And define some utulity functions to adjust the p-values
adj_p_vals <- function(p_table)
{
	sapply(p.adjust.methods,
		   function(meth) p.adjust(p_table$p.value, meth))
}


add_p_vals <- function(p_table)
{
	cbind(p_table, adj_p_vals(p_table))
}

adjust_wilcoxon <- function(column, hypothesis)
{
	map_df(c("B", "BS", "C", "L", "N", "P", "S"),
		   ~wilcoxon_test(., column, hypothesis)) %>%
		add_p_vals
}


## Perform the Wilcoxon tests for each treatment against the ancestors
bind_rows(adjust_wilcoxon("CtoN_Molar", "less"),
		  adjust_wilcoxon("CtoP_Molar", "less"),
		  adjust_wilcoxon("Respiration", "greater"),
		  adjust_wilcoxon("Cass", "greater"),
		  adjust_wilcoxon("Ratio", "two.sided"),
		  adjust_wilcoxon("Biomass", "two.sided"),
		  adjust_wilcoxon("NO3toBiomass", "two.sided"),
		  adjust_wilcoxon("PO4toBiomass", "two.sided")) %>%
	filter(fdr < 0.1) %>%
	dplyr::select(Treatment, Stoich, fdr)


# Trade-offs in umax and ks -----------------------------------------------

source(here("R-scripts", "predict_monod.R"))

all_stars3 %>% 
	ggplot(aes(x = p_ks, y = p_umax)) + geom_point() +
	geom_smooth(method = "lm", color = "grey") + ylab("Max growth rate (umax)") +
	xlab("Half saturation constant (uM P)") 
ggsave("figures/phosphate-umax-ks-tradeoff.png", width = 8, height = 6)

all_stars3 %>% 
	ggplot(aes(x = p_star, y = p_umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	 ylab("Max growth rate (umax)") +
	xlab("P*") 
ggsave("figures/phosphate-umax-rstar-tradeoff.png", width = 8, height = 6)

all_stars3 %>% 
	ggplot(aes(x = n_star, y = n_umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	ylab("Max growth rate (umax)") +
	xlab("N*") 
ggsave("figures/nitrate-umax-rstar-tradeoff.png", width = 8, height = 6)

all_stars3 %>% 
	ggplot(aes(x = n_ks, y = n_umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	ylab("Max growth rate (umax)") +
	xlab("Half saturation constant (um N)") 
ggsave("figures/nitrate-umax-ks-tradeoff.png", width = 8, height = 6)

## maybe try this with the Monod curve instead of the EP curve?
all_stars3 %>% 
	ggplot(aes(x = i_star, y = light_umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	ylab("Max growth rate (umax)") +
	xlab("I*") 
ggsave("figures/light-umax-rstar-tradeoff.png", width = 8, height = 6)

light_monod <- read_csv("data-processed/light-monod-params-direct.csv") 
m <- 0.56
lm1 <- light_monod %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(i_star = ks*m/(umax-m)) 

lm1 %>% 
	ggplot(aes(x = i_star, y = umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	ylab("Max growth rate (umax)") +
	xlab("I*") 
ggsave("figures/light-umax-rstar-tradeoff-monod.png", width = 8, height = 6)


lm1 %>% 
	lm(umax ~ i_star, data = .) %>% summary()

all_stars3 %>% 
	lm(n_umax ~ n_star, data = .) %>% summary()
all_stars3 %>% 
	lm(p_umax ~ p_star, data = .) %>% summary()

lm2 <- lm1 %>% 
	select(population, i_star) %>%
	rename(i_star_m = i_star)

r3 <- left_join(all_stars3, lm2) 


# umax trade-off plots ----------------------------------------------------


p_params <- all_stars3 %>% 
	mutate(p_star = p_ks*m/(p_umax-m))  %>% 
	select(population, p_star, p_umax, treatment) %>% 
	dplyr::rename(r_star = p_star,
		   umax = p_umax) %>% 
	mutate(resource = "Phosphorus") %>% 
	left_join(., treatments)  
n_params <- all_stars3 %>% 
	mutate(n_star = n_ks*m/(n_umax-m))  %>% 
	select(population, n_star, n_umax, treatment) %>% 
	dplyr::rename(r_star = n_star,
		   umax = n_umax) %>% 
	mutate(resource = "Nitrogen") %>% 
	left_join(., treatments) 
i_params <- lm1 %>% 
	select(population, i_star, umax, treatment) %>% 
	dplyr::rename(r_star = i_star,
		   umax = umax) %>% 
	mutate(resource = "Light") %>% 
	left_join(., treatments) 

library(broom)
modi <- lm(umax ~ r_star + ancestor_id, data = i_params)
tidy(modi, conf.int = TRUE)
summary(modi)

modn <- lm(umax ~ r_star + ancestor_id, data = n_params)
tidy(modn, conf.int = TRUE)
summary(modn)

modp <- lm(umax ~ r_star + ancestor_id, data = p_params)
tidy(modp, conf.int = TRUE)
summary(modp)
cor(p_params$umax, p_params$r_star)
cor(n_params$umax, n_params$r_star)
cor(i_params$umax, i_params$r_star)

iplot <- i_params %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	coord_cartesian() +
	# scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	geom_point(size = 3) +
	xlab(expression("I*" ~ (mu * mol ~ m^{-2} * s^{-1}))) + 
	# ylab("")
	ylab(expression("" ~ mu * max ~ (day^{-1}))) 
ggsave("figures/i-umax.png", width = 4, height = 3)

nplot <- n_params %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	coord_cartesian() +
	# scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	geom_point(size = 3) +
	xlab("N* (uM N)") + 
	# ylab("")
	ylab(expression("" ~ mu * max ~ (day^{-1}))) 
ggsave("figures/n-umax.png", width = 4, height = 3)

pplot <- p_params %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	coord_cartesian() +
	# scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	geom_point(size = 3) +
	xlab("P* (uM P)") + 
	ylab(expression("" ~ mu * max ~ (day^{-1}))) 
ggsave("figures/p-umax.png", width = 4, height = 3)

multi_tradeoffs <- plot_grid(pplot, nplot, iplot, pca_plot, labels = c("A", "B", "C", "D"), align = "h", ncol = 2, nrow = 2)
save_plot("figures/trade_offs2_0.56.png", multi_tradeoffs,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 2, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.1
)


# multipanel trade-offs ---------------------------------------------------


first_row_p <-  plot_grid(pplot, nplot, iplot, labels = c("A", "B", "C"), align = "v", nrow = 1, ncol = 3, label_x = 0.34, label_y = 0.95) +
	theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.spacing = unit(c(0, 0), "cm"))
second_row_p <-  plot_grid(pca_plot1, labels = c(''), nrow = 1) +
	theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

gg_allp <-  plot_grid(first_row_p, second_row_p, labels=c('', 'D'), ncol=1, nrow = 2, 
					  rel_heights = c(7, 8), label_x = 0.08, label_y = 0.98, rel_widths = c(0.3, 1), align = "v")


save_plot("figures/umax-tradeoffs-multi.png", gg_allp,
		  # each individual subplot should have an aspect ratio of 1.3
		  # base_aspect_ratio = 1.2,
		  base_height= 3, ncol=1, nrow=2, base_width = 8
)



all_params <- bind_rows(p_params, n_params, i_params)

all_params %>% 
	left_join(., treatments) %>% 
	# filter(resource == "Phosphorus") %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	geom_point(alpha = 0.7) +
	facet_wrap(~ resource, scales = "free") +
coord_cartesian() +
	# scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=14, family = "Arial", color = "black"),
		axis.title=element_text(size=14, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1) +
	xlab("Minimum resource requirement") + ylab("Maximum population growth rate")
# ggsave("figures/poster-umax-rstar-tradeoffs.pdf", width = 10, height = 3.5)
ggsave("figures/umax-rstar-tradeoffs-grey.png", width = 10, height = 4)


umax_tradeoffs_plot <- all_params %>% 
	left_join(., treatments) %>% 
	# filter(resource == "Phosphorus") %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	geom_point(alpha = 0.7) +
	facet_wrap(~ resource, scales = "free") +
	coord_cartesian() +
	# scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=14, family = "Arial", color = "black"),
		axis.title=element_text(size=14, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1) +
	xlab("Minimum resource requirement") + ylab("Maximum population growth rate")


multi_plot_to <- plot_grid(umax_tradeoffs_plot, pca_plot, labels = c("A", "B"), align = "h")

save_plot("figures/umax-trade-offs-pca.png", multi_plot_to,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 0.9
)


p_params <- all_params %>% 
	filter(resource == "Phosphorus") 

all_params %>% 
	filter(resource == "Nitrogen") %>% 
	lm(umax ~ r_star, data = .) %>% summary()

all_params %>% 
	filter(resource == "Light") %>% 
	lm(umax ~ r_star, data = .) %>% summary()


library(modelr)
library(lmodel2)
rma <- lmodel2(umax ~ r_star, data = i_params, range.y = "interval", range.x = "interval")
rma$regression.results
rma$confidence.intervals
rma$rsquare
rma$P.param
summary(rma)


all_params %>% 
	# filter(resource == "Phosphorus") %>% 
	ggplot(aes(x = r_star, y = umax)) +  geom_smooth(method = "lm", color = "grey") +
	geom_point(alpha = 0.7) +
	facet_wrap(~ resource, scales = "free") +
	coord_cartesian() +
	scale_x_reverse() +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank(),
		  strip.text.x = element_text(size=14, family = "Arial", color = "black")) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=14, family = "Arial", color = "black"),
		axis.title=element_text(size=14, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1) +
	xlab("Minimum resource requirement") + ylab("Maximum population growth rate")
ggsave("figures/umax-rstar-tradeoffs.png", width = 9, height = 3.4)


all_params %>% 
	group_by(resource) %>% 
	do(tidy(lm(umax ~ r_star, data = .), conf.int = TRUE)) %>% View

all_stars3 %>% 
	ggplot(aes(x = n_star, y = n_umax)) +
	geom_smooth(method = "lm", color = "grey") +
	geom_point() +
	ylab("Max growth rate (umax)") +
	xlab("N*") 




all_preds_p <- all_stars3 %>% 
	rename(ks = p_ks,
		   umax = p_umax) %>% 
	split(.$population) %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population") %>% 
	rename(phosphate_concentration = nitrate_concentration.x) %>% 
	filter(phosphate_concentration < 50)

all_preds_n <- all_stars3 %>% 
	rename(ks = n_ks,
		   umax = n_umax) %>% 
	split(.$population) %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population") %>% 
	rename(nitrate_concentration = nitrate_concentration.x) 


library(wesanderson)
library(beyonce)
pal <- wes_palette("FantasticFox1", 6, type = "continuous")
all_preds_p2 <- all_preds_p %>% 
	left_join(., treatments, by = "population")

all_preds_n2 <- all_preds_n %>% 
	left_join(., treatments, by = "population")

all_preds_p2 %>% 
	filter(treatment != "none") %>% 
	ggplot(aes(x = phosphate_concentration, y = growth_rate, color = ancestor_id, group = population)) + geom_line(size = 1) +
	# facet_wrap( ~ treatment) +
	geom_line(aes(x = phosphate_concentration, y = growth_rate), 
			  data = filter(all_preds_p2, treatment == "none"), color = 'black', size = 1) +
	# xlim(0, 5) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	# scale_color_manual(values = pal, name =  "Ancestor") +
	ylab("Growth rate (/day)") + xlab("Phosphate concentration (um P)")
ggsave("figures/phosphate-monod-color.png", width = 6, height= 4)

all_preds_n2 %>% 
	filter(treatment != "none") %>% 
	ggplot(aes(x = nitrate_concentration, y = growth_rate, color = ancestor_id, group = population)) + geom_line(size = 1) +
	# facet_wrap( ~ treatment) +
	geom_line(aes(x = nitrate_concentration, y = growth_rate), 
			  data = filter(all_preds_n2, treatment == "none"), color = 'black', size = 1) +
	# xlim(0, 5) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	# scale_color_manual(values = pal, name =  "Ancestor") +
	ylab("Growth rate (/day)") + xlab("Nitrate concentration (um P)")
ggsave("figures/nitrate-monod-color.png", width = 6, height= 4)


all_preds_n2 %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>%
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	group_by(treatment, nitrate_concentration) %>% 
	summarise_each(funs(mean, std.error), growth_rate) %>% 
	ggplot(aes(x = nitrate_concentration, y = mean, group = treatment, color = treatment)) + geom_line(size = 1) +
	geom_ribbon(aes(ymin = mean-std.error, ymax = mean + std.error, fill = treatment), alpha = 0.4, linetype = "blank") +
	# facet_wrap(~ treatment, nrow = 2) +
	ylab("Growth rate (per day)") + xlab("Nitrate (uM N)") + 
	scale_color_manual(values = cols, name = "Selection treatment") +
	scale_fill_manual(values = cols, name = "Selection treatment") +
	xlim(0, 25) + ylim(0, 1)
ggsave("figures/nitrate-monod-treatment-average.png", width = 8, height= 4)


all_preds_p2 %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>%
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	group_by(treatment, phosphate_concentration) %>% 
	summarise_each(funs(mean, std.error), growth_rate) %>% 
	ggplot(aes(x = phosphate_concentration, y = mean, group = treatment, color = treatment)) + geom_line(size = 1) +
	# geom_ribbon(aes(ymin = mean-std.error, ymax = mean + std.error, fill = treatment), alpha = 0.4, linetype = "blank") +
	# facet_wrap(~ treatment, nrow = 2) +
	ylab("Growth rate (per day)") + xlab("Phosphate (uM P)") + 
	scale_color_manual(values = cols, name = "Selection treatment") +
	scale_fill_manual(values = cols, name = "Selection treatment") +
	xlim(0, 2.5) + ylim(0, 1)
ggsave("figures/phosphate-monod-treatment-average-no-error-lowP.png", width = 8, height= 4)




growth_16p <- all_preds_p %>% 
	select(population, phosphate_concentration, growth_rate) %>% 
	# group_by(population, nitrate_concentration.x) %>% 
	# summarise(growth_rate_mean = mean(growth_rate_mean)) %>% 
	# filter(temperature %in% c(10, 16, 22, 28, 34, 40)) %>% 
	spread(key = phosphate_concentration, value = growth_rate, 2:3) %>% 
	ungroup() %>% 
	select(-population)

unique(predictions_summary$population)


pca16 <- rda(growth_16p)
summary(pca16)


pca16 <- rda(growth_16p, scale=TRUE)

loadings16 <- scores(pca16,choices=c(1,2))
summary(eigenvals(pca16))

phosphate_conc <- unique(all_preds_p$phosphate_concentration)
pcs16 <- as_data_frame((loadings16[[1]]))
pc1_16 <- pcs16 %>% 
	mutate(phosphate = phosphate_conc) 
pc1_16 %>% 
	filter(phosphate != 0) %>% 
	ggplot(aes(x = phosphate, y = PC2)) + geom_line() +
	xlim(0, 50) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Phosphate concentration (uM N)") +xlim(0, 50)
ggsave("figures/phosphate-pca.png", width = 8, height = 6)



#### plots for Anita
cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")


all_stars4 <- all_stars3 %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
		   							 levels=c("A", "C", "L", "N",
		   							 		 "P", "B", "S", "BS")))
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

all_stars4_sum <- all_stars4 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), p_star, n_star, i_star)
