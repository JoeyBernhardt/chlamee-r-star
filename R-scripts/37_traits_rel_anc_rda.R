


### RDA plots with the trait changes relative to the ancestor


library(tidyverse)
library(cowplot)
library(vegan)
library(here)
library(broom)

cbbPalette <- c("black", "#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")

all_sizes3 <- read_csv("data-processed/all_cell_sizes_low_resources.csv")

light_monod <- read_csv("data-processed/light-monod-params-direct.csv")

salt_tol_CI <- read_csv("data-processed/salt-tol-CI-boot-12.csv") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	select(treatment, ancestor_id, population, change_salt_tol_mean) %>% 
	dplyr::rename(salt_tol = change_salt_tol_mean)
library(dplyr)

ancestors <- all_sizes3 %>% 
	filter(treatment %in% c("Ancestors")) %>% 
	dplyr::rename(ancestor_pstar = p_star) %>% 
	dplyr::rename(ancestor_istar = i_star) %>% 
	dplyr::rename(ancestor_nstar = n_star) %>% 
	dplyr::rename(ancestor_biovolume = phosphate_biovolume) %>% 
	dplyr::rename(ancestor_pc = pc) %>% 
	dplyr::rename(ancestor_nc = nc) %>% 
	select(ancestor_id, contains("ancestor"))

### how what was the % change relative to ancestors?
prot_table <- left_join(all_sizes3, ancestors) %>% 
	left_join(., salt_tol_CI) %>% 
	mutate(change_salt_tol = salt_tol) %>% 
	mutate(change_pstar = (((ancestor_pstar - p_star)/ancestor_pstar)*-100)) %>% 
	mutate(change_istar = (ancestor_istar - i_star)/ancestor_istar*-100) %>% 
	mutate(change_nstar = (ancestor_nstar - n_star)/ancestor_nstar*-100) %>%
	mutate(change_pc = scale(ancestor_pc - pc, center = FALSE)) %>% 
	mutate(change_nc = scale(ancestor_nc - nc, center = FALSE)) %>% 
	mutate(change_biovolume = scale(ancestor_biovolume - phosphate_biovolume, center = FALSE)) %>% 
	# filter(treatment != "Ancestors") %>% 
	select(contains("change"), ancestor_id, treatment, -contains("fold")) %>%
	dplyr::rename(Strain = ancestor_id) %>% 
	dplyr::rename(Treatment = treatment) %>% 
	mutate(Treatment = str_replace(Treatment, "C", "control")) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("Ancestors", "control", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(change_salt_tol = ifelse(Treatment == "Ancestors", 0, change_salt_tol))

## check for % change in salt tolerance

salts <- read_csv("data-processed/salt-tolerances-sigmoid.csv") 

salt_anc <- salts %>% 
	filter(treatment == "A") %>% 
	dplyr::rename(ancestor_salt_tol = salt_tolerance) %>% 
	select(ancestor_id, ancestor_salt_tol)

salts2 <- left_join(salts, salt_anc)

salts2 %>% 
	mutate(change_salt_tol = (ancestor_salt_tol - salt_tolerance)/ancestor_salt_tol*-100) %>% 
	filter(treatment %in% c("S", "BS")) %>% View




prot_table <- left_join(all_sizes3, ancestors) %>% 
	left_join(., salt_tol_CI) %>% 
	mutate(change_salt_tol = scale(salt_tol, center = FALSE)) %>% 
	mutate(change_pstar = scale(ancestor_pstar - p_star, center = FALSE)) %>% 
	mutate(change_istar = scale(ancestor_istar - i_star, center = FALSE)) %>% 
	mutate(change_nstar = scale(ancestor_nstar - n_star, center = FALSE)) %>% 
	mutate(change_pc = scale(ancestor_pc - pc, center = FALSE)) %>% 
	mutate(change_nc = scale(ancestor_nc - nc, center = FALSE)) %>% 
	mutate(change_biovolume = scale(ancestor_biovolume - phosphate_biovolume, center = FALSE)) %>% 
	# filter(treatment != "Ancestors") %>% 
	select(contains("change"), ancestor_id, treatment, -contains("fold")) %>% 
	dplyr::rename(Strain = ancestor_id) %>% 
	dplyr::rename(Treatment = treatment) %>% 
	mutate(Treatment = str_replace(Treatment, "C", "control")) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("Ancestors", "control", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(change_salt_tol = ifelse(Treatment == "Ancestors", 0, change_salt_tol))

allsizes5 <- left_join(all_sizes3, ancestors) %>% 
	mutate(change_pstar = ancestor_pstar - p_star) %>% 
	mutate(change_istar = ancestor_istar - i_star) %>% 
	mutate(change_nstar = ancestor_nstar - n_star) %>% 
	mutate(change_pc = ancestor_pc - pc) %>% 
	mutate(change_nc = ancestor_nc - nc) %>% 
	mutate(change_biovolume = ancestor_biovolume - phosphate_biovolume) %>% 
	select(population, treatment, ancestor_id, p_star, ancestor_pstar, change_pstar, everything()) %>% 
	select(population, treatment, ancestor_id, phosphate_biovolume, ancestor_biovolume, change_biovolume, everything()) 

spe <- dplyr::select(prot_table, -Treatment, -Strain)
FULL.cap <- capscale(spe ~ Treatment + Condition(Strain), data=prot_table)
FULL.cap_str <- capscale(spe ~ Strain + Condition(Treatment), data=prot_table)

library(vegan)
anova.cca(FULL.cap, by = "axis", step = 1000)
anova.cca(FULL.cap, step = 1000)

coef(FULL.cap) ### gives the the equivalent of regression coefficients for each explanatory variable
# on each canonical axis

RsquareAdj(FULL.cap)$adj.r.squared
RsquareAdj(FULL.cap)$r.squared

summary(FULL.cap)

plot(FULL.cap, scaling = 1)
full.sc <- scores(FULL.cap, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, full.sc[,1], full.sc[,2], length = 0, lty = 1, col = "green")

plot(FULL.cap, scaling = 2)
full2.sc <- scores(FULL.cap, choices = 1:2, display = "sp")
arrows(0, 0, full2.sc[,1], full2.sc[,2], length = 0, lty = 1, col = "red")

plot(FULL.cap, scaling = 1, display = c("sp", "lc", "cn"))
arrows(0, 0, full.sc[,1], full.sc[,2], length = 0, lty = 1, col = "red")

plot(FULL.cap, scaling = 2, display = c("sp", "lc", "cn"))
arrows(0, 0, full2.sc[,1], full2.sc[,2], length = 0, lty = 1, col = "red")

color_column <- "Treatment"
CAP1 <- scores(FULL.cap, display="wa", scaling=3)[,1]
CAP2 <- scores(FULL.cap, display="wa", scaling=3)[,2]
Res.dim <- as.data.frame(scores(FULL.cap, display="wa", scaling=3)[,1:2])
Res.dim$Strain <- prot_table$Strain
Res.dim$Treatment <- prot_table$Treatment
Res.dim
write_csv(Res.dim, "data-processed/RDA-scores-trait-change.csv")

dims <- Res.dim

dims_summary <- dims %>% dplyr::group_by(Treatment) %>% 
	summarise(CAP1mean=mean(CAP1), CAP2mean=mean(CAP2), CAP1sd=sd(CAP1), CAP2sd=sd(CAP2))


length(unique(dims_summary$Treatment))
length(unique(dims$Treatment))

		ggplot() + 
		geom_errorbar(aes(x = CAP1mean, ymin = CAP2mean - CAP2sd, ymax = CAP2mean + CAP2sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
		geom_errorbarh(aes(y = CAP2mean, xmin = CAP1mean - CAP1sd, xmax = CAP1mean + CAP1sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
			geom_point(aes(x = CAP1, y = CAP2, color = Treatment), data = dims, size = 2) +
		scale_colour_manual(values=cbbPalette) 
	
	
	
	# geom_point(data=Res.dim, aes(x=CAP1, y= CAP2, color=color_column)) +
	geom_errorbarh(aes(y = CAP2, xmin = CAP1mean - CAP1sd, xmax = CAP1mean + CAP1sd)) +
	geom_errorbar(aes(ymin = CAP2mean - CAP2sd, ymax = CAP2mean + CAP2sd), width = 0) +
	scale_colour_manual(values=cbbPalette) 
ggsave("figures/rda-trait-change-rel-ancestors-salt.png", width = 8, height = 6)


## ok now let's get the angles and distances of each treatment centroid from the Ancestor centroid

View(dims_summary)
ancestor_y <- dims_summary$CAP2mean[dims_summary$Treatment == "Ancestors"]
ancestor_x <- dims_summary$CAP1mean[dims_summary$Treatment == "Ancestors"]

dims_anc <- dims_summary %>% 
	mutate(new_x = CAP1mean - ancestor_x) %>% 
	mutate(new_y = CAP2mean - ancestor_y)

View(dims_anc)

angs_dist <- dims_anc %>% 
	mutate(angle = atan2(new_y, new_x)*(180/pi)) %>% 
	mutate(distance = sqrt((new_x^2) + (new_y^2)))

angs_dist %>% 
	ggplot(aes(x = angle, y = distance, color = Treatment)) + 
	geom_point(size = 2) + 
	geom_segment(aes(x = 0, y = 0, xend = angle, yend = distance, color = Treatment)) + 
	# coord_polar(theta="x", start= -135) +
	xlim(-180, 180) + 
	# scale_x_continuous(breaks = c(0, 90, 180)) +
	scale_colour_manual(values=cbbPalette) + theme_bw() + xlab("Angle from ancestor centroid") +
	ylab("Distance from ancestor centroid")
ggsave("figures/angle-distance-plot-rda-non-polar.png", width = 8, height = 6)
	

ggplot() +
	geom_point(data=rda_norm, aes(x=Angle, y=Dist, color=Treatment)) + 
	geom_point(data=rda_means,
			   aes(x=Mean_angle, y=Mean_dist, color=Treatment), size = 5) +
	coord_polar(theta="x", start=0) +
	scale_colour_manual(values=cbbPalette)


### joey's significance tests
dist_test_cap2 <- aov(lm(CAP2~Treatment, data=Res.dim))
TukeyHSD(dist_test_cap2)


dist_test_cap1 <- aov(lm(CAP1~Treatment, data=Res.dim))
TukeyHSD(dist_test_cap1)

# Test for significant differences between Strains and Treatments
anova(FULL.cap, data=prot_table) # Treatments: p < 0.001
anova(FULL.cap_str, data=prot_table) # Strains: p = 0.762

cos(150*pi/180)


### now the raw trait associations
library(vegan)
library(tidyverse)
all_sizes3 <- read_csv("data-processed/all_cell_sizes_low_resources.csv") ### ok this has the old R*s
salt_tolerances <- read_csv("data-processed/salt-tolerances-sigmoid.csv") %>% 
	select(population, salt_tolerance)
m <- 0.56

nmonod <- read.csv("data-processed/nitrate-monod-params-bootstrap.csv") %>% 
	mutate(population = as.character(population_index)) %>% 
	mutate(n_star = ks*m/(umax-m)) %>% 
	dplyr::rename(n_umax = umax,
		   n_ks = ks) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), n_star, n_umax, n_ks)
pmonod <- read.csv("data-processed/phosphate-monod-params-bootstrap.csv") %>% 
	mutate(population = as.character(population_index)) %>% 
	mutate(p_star = ks*m/(umax-m)) %>% 
	dplyr::rename(p_umax = umax,
		   p_ks = ks) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), p_star, p_umax, p_ks)
imonod <- read.csv("data-processed/light_monod_params_bootstrapped.csv") %>% 
	mutate(population = as.character(population_index)) %>% 
	mutate(i_star = ks*m/(umax-m)) %>% 
	dplyr::rename(i_umax = umax,
		   i_ks = ks) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), i_star, i_umax, i_ks)

all_sizes4 <- all_sizes3 %>% 
	select(-p_ks, -p_umax, -p_star, -n_ks, -n_umax, -n_star, -light_umax, -i_star)

### note to self -- find out why the bootstrapped R*s are different than the estimated ones.
traits <- left_join(all_sizes4, salt_tolerances) %>% 
	left_join(imonod) %>% 
	left_join(pmonod) %>% 
	left_join(nmonod)


prot_table <- traits %>% 
	mutate(salt_tolerance = scale(salt_tolerance, center = FALSE)) %>% 
	mutate(pstar = scale(p_star, center = FALSE)) %>% 
	mutate(istar = scale(i_star, center = FALSE)) %>% 
	mutate(nstar = scale(n_star, center = FALSE)) %>% 
	mutate(pc = scale(pc, center = FALSE)) %>% 
	mutate(nc = scale(nc, center = FALSE)) %>% 
	mutate(cons_vector_new = scale(cons_vector_new, center = FALSE)) %>% 
	mutate(pbiovolume = scale(phosphate_biovolume, center = FALSE)) %>% 
	mutate(nbiovolume = scale(nitrate_biovolume, center = FALSE)) %>% 
	mutate(lbiovolume = scale(light_biovolume, center = FALSE)) %>% 
	# select(ancestor_id, treatment, salt_tolerance, pstar, nstar, istar, pc, nc, biovolume, cons_vector_new) %>%
	select(ancestor_id, treatment, salt_tolerance, pstar, nstar, istar, cons_vector_new, pbiovolume) %>%
	dplyr::rename(Strain = ancestor_id) %>% 
	dplyr::rename(Treatment = treatment) %>% 
	mutate(Treatment = str_replace(Treatment, "C", "control")) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("Ancestors", "control", "L", "N",
							  		 "P", "B", "S", "BS"))) 


spe <- dplyr::select(prot_table, -Treatment, -Strain)
FULL.cap <- capscale(spe ~ Treatment + Condition(Strain), data=prot_table)
# FULL.cap <- capscale(spe ~ Treatment, data=prot_table)
FULL.cap_str <- capscale(spe ~ Strain + Condition(Treatment), data=prot_table)

anova.cca(FULL.cap, by = "axis", step = 1000)
anova.cca(FULL.cap)

coef(FULL.cap) ### gives the the equivalent of regression coefficients for each explanatory variable
# on each canonical axis

RsquareAdj(FULL.cap)$adj.r.squared
RsquareAdj(FULL.cap)$r.squared

summary(FULL.cap)


# Test for significant differences between Strains and Treatments
anova(FULL.cap, data=prot_table) # Treatments: p < 0.001
anova(FULL.cap_str, data=prot_table) # Strains: p = 0.762

dist_test_cap2 <- aov(lm(CAP2~Treatment, data=Res.dim))
summary(dist_test_cap2)
TukeyHSD(dist_test_cap2)

dist_test_cap1 <- aov(lm(CAP1~Treatment, data=Res.dim))
summary(dist_test_cap1)
TukeyHSD(dist_test_cap1)


color_column <- "Treatment"
CAP1 <- scores(FULL.cap, display="wa", scaling=3)[,1]
CAP2 <- scores(FULL.cap, display="wa", scaling=3)[,2]
Res.dim <- as.data.frame(scores(FULL.cap, display="wa", scaling=3)[,1:2])
Res.dim$Strain <- prot_table$Strain
Res.dim$Treatment <- prot_table$Treatment
Res.dim





dims <- Res.dim

library(tidyverse)
dims_summary <- dims %>% 
	dplyr::group_by(Treatment) %>% 
	summarise(CAP1mean=mean(CAP1), CAP2mean=mean(CAP2), CAP1sd=sd(CAP1), CAP2sd=sd(CAP2)) %>% 
	mutate(Treatment = str_replace(Treatment, "control", "C")) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) 

dims <- dims %>% 
	mutate(Treatment = str_replace(Treatment, "control", "C")) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) 

trait_vectors <- as.data.frame(scores(FULL.cap, choices = 1:2, display = "sp", scaling = 2)) %>% 
	mutate(trait = rownames(.)) %>%
	mutate(trait = str_replace(trait, "cons_vector_new", "P:N")) %>% 
	mutate(trait = str_replace(trait, "pbiovolume", "Biovolume")) %>% 
	mutate(trait = str_replace(trait, "istar", "I*")) %>% 
	mutate(trait = str_replace(trait, "pstar", "P*")) %>% 
	mutate(trait = str_replace(trait, "nstar", "N*")) %>% 
	mutate(trait = str_replace(trait, "salt_tolerance", "Salt \n tolerance")) 
	
# arrows(0, 0, full2.sc[,1], full2.sc[,2], length = 0, lty = 1, col = "red")

library(cowplot)
rda_plot <- ggplot() + 
	geom_segment(aes(x = 0, y = 0, xend = CAP1, yend = CAP2, text =  trait), data = trait_vectors, color = "black",
				 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
	geom_hline(yintercept = 0, color = "grey") + geom_vline(xintercept = 0, color = "grey") +
	geom_errorbar(aes(x = CAP1mean, ymin = CAP2mean - CAP2sd, ymax = CAP2mean + CAP2sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
	geom_errorbarh(aes(y = CAP2mean, xmin = CAP1mean - CAP1sd, xmax = CAP1mean + CAP1sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
	geom_point(aes(x = CAP1, y = CAP2, color = Treatment), data = dims, size = 4) +
	geom_point(aes(x = CAP1, y = CAP2), data = dims, size = 4, shape = 1, color = "black") +
	scale_colour_manual(values= cbbPalette, name = "Selection environment") +
	xlab("RDA1 (56.61% of constrained variation)") + ylab("RDA2 (28.34% of constrained variation)") +
	# geom_text(aes(x = CAP1, y = CAP2, label = trait), data = trait_vectors, color = "black", size = 5.5,
	# 		  hjust = 0, nudge_x = -0.29, nudge_y = 0.27) + 
	theme(legend.position = c(0.01, 0.2),
		  axis.text = element_text(size=18),
		  axis.title=element_text(size=18),
		  legend.title = element_text(colour="black", size=16),
		  legend.text = element_text(colour="black", size=16))
ggsave("figures/rda-traits-all-traits.png", width = 8, height = 6)

dims_anc <- dims %>% 
	filter(Treatment == "Ancestors") %>% 
	dplyr::rename(anc_cap1 = CAP1) %>% 
	dplyr::rename(anc_cap2 = CAP2) %>% 
	select(-Treatment)

dims_distances <- left_join(dims, dims_anc) %>% 
	mutate(change_x = CAP1 - anc_cap1) %>% 
	mutate(change_y = CAP2 - anc_cap2) %>% 
	mutate(angle = atan2(change_y, change_x)*(180/pi)) %>% 
	mutate(distance = sqrt((change_x^2) + (change_y^2))) %>% 
	mutate(angle_360 = ifelse(angle < 0, 360 + angle, angle))

dims_dist_sum <- dims_distances %>% 
	group_by(Treatment) %>% 
	summarise_each(funs(mean), angle, distance, angle_360)

ggplot() + 
	# geom_errorbar(aes(x = CAP1mean, ymin = CAP2mean - CAP2sd, ymax = CAP2mean + CAP2sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
	# geom_errorbarh(aes(y = CAP2mean, xmin = CAP1mean - CAP1sd, xmax = CAP1mean + CAP1sd, color = Treatment), width = 0, data = dims_summary, size = 1) +
	geom_point(aes(x = CAP1, y = CAP2, color = Treatment), data = dims, size = 2, alpha = 0.3) +
	geom_point(aes(x = CAP1, y = CAP2), color = "black", data = filter(dims, Treatment == "Ancestors"), size = 2, alpha = 1) +
	geom_segment(aes(x = anc_cap1, y = anc_cap2, xend = CAP1, yend = CAP2, color = Treatment),
				 alpha = 1, data = filter(dims_distances, Treatment != "Ancestors"),
				 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
	geom_point(aes(x = CAP1, y = CAP2), color = "black", data = filter(dims, Treatment == "Ancestors"), size = 2, alpha = 1) +
	scale_colour_manual(values=cbbPalette) +ylab("RDA2") + xlab("RDA1")
ggsave("figures/rda-traits-salt-vectors.png", width = 6, height = 4)


## ok now let's get the angles and distances of each treatment centroid from the Ancestor centroid




# ggplot() +
# 	geom_point(aes(x = CAP1, y = CAP2, color = Treatment), data = filter(dims, Treatment == "Ancestors")) +
# 	geom_segment(aes(x = anc_cap1, y = anc_cap2, xend = CAP1, yend = CAP2, color = Treatment),
# 				 alpha = 1, data = filter(dims_distances, Treatment != "Ancestors"),
# 				 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
# 	scale_colour_manual(values=cbbPalette) 



dims_distances %>% 
	ggplot(aes(x = angle_360, y = distance, color = Treatment)) + 
	geom_point(size = 2) + 
	# geom_segment(aes(x = 0, y = 0, xend = angle, yend = distance, color = Treatment), alpha = 0.5) + 
	
	# geom_segment(aes(x = 0, y = 0, xend = angle, yend = distance, color = Treatment), alpha = 1, data = dims_dist_sum, size = 1.5) + 
	coord_polar(theta="x", start= 0) +
	xlim(0, 360) +
	geom_point(aes(x = angle_360, y = distance, color = Treatment), data = dims_dist_sum, size = 4) +
	geom_point(aes(x = angle_360, y = distance), data = dims_dist_sum, size = 4, color = "black", shape = 1) +
	# xlim(-180, 180) + 
	# scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), limits = c(-180, 180)) +
	scale_colour_manual(values=cbbPalette) +
	theme_bw() +
	xlab("Angle from ancestor") +
	ylab("Distance from ancestor")
ggsave("figures/angle-distance-plot-traits-raw-rda-polar.png", width = 8, height = 5)


dims_distances %>% 
	ggplot(aes(x = angle, y = distance, color = Treatment)) + 
	geom_point(size = 2) + 

	geom_segment(aes(x = 0, y = 0, xend = angle, yend = distance, color = Treatment), alpha = 0.5) + 
	
	geom_segment(aes(x = 0, y = 0, xend = angle, yend = distance, color = Treatment), alpha = 1, data = dims_dist_sum, size = 1) + 
	# coord_polar(theta="x", start= 0) +
	geom_point(aes(x = angle, y = distance, color = Treatment), data = dims_dist_sum, size = 4) +

	# xlim(0, 360) +

	# xlim(-180, 180) + 
	scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), limits = c(-180, 180)) +
	scale_colour_manual(values=cbbPalette) +
	geom_point(size = 2) + 
	geom_point(aes(x = angle, y = distance), data = dims_dist_sum, size = 4, color = "black", shape = 1) +
	# theme_bw() +
	xlab("Angle from ancestor") +
	ylab("Distance from ancestor")
ggsave("figures/angle-distance-plot-traits-raw-rda-non-polar.png", width = 8, height = 6)







summary(FULL.cap)

plot(FULL.cap, scaling = 1) ## distance biplot, (2) The angles between response and explanatory
# variables in the biplot reflect their correlations (but not the angles among
# 													response variables).
full.sc <- scores(FULL.cap, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, full.sc[,1], full.sc[,2], length = 0, lty = 1, col = "purple")

plot(FULL.cap, scaling = 2) ## correlation biplot, The angles in the biplot between response and
# explanatory variables, and between response variables themselves or explanatory
# variables themselves, reflect their correlations. The relationship
# between the centroid of a qualitative explanatory variable and a response variable
# (species) is found by projecting the centroid at right angle on the species. 
# variable (as for individual objects). (4) Distances among centroids, and between
# centroids and individual objects, do not approximate their Euclidean
# distances.

full2.sc <- scores(FULL.cap, choices = 1:2, display = "sp", scaling = 2)
arrows(0, 0, full2.sc[,1], full2.sc[,2], length = 0, lty = 1, col = "red")

plot(FULL.cap, scaling = 1, display = c("cn", "lc", "wa", "sp"))
arrows(0, 0, full.sc[,1], full.sc[,2], length = 0, lty = 1, col = "red")

plot(FULL.cap, scaling = 2, display = c("sp"))
arrows(0, 0, full2.sc[,1], full2.sc[,2], length = 0, lty = 1, col = "red")
