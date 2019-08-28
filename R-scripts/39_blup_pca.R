

library(tidyverse)

nitrate_growth <- read_csv("data-processed/nitrate-sequential-fits.csv") %>% 
	rename(growth = estimate) %>% 
	mutate(experiment = "nitrate") %>% 
	mutate(environment = paste(experiment, nitrate_concentration, sep = "_")) %>% 
	separate(well_plate, into = c("well", "plate")) %>% 
	select(experiment, environment, population, growth, well)
light_growth <- read_csv("data-processed/light-sequential-fits.csv") %>% 
	rename(growth = estimate) %>% 
	mutate(experiment = "light") %>% 
	mutate(environment = paste(experiment, light, sep = "_")) %>% 
	separate(well_plate, into = c("well", "plate")) %>% 
	select(experiment, environment, population, growth, well)

phosphate_growth <- read_csv("data-processed/phosphate-sequential-fits.csv") %>% 
	rename(growth = estimate) %>% 
	mutate(experiment = "phosphate") %>% 
	mutate(environment = paste(experiment, phosphate_concentration, sep = "_")) %>% 
	separate(well_plate, into = c("well", "plate")) %>% 
	select(experiment, environment, population, growth, well)

salt_growth <- read.csv("data-processed/salt_growth_rates.csv") %>% 
	separate(salt_level, into = c("S", "concentration"), sep = 1) %>% 
	mutate(concentration = as.numeric(concentration)) %>% 
	group_by(unique_id) %>% 
	distinct(mu, .keep_all = TRUE)


dat <- bind_rows(nitrate_growth, phosphate_growth, light_growth) %>% 
	mutate(environment = as.factor(environment)) %>% 
	mutate(population = as.factor(population)) %>% 
	mutate(well = as.factor(well))


library(lme4)
library(plyr)
library(dplyr)

fit.1 <- lmer(dat$growth ~ 1 + (1|dat$population))
fit.2 <- lmer(dat$growth ~ 1 + dat$environment + (1|dat$population))
fit.3 <- lmer(dat$growth ~ 1 + dat$environment + (1|dat$population) + (1|dat$population:dat$environment))
fit_4 <- lmer(dat$growth ~ 1 + dat$environment  + (1|dat$population) + (1|dat$well) +
			  	(1|dat$population:dat$environment))
resid_4 <- resid(fit_4)
hist(resid_4)

fit_aov <- anova(fit.1, fit.2, fit.3, fit_4)
fit_aov
summary(fit.3)
summary(fit_4)


ran.3 <- ranef(fit.3)
tab.gxe <- ran.3$`dat$population:dat$environment`

fix_4 <- fixef(fit_4)
ran_4 <- ranef(fit_4)
new_gxe <- ran_4$`dat$population:dat$environment`

plot(new_gxe$`(Intercept)`, tab.gxe$`(Intercept)`)
cor.test(new_gxe$`(Intercept)`, tab.gxe$`(Intercept)`)


library(readxl)
library(janitor)
treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "Ancestors", treatment))
newg1 <- new_gxe %>% 
	mutate(pop_env = rownames(.)) %>% 
	separate(pop_env, into =c("population", "experiment"), sep = c(":")) %>% 
	# separate(experiment, into = c("experiment", "resource_level"), sep = "_") %>% 
	dplyr::rename(blup = "(Intercept)") %>% 
	spread(key = experiment, value = blup) %>% 
	left_join(., treatments)

write_csv(newg1, "data-processed/gxe_blups.csv")


newg <- new_gxe %>% 
	mutate(pop_env = rownames(.)) %>% 
	separate(pop_env, into =c("population", "experiment"), sep = c(":")) %>% 
	# separate(experiment, into = c("experiment1", "resource_level"), sep = "_", remove = FALSE) %>% 
	dplyr::rename(blup = "(Intercept)") %>% 
	spread(key = experiment, value = blup) %>% 
	select(-population)

mean_cent <- scale(t(newg))

pca_blup <- prcomp(newg)
pca_blup <- prcomp(mean_cent)
pca_blup1 <- rda(mean_cent, scale = FALSE)

summary(pca_blup1)
summary(pca_blup)



biplot(pca_blup1, 
	   col=c("black", "darkgray"), 
	   cex=0.8, 
	   xlim=c(-2,2),
	   ylim=c(-2, 2),
	   xlab="PC1 (21.93% explained var.)",
	   ylab="PC2 (16.05% explained var.)")

biplot(pca_blup, 
	   col=c("black", "darkgray"), 
	   cex=0.8, 
	   xlim=c(-0.5,0.5),
	   ylim=c(-0.5, 0.5),
	   xlab="PC1 (21.93% explained var.)",
	   ylab="PC2 (16.05% explained var.)")

ggbiplot(pca_blup, ellipse = FALSE, size = 6, circle = FALSE, scale = 1)

biplot(pca_blup1, 
	   col=c("black", "darkgray"), 
	   cex=0.8, 
	   xlim=c(-1,1),
	   ylim=c(-1, 1),
	   xlab="PC1",
	   ylab="PC2")



scores(pca_blup, choices = 1:2)
library(vegan)

all_sizes3 <- read_csv("data-processed/all_cell_sizes_low_resources.csv") %>% 
	select(population, p_star, n_star, i_star, cons_vector_new, phosphate_biovolume, nitrate_biovolume, light_biovolume) %>% 
	arrange((population)) %>% 
	select(-population)
(fit <- envfit(pca_blup1, all_sizes3, perm = 999))
trait_vectors1 <- as.data.frame(fit$vectors$arrows) %>% 
	mutate(trait = rownames(.))

fit

plot(pca_blup1)
plot(fit)
plot(fit, p.max = 0.05, col = "red")

cols_anc  <- c("black", "#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")

library(ggbiplot)
ggbiplot(pca_blup, group = newg1$treatment,
				   ellipse = FALSE, size = 6, circle = FALSE, scale = 1) + 
	scale_color_manual(values = cols_anc, name = "Selection treatment") +
	geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = trait_vectors1, color = "black",
				 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
	geom_text(aes(x = PC1, y = PC2, label = trait), data = trait_vectors1) +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) +xlim(-4, 4) + ylim(-3, 3)
ggsave("figures/blup-pca.png", width = 8, height = 6)




mean_cent <- scale(t(newg))

pca_blup <- prcomp(newg)
pca_blup <- prcomp(mean_cent)
pca_blup1 <- rda(mean_cent)

summary(pca_blup1)


loadings <- scores(pca_blup1, choices = c(1,2)) 


	pcs <- as.data.frame(loadings$species) %>% 
		mutate(population = newg1$population) %>% 
		left_join(., treatments) %>% 
		mutate(treatment = factor(treatment,
								  levels=c("Ancestors", "C", "L", "N",
								  		 "P", "B", "S", "BS"))) 
		
	
	pcs_envs <- as.data.frame(loadings$sites) %>% 
		mutate(environment = rownames(.))
	
	(fit <- envfit(pca_blup1, all_sizes3, perm = 999))
	trait_vectors1 <- as.data.frame(fit$vectors$arrows) %>% 
		mutate(trait = rownames(.))
	
ggplot() +
	 
		geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = pcs_envs, color = "grey",
					 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
		geom_text(aes(x = PC1, y = PC2, label = environment), data = pcs_envs, color = "black") +
	geom_point(aes(x = PC1, y = PC2, color = treatment), data = pcs, size = 2) +
		geom_hline(yintercept = 0) +
		geom_vline(xintercept = 0) +
	xlim(-3, 3) + ylim(-3.5, 2) +
		scale_color_manual(values = cols_anc, name = "Selection treatment") +
	xlab("PC1 (21.93% explained var.)") + 
	ylab("PC2 (16.05% explained var.)")
ggsave("figures/blup_pca_environment_vectors.png", width = 8, height = 6)
