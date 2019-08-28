

library(tidyverse)
library(readxl)
library(janitor)
library(cowplot)

biomass <- read_excel("data-raw/stoichiometry/Biomass_calculations_pre_post_weight_47mm_GFF.xlsx") %>% 
	clean_names()

treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


bio2 <- left_join(biomass, treatments, by = "population") %>% 
	filter(filter_number != "34") %>% View
	mutate(biomass_mg = ifelse(filter_number %in% c(33, 24), biomass_mg*4, biomass_mg)) %>% 
	filter(!is.na(population))

bio2 %>% 
	ggplot(aes(x = reorder(population, biomass_mg), y = biomass_mg, fill = treatment)) + geom_bar(stat = "identity")


bio2 %>% 
	ggplot(aes(x = treatment, y= biomass_mg, color = ancestor_id)) + geom_point()

# light <- read_csv("data-processed/all_light_consumption_absorbance.csv")
light <- read_csv("data-processed/light_consumption.csv")
biomasses <- read_csv("data-processed/biomasses-stoichiometry.csv")

light2 <- light %>% 
	left_join(., biomasses) %>% 
	mutate(biomass_l = biomass_mg*(1000/volume_filtered_47mm)) %>%
	mutate(biomass_400ml = biomass_l/1000*400) %>% 
	mutate(k = mean_consumption/biomass_400ml) %>% 
	left_join(., treatments)

write_csv(light2, "data-processed/light_consumption_k.csv")


light2 %>% 
	ggplot(aes(x = treatment, y= k, color = ancestor_id)) + geom_point() +
	ylab("Per biomass light consumption (umols/m2/s/mg)")

light_anc <- light2 %>% 
	filter(treatment == "none") %>% 
	rename(k_anc = k) %>% 
	select(ancestor_id, k_anc)

light3 <- left_join(light2, light_anc) %>% 
	mutate(change_k = k - k_anc) %>% 
	group_by(treatment) %>% 
	mutate(mean_change_k = mean(change_k)) %>% 
	mutate(se_change_k = std.error(change_k)) %>% 
	mutate(Ft = (k*biomass_400ml)/((k*biomass_400ml) +0.1298707)) %>% ## kbg = 0.1298707
	mutate(mean_ft = mean(Ft)) %>% 
	mutate(se_ft = std.error(Ft))

light3 %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_k)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_change_k), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_change_k - se_change_k, ymax = mean_change_k + se_change_k), color = "purple", width = 0.1) +
	geom_hline(yintercept = 0) + ylab("Change in k") + xlab("")
ggsave("figures/change_k_light.png", width = 6, height = 4)

ft_anc <- light3 %>% 
	filter(treatment == "none") %>%
	rename(ft_anc = Ft) %>% 
	ungroup() %>% 
	select(ancestor_id, ft_anc)

light4 <- light3 %>% 
	left_join(., ft_anc) %>% 
	mutate(change_ft = Ft - ft_anc) %>% 
	mutate(mean_change_ft = mean(change_ft)) %>% 
	mutate(se_change_ft = std.error(change_ft))
	

light4 %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_ft)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_change_ft), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_change_ft - se_change_ft, ymax = mean_change_ft + se_change_ft),
				  color = "purple", width = 0.1)  +
	geom_hline(yintercept = 0)
ggsave("figures/change-ft-light.png", width = 6, height = 4)
	