

### stoichiometry samples

library(readxl)
library(janitor)
library(tidyverse)
library(cowplot)

all_rstars <- read_csv("data-processed/all-rstars.csv") %>% 
	mutate(cons_vector = pc/nc) ## ## pc is the P:C molar ratio, nc is the N:C molar ratio

pn <- read_excel("data-raw/stoichiometry/ChlamEE-24062019-PNC.xls", sheet = "Nutriens P-N", skip = 4) %>% 
	clean_names() %>% 
	filter(!is.na(name)) %>% 
	filter(!grepl("A", name)) %>% 
	rename(sample = sample_identity) %>% 
	mutate(sample = as.character(sample)) %>% 
	filter(sample != "Drift", sample != "Wash")

pn_ids <- read_csv("data-raw/stoichiometry/ChlamEE-PN-filter-IDs.csv") %>% 
	mutate(sample = as.character(sample))## this is the id key that converts the PN sample ids from Dany's data to our populations

cn <- read_excel("data-raw/stoichiometry/ChlamEE-24062019-PNC.xls", sheet = "Nutiens C-N", skip = 7) %>% 
	clean_names() %>% 
	filter(!grepl("n", name)) %>% 
	rename(filter_number = name) %>% 
	mutate(filter_number = as.numeric(filter_number))

## population 34 does not exist

volumes <- read_excel("data-raw/stoichiometry/ClamEE_stoichiometry_15.05.xlsx") %>% 
	clean_names()


weights <- read_excel("data-raw/stoichiometry/Biomass_calculations_pre_post_weight_47mm_GFF.xlsx") %>% 
	clean_names() %>% 
	filter(!is.na(filter_number))

cn2 <- left_join(cn, weights) %>% 
	filter(filter_number != "34")



pn2 <- left_join(pn, pn_ids) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), n_µg_filter, p_µg_filter) %>% 
	rename(p_ug_filter = p_µg_filter,
		   n_ug_filter = n_µg_filter)

all_stoich <- left_join(cn2, pn2) %>% 
	mutate(volume_filtered_25mm = 120) %>% 
	mutate(volume_filtered_47mm = 200) %>% 
	mutate(volume_filtered_47mm = ifelse(population == "3", 50, volume_filtered_47mm)) %>% 
	mutate(volume_filtered_25mm = ifelse(population == "19", 180, volume_filtered_25mm)) 

biomasses <- all_stoich %>% 
	select(population, filter_number, biomass_mg, volume_filtered_47mm)
write_csv(biomasses, "data-processed/biomasses-stoichiometry.csv")


area_big <- pi*(3.5/2)^2
area_small <- pi*(1.8/2)^2	

1/(1/(area_big/area_small)/0.71)


all_stoich2 <- all_stoich %>% 
	mutate(n_per_l_g = n_ug_filter*1000/volume_filtered_25mm/1000/1000) %>% 
	mutate(p_per_l_g = p_ug_filter*1000/volume_filtered_25mm/1000/1000) %>% 
	mutate(n_per_l_big_g = n_mg*1000/volume_filtered_47mm*(area_big/area_small)/1000) %>% 
	mutate(c_per_l_big_g = c_mg*1000/volume_filtered_47mm*(area_big/area_small)/1000) %>% 
	mutate(n_moles_per_l = n_per_l_big_g /14.0067) %>% 
	mutate(c_moles_per_l = c_per_l_big_g /12.0107) %>% 
	mutate(p_moles_per_l = p_per_l_g /30.973762) %>% ### should this be multiplied by 2??
	mutate(nc = n_moles_per_l/c_moles_per_l) %>% 
	mutate(pc = p_moles_per_l / c_moles_per_l) %>% 
	mutate(cn = c_moles_per_l/n_moles_per_l) %>% 
	mutate(cp = c_moles_per_l/p_moles_per_l)

all_stoich3 <- all_stoich2 %>% 
	select(population, nc, pc) %>% 
	mutate(type = "new") %>% 
	rename(nc_new = nc,
		   pc_new = pc)

all_stoich4 <- left_join(all_rstars, all_stoich3) 

all_stoich4 %>% 
	ggplot(aes(x = pc, y = pc_new, color = population)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of P:C") + xlab("Old round of P:C")

all_stoich4 %>% 
	ggplot(aes(x = 1/pc, y = 1/pc_new)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of C:P") + xlab("Old round of C:P")

all_stoich4 %>% 
	ggplot(aes(x = 1/nc, y = 1/nc_new)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of C:N") + xlab("Old round of C:N")

all_stoich5 <- all_stoich4 %>% 
	rename(nc_old = nc,
		   pc_old = pc,
		   nc = nc_new,
		   pc = pc_new) %>% 
	mutate(cons_vector_new = pc/nc) %>% 
	select(-type)
	

all_stoich5 %>% 
	ggplot(aes(x = cons_vector, y = cons_vector_new)) + geom_point() +
	geom_abline(slope = 1, intercept = 0)

write_csv(all_stoich5, "data-processed/all_rstars_new_stoich.csv")

stoich_corr <- read_excel("data-raw/stoichiometry/ChlamEE_140218_stoich_growthrates_corrections.xlsx") %>% 
	clean_names() %>% 
	mutate(nc_corr = 1/cto_n_molar_corr) %>% 
	mutate(nc_uncorr = 1/cto_n_molar) %>% 
	mutate(pc_corr = 1/cto_p_molar_corr) %>% 
	mutate(pc_uncorr = 1/cto_p_molar) %>% 
	select(strain_id, starts_with("nc"), starts_with("pc")) %>% 
	rename(population = strain_id) %>% 
	mutate(population = ifelse(population == "Anc 2", "anc2", population)) %>% 
	mutate(population = ifelse(population == "Anc 3", "anc3", population)) %>% 
	mutate(population = ifelse(population == "Anc 4", "anc4", population)) %>% 
	mutate(population = ifelse(population == "Anc 5", "anc5", population))
	

stoich <- left_join(all_stoich4, stoich_corr) %>% 
	mutate(cons_vectors_corr = pc_corr/nc_corr) %>% 
	mutate(cons_vectors_uncorr = pc_uncorr/nc_uncorr) %>% 
	mutate(cons_vectors_new = pc_new/nc_new)

stoich %>% 
	ggplot(aes(x = pc_corr, y = pc_new)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of P:C") + xlab("Old round of P:C")

stoich %>% 
	ggplot(aes(x = pc_corr, y = pc_uncorr)) + geom_point() + 
	geom_abline(slope = 2.64, intercept = 0) +
	ylab("New round of P:C") + xlab("Old round of P:C")

stoich %>% 
	ggplot(aes(x = nc_corr, y = nc_uncorr)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of N:C") + xlab("Old round of N:C")

stoich %>% 
	ggplot(aes(x = cons_vectors_corr, y = cons_vectors_new)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of CVs") + xlab("Old round of CVs (corrected)")

stoich %>% 
	ggplot(aes(x = cons_vectors_uncorr, y = cons_vectors_new)) + geom_point() + 
	geom_abline(slope = 1, intercept = 0) +
	ylab("New round of CVs") + xlab("Old round of CVs (uncorrected)")

lm(pc_uncorr ~ pc_corr, data = stoich) %>% summary()

stoich %>% 
	mutate(diff = pc_uncorr/pc_corr) %>% View

all_stoich2 %>% 
	ggplot(aes(x = population, y = p_moles_per_l)) + geom_point()


stoich %>% 
	ggplot(aes(x = treatment, y = nc_new, color = treatment)) + geom_boxplot() + ylab("N:C")

stoich %>% 
	ggplot(aes(x = treatment, y = pc_new, color = treatment)) + geom_boxplot() + ylab("P:C")
