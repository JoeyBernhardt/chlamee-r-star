


library(tidyverse)
library(cowplot)
library(here)

treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")
m <- 0.56
all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	# read_csv("data-processed/all-rstars.csv") %>% 
	mutate(cons_vector = pc/nc) %>% 
	mutate(n_star = n_ks*m/(n_umax-m)) %>%  
	mutate(p_star = p_ks*m/(p_umax-m))
		   	## ## pc is the P:C molar ratio, nc is the N:C molar ratio
all_combos_sympatry <- read_csv("data-processed/chlamee-unique-combos2.csv") %>% 
	mutate(competition_type = "sympatry") %>% 
	mutate(comp1 = as.character(comp1))%>% 
	mutate(combination = as.character(combination))
## all the sympatric combos
all_combos_allopatry <- read_csv("data-processed/all_combos_allopatry.csv") %>% 
	mutate(competition_type = "allopatry") %>% 
	mutate(comp1 = as.character(comp1)) %>% 
	mutate(combination = as.character(combination))
## all the allopatric combos
all_combos_before_after <- read_csv("data-processed/all-before-after-combos.csv") %>% 
	mutate(competition_type = "same_ancestor_before_after") %>% 
	mutate(combination = rownames(.)) %>% 
	mutate(comp1 = as.character(comp1)) %>% 
	mutate(combination = as.character(combination))

names(all_combos_sympatry)
names(all_combos_allopatry)
names(all_combos_before_after)

all_comparisons <- bind_rows(all_combos_allopatry, all_combos_sympatry, all_combos_before_after) %>% 
	mutate(unique_combination = rownames(.))

View(all_comparisons)
### conditions for coexistence
## 1. zngis cross
## 2. populations consume more of the resource that limits them most

### this function tells us if there is a possibility for coexistence at at least one supply point

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

## find the coexistence outcomes for all 
outcomes_all_together <- all_comparisons %>% 
	split(.$unique_combination) %>% 
	map_df(coexist_or_not, .id = "unique_combination") %>%
	left_join(., all_comparisons) %>% 
	mutate(competition_outcome = case_when(zngis_cross == TRUE & steeper_zngi == TRUE ~ "stable coexistence",
										   zngis_cross == TRUE & steeper_zngi == FALSE ~ "unstable coexistence",
										   zngis_cross == FALSE ~ "exclusion"))

### ok now determine the outcomes!

before_after <- outcomes_all_together %>% 
	filter(competition_type == "same_ancestor_before_after") %>% 
	mutate(ancestral1 = ifelse(grepl("c", comp1), "ancestral", "descendant")) %>% 
	mutate(ancestral2 = ifelse(grepl("c", comp2), "ancestral", "descendant")) %>% 
	mutate(ancestry = case_when(ancestral1 == "ancestral" & ancestral2 == "ancestral" ~ "1. Before selection (ancestors)",
								ancestral1 == "descendant" & ancestral2 == "descendant" ~ "2. After selection in different envts",
								TRUE ~ "mixed")) %>% 
	filter(ancestry != "mixed")

before_after %>% 
	group_by(ancestry, competition_outcome) %>% 
	tally() %>% 
	group_by(ancestry) %>% 
	mutate(percentage = 100*n/sum(n)) %>% 
	ggplot(aes(x = ancestry, y = percentage, fill = competition_outcome, label = n)) + geom_bar(stat = "identity") +
	scale_fill_viridis_d(name = "Competition outcome") + geom_text() + xlab("")
ggsave("figures/before_after_selection_competition_outcomes.png", width = 10, height = 6)
	

before_after %>% 
	group_by(ancestry, ancestor_combo, competition_outcome) %>% 
	tally() %>% View


before_after %>% 
	filter(ancestor_combo == "anc3_4") %>% View
	ggplot(aes(x = ancestry, y = competition_outcome)) + geom_point() + 
	facet_wrap( ~ ancestor_combo)

	
pop1 <- before_after %>% 
	select(unique_combination, pop1) %>% 
	left_join(., treatments, by = c("pop1" = "population")) %>% 
	dplyr::rename(treatment_pop1 = treatment) %>% 
	dplyr::rename(ancestor_id1 = ancestor_id) %>% 
	select(-note)
pop2 <- before_after %>% 
	select(unique_combination, pop2) %>% 
	left_join(., treatments, by = c("pop2" = "population")) %>% 
	dplyr::rename(treatment_pop2 = treatment) %>% 
	dplyr::rename(ancestor_id2 = ancestor_id) %>% 
	select(-note)

before_after %>% 
	left_join(., pop1) %>% 
	left_join(., pop2) %>% 
	mutate(both_treatments = paste(treatment_pop1, treatment_pop2, sep = "_")) %>% 
ggplot(aes(x = ancestry, y = competition_outcome, color = both_treatments)) + geom_jitter() + 
	facet_wrap( ~ ancestor_combo, nrow = 2)
ggsave("figures/before_after_selection_treatments.png", width = 20, height = 10)

before_after %>% 
	left_join(., pop1) %>% 
	left_join(., pop2) %>% 
	mutate(both_treatments = paste(treatment_pop1, treatment_pop2, sep = "_")) %>% 
	ggplot(aes(x = ancestry, y = competition_outcome, label = unique_combination, color = competition_outcome)) + geom_text(position=position_jitter(width=.5,height=.5)) + 
	facet_wrap( ~ ancestor_combo, nrow = 2)
ggsave("figures/before_after_selection_treatments_numbers.png", width = 20, height = 10)

ba_treatments <- before_after %>% 
	left_join(., pop1) %>% 
	left_join(., pop2) %>% 
	mutate(both_treatments = paste(treatment_pop1, treatment_pop2, sep = "_")) 

ba_anc <- before_after %>% 
	left_join(., pop1) %>% 
	left_join(., pop2) %>% 
	mutate(both_treatments = paste(treatment_pop1, treatment_pop2, sep = "_")) %>% 
	filter(ancestry == "1. Before selection (ancestors)")

ba_treatments %>% 
	ggplot(aes(x = ancestry, y = competition_outcome, color = treatment_pop1, fill = treatment_pop2)) + geom_point(shape = 21, position=position_jitter(width=.5,height=.3), stroke = 1.5) + 
	# geom_point(aes(x = ancestry, y = competition_outcome, color = treatment_pop1, fill = treatment_pop2),data = ba_anc) +
	facet_wrap( ~ ancestor_combo, nrow = 2)
ggsave("figures/before_after_selection_treatments_colours_new.png", width = 30, height = 15)


ba_treatments2 <- ba_treatments %>% 
	mutate(same_envt = ifelse(treatment_pop1 == treatment_pop2, "same_envt", "diff_envt"))

ba_treatments2 %>% 
	group_by(ancestry, competition_outcome, same_envt) %>% 
	tally() %>% 
	group_by(ancestry) %>% 
	mutate(percentage = 100*n/sum(n)) %>% View
	ggplot(aes(x = ancestry, y = percentage, fill = competition_outcome, label = n)) + geom_bar(stat = "identity") +
	scale_fill_viridis_d(name = "Competition outcome") + geom_text() + xlab("")
		   	
	
percentages <-	ba_treatments2 %>% 
		group_by(ancestry, competition_outcome, same_envt) %>% 
		tally() %>% 
		group_by(ancestry) %>% 
		mutate(percentage = 100*n/sum(n)) 

sum(percentages$n)

### does the breakdown of competitive outcomes differ depending on whether they were selected in the same or different envts?

total_diff <- 374+193+131
total_same <- 29+19+11

11/total_same*100
131/total_diff*100

ba_treatments3 <- ba_treatments %>% 
	mutate(same_envt = ifelse(treatment_pop1 == treatment_pop2, "same_envt", "diff_envt")) %>% 
	mutate(resource_limitation = ifelse(both_treatments %in% c("N_N", "P_P", "L_L"), "both_resource_limitation", "other")) %>% 
	group_by(ancestry, competition_outcome, same_envt, resource_limitation) %>% 
	tally() %>% 
	group_by(ancestry) %>% 
	mutate(percentage = 100*n/sum(n)) %>% View


treatment_vec <- c("N", "C", "BS", "S", "C", "L", "none", "P")


jitters_h <- runif(min = 0.001, max = 0.5, n = nrow(ba_treatments))
jitters_v <- runif(min = 0.001, max = 0.5, n = nrow(ba_treatments))

pies <- ba_treatments %>% 
	select(unique_combination, treatment_pop1, treatment_pop2, ancestry, competition_outcome, ancestor_combo) %>% 
	mutate(N = ifelse(grepl("N",treatment_pop2) | grepl("N",treatment_pop1), .50, 0)) %>% 
	mutate(P = ifelse(grepl("P",treatment_pop2) | grepl("P",treatment_pop1), .50, 0)) %>% 
	mutate(B = ifelse(treatment_pop2 == "B" | treatment_pop1 == "B", .50, 0)) %>% 
	mutate(BS = ifelse(treatment_pop2 == "BS" | treatment_pop1 == "BS", .50, 0)) %>% 
	mutate(C = ifelse(treatment_pop2 == "C" | treatment_pop1 == "C", .50, 0)) %>% 
	mutate(L = ifelse(treatment_pop2 == "L" | treatment_pop1 == "L", .50, 0)) %>%
	mutate(none = ifelse(treatment_pop2 == "none" | treatment_pop1 == "none", .50, 0)) %>% 
	mutate(S = ifelse(treatment_pop2 == "S" | treatment_pop1 == "S", .50, 0)) %>% 
	mutate(comp_out = as.numeric(as.factor(competition_outcome))) %>%
	mutate(anc_out = as.numeric(as.factor(ancestry))) %>% 
	mutate(jitters_h = jitters_h) %>% 
	mutate(jitters_v = jitters_v) %>% 
	mutate(anc_out2 = anc_out + jitters_h) %>% 
	mutate(comp_out2 = comp_out + jitters_v) %>% 
	select(6:20, unique_combination) %>% 
	select(anc_out, comp_out, everything()) 
	

colors <- names(pies)[4:11]
library(scatterpie)
str(pies)

cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")
	ggplot() +
	geom_scatterpie(aes(x = anc_out2, y = comp_out2), cols = colors, data = pies, color = NA, size = 3)  +
		coord_equal() + facet_wrap( ~ ancestor_combo, nrow  =2) + 
		scale_color_manual(values = cols_no_anc, name = "Selection environment") +
		# scale_fill_brewer(name = "Treatment", type = "qual") +
		xlab("Ancestry") + ylab("Competition outcome") +
		theme(axis.title.x=element_blank(),
			  axis.text.x=element_blank(),
			  axis.ticks.x=element_blank()) +
		theme(axis.title.y=element_blank(),
			  axis.text.y=element_blank(),
			  axis.ticks.y=element_blank())
		
ggsave("figures/before_after_selection_treatments_colours_mini_pies.png", width = 15, height = 8)


	


### all the sympatry cases

outcomes_all_together %>% 
	filter(competition_type == "sympatry") %>% 
	group_by(treatment, competition_outcome) %>%
	tally() %>% 
	mutate(percentage = 100*n/sum(n)) %>% 
	ggplot(aes(x = treatment, y = percentage, fill = competition_outcome, label = n)) + geom_bar(stat = "identity") +
	scale_fill_viridis_d() + geom_text()
ggsave("figures/sympatry_competition_outcomes.png", width = 10, height = 6)

symp_outcomes <- outcomes_all_together %>% 
	filter(competition_type == "sympatry") %>% 
	filter(treatment != "none") %>% 
	group_by(treatment, competition_outcome) %>%
	tally() %>% 
	mutate(percentage = 100*n/sum(n)) %>% 
	group_by(competition_outcome) %>% 
	summarise(percentage = mean(percentage)) %>% 
	mutate(ancestry = "3. After selection in the same environment")


bf_outcomes <- before_after %>% 
	group_by(ancestry, competition_outcome) %>% 
	tally() %>% 
	group_by(ancestry) %>% 
	mutate(percentage = 100*n/sum(n))

all_comp_outcomes <- bind_rows(bf_outcomes, symp_outcomes)

all_comp_outcomes %>% 
	ggplot(aes(x = ancestry, y = percentage, group = competition_outcome, fill = competition_outcome, label = n)) + geom_bar(stat = "identity") +
	scale_fill_viridis_d()
ggsave("figures/all_competition_outcomes_barchart.png", width = 14, height = 6)
	


# ## find all the coexistence outcomes for the allopatric cases
# outcomes_allopatry <- all_combos_allopatry %>% 
# 	split(.$combination) %>% 
# 	map_df(coexist_or_not, .id = "combination") ## go through all the combos, apply the coexist_or_not function
# 
# ## for the sympatry cases, find the combos that are just ancestors
# ancestral_combos <- all_combos_sympatry %>% 
# 	filter(treatment == "none")
# 
# ## find all the coexistence outcomes for the sympatric cases cases
# outcomes_sympatry <- all_combos_sympatry %>% 
# 	split(.$combination) %>% 
# 	map_df(coexist_or_not, .id = "combination") %>% 
# 	mutate(state = ifelse(combination %in% c(ancestral_combos$combination), "ancestral", "descendant")) %>%
# 	mutate(combination = as.numeric(combination)) %>% 
# 	left_join(all_combos_sympatry, by = "combination")
# 
# 
# outcomes_symp <- outcomes_sympatry %>% 
# 	group_by(state, coexist) %>% 
# 	tally()
# 	
# symp_percent_coexist_descendent <- 100*(outcomes_symp$n[outcomes_symp$state == "descendant" & outcomes_symp$coexist == TRUE]) / (outcomes_symp$n[outcomes_symp$state == "descendant" & outcomes_symp$coexist == TRUE] +
# 	outcomes_symp$n[outcomes_symp$state == "descendant" & outcomes_symp$coexist == FALSE])
# 
# symp_percent_coexist_ancestral <- 100*(outcomes_symp$n[outcomes_symp$state == "ancestral" & outcomes_symp$coexist == TRUE]) / (outcomes_symp$n[outcomes_symp$state == "ancestral" & outcomes_symp$coexist == TRUE] +
# 																																outcomes_symp$n[outcomes_symp$state == "ancestral" & outcomes_symp$coexist == FALSE])
# allopatry_percent_coexist <- 100*(length(outcomes_allopatry$coexist[outcomes_allopatry$coexist == TRUE]))/(length(outcomes_allopatry$coexist))
# 																															   	
# allo_long <- outcomes_allopatry %>% 
# 	 	gather(key = competitor, value = population, pop1, pop2) %>% 
# 	left_join(., all_rstars)

write_csv(outcomes_allopatry, "data-processed/tilman-competition-outcomes-allopatry.csv")

allo_long <- outcomes_all_together %>% 
	filter(competition_type == "allopatry") %>% 
	gather(key = competitor, value = population, pop1, pop2) %>% 
	left_join(., all_rstars, by = "population")

symp_long <- outcomes_all_together %>% 
	filter(competition_type == "sympatry") %>% 
	gather(key = competitor, value = population, pop1, pop2) %>% 
	left_join(., all_rstars, by = "population")

## plot it! 
allo_long %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	mutate(coexist2 = ifelse(coexist == TRUE, "coexistence", "exclusion")) %>%
	ggplot(aes(x = n_star, y = p_star, color =competitor, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = competitor, x = max_n_star, y = max_p_star, xend = max_n_star + nc*9, yend = max_p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_text(aes(x = 3.5, y = 0.75),
			  size = 3, 
			  fontface = "bold", color = "black") +
	facet_wrap( ~ unique_combination) + 
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")  
ggsave("figures/comp_outcome_tilman_new_stoich_unique.png", width = 20, height = 20)



allo_long %>% 
	filter(competition_outcome == "stable coexistence") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	mutate(coexist2 = ifelse(coexist == TRUE, "coexistence", "exclusion")) %>%
	ggplot(aes(x = n_star, y = p_star, color =competitor, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = competitor, x = max_n_star, y = max_p_star, xend = max_n_star + nc*9, yend = max_p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_text(aes(x = 3.5, y = 0.75),
			  size = 3, 
			  fontface = "bold", color = "black") +
	facet_wrap( ~ unique_combination) + 
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")  
ggsave("figures/comp_outcome_tilman_new_stoich_unique_coexist.png", width = 20, height = 20)



### poster figures

allo_long %>% 
	filter(competition_outcome == "stable coexistence") %>% 
	filter(unique_combination %in% c("76", "83")) %>% View


	
ba2 <- ba_treatments %>% 
	gather(key = competitor, value = population, pop1, pop2) %>% 
	left_join(., all_rstars, by = "population") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star))


(ba2$unique_combination[ba2$combination == "76"])
symp_snip <- symp_long %>% 
	filter(treatment.x %in% c("none")) %>% 
	filter(ancestor_id.y == c("anc3")) %>% 
	filter(unique_combination == 120) 

cols <- c("#CC5A9F", "#FF0000", "#2F8C55", "#E9EF28", "#26549C", "#333333", "#DB3C01", "#01A0C6")


# allo_plot <-
	allo_long %>% 
	filter(unique_combination == 29) %>% 
	filter(competition_outcome == "stable coexistence") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment.y, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = treatment.y, x = max_n_star, y = max_p_star, xend = max_n_star + nc, yend = max_p_star + pc),
				 size = 1, linetype = 2) +
		geom_segment(aes(color = treatment.y, x = max_n_star, y = max_p_star, xend = max_n_star - nc, yend = max_p_star - pc),
					 size = 1, linetype = 2,
					 arrow = arrow(type = "closed", length = unit(0.09, "inches"))) +
	geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 2),
				 size = 1.5, linetype = 1) +
	geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = 14, yend = p_star),
				 size = 1.5, linetype = 1) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star + nc, yend = p_star + pc),
				 size = 1, linetype = 2) +
		geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star - nc, yend = p_star - pc),
					 size = 1, linetype = 2,
					 arrow = arrow(type = "closed", length = unit(0.09, "inches"))) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 2),
				 size = 1.5, linetype = 1) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = 14, yend = p_star),
				 size = 1.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	# geom_point(aes(x = 1.75, y = 0.096), color = "black") + 
	# facet_wrap( ~ unique_combination) + 
	scale_color_manual(values = c(cols[2], "grey", cols[7])) +
 	# ylim(0.06, 0.15) + xlim(1.2, 1.8) +
	coord_cartesian() +
	theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 1.5)  
	ggsave("figures/tilman-combination-colors-unique-combination29.pdf", width = 6, height = 5)
	
	
	cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")
	
	allo_long %>% 
		filter(unique_combination == 29) %>% 
		filter(competition_outcome == "stable coexistence") %>% 
		group_by(combination) %>% 
		mutate(max_n_star = max(n_star)) %>% 
		mutate(max_p_star = max(p_star)) %>% 
		ggplot(aes(x = n_star, y = p_star, color = treatment.y)) + geom_point() +
		geom_segment(aes(color = treatment.y, x = max_n_star, y = max_p_star, xend = max_n_star + nc*7, yend = max_p_star + pc*7),
					 size = 1, linetype = 2) +
		geom_segment(aes(color = treatment.y, x = max_n_star, y = max_p_star, xend = max_n_star - nc*7, yend = max_p_star - pc*7),
					 size = 1, linetype = 2,
					 arrow = arrow(type = "closed", length = unit(0.09, "inches"))) +
		geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 1.2),
					 size = 1, linetype = 1) +
		geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = 13, yend = p_star),
					 size = 1, linetype = 1) +
		geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star + nc*7, yend = p_star + pc*7),
					 size = 1, linetype = 2) +
		geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star - nc*7, yend = p_star - pc*7),
					 size = 1, linetype = 2,
					 arrow = arrow(type = "closed", length = unit(0.09, "inches"))) +
		geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 1.2),
					 size = 1, linetype = 1) +
		geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = 13, yend = p_star),
					 size = 1, linetype = 1) +
		ylab("P (uM)") + xlab("N (uM)") +
		scale_color_manual(values = c(cols_no_anc[2], "black", cols_no_anc[6])) +
		# ylim(0.5, 2) + xlim(1.2, 1.8) +
		coord_cartesian() +
		theme( 
			plot.margin = unit(c(1,1,1,1), "lines"),
			axis.text = element_text(size=14),
			axis.title=element_text(size=14),
			rect = element_rect(fill = "transparent"),
			legend.position = "none") +
		panel_border(colour = "black", linetype = "solid", size = 0.75)
	ggsave("figures/tilman-combination-colors-unique-combination29.png", width = 6, height = 5)
	
#### chesson conversions
	
	allo_snip <- allo_long %>% 
		filter(unique_combination == 29) 
	
	symp_snip <- symp_long %>% 
		filter(treatment.x %in% c("none")) %>% 
		filter(ancestor_id.y == c("anc3")) %>% 
		filter(unique_combination == 120) 
	
	all_combos_poster <- bind_rows(allo_snip, symp_snip)
	write_csv(all_combos_poster, "data-processed/all_combos_poster.csv")
	
	poster_combos <- combn(all_combos_poster$population, 2) %>% 
		t() %>% 
		as.data.frame() %>% 
		rename(comp1 = V1,
			   comp2 = V2) %>% 
		mutate(combo = rownames(.))
	
	
symp_plot <- symp_long %>% 
		filter(treatment.x %in% c("none")) %>% 
		filter(ancestor_id.y == c("anc3")) %>% 
		filter(unique_combination == 120) %>% 
		# filter(competition_outcome == "stable coexistence") %>% 
		group_by(combination) %>% 
		mutate(max_n_star = max(n_star)) %>% 
		mutate(max_p_star = max(p_star)) %>% 
		mutate(coexist2 = ifelse(coexist == TRUE, "coexistence", "exclusion")) %>%
		ggplot(aes(x = n_star, y = p_star, color = ancestor_id.y, label = competition_outcome)) + geom_point() +
		geom_segment(aes(color = ancestor_id.y, x = n_star, y = p_star, xend = n_star + nc, yend = p_star + pc),
					 size = 0.5, linetype = 2) +
		geom_segment(aes(color = ancestor_id.y, x = n_star, y = p_star, xend = n_star, yend = 0.15),
					 size = 0.5, linetype = 1) +
		geom_segment(aes(color = ancestor_id.y, x = n_star, y = p_star, xend = 1.8, yend = p_star),
					 size = 0.5, linetype = 1) +
		ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = 1.76, y = 0.097), color = "black") + 
		# facet_wrap( ~ unique_combination) + 
		#ylim(0.05, 0.15) + xlim(1.2, 2)
	ylim(0.06, 0.15) + xlim(1.2, 1.8) +
		coord_cartesian() +
		theme( 
			plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
			axis.text = element_text(size=13),
			axis.title=element_text(size=14)) +
		panel_border(colour = "black") 
	ggsave("figures/Tilman-plot-unique-combination-120.pdf", width = 8, height = 6)
	
	
## ok here let's just plot the ancestors and descendants	
	
	multi_plot <- plot_grid(allo_plot, symp_plot, align = "h", nrow = 1)
	save_plot("figures/tilman-plots-combo.pdf", multi_plot,
			  ncol = 2, # we're saving a grid plot of 2 columns
			  nrow = 1, # and 2 rows
			  # each individual subplot should have an aspect ratio of 1.3
			  base_aspect_ratio = 1.2
	)
	
### anc 3 and 4 could coexist in their ancestral states

anc34 <- ba_treatments %>% 
	filter(ancestor_combo == "anc3_4") %>%
	gather(key = competitor, value = population, pop1, pop2) %>% 
	left_join(., all_rstars, by = "population")

View(anc34)

anc34 %>% 
	filter(both_treatments == "none_none") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, label = competition_outcome)) + geom_point() +
	geom_segment(aes(linetype = ancestor_id.y, x = max_n_star, y = max_p_star, xend = max_n_star + nc*9, yend = max_p_star + pc*9),
				 size = 0.5) +
	geom_segment(aes(linetype = ancestor_id.y, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5) +
	geom_segment(aes(linetype = ancestor_id.y, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5) +
	ylab("P (uM)") + xlab("N (uM)") +
	scale_linetype_manual(values = c("solid", "dotted"), name = "Ancestor") +
	coord_cartesian() +
	# scale_color_manual(values = c("orange", "deepskyblue4")) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")

ggsave("figures/anc34-before-after-black.png", width = 5, height = 3)


anc34 %>% 
	filter(both_treatments == "P_P") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color =competitor, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	coord_cartesian() +
	scale_color_manual(values = c("orange", "deepskyblue4")) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")
ggsave("figures/anc34-before-after-lowP.png", width = 5, height = 3)

anc34 %>% 
	filter(both_treatments == "N_N") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color =competitor, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	coord_cartesian() +
	scale_color_manual(values = c("orange", "deepskyblue4")) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")
ggsave("figures/anc34-before-after-lowN.png", width = 5, height = 3)

anc34 %>% 
	filter(both_treatments == "N_N") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color =competitor, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star + nc*9, yend = p_star + pc*9),
				 size = 0.5, linetype = 2) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5, linetype = 1) +
	geom_segment(aes(color = competitor, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	coord_cartesian() +
	scale_color_manual(values = c("orange", "deepskyblue4")) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")
ggsave("figures/anc34-before-after-lowN.png", width = 5, height = 3)

anc34 %>% 
	filter(both_treatments %in% c("N_P", "P_N", "P_P", "N_N", "none_none")) %>% 
	# filter(unique_combination %in% c(756, 774, ))
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment.y, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = n_star + nc*9, yend = p_star + pc*9, linetype = ancestor_id.y),
				 size = 0.5) +
	geom_segment(aes(linetype = ancestor_id.y, color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 1),
				 size = 0.5) +
	geom_segment(aes(linetype = ancestor_id.y, color =treatment.y, x = n_star, y = p_star, xend = 4, yend = p_star),
				 size = 0.5) +
	ylab("P (uM)") + xlab("N (uM)") +
	facet_wrap( ~ unique_combination) +
	coord_cartesian() +
	scale_color_manual(values = c("orange", "black", "deepskyblue4"), name = "Selection environment") +
	scale_linetype_manual(values = c("solid", "dotted"), name = "Ancestor") +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black")  +
	theme(strip.background = element_blank(),
		  strip.text.y = element_blank())
ggsave("figures/anc34-before-after-lowNP_PP_NN.png", width = 8, height = 6)

