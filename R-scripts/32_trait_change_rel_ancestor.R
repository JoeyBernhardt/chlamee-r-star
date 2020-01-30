


library(tidyverse)
library(cowplot)
library(plotrix)

theme_set(theme_cowplot())


cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")
cols_no_anc2  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")

### Figure 2 trait changes relative to ancestors

salt_tol_CI <- read_csv("data-processed/salt-tol-CI-boot-12.csv") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal"))


salt_tol_CI_sum <- salt_tol_CI %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_salt_tol_mean) %>% 
	mutate(change_salt_tol_upper = mean + std.error) %>% 
	mutate(change_salt_tol_lower = mean - std.error) %>% 
	dplyr::rename(change_salt_tol_mean = mean) %>% 
	mutate(diversity = "treatment_average")

salt_tol_CI2 <- bind_rows(salt_tol_CI, salt_tol_CI_sum)

rstar_CI_p <- read_csv("data-processed/change-pstar-monod-boot-0.56.csv")%>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal"))

rstar_CI_p_sum <- rstar_CI_p %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_rstar_mean) %>% 
	mutate(change_rstar_upper = mean + std.error) %>% 
	mutate(change_rstar_lower = mean - std.error) %>% 
	dplyr::rename(change_rstar_mean = mean) %>% 
	mutate(diversity = "treatment_average")

rstar_CI_p2 <- bind_rows(rstar_CI_p, rstar_CI_p_sum)


rstar_CI_n <- read_csv("data-processed/change-nstar-monod-boot-0.56.csv")%>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal"))


rstar_CI_n_sum <- rstar_CI_n %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_rstar_mean) %>% 
	mutate(change_rstar_upper = mean + std.error) %>% 
	mutate(change_rstar_lower = mean - std.error) %>% 
	dplyr::rename(change_rstar_mean = mean) %>% 
	mutate(diversity = "treatment_average")

rstar_CI_n2 <- bind_rows(rstar_CI_n, rstar_CI_n_sum)


rstar_CI_i <- read_csv("data-processed/change-istar-monod-boot-0.56.csv")%>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) 

rstar_CI_i_sum <- rstar_CI_i %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_rstar_mean) %>% 
	mutate(change_rstar_upper = mean + std.error) %>% 
	mutate(change_rstar_lower = mean - std.error) %>% 
	dplyr::rename(change_rstar_mean = mean) %>% 
	mutate(diversity = "treatment_average")

rstar_CI_i2 <- bind_rows(rstar_CI_i, rstar_CI_i_sum)


### now ask about changes in consumption vectors
all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	mutate(np = 1/cons_vector_new)


pc_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, pc) %>% 
	dplyr::rename(pc_anc = pc)

pcs <- all_rstars %>% 
	left_join(., pc_anc) %>% 
	mutate(change_pc = pc - pc_anc)%>% 
	group_by(treatment) %>% 
	mutate(mean_change_pc = mean(change_pc)) %>% 
	mutate(se_change_pc = std.error(change_pc)) %>% 
	mutate(mean_pc = mean(pc)) %>% 
	mutate(se_pc = std.error(pc)) 


pcs2 <- pcs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	dplyr::rename(change_pc_mean = change_pc)

pcs_sum <- pcs2 %>% 
	group_by(treatment) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_pc_mean) %>% 
	mutate(change_pc_upper = mean + std.error) %>% 
	mutate(change_pc_lower = mean - std.error) %>% 
	dplyr::rename(change_pc_mean = mean) %>% 
	mutate(diversity = "treatment_average")

pcs3 <- bind_rows(pcs2, pcs_sum) %>% 
	mutate(change_pc_lower = ifelse(is.na(change_pc_lower), 0, change_pc_lower)) %>% 
	mutate(change_pc_upper = ifelse(is.na(change_pc_upper), 0, change_pc_upper))


consvec_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, cons_vector_new) %>% 
	dplyr::rename(consvec_anc = cons_vector_new)

consvecs <- all_rstars %>% 
	left_join(., consvec_anc) %>% 
	mutate(change_consvec = cons_vector_new - consvec_anc) %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	dplyr::rename(change_consvec_mean = change_consvec)

consvecs_sum <- consvecs %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_consvec_mean) %>% 
	mutate(change_consvec_upper = mean + std.error) %>% 
	mutate(change_consvec_lower = mean - std.error) %>% 
	mutate(change_consvec_mean = mean) %>% 
	mutate(diversity = "treatment_average")

consvecs2 <- bind_rows(consvecs, consvecs_sum) %>% 
	mutate(change_consvec_lower = ifelse(is.na(change_consvec_lower), 0, change_consvec_lower)) %>% 
	mutate(change_consvec_upper = ifelse(is.na(change_consvec_upper), 0, change_consvec_upper))

# start plotting here -----------------------------------------------------



plot_pn <- consvecs2 %>%
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	ggplot(aes(x = treatment, y = change_consvec_mean, color = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_consvec_mean, ymin = change_consvec_lower, ymax = change_consvec_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0, seed = NA)) +
	geom_hline(yintercept = 0) +
	ylab("Change in P:N") +
	xlab("") + xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none") +
	scale_size(range = c(0.2, 0.7)) +
	scale_shape_manual(values = c(19, 25, 19))
ggsave("figures/change-pn.png", width = 6, height = 4)

plot5 <- pcs3 %>%
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	ggplot(aes(x = treatment, y = change_pc_mean, color = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_pc_mean, ymin = change_pc_lower, ymax = change_pc_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0, seed = NA)) +
	geom_hline(yintercept = 0) +
	ylab("Change in P:C") +
	# xlab("") +
	xlab("") + 
	# scale_color_manual(values = c("orange", "deepskyblue4", "lightblue")) +
	# theme(legend.position="none") + scale_size(range = c(0.2, 0.7))
scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none") + scale_size(range = c(0.2, 0.7)) + scale_shape_manual(values = c(19, 25, 19))
	ggsave("figures/change-pc.png", width = 6, height = 4)

pcs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = pc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_pc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_pc - se_pc, ymax = mean_pc + se_pc), color = "purple", width = 0.1) +
	ylab("P:C") + xlab("")
ggsave("figures/pc.png", width = 6, height = 4)

### now N
nc_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, nc) %>% 
	dplyr::rename(nc_anc = nc)

ncs <- all_rstars %>% 
	left_join(., nc_anc) %>% 
	mutate(change_nc = nc - nc_anc)%>% 
	group_by(treatment) %>% 
	mutate(mean_change_nc = mean(change_nc)) %>% 
	mutate(se_change_nc = std.error(change_nc)) %>% 
	mutate(mean_nc = mean(nc)) %>% 
	mutate(se_nc = std.error(nc)) 



ncs2 <- ncs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	dplyr::rename(change_nc_mean = change_nc)

ncs_sum <- ncs2 %>% 
	group_by(treatment) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_nc_mean) %>% 
	mutate(change_nc_upper = mean + std.error) %>% 
	mutate(change_nc_lower = mean - std.error) %>% 
	dplyr::rename(change_nc_mean = mean) %>% 
	mutate(diversity = "treatment_average")

ncs3 <- bind_rows(ncs2, ncs_sum) %>% 
	mutate(change_nc_lower = ifelse(is.na(change_nc_lower), 0, change_nc_lower)) %>% 
	mutate(change_nc_upper = ifelse(is.na(change_nc_upper), 0, change_nc_upper))


plot6 <- ncs3 %>%
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	ggplot(aes(x = treatment, y = change_nc_mean, color = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_nc_mean, ymin = change_nc_lower, ymax = change_nc_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0, seed = NA)) +
	geom_hline(yintercept = 0) +
	ylab("Change in N:C") +
	xlab("") + xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none") + scale_size(range = c(0.2, 0.7)) + scale_shape_manual(values = c(19, 25, 19))
ggsave("figures/change-nc.png", width = 6, height = 4)


nc_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, nc) %>% 
	dplyr::rename(nc_anc = nc)

ncs <- all_rstars %>% 
	left_join(., nc_anc) %>% 
	mutate(change_nc = nc - nc_anc) %>% 
	group_by(treatment) %>% 
	mutate(mean_change_nc = mean(change_nc)) %>% 
	mutate(se_change_nc = std.error(change_nc)) %>% 
	mutate(mean_nc = mean(nc)) %>% 
	mutate(se_nc = std.error(nc)) 



ncs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_nc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_change_nc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_change_nc - se_change_nc, ymax = mean_change_nc + se_change_nc), color = "purple", width = 0.1) + 
	geom_hline(yintercept = 0) + ylab("Change in N:C") + xlab("")
ggsave("figures/change-nc.png", width = 6, height = 4)

ncs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = nc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_nc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_nc - se_nc, ymax = mean_nc + se_nc), color = "purple", width = 0.1) +
	ylab("N:C") + xlab("")
ggsave("figures/nc.png", width = 6, height = 4)

consvec_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, cons_vector_new) %>% 
	dplyr::rename(consvec_anc = cons_vector_new)

consvecs <- all_rstars %>% 
	left_join(., consvec_anc) %>% 
	mutate(change_consvec = cons_vector_new - consvec_anc)

consvecs %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_consvec, color = ancestor_id)) + geom_point(alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	ylab("Change in consumption vector \n (P:N molar ratio)") + xlab("") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "")
ggsave("figures/change-cons-vec.png", width = 6, height = 4)





salt_tol_CI3 <- salt_tol_CI2 %>%
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2))

plot_pn <- consvecs2 %>%
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	ggplot(aes(x = treatment, y = change_consvec_mean, color = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_consvec_mean, ymin = change_consvec_lower, ymax = change_consvec_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0, seed = NA)) +
	geom_hline(yintercept = 0) +
	ylab("Change in P:N") +
	xlab("") + xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none",
		  axis.text = element_text(size=16),
		  axis.title=element_text(size=16)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(19, 25, 19))
ggsave("figures/change-cons-vecs.png", width = 6, height = 4)



# Salt plot ---------------------------------------------------------------

salt <- salt_tol_CI2 %>% 
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>%
	mutate(treatment2 = as.character(treatment)) %>% 
	group_by(treatment) %>% 
	mutate(rank = dense_rank((change_salt_tol_mean))) %>% 
	mutate(diversity = str_replace(diversity, "Genotypically", "AGenotypically")) %>% 
	mutate(treatment_order = case_when(treatment == "A" ~ "1A",
									   treatment == "C" ~ "2C",
									   treatment == "L" ~ "3L",
									   treatment == "N" ~ "4N",
									   treatment == "P" ~ "5P",
									   treatment == "B" ~ "6B",
									   treatment == "S" ~ "7S",
									   treatment == "BS" ~ "8BS")) %>% 
	mutate(unique_treatment = paste(treatment_order, rank, sep = "_")) %>% 
	mutate(unique_treatment = str_replace(unique_treatment, "anc", "danc")) 


splot <- salt %>% 
	ggplot(aes(x = treatment, y = change_salt_tol_mean, color = treatment, fill = treatment)) + 
	geom_pointrange(alpha = 1, aes(shape = diversity, size = diversity2, x = unique_treatment,
								   y = change_salt_tol_mean, ymin = change_salt_tol_lower, ymax = change_salt_tol_upper),
					position=position_dodge2(width = 0.01), data = salt) +
	geom_segment(alpha = 0.5,  aes(x=unique_treatment, xend=unique_treatment, y=0, yend=change_salt_tol_mean, color = treatment)) +
	# geom_bar(width = 0.4, alpha = 0.7, aes(shape = diversity, size = diversity2, x = unique_treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
	# 				position=position_dodge2(width = 0.5), data = ri, stat = "identity") +
	geom_hline(yintercept = 0) +
	ylab(expression("Change in salt tolerance" ~ (g ~ L^{-1}))) +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	scale_fill_manual(values = c("black", cols_no_anc)) +
	# theme(legend.position="none",
	# 	  axis.text = element_text(size=19),
	# 	  axis.title=element_text(size=19)) +
	
	theme(
		  axis.text = element_text(size=19),
		  axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
	theme(axis.title.x=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	scale_y_continuous(labels = scales::number_format(accuracy = 1)) 
ggsave("figures/salt-trait-change-95CI-sub-S-order-wide.pdf", width = 6, height = 4.5)
ggsave("figures/salt-trait-change-95CI-sub-S-order.pdf", width = 7, height = 3.95)
ggsave("figures/salt-trait-change-95CI-sub-S-order-hor.pdf", width = 4.5, height = 5.6)


salt %>% 
	filter(treatment == "A") %>% 
	ggplot(aes(x = treatment, y = change_salt_tol_mean, color = treatment, fill = treatment)) + 
	geom_point(alpha = 1, aes(shape = diversity, size = diversity2)) +
	# geom_hline(yintercept = 0) +
	ylab(expression("Change in salt tolerance" ~ (g ~ L^{-1}))) +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	scale_fill_manual(values = c("black", cols_no_anc)) +
	theme(
		axis.text = element_text(size=19),
		axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
	theme(axis.title.x=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	scale_y_continuous(labels = scales::number_format(accuracy = 1)) 
ggsave("figures/change-fig-legend.pdf", width = 4, height = 4)



plot_salt <- salt_tol_CI2 %>%
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	# filter(treatment %in% c("A", "S")) %>% 
	# filter(diversity == "Genotypically diverse") %>% 
	ggplot(aes(x = treatment, y = change_salt_tol_mean, color = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_salt_tol_mean, ymin = change_salt_tol_lower, ymax = change_salt_tol_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0.5, seed = 1)) +
	geom_hline(yintercept = 0) +
	ylab(expression("Change in salt tolerance" ~ (g ~ L^{-1}))) +
	# ylab("Change in salt tolerance (g/L)") +
	xlab("") + xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none",
		  axis.text = element_text(size=16),
		  axis.title=element_text(size=16)) + 
	scale_size(range = c(0.7, 1.2)) +
	# scale_size(range = c(0.4)) +
	# scale_shape_manual(values = c(25)) +
	scale_shape_manual(values = c(19, 25, 19)) 
ggsave("figures/salt-tol-trait-change-95CI.png", width = 6, height = 4)
ggsave("figures/salt-tol-trait-change-95CI-sub-S.png", width = 2, height = 5)

library(viridis)
library(scales)
colv <- viridis_pal()(10)
show_col(viridis_pal()(10))
show_col(cola)
cola <- cols_no_anc2

# plot_I <- 
	ri <- rstar_CI_i2 %>% 
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>%
		mutate(treatment2 = as.character(treatment)) %>% 
		group_by(treatment) %>% 
		mutate(rank = dense_rank((change_rstar_mean))) %>% 
		mutate(diversity = str_replace(diversity, "Genotypically", "AGenotypically")) %>% 
		mutate(treatment_order = case_when(treatment == "A" ~ "1A",
										   treatment == "C" ~ "2C",
										   treatment == "L" ~ "3L",
										   treatment == "N" ~ "4N",
										   treatment == "P" ~ "5P",
										   treatment == "B" ~ "6B",
										   treatment == "S" ~ "7S",
										   treatment == "BS" ~ "8BS")) %>% 
		mutate(unique_treatment = paste(treatment_order, rank, sep = "_")) %>% 
		mutate(unique_treatment = str_replace(unique_treatment, "anc", "danc")) 
	

	
	# plot_I <- 
	
# light plot --------------------------------------------------------------		
	iplot <- 	ri %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment, fill = treatment)) + 
		geom_pointrange(alpha = 1, aes(shape = diversity, size = diversity2, x = unique_treatment,
										 y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
						position=position_dodge2(width = 0.01), data = ri) +
			geom_segment(alpha = 0.5,  aes(x=unique_treatment, xend=unique_treatment, y=0, yend=change_rstar_mean, color = treatment)) +
			# geom_bar(width = 0.4, alpha = 0.7, aes(shape = diversity, size = diversity2, x = unique_treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
			# 				position=position_dodge2(width = 0.5), data = ri, stat = "identity") +
		geom_hline(yintercept = 0) +
	ylab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1}))) +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
			scale_fill_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none",
		  axis.text = element_text(size=19),
		  axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
		theme(axis.title.x=element_blank(),
			  axis.text.x =element_blank(),
			  axis.ticks.x=element_blank())
# ggsave("figures/i-star-trait-change-95CI-sub-L-order.png", width = 10, height = 4)
ggsave("figures/i-star-trait-change-95CI-sub-L-order.pdf", width = 7, height = 3.92)
ggsave("figures/i-star-trait-change-95CI-sub-L-order-hor.pdf", width = 4.5, height = 5.6)




ggplot() +
geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = reorder(unique_treatment, change_rstar_mean), y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
				position=position_dodge2(width = 0.1), data = filter(ri, treatment == "N"))





rn <- rstar_CI_n2 %>% 
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>%
	mutate(treatment2 = as.character(treatment)) %>% 
	group_by(treatment) %>% 
	mutate(rank = dense_rank((change_rstar_mean))) %>% 
	mutate(diversity = str_replace(diversity, "Genotypically", "AGenotypically")) %>% 
	mutate(treatment_order = case_when(treatment == "A" ~ "1A",
									   treatment == "C" ~ "2C",
									   treatment == "L" ~ "3L",
									   treatment == "N" ~ "4N",
									   treatment == "P" ~ "5P",
									   treatment == "B" ~ "6B",
									   treatment == "S" ~ "7S",
									   treatment == "BS" ~ "8BS")) %>% 
	mutate(unique_treatment = paste(treatment_order, rank, sep = "_")) %>% 
	mutate(unique_treatment = str_replace(unique_treatment, "anc", "danc")) 


# Nitrogen plot -----------------------------------------------------------


nplot <- rn %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment, fill = treatment)) + 
	geom_pointrange(alpha = 1, aes(shape = diversity, size = diversity2, x = unique_treatment,
								   y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
					position=position_dodge2(width = 0.01), data = rn) +
	geom_segment(alpha = 0.5,  aes(x=unique_treatment, xend=unique_treatment, y=0, yend=change_rstar_mean, color = treatment)) +
	# geom_bar(width = 0.4, alpha = 0.7, aes(shape = diversity, size = diversity2, x = unique_treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
	# 				position=position_dodge2(width = 0.5), data = ri, stat = "identity") +
	geom_hline(yintercept = 0) +
	ylab("Change in N* (uM N)") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	scale_fill_manual(values = c("black", cols_no_anc)) +
	theme(legend.position="none",
		  axis.text = element_text(size=19),
		  axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
	theme(axis.title.x = element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) 
ggsave("figures/n-star-trait-change-95CI-sub-N-order.pdf", width = 6, height = 4.5)
ggsave("figures/n-star-trait-change-95CI-sub-N-order.pdf", width = 7, height = 3.92)
ggsave("figures/n-star-trait-change-95CI-sub-N-order-hor.pdf", width = 4.5, height = 5.6)


# P plot ------------------------------------------------------------------

rp <- rstar_CI_p2 %>% 
	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>%
	mutate(treatment2 = as.character(treatment)) %>% 
	group_by(treatment) %>% 
	mutate(rank = dense_rank((change_rstar_mean))) %>% 
	mutate(diversity = str_replace(diversity, "Genotypically", "AGenotypically")) %>% 
	mutate(treatment_order = case_when(treatment == "A" ~ "1A",
									   treatment == "C" ~ "2C",
									   treatment == "L" ~ "3L",
									   treatment == "N" ~ "4N",
									   treatment == "P" ~ "5P",
									   treatment == "B" ~ "6B",
									   treatment == "S" ~ "7S",
									   treatment == "BS" ~ "8BS")) %>% 
	mutate(unique_treatment = paste(treatment_order, rank, sep = "_")) %>% 
	mutate(unique_treatment = str_replace(unique_treatment, "anc", "danc")) 


pplot <- rp %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment, fill = treatment)) + 
	geom_pointrange(alpha = 1, aes(shape = diversity, size = diversity2, x = unique_treatment,
								   y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
					position=position_dodge2(width = 0.01), data = rp) +
	geom_segment(alpha = 0.5,  aes(x=unique_treatment, xend=unique_treatment, y=0, yend=change_rstar_mean, color = treatment)) +
	# geom_bar(width = 0.4, alpha = 0.7, aes(shape = diversity, size = diversity2, x = unique_treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
	# 				position=position_dodge2(width = 0.5), data = ri, stat = "identity") +
	geom_hline(yintercept = 0) +
	ylab("Change in P* (uM P)") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	scale_fill_manual(values = c("black", cols_no_anc)) +
	theme(legend.position= "none",
		  axis.text = element_text(size=19),
		  axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
	theme(axis.title.x=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) 
ggsave("figures/p-star-trait-change-95CI-sub-P-order.pdf", width = 7, height = 3.92)
ggsave("figures/p-star-trait-change-95CI-sub-P-order-wide.pdf", width = 6, height = 4.5)
ggsave("figures/p-star-trait-change-95CI-sub-P-order-hor.pdf", width = 4.5, height = 5.6)




# all plots ---------------------------------------------------------------



allplots <- plot_grid(pplot, nplot, iplot, splot, labels = c("A", "B", "C", "D"),
					  align = "v", nrow = 2, ncol = 2, label_size = 18, label_x = 0.08, label_y = 1.02)


save_plot("figures/all-change-plots.pdf", allplots,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 2, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.4
)



# rn %>% 
# 	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
# 	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
# 					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
# 												  dodge.width = 0, seed = 1)) +
# 	geom_hline(yintercept = 0) +
# 	ylab("Change in N* (uM N)") +
# 	xlab("") +  scale_color_manual(values = c("black", cols_no_anc)) +
# 	theme(legend.position="none",
# 		  axis.text = element_text(size=16),
# 		  axis.title=element_text(size=16)) +
# 	scale_size(range = c(0.7, 1.2)) +
# 	# scale_size(range = c(0.4)) +
# 	scale_shape_manual(values = c(19, 25, 19)) 
# ggsave("figures/n-star-trait-change-95CI-sub-N.png", width = 2, height = 5)


# plot_P <- rstar_CI_p2 %>% 
# 	mutate(treatment = str_replace(treatment, "Ancestors", "A")) %>% 
# 	mutate(treatment = factor(treatment,
# 							  levels=c("A", "C", "L", "N",
# 							  		 "P", "B", "S", "BS"))) %>% 
# 	mutate(diversity = factor(diversity,
# 							  levels=c("treatment_average", "Genotypically diverse", "Isoclonal"))) %>% 
# 	mutate(diversity2 = ifelse(diversity == "treatment_average", 0.2, 0.1)) %>% 
# 	mutate(diversity2 = as.numeric(diversity2)) %>% 
# 	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
# 	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
# 					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
# 												  dodge.width = 0, seed = NA)) +
# 	geom_hline(yintercept = 0) +
# 	ylab("Change in P* (uM P)") +
# 	xlab("") +
# 	scale_color_manual(values = c("black", cols_no_anc)) +
# 	theme(legend.position="none",
# 		  axis.text = element_text(size=16),
# 		  axis.title=element_text(size=16)) + 
# 	# scale_size(range = c(0.4)) +
# 	scale_size(range = c(0.7, 1.2)) +
# 	# scale_shape_manual(values = c(25)) +
# 	scale_shape_manual(values = c(19, 25, 19)) 
# ggsave("figures/p-star-trait-change-95CI.png", width = 6, height = 4)
ggsave("figures/p-star-trait-change-95CI-sub.png", width = 2, height = 5)

plot4leg <- rstar_CI_p2 %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	mutate(diversity = ifelse(diversity == "treatment_average", "Treatment average", diversity)) %>% 
	mutate(diversity = factor(diversity,
							  levels=c("Treatment average", "Genotypically diverse", "Isoclonal"))) %>% 
	mutate(diversity2 = ifelse(diversity == "Treatment average", 0.2, 0.1)) %>% 
	mutate(diversity2 = as.numeric(diversity2)) %>% 
	filter(diversity %in% c("Genotypically diverse", "Isoclonal")) %>% 
	mutate(diversity = factor(diversity, levels = c("Isoclonal", "Genotypically diverse"))) %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
	geom_pointrange(alpha = 1, aes(shape = diversity, x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper),
					position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
												  dodge.width = 0, seed = 1)) +
	geom_hline(yintercept = 0) +
	ylab("Change in P* (uM N)") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(19, 25, 19)) +
	guides(col = guide_legend(override.aes = list(size = 1, alpha = 1), 
								nrow = 2, title.position = "left"))

legend <- get_legend(plot4leg)

# multi_plot_change <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, legend, labels = c("A", "B", "C", "D", "E", "F", "G"), align = "h")
multi_plot_change <- plot_grid(plot_P, plot_N, plot_I, plot_salt, plot_pn, rda_plot, labels = c("A", "B", "C", "D", "E", "F"), align = "h", nrow = 3, ncol = 2)

?plot_grid
save_plot("figures/trait_change_rel_anc-colors-0.56-rda-bigger.png", multi_plot_change,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.3
)

### new version -- take out pn plot.

first_row <-  plot_grid(plot_P, plot_N, plot_I, plot_salt, labels = c("A", "B", "C", "D"),
						align = "v", nrow = 2, ncol = 2, label_x = 0.23, label_y = 0.97, label_size = 16)
second_row <-  plot_grid(rda_plot, labels = c('E'), nrow = 1, label_x = 0.08, label_y = 0.98, label_size = 16)

gg_all <-  plot_grid(first_row, second_row, labels=c('', ''), ncol=1)
save_plot("figures/trait_change_rel_anc-colors-0.56-rda-even-bigger.png", gg_all,
		  # each individual subplot should have an aspect ratio of 1.3
		  # base_aspect_ratio = 1.2,
		  base_height=6.3, ncol=1, nrow=2, base_width = 9.2)


### now figure 3

rstar_CI_p2 <- rstar_CI_p %>% 
	mutate(change_pstar_lower = scale(change_rstar_lower, center = FALSE)) %>% 
	mutate(change_pstar_upper = scale(change_rstar_upper, center = FALSE)) %>% 
	mutate(change_pstar_mean = scale(change_rstar_mean, center = FALSE)) %>% 
	select(-contains("rstar"))
rstar_CI_n2 <- rstar_CI_n %>% 
	mutate(change_nstar_lower = scale(change_rstar_lower, center = FALSE)) %>% 
	mutate(change_nstar_upper = scale(change_rstar_upper, center = FALSE)) %>% 
	mutate(change_nstar_mean = scale(change_rstar_mean, center = FALSE))%>% 
	select(-contains("rstar"))

rstar_CI_i2 <- rstar_CI_i %>% 
	mutate(change_istar_lower = scale(change_rstar_lower, center = FALSE)) %>% 
	mutate(change_istar_upper = scale(change_rstar_upper, center = FALSE)) %>% 
	mutate(change_istar_mean = scale(change_rstar_mean, center = FALSE))%>% 
	select(-contains("rstar"))
salt_tol_CI2 <- salt_tol_CI %>% 
	mutate(change_salt_tol_lower = scale(change_salt_tol_lower, center = FALSE)) %>% 
	mutate(change_salt_tol_upper = scale(change_salt_tol_upper, center = FALSE)) %>% 
	mutate(change_salt_tol_mean = scale(change_salt_tol_mean, center = FALSE))
	

all_changes <- left_join(rstar_CI_p2, rstar_CI_n2) %>% 
	left_join(., rstar_CI_i2) %>% 
	left_join(., salt_tol_CI2) %>% 
	filter(treatment != "Ancestors")
View(all_changes)

write_csv(all_changes, "data-processed/all_changes-0.56.csv")
all_changes_raw <- read_csv("data-processed/all_changes-0.56.csv") 
size_changes <- read_csv("data-processed/all-sizes-ancestors-desc.csv") %>% 
	select(ancestor_id, treatment, 38:40) %>% 
	mutate(change_phosphate_size = as.numeric(scale(change_phosphate_size, center = FALSE))) %>% 
	mutate(change_nitrate_size = scale(change_nitrate_size, center = FALSE)) %>% 
	mutate(change_light_size = scale(change_light_size, center = FALSE)) %>% 
	as.data.frame()

size_changes_meta <- read_csv("data-processed/all-sizes-ancestors-desc.csv") %>% 
	select(1:2)

size_changes <- read_csv("data-processed/all-sizes-ancestors-desc.csv") %>% 
	select(38:40) %>% 
	scale(., center = FALSE) %>% 
	as.data.frame(.) %>% 
	bind_cols(., size_changes_meta)


all_changes <- all_changes_raw %>% 
	left_join(., size_changes)

all_changes %>% 
	ggplot(aes(x = change_salt_tol_mean, y = change_pstar_mean)) + geom_point()

all_changes %>% 
	ggplot(aes(x = change_salt_tol_mean, y = change_phosphate_size)) + geom_point()


mod_change <- lm(change_pstar_mean ~  change_nstar_mean + change_istar_mean + change_salt_tol_mean + change_phosphate_size, data = all_changes)
mod_change1a <- lm(change_pstar_mean ~  change_istar_mean + change_phosphate_size + ancestor_id, data = all_changes)
mod_change1b <- lm(change_pstar_mean ~  change_nstar_mean + change_phosphate_size + ancestor_id, data = all_changes)
mod_change1c <- lm(change_pstar_mean ~  change_salt_tol_mean + change_phosphate_size + ancestor_id, data = all_changes)
# mod_change1d <- lm(change_pstar_mean ~  change_salt_tol_mean, data = all_changes)
summary(mod_change1d)

library(visreg)
## pstar vs. salt tolerance
mod_salt <- visreg(mod_change1c, "change_salt_tol_mean", plot = FALSE)
resids_salt <- data.frame(mod_salt$res) %>% 
	left_join(all_changes, by = "change_salt_tol_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_salt <- data.frame(mod_salt$fit)

plot1 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	scale_x_reverse() + 
	geom_line(aes(x = change_salt_tol_mean, y = visregFit), data = fits_salt, color = "black", size = 1) +
	geom_ribbon(aes(x = change_salt_tol_mean, ymin = visregLwr, ymax = visregUpr), data = fits_salt, alpha = 0.1) +
	ylab("Change in P* (uM P)") + xlab("Change in salt tol. (g/L)") +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes, color = treatment), data = resids_salt, size = 6) +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes), data = resids_salt, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  	  axis.title = element_text(size = 24),
		  	  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24))
ggsave("figures/change_pstar_salt_tol_RMA-OLS.png", width = 4.5, height = 3)
ggsave("figures/change_pstar_salt_tol_OLS.png", width = 4.5, height = 3)

### pstar vs. nstar
mod_change1b <- lm(change_pstar_mean ~  change_nstar_mean + change_phosphate_size + ancestor_id, data = all_changes)
summary(mod_change1b)
mod_pn <- visreg(mod_change1b, "change_nstar_mean", plot = FALSE)
resids_pn <- data.frame(mod_pn$res) %>% 
	left_join(all_changes, by = "change_nstar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_pn <- data.frame(mod_pn$fit)

plot2 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_nstar_mean, y = visregFit), data = fits_pn, color = "black", size = 1) +
	geom_ribbon(aes(x = change_nstar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_pn, alpha = 0.1) +
	ylab("Change in P* (uM P)") + xlab("Change in N* (uM N)") +
	geom_point(aes(x = change_nstar_mean, y = visregRes, color = treatment), data = resids_pn, size = 4) +
	geom_point(aes(x = change_nstar_mean, y = visregRes), data = resids_pn, shape = 1, color = "black", size = 4) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24)) +
	ylim(-3, 2) +
	xlim(-2, 2)
ggsave("figures/change_pstar_salt_tol_OLS.png", width = 4.5, height = 3)

## pstar vs. istar
mod_change1a <- lm(change_pstar_mean ~  change_istar_mean + change_phosphate_size + ancestor_id, data = all_changes)
mod_istar <- visreg(mod_change1a, "change_istar_mean", plot = FALSE)
resids_istar <- data.frame(mod_istar$res) %>% 
	left_join(all_changes, by = "change_istar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_istar <- data.frame(mod_istar$fit)

plot3 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_istar_mean, y = visregFit), data = fits_istar, color = "black", size = 1) +
	# geom_abline(intercept = fit$regression.results$Intercept[n], slope = fit$regression.results$Slope[n], col = "grey", size = 1) +
	geom_ribbon(aes(x = change_istar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_istar, alpha = 0.1) +
	ylab("Change in P* (uM P)") +
	xlab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))+
	geom_point(aes(x = change_istar_mean, y = visregRes, color = treatment), data = resids_istar, size = 4) +
	geom_point(aes(x = change_istar_mean, y = visregRes), data = resids_istar, shape = 1, color = "black", size = 4) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24)) +
	ylim(-3, 2) +
	xlim(-2, 2)
ggsave("figures/change_pstar_istar_OLS.png", width = 4, height = 3)


### nstar vs. pstar
mod_np <- lm(change_nstar_mean ~  change_pstar_mean + change_nitrate_size + ancestor_id, data = all_changes)
summary(mod_np)
mod_np <- visreg(mod_np, "change_pstar_mean", plot = FALSE)
resids_np <- data.frame(mod_np$res) %>% 
	left_join(all_changes, by = "change_pstar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_np <- data.frame(mod_np$fit)

plot4 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_pstar_mean, y = visregFit), data = fits_np, color = "black", size = 1) +
	geom_ribbon(aes(x = change_pstar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_np, alpha = 0.1) +
	ylab("Change in N* (uM N)") + xlab("Change in P* (uM P)") +
	geom_point(aes(x = change_pstar_mean, y = visregRes, color = treatment), data = resids_np, size = 6) +
	geom_point(aes(x = change_pstar_mean, y = visregRes), data = resids_np, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24))
ggsave("figures/change_nstar_pstar_OLS.png", width = 4.5, height = 3)


### nstar vs. istar
mod_ni <- lm(change_nstar_mean ~  change_istar_mean + change_nitrate_size + ancestor_id, data = all_changes)
summary(mod_ni)
mod_ni <- visreg(mod_ni, "change_istar_mean", plot = FALSE)
resids_ni <- data.frame(mod_ni$res) %>% 
	left_join(all_changes, by = "change_istar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_ni <- data.frame(mod_ni$fit)

plot5 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_istar_mean, y = visregFit), data = fits_ni, color = "black", size = 1) +
	geom_ribbon(aes(x = change_istar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_ni, alpha = 0.1) +
	ylab("Change in N* (uM N)") +
	xlab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))+
	geom_point(aes(x = change_istar_mean, y = visregRes, color = treatment), data = resids_ni, size = 4) +
	geom_point(aes(x = change_istar_mean, y = visregRes), data = resids_ni, shape = 1, color = "black", size = 4) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24)) + 
	ylim(-3, 2) +
	xlim(-2, 2)
ggsave("figures/change_nstar_istar_OLS.png", width = 4.5, height = 3)

### nstar vs. salt tolerance
mod_ns <- lm(change_nstar_mean ~  change_salt_tol_mean + change_nitrate_size + ancestor_id, data = all_changes)
summary(mod_ns)
mod_ns <- visreg(mod_ns, "change_salt_tol_mean", plot = FALSE)
resids_ns <- data.frame(mod_ns$res) %>% 
	left_join(all_changes, by = "change_salt_tol_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_ns <- data.frame(mod_ns$fit)
plot6 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	scale_x_reverse() +
	geom_line(aes(x = change_salt_tol_mean, y = visregFit), data = fits_ns, color = "black", size = 1) +
	geom_ribbon(aes(x = change_salt_tol_mean, ymin = visregLwr, ymax = visregUpr), data = fits_ns, alpha = 0.1) +
	ylab("Change in N* (uM N)") + xlab("Change in salt tol. (g/L)") +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes, color = treatment), data = resids_ns, size = 6) +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes), data = resids_ns, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24))
ggsave("figures/change_nstar_salt_tol_OLS.png", width = 4.5, height = 3)


### I* models
### istar vs. pstar
mod_ip <- lm(change_istar_mean ~  change_pstar_mean + change_light_size + ancestor_id, data = all_changes)
summary(mod_ip)
mod_ip <- visreg(mod_ip, "change_pstar_mean", plot = FALSE)
resids_ip <- data.frame(mod_ip$res) %>% 
	left_join(all_changes, by = "change_pstar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_ip <- data.frame(mod_ip$fit)

plot7 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_pstar_mean, y = visregFit), data = fits_ip, color = "black", size = 1) +
	geom_ribbon(aes(x = change_pstar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_ip, alpha = 0.1) +
	ylab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))+
	xlab("Change in P* (uM P)") +
	geom_point(aes(x = change_pstar_mean, y = visregRes, color = treatment), data = resids_ip, size = 6) +
	geom_point(aes(x = change_pstar_mean, y = visregRes), data = resids_ip, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24))
ggsave("figures/change_istar_pstar_OLS.png", width = 4.5, height = 3)


### istar vs. nstar
mod_in <- lm(change_istar_mean ~  change_nstar_mean + change_light_size + ancestor_id, data = all_changes)
summary(mod_in)
mod_in <- visreg(mod_in, "change_nstar_mean", plot = FALSE)
resids_in <- data.frame(mod_in$res) %>% 
	left_join(all_changes, by = "change_nstar_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_in <- data.frame(mod_in$fit)

plot8 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	geom_line(aes(x = change_nstar_mean, y = visregFit), data = fits_in, color = "black", size = 1) +
	geom_ribbon(aes(x = change_nstar_mean, ymin = visregLwr, ymax = visregUpr), data = fits_in, alpha = 0.1) +
	ylab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))+
	xlab("Change in N* (um N)") +
	geom_point(aes(x = change_nstar_mean, y = visregRes, color = treatment), data = resids_in, size = 6) +
	geom_point(aes(x = change_nstar_mean, y = visregRes), data = resids_in, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24))
ggsave("figures/change_istar_nstar_OLS.png", width = 4.5, height = 3)

### istar vs. salt tolerance
mod_is <- lm(change_istar_mean ~  change_salt_tol_mean + change_light_size + ancestor_id, data = all_changes)
summary(mod_is)

library(stargazer)

mod_pn <- lm(change_pstar_mean ~  change_nstar_mean + change_phosphate_size + ancestor_id, data = all_changes)
mod_ps <- lm(change_pstar_mean ~  change_salt_tol_mean + change_phosphate_size + ancestor_id, data = all_changes)
mod_pi <- lm(change_pstar_mean ~  change_istar_mean + change_phosphate_size + ancestor_id, data = all_changes)


# latex tables with OLS outputs -------------------------------------------
library(stargazer)
stargazer(mod_is, mod_in, title="Results", align=TRUE)
 
stargazer(mod_ps, mod_pi, mod_pn, title="",
		  align=TRUE, dep.var.labels= "Change in P*",
		  covariate.labels=c("Change in salt tol","Change in I*",
		  				   "Change in N*","Change in size","Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "")

stargazer(mod_ps, mod_pi, mod_pn, title="", type="html",
		  align=TRUE, dep.var.labels= "Change in P*",
		  covariate.labels=c("Change in salt tol","Change in I*",
		  				   "Change in N*","Change in size","Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "",  out="tables/p-models.htm")
mod_is <- lm(change_istar_mean ~  change_salt_tol_mean + change_light_size + ancestor_id, data = all_changes)
mod_in <- lm(change_istar_mean ~  change_nstar_mean + change_light_size + ancestor_id, data = all_changes)
mod_ip <- lm(change_istar_mean ~  change_pstar_mean + change_light_size + ancestor_id, data = all_changes)

stargazer(mod_is, mod_in, mod_ip, title="", type="html",
		  align=TRUE, dep.var.labels= "Change in I*",
		  covariate.labels=c("Change in salt tol","Change in N*",
		  				   "Change in P*","Change in size","Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "", out="tables/i-models.htm")

mod_ns <- lm(change_nstar_mean ~  change_salt_tol_mean + change_nitrate_size + ancestor_id, data = all_changes)
mod_ni <- lm(change_nstar_mean ~  change_istar_mean + change_nitrate_size + ancestor_id, data = all_changes)
mod_np <- lm(change_nstar_mean ~  change_pstar_mean + change_nitrate_size + ancestor_id, data = all_changes)
stargazer(mod_ns, mod_ni, mod_np, title="", type="html",
		  align=TRUE, dep.var.labels= "Change in N*",
		  covariate.labels=c("Change in salt tol","Change in I*",
		  				   "Change in P*","Change in size","Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "", out="tables/n-models.htm")


mod_is <- lm(change_istar_mean ~  change_salt_tol_mean + change_light_size + ancestor_id, data = all_changes)
mod_is <- visreg(mod_is, "change_salt_tol_mean", plot = FALSE)
resids_is <- data.frame(mod_is$res) %>% 
	left_join(all_changes, by = "change_salt_tol_mean") %>% 
	mutate(treatment = factor(treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
fits_is <- data.frame(mod_is$fit)
plot9 <- ggplot() +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)  +
	scale_x_reverse() +
	geom_line(aes(x = change_salt_tol_mean, y = visregFit), data = fits_is, color = "black", size = 1) +
	geom_ribbon(aes(x = change_salt_tol_mean, ymin = visregLwr, ymax = visregUpr), data = fits_is, alpha = 0.1) +
	ylab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))+
	xlab("Change in salt tol. (g/L)") +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes, color = treatment), data = resids_is, size = 6) +
	geom_point(aes(x = change_salt_tol_mean, y = visregRes), data = resids_is, shape = 1, color = "black", size = 6) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 24, family = "Lato"),
		  axis.text.x = element_text(size = 24),
		  axis.text.y = element_text(size = 24)) 
ggsave("figures/change_istar_salt_tol_OLS.png", width = 4.5, height = 3)


multi_plot_changes <- plot_grid(plot3, plot2, plot1, plot4, plot5, plot6, plot7,plot8, plot9,  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), align = "v",
								label_size = 20, label_x = 0.22)

save_plot("figures/changes-OLS-tradeoffs-big.png", multi_plot_changes,
		  ncol = 3, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.2
)

### new subset

multi_plot_changes_sub <- plot_grid(plot3, plot5, plot2,  labels = c("A", "B", "C"), align = "h",
								label_size = 20, label_x = 0.22, nrow = 1, ncol = 3)

save_plot("figures/changes-OLS-tradeoffs-big-sub.png", multi_plot_changes_sub,
		  ncol = 3, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1, base_height = 4.2, base_width = 4.2
)


summary(mod_change)
visreg(mod_change)
summary(mod_change1a)
summary(mod_change1b)
summary(mod_change1c)
visreg(mod_change1a, "change_istar_mean", gg= TRUE) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_pstar_istar.png", width = 4.5, height = 4)
visreg(mod_change1a, "change_phosphate_size", gg= TRUE) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_pstar_size.png", width = 4.5, height = 4)

visreg(mod_change1b, "change_nstar_mean", gg= TRUE) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_pstar_nstar.png", width = 4.5, height = 4)

mod_change1c <- lm(change_pstar_mean ~  change_salt_tol_mean + change_phosphate_size + ancestor_id, data = all_changes)


summary(mod_change1d)
# mc <- lmodel2(change_pstar_mean ~  change_salt_tol_mean + change_phosphate_size + ancestor_id, data = all_changes, range.y = "interval", range.x = "interval")
# mc1 <- lmodel2(change_pstar_mean ~  -change_salt_tol_mean, data = all_changes, range.y = "interval", range.x = "interval")
# summary(mc1)
# mc$regression.results
# mc$confidence.intervals
# mc1
# visreg(mc)

ggplot() +
	geom_point(data = all_changes, aes(x = change_salt_tol_mean, y = change_pstar_mean)) +
	geom_abline(slope = -0.9528606, intercept = -0.0818887) 




all_changes %>% 
ggplot(aes(x = change_phosphate_size, y = change_pstar_mean)) + geom_point() +
	geom_smooth(method = "lm")

str(all_changes)

mod_change2 <- lm(change_nstar_mean ~ change_pstar_mean + change_istar_mean + change_salt_tol_mean + change_nitrate_size, data = all_changes)
mod_change2a <- lm(change_nstar_mean ~  change_pstar_mean + change_nitrate_size + ancestor_id, data = all_changes)
mod_change2b <- lm(change_nstar_mean ~  change_istar_mean + change_nitrate_size + ancestor_id, data = all_changes)
summary(mod_change2a)
summary(mod_change2b)

visreg(mod_change2a, "change_pstar_mean", gg = TRUE) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_nstar_pstar.png", width = 4.5, height = 4)

visreg(mod_change2b, "change_istar_mean", gg = TRUE) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_nstar_istar.png", width = 4.5, height = 4)


mod_change3a <- lm(change_istar_mean ~  change_pstar_mean + change_light_size + ancestor_id, data = all_changes)
mod_change3b <- lm(change_istar_mean ~  change_nstar_mean + change_light_size + ancestor_id, data = all_changes)
summary(mod_change3a)
summary(mod_change3b)
visreg(mod_change2)

visreg(mod_change3a, "change_pstar_mean", gg = TRUE) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_istar_pstar.png", width = 4.5, height = 4)

visreg(mod_change3b, "change_nstar_mean", gg = TRUE) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("figures/change_istar_nstar.png", width = 4.5, height = 4)

# 
# set.seed(101) ## for reproducibility
# nsim <- 1000
# res <- data.frame() ## set aside space for results
# for (i in 1:nsim) {
# 	## scramble both changes
# 	perm1 <- sample(nrow(all_changes))
# 	perm2 <- sample(nrow(all_changes))
# 	hold <- data.frame(change_istar_mean = all_changes$change_istar_mean[perm1],
# 					 change_pstar_mean = all_changes$change_pstar_mean[perm2]) %>% 
# 		mutate(trade_off = case_when(change_istar_mean < 0 & change_pstar_mean >0 ~ "tradeoff",
# 									 change_istar_mean > 0 & change_pstar_mean < 0 ~ "tradeoff",
# 									 change_istar_mean < 0 & change_pstar_mean < 0 ~ "bothbetter",
# 									 change_istar_mean > 0 & change_pstar_mean > 0 ~ "bothworse")) %>% 
# 		group_by(trade_off) %>% 
# 		tally()
# 	hold$replicate <- i
# 	res <- bind_rows(res, hold)
# 	
# }

obs1 <- all_changes %>% 
	mutate(trade_off = case_when(change_istar_mean < 0 & change_pstar_mean >0 ~ "tradeoff1",
								 change_istar_mean > 0 & change_pstar_mean < 0 ~ "tradeoff2",
								 change_istar_mean < 0 & change_pstar_mean < 0 ~ "bothbetter",
								 change_istar_mean > 0 & change_pstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

obs1 <- all_changes %>% 
	mutate(trade_off = case_when(change_istar_mean < 0 & change_pstar_mean >0 ~ "tradeoff1",
								 change_istar_mean > 0 & change_pstar_mean < 0 ~ "tradeoff2",
								 change_istar_mean < 0 & change_pstar_mean < 0 ~ "bothbetter",
								 change_istar_mean > 0 & change_pstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

binom.test(14,32,(1/4),alternative="greater")

obs2 <- all_changes %>% 
	mutate(trade_off = case_when(change_nstar_mean < 0 & change_pstar_mean >0 ~ "tradeoff1",
								 change_nstar_mean > 0 & change_pstar_mean < 0 ~ "tradeoff2",
								 change_nstar_mean < 0 & change_pstar_mean < 0 ~ "bothbetter",
								 change_nstar_mean > 0 & change_pstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

obs3 <- all_changes %>% 
	mutate(trade_off = case_when(change_istar_mean < 0 & change_nstar_mean >0 ~ "tradeoff1",
								 change_istar_mean > 0 & change_nstar_mean < 0 ~ "tradeoff2",
								 change_istar_mean < 0 & change_nstar_mean < 0 ~ "bothbetter",
								 change_istar_mean > 0 & change_nstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

obs4 <- all_changes %>% 
	mutate(trade_off = case_when(change_salt_tol_mean > 0 & change_istar_mean > 0 ~ "tradeoff1",
								 change_salt_tol_mean < 0 & change_istar_mean < 0 ~ "tradeoff2",
								 change_salt_tol_mean > 0 & change_istar_mean < 0 ~ "bothbetter",
								 change_salt_tol_mean < 0 & change_istar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()



obs5 <- all_changes %>% 
	mutate(trade_off = case_when(change_salt_tol_mean > 0 & change_pstar_mean > 0 ~ "tradeoff1",
								 change_salt_tol_mean < 0 & change_pstar_mean < 0 ~ "tradeoff2",
								 change_salt_tol_mean > 0 & change_pstar_mean < 0 ~ "bothbetter",
								 change_salt_tol_mean < 0 & change_pstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

obs6 <- all_changes %>% 
	mutate(trade_off = case_when(change_salt_tol_mean > 0 & change_nstar_mean > 0 ~ "tradeoff1",
								 change_salt_tol_mean < 0 & change_nstar_mean < 0 ~ "tradeoff2",
								 change_salt_tol_mean > 0 & change_nstar_mean < 0 ~ "bothbetter",
								 change_salt_tol_mean < 0 & change_nstar_mean > 0 ~ "bothworse")) %>% 
	group_by(trade_off) %>% 
	tally()

# res %>% 
# 	filter(trade_off == "bothbetter") %>% 
# 	ggplot(aes(x = n)) + geom_histogram(bins = 9) + geom_vline(xintercept = 19, color = "blue")
# 
# both_better <- res %>% 
# 	filter(trade_off == "bothbetter")
# 
# sum(both_better$n > 19)/ length(both_better$n > 19)
# thing <- ecdf(both_better$n)
# 1- thing(19)


all_changes_anc_sum <- all_changes %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_istar_mean, change_nstar_mean, change_pstar_mean, change_salt_tol_mean)

plota <- ggplot() +
	geom_point(aes(x = change_istar_mean, y = change_pstar_mean, color = treatment), data = filter(all_changes, treatment != "Ancestors"), size = 4) +
	theme(legend.position = "top", legend.direction = "horizontal") +
	ylab("Change in P* (um P)") + 
	scale_color_manual(values = c(cols_no_anc), name = "") +
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')')) + guides(col = guide_legend(override.aes = list(size = 5, alpha = 1), 
																					   nrow = 1, title.position = "left"))

plot1a <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_pstar_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_pstar_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_pstar_mean_mean + change_pstar_mean_std.error,
						ymin = change_pstar_mean_mean - change_pstar_mean_std.error, color = treatment),
					data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_pstar_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = treatment),
				  data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in P* (um P)") + 
	xlab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1}))) 
	

plot1b <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_nstar_mean, y = change_pstar_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_nstar_mean_mean, y = change_pstar_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_nstar_mean_mean, ymax = change_pstar_mean_mean + change_pstar_mean_std.error,
					  ymin = change_pstar_mean_mean - change_pstar_mean_std.error, color = treatment),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_nstar_mean_mean - change_nstar_mean_std.error, y = change_pstar_mean_mean,
					   xmax = change_nstar_mean_mean + change_nstar_mean_std.error, color = treatment),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in P* (um P)") + 
	xlab("Change in N* (um N)") 

plot1c <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_nstar_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_nstar_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_nstar_mean_mean + change_nstar_mean_std.error,
					  ymin = change_nstar_mean_mean - change_nstar_mean_std.error, color = treatment),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_nstar_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = treatment),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in N* (um N)") + 
	xlab(expression("Change in I*" ~ (mu * mol ~ m^{-2} * s^{-1})))

plot1d <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_salt_tol_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_salt_tol_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = treatment),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = treatment),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')'))

plot1e <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_pstar_mean, y = change_salt_tol_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_pstar_mean_mean, y = change_salt_tol_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_pstar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = treatment),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_pstar_mean_mean - change_pstar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_pstar_mean_mean + change_pstar_mean_std.error, color = treatment),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab("Change in P* (uM P)")

plot1f <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_nstar_mean, y = change_salt_tol_mean, color = treatment), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_nstar_mean_mean, y = change_salt_tol_mean_mean, color = treatment), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_nstar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = treatment),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_nstar_mean_mean - change_nstar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_nstar_mean_mean + change_nstar_mean_std.error, color = treatment),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab("Change in N* (uM N)")

legend <- get_legend(plota)

multi_plot_trade_offs <- plot_grid(plot1a, plot1d, plot1b, plot1e, plot1c, plot1f, labels = c("A", "B", "C", "D", "E", "F"), align = "h", ncol = 2, nrow = 3)
save_plot("figures/trait_trade_offs2_0.56.png", multi_plot_trade_offs,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.1
)

multi_plot_trade_offs_legend <- plot_grid(legend, align = "h", ncol = 1, nrow = 1)
save_plot("figures/trait_trade_offs2_legend_0.56.png", multi_plot_trade_offs_legend,
		  ncol = 1, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.2
)


### plot the ancestors and the descendants side by side

View(all_changes)


nitrate_boot <- read.csv("data-processed/nitrate-monod-params-bootstrap.csv") %>% 
	dplyr::rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = str_replace(treatment, "none", "Ancestors")) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	summarise(nstar_lower =quantile(rstar, probs=0.025),
			  nstar_upper =quantile(rstar, probs=0.975),
			  nstar_mean = mean(rstar))

low_n <- nitrate_boot %>% 
	filter(treatment %in% c("Ancestors", "N")) 



	ggplot() +
	geom_pointrange(aes(x = treatment, y = nstar_mean, ymax = nstar_upper, ymin = nstar_lower,
						group = ancestor_id, color = ancestor_id), data = low_n) +
		geom_line(aes(x = treatment, y = nstar_mean, group = ancestor_id, color = ancestor_id), data = low_n) +
		scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "") 
	
light_boot <- read.csv("data-processed/light_monod_params_bootstrapped.csv") %>% 
	dplyr::rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = str_replace(treatment, "none", "Ancestors")) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	summarise(istar_lower =quantile(rstar, probs=0.025),
			  istar_upper =quantile(rstar, probs=0.975),
			  istar_mean = mean(rstar))


low_i <- light_boot %>% 
	filter(treatment %in% c("Ancestors", "L")) 



ggplot() +
	geom_pointrange(aes(x = treatment, y = istar_mean, ymax = istar_upper, ymin = istar_lower,
						group = ancestor_id, color = ancestor_id), data = low_i) +
	geom_line(aes(x = treatment, y = istar_mean, group = ancestor_id, color = ancestor_id), data = low_i) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "") +
	ylab(bquote('I* ('*mu~'mol'~ m^-2~s^-1*')')) + xlab("")
ggsave("figures/istars-ancestors.png", width = 6, height = 4)

### now ask about changes in consumption vectors
all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	mutate(np = 1/cons_vector_new)


all_rstars %>% 
	ggplot(aes(x = treatment, y = pc)) + geom_point()


pc_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, pc) %>% 
	dplyr::rename(pc_anc = pc)

pcs <- all_rstars %>% 
	left_join(., pc_anc) %>% 
	mutate(change_pc = pc - pc_anc)%>% 
	group_by(treatment) %>% 
	mutate(mean_change_pc = mean(change_pc)) %>% 
	mutate(se_change_pc = std.error(change_pc)) %>% 
	mutate(mean_pc = mean(pc)) %>% 
	mutate(se_pc = std.error(pc)) 


pcs %>% 
	ggplot(aes(x = treatment, y = change_pc)) + geom_point() +
	geom_hline(yintercept = 0)

pcs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_pc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_change_pc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_change_pc - se_change_pc, ymax = mean_change_pc + se_change_pc), color = "purple", width = 0.1) + 
	geom_hline(yintercept = 0) + ylab("Change in P:C") + xlab("")
ggsave("figures/change-pc.png", width = 6, height = 4)

pcs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = pc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_pc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_pc - se_pc, ymax = mean_pc + se_pc), color = "purple", width = 0.1) +
	ylab("P:C") + xlab("")
ggsave("figures/pc.png", width = 6, height = 4)


nc_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, nc) %>% 
	dplyr::rename(nc_anc = nc)

ncs <- all_rstars %>% 
	left_join(., nc_anc) %>% 
	mutate(change_nc = nc - nc_anc) %>% 
	group_by(treatment) %>% 
	mutate(mean_change_nc = mean(change_nc)) %>% 
	mutate(se_change_nc = std.error(change_nc)) %>% 
	mutate(mean_nc = mean(nc)) %>% 
	mutate(se_nc = std.error(nc)) 
	


ncs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_nc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_change_nc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_change_nc - se_change_nc, ymax = mean_change_nc + se_change_nc), color = "purple", width = 0.1) + 
	geom_hline(yintercept = 0) + ylab("Change in N:C") + xlab("")
ggsave("figures/change-nc.png", width = 6, height = 4)

ncs %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "L", "C", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = nc)) + geom_point() +
	geom_point(aes(x = treatment, y = mean_nc), size = 3, color = "purple") +
	geom_errorbar(aes(x = treatment, ymin = mean_nc - se_nc, ymax = mean_nc + se_nc), color = "purple", width = 0.1) +
	ylab("N:C") + xlab("")
ggsave("figures/nc.png", width = 6, height = 4)

consvec_anc <- all_rstars %>% 
	filter(treatment == "none") %>% 
	select(ancestor_id, cons_vector_new) %>% 
	dplyr::rename(consvec_anc = cons_vector_new)

consvecs <- all_rstars %>% 
	left_join(., consvec_anc) %>% 
	mutate(change_consvec = cons_vector_new - consvec_anc)

consvecs_sum <- consvecs %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), change_consvec) %>% 
	mutate(change_consvec_upper = mean + std.error) %>% 
	mutate(change_consvec_lower = mean - std.error) %>% 
	mutate(change_consvec_mean = mean) %>% 
	mutate(diversity = "treatment average")

consvecs %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = change_consvec, color = ancestor_id)) + geom_point(alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	ylab("Change in consumption vector \n (P:N molar ratio)") + xlab("") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "")
ggsave("figures/change-cons-vec.png", width = 6, height = 4)


all_rstars %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = cons_vector_new, color = ancestor_id)) + geom_point() +
	ylab("Consumption vector \n (P:N molar ratio)") + xlab("") +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "")
ggsave("figures/cons-vec.png", width = 6, height = 4)
	


## chi square test for trade-offs

td <- read_csv("data-processed/istar-pstar-trade-off.csv")
td2 <- table(td)
td2

df <- read.csv("https://goo.gl/j6lRXD")


td1 <- all_changes %>% 
	mutate(change_istar = ifelse(change_istar_mean < 0, "better", "worse")) %>% 
	mutate(change_pstar = ifelse(change_pstar_mean < 0, "better", "worse"))

fisher.test(table(td1$change_istar, td1$change_pstar), alternative = "greater")
chisq.test(td1$change_istar, td1$change_pstar, correct=FALSE)
?fisher.test

fisher.test(rbind(c(19,3),c(8,2)), alternative="two.sided")$p.value
fisher.test(rbind(c(1,9),c(11,3)), alternative="less")$p.value

binom.test(12,32,(1/4),alternative="greater")
