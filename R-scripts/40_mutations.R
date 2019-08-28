
library(tidyverse)
library(cowplot)
library(purrr)


cbbPalette <- c("#6b6b6b",
				"#f9a729",
				"#97cfd0",
				"#00a2b3",
				"#f1788d",
				"#cf3e53",
				"#b9ca5d")

long_variant_table <-
	read_csv("genomics/long_variant_table.csv") %>% 
	mutate(Variant_type = ifelse(Proportion > 70, 1, 0)) %>% 
	filter(Depth > 20,
		   Depth < 250)

mut_table <- 
	long_variant_table %>% 
	unite(Loc, Chromosome, Locus, sep="_") %>%
	dplyr::select(Loc, Variant_type, Strain) %>%
	spread(Strain, Variant_type) %>%
	mutate(rsum = rowSums(.[-1], na.rm = TRUE)) %>%
	filter(rsum != 37) %>%
	dplyr::select(-rsum) %>% 
	gather(Str, Val, -Loc) %>%
	separate(Str, c("Strain", "Treatment"), "_")


########################################################################################
## Prepare pairwise comparisons of each ancestor to all of its respective descendants ##
########################################################################################

strain  <- "Anc2"

mut_count <- function(strain)
{
	exp_subset <- 
		filter(mut_table,
			   Strain == strain,
			   !is.na(Val))
	left_join(filter(exp_subset, Treatment == "A"), filter(exp_subset, Treatment != "A"), by="Loc") %>% 
		filter(Val.x != Val.y) %>% 
		group_by(Treatment.y) %>% 
		tally() %>% 
		mutate(Strain = strain) %>%
		dplyr::rename(Treatment = Treatment.y)
}

cbbPalette <- c("#6b6b6b",
				"#f9a729",
				"#97cfd0",
				"#00a2b3",
				"#f1788d",
				"#cf3e53",
				"#b9ca5d")

map_dfr(
	c("CC1690", "Anc2", "Anc3", "Anc4", "Anc5"),
	mut_count) %>%
	spread(Treatment, n)

unique(mut_table$Strain)

## 111100715 bases in Chlamy genome
## 285 generations

mutations <- map_dfr(
	c("CC1690", "Anc2", "Anc3", "Anc4", "Anc5"),
	mut_count) %>%
	mutate(mu = n) %>% 
	mutate(mu_per_gen = n / 285 / 111100715) %>% 
	mutate(Treatment = factor(Treatment,
							  levels=c("C", "L", "N",
							  		 "P", "B", "S", "BS")))
write_csv(mutations, "data-processed/mutations.csv")

plota <- mutations %>% 
	filter(Strain != "CC1690") %>% 
	ggplot(aes(x=Treatment, y=mu)) +
	geom_boxplot() +
	geom_point(aes(color=Strain), size = 3) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18, n = 5), name = "") +
	ylab("Number of mutations") + theme(legend.position = "top") + xlab("Selection environment")



plotb <- mutations %>% 
	filter(Strain != "CC1690") %>% 
	ggplot(aes(x=Strain, y=mu)) +
	geom_boxplot() +
	geom_point(aes(color=Treatment), size = 3) +
	scale_color_manual(values = cbbPalette, name = "Selection environment") +
	ylab("Number of mutations") + xlab("Ancestor") +
	theme(legend.position = "top")

multi_plot_muts <- plot_grid(plotb, plota, labels = c("A", "B"), align = "h", nrow = 1, ncol = 2, rel_widths = c(1, 0.8))
save_plot("figures/plot_muts.png", multi_plot_muts,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 1, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1, 
		  base_width = 4.3, base_height = 4
)



map_dfr(
	c("CC1690", "Anc2", "Anc3", "Anc4", "Anc5"),
	mut_count) %>% View
	mutate(mu = n / 285 / 111100715) %>% 
	ggplot(aes(x=Treatment, y=mu)) +
	geom_boxplot() +
	geom_point(aes(color=Strain)) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "")
ggsave("treatment_effect_on_mu.pdf", last_plot(), useDingbats=FALSE)

map_dfr(
	c("CC1691", "Anc2", "Anc3", "Anc4", "Anc5"),
	mut_count) %>%
	mutate(mu = n / 285 / 111100715) %>% 
	ggplot(aes(x=Strain, y=mu)) +
	geom_boxplot() +
	geom_point(aes(color=Treatment)) +
	scale_color_manual(values = cbbPalette)