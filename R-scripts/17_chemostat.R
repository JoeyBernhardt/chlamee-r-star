### plot the RFUs over time


library(tidyverse)
library(beyonce)
library(janitor)
library(lubridate)


rfus <- read_csv("data-raw/chemostats/RFU-over-time.csv") %>% 
	clean_names()


rfus2 <- rfus %>% 
	mutate(ancestor_id = case_when(id == "2" ~"anc2",
								   id == "3" ~ "anc3",
								   id == "4" ~ "anc4",
								   id == "5" ~ "anc5",
								   id == "pop" ~ "cc1690")) %>% 
	filter(!is.na(treatment)) 


rfus2 %>% 
	filter(is.na(treatment)) %>% View


beyonce_palette(type = "discrete", 18)
beyonce_palette()
c(beyonce_palette(type = "discrete", 18))

rfus2 %>% 
	ggplot(aes(x = date, y = rfu, color = ancestor_id)) + geom_point() +
	geom_line() +
	scale_color_manual(values = c("#424395", "#018AC4", "#5EC2DA", "#EBC915", "#550133"), name = "Ancestor") +
	# scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	facet_wrap(~ treatment) + ylab("RFU") + xlab("Date") + scale_y_log10()
ggsave("figures/chemostat-rfus-time.png", width = 8, height = 6)


all_temps <- read_csv("data-processed/all_temps_cells_RFU.csv") %>% 
	filter(cells_per_ml < 10000, temperature == 20)

all_temps %>% 
	filter(cells_per_ml < 10000, temperature == 20) %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, shape = factor(temperature))) +
	geom_point(size = 3) + geom_smooth(method = "lm")

mod1 <- all_temps %>% 
	filter(cells_per_ml < 10000, temperature == 20) %>% 
	lm(cells_per_ml ~ RFU, data = .) 

library(modelr)
rfus2 %>% 
	rename(RFU = rfu) %>% 
	add_predictions(mod1) %>%
	mutate(cell_count = pred*27) %>% 
	ggplot(aes(x = date, y = cell_count, color = ancestor_id)) + geom_point() +
	geom_line() +
	scale_color_manual(values = c("#424395", "#018AC4", "#5EC2DA", "#EBC915", "#550133"), name = "Ancestor") +
	# scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	facet_wrap(~ treatment) + ylab("Population abundance (cell count)") + xlab("Date") + scale_y_log10() +
	geom_hline(yintercept = 10000, color = "grey") 
	ggsave("figures/chemostat-cell-counts-time.png", width = 8, height = 6)
	
	cell_counts_chemostat <- rfus2 %>% 
		rename(RFU = rfu) %>% 
		add_predictions(mod1) %>%
		mutate(cell_count = pred*27)

	write_csv(cell_counts_chemostat, "data-processed/chemostat-cell-counts.csv")
	
	
library(plotrix)	
	cell_counts_chemostat %>% 
		filter(week > 29) %>% 
		group_by(treatment, ancestor_id) %>% 
		summarise_each(funs(mean, std.error), cell_count) %>% View
	