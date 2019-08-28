

library(tidyverse)
library(cowplot)


tiff_by_hand <- read_csv("imaging/Fluo_per_cell/results/B2-by-hand-tif.csv") %>% 
	mutate(type = "tiff by hand")
tiff_by_macro <- read_csv("imaging/Fluo_per_cell/results/B2_01_1_1_Chlorophyll_A_001_plate18_results.csv") %>% 
	mutate(type = "tiff by macro")
png_by_macro <- read_csv("imaging/Fluo_per_cell/results/B2_01_1_1_ChlorophyllA_001_plate18.png_results.csv") %>% 
	mutate(type = "png by macro")


all_cells <- bind_rows(tiff_by_hand, tiff_by_macro, png_by_macro)


all_cells %>% 
	filter(type != "png by macro") %>% 
	ggplot(aes(x = Mean, fill = type)) + geom_histogram() +
	facet_wrap( ~ type)


setdiff(tiff_by_macro$Mean, png_by_macro$Mean)
identical(tiff_by_macro$Mean, png_by_macro$Mean)
