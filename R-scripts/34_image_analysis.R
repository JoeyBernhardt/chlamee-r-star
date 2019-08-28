

### imaging
library(tidyverse)

stoich_images_plate1 <- read_csv("imaging/Stoichiometry_imaging_15+16_05_19_manual_BF/
								 Stoich_abs_ChlamEE_1_BF/results/stoich-plate1-results/B2_Stoich_abs_BF_plate_1.png_results.csv")
s1 <- read_csv("imaging/Stoichiometry_imaging_15+16_05_19_manual_BF/Stoich_abs_ChlamEE_1_BF/results/stoich-plate1-results/B2_Stoich_abs_BF_plate_1.png_results.csv")


s1 %>% 
	filter(Area > 200) %>% View
	ggplot(aes(x = Area)) + geom_density()
