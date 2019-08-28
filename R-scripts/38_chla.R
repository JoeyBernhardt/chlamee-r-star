


library(emoGG)
library(readxl)
library(janitor)
library(cowplot)

treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629")

chl <- read_excel("data-raw/stoichiometry/chla_analysis_joey_18.07.19.xlsx", skip = 1) %>% 
	clean_names() %>% 
	filter(population != "ethanol") %>% 
	filter(! number %in% c(27)) %>% 
	mutate(population = str_replace(population, "11 LR", "14")) %>% 
	mutate(population = str_replace(population, " ", "")) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = str_replace(treatment, "none", "Ancestors")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	rename(chla = chla_µg_l)
names(chl)

emoji_search("car")

chl %>% 
	ggplot(aes(x = treatment, y = chla)) + geom_emoji(emoji="1f697") +
	ylab("Chla (ug/L)") + xlab("")
ggsave("figures/chla-cars.png", width = 10, height = 8)
