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

?beyonce
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
