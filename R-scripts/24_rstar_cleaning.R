

### data processing for Jurg

all_rstars <- read_csv("data-processed/all-rstars.csv")
combos <- read_csv("data-processed/chlamee-unique-combos2.csv") %>% 
	gather(key = competitor, value = identity, comp1, comp2)

all_combos_all <- read_csv("data-processed/all_combos_allopatry.csv")


all_rstars_clean <- all_rstars %>% 
	select(1:9, i_star, light_umax, nc, pc) %>% 
	rename(i_umax = light_umax)

write_csv(all_rstars_clean, "data-processed/rstars_clean.csv")
