

## load packages
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(plotrix)


i_ins <- read_csv("data-raw/stoichiometry/I_ins_i_outs_2019-05-15.csv") %>% 
	clean_names() %>% 
	mutate(position = case_when(in_out == "out" & position == 1 ~ 2,
								in_out == "out" & position == 2 ~ 5,
								in_out == "out" & position == 3 ~ 8,
								in_out == "out" & position == 4 ~ 10,
								TRUE ~ position))

plate_layout <- read_excel("data-raw/stoichiometry/absorbances-plate-layout.xlsx") %>% 
	select(1:3)


### calculate kbg
## kbg = (ln(Iin) - ln(Iout0))/z

## ok so here the Iout0 is the COMBO, and the Iin is the no bottle
i_ins %>% 
	filter(type %in% c("COMBO", "no bottle")) %>% 
	group_by(position, type) %>% 
	summarise_each(funs(mean, std.error), irradiance) %>% 
	ggplot(aes(x = position, y = mean, color = type)) + geom_point()


kbg <- i_ins %>% 
	filter(type %in% c("COMBO", "no bottle")) %>% 
	filter(position %in% c(2, 5, 8)) %>% 
	group_by(type) %>% 
	summarise_each(funs(mean), irradiance) %>% 
	spread(key = type, value = irradiance) %>% 
	rename(Iin = `no bottle`,
		   Iout0 = COMBO) %>% 
	mutate(kbg = (log(Iin) - log(Iout0))/1)
kbg$kbg[[1]]


absorbances1 <- read_excel("data-raw/stoichiometry/chlmaEE_Absorbance_600_1.xlsx", range = "B25:N33")
absorbances2 <- read_excel("data-raw/stoichiometry/chlmaEE_Absorbance_600_2.xlsx", range = "B25:N33")
absorbances3 <- read_excel("data-raw/stoichiometry/chlmaEE_Absorbance_600_3.xlsx", range = "B25:N33")


abs1 <- absorbances1 %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	rename(absorbance = RFU) %>% 
	mutate(plate = 1)

abs2 <- absorbances2 %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	rename(absorbance = RFU) %>% 
	mutate(plate = 2)

abs3 <- absorbances3 %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	rename(absorbance = RFU) %>% 
	mutate(plate = 3)


all_abs <- bind_rows(abs1, abs2, abs3) %>% 
	left_join(., plate_layout)


i_ins %>%
	filter(is.na(population)) %>% 
	ggplot(aes(x = irradiance, fill = type)) + geom_density()


in2 <- i_ins %>% 
	filter(type == "no bottle") %>% 
	group_by(position) %>% 
	summarise_each(funs(mean, std.error), irradiance) %>% 
	filter(position %in% c(2, 5, 8, 10)) %>% 
	rename(i_in = mean)

in3 <- i_ins %>% 
	filter(!is.na(population)) %>% 
	left_join(., in2, by = "position") 

in4 <- in3 %>% 
	mutate(light_consumption = i_in - irradiance) %>% 
	select(population, position, light_consumption, everything()) %>% 
	mutate(population = str_replace(population, "AC", "anc"))

in4b <- in3 %>% 
	mutate(light_consumption = ((log(i_in) - log(irradiance))/1) - kbg$kbg[[1]]) %>% 
	select(population, position, light_consumption, everything()) %>% 
	mutate(population = str_replace(population, "AC", "anc"))


in4 %>% 
	ggplot(aes(x = position, y = light_consumption)) + geom_bar(stat = "identity") +
	facet_wrap( ~ population)


in5 <- in4 %>% 
	filter(position %in% c(2, 5, 8)) %>% 
	group_by(population) %>% 
	summarise(mean_consumption = mean(light_consumption))

in5b <- in4b %>% 
	filter(position %in% c(2, 5, 8)) %>% 
	group_by(population) %>% 
	summarise(mean_consumption = mean(light_consumption))

write_csv(in5b, "data-processed/light_consumption.csv")

View(in5)

all_abs2 <- all_abs %>% 
	left_join(., in5)


all_abs2 %>% 
	filter(absorbance < 0.110) %>% 
	ggplot(aes(x = mean_consumption, y = absorbance)) + geom_point()



# absorbance is poopy, no bueno, moving on to RFUs -------------------------------------------------------

rfu1 <- read_excel("data-raw/stoichiometry/ChlamEE-Stoich-ChlaFluorescence-1_190515_111323_0001_.xlsx", range = "B56:N64")
rfu2 <- read_excel("data-raw/stoichiometry/ChlamEE-Stoich-ChlaFluorescence-2_190515_112045_0001_.xlsx", range = "B56:N64")
rfu3 <- read_excel("data-raw/stoichiometry/ChlamEE-Stoich-ChlaFluorescence-3_190515_112746_0001_.xlsx", range = "B56:N64")


rfu1b <- rfu1 %>%
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = 1)

rfu2b <- rfu2 %>%
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = 2)

rfu3b <- rfu3 %>%
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = 3)


all_rfu <- bind_rows(rfu1b, rfu2b, rfu3b) %>% 
	left_join(., plate_layout) %>% 
	filter(!is.na(population))

all_abs3 <- all_rfu %>% 
	left_join(., in5) %>% 
	left_join(., treatments) %>% 
	left_join(., all_abs2)

View(all_abs3)

all_abs3 %>% 
	ggplot(aes(x = mean_consumption, y = RFU, color = treatment)) + geom_point(size = 3) +
	geom_point(size = 3, shape = 1, color = "black") +
	geom_smooth(method = "lm", color ="black") +
	scale_color_brewer(type = "qual") + 
	facet_wrap( ~ treatment) + xlab("Light consumption (umols/m2/s)")

write_csv(all_abs3, "data-processed/all_light_consumption_absorbance.csv")

library(here)
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")
