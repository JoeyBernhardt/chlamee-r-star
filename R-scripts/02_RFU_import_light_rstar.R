
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 

plate_key <- read_excel("data-general/ChlamEElightexpt_Plate_key_plate_number_20181217.xlsx") %>% 
	clean_names()

plate_info <- left_join(plate_layout, plate_key, by = c("plate_key")) %>% 
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population)) %>% 
	rename(plate = plate_number)


write_csv(plate_info, "data-processed/chlamee-light-rstar-plate-info.csv")


RFU_files <- c(list.files("data-raw/light-rstar-rfu", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")

# RFU_files[grepl("104", RFU_files)]

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>%
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	mutate(file_name = str_replace(file_name, " ", ""))


all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	spread(key = plate_number, value = plate_1) %>%  
	separate(Time, into = c("crap", "time"), sep = " ") %>% 
	select(-crap) %>% 
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE) %>% 
	separate(plate, into = c("plate", "other"), sep = "_") %>% 
	select(-other)


all_plates2 <- dplyr::left_join(all_plates, all_times, by = "file_name")

all_temp_RFU <- all_plates2 %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))





all_rfus_raw <- left_join(all_temp_RFU, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate)) %>% 
	filter(!is.na(plate_key))



all_rfus2 <- all_rfus_raw %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time)) %>% 
	mutate(population = ifelse(population == "cc1629", "COMBO", population))


all_rfus3 <- all_rfus2 %>% 
	group_by(expt_light_level) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 

# library(plotrix)
# all_rfus3 %>%
# 	filter(population != "COMBO") %>% 
# 	summarise_each(funs(mean, std.error), RFU) %>% View
# 

all_rfus2 %>% 
	ggplot(aes(x = population, y = RFU, color = acclimation_level)) + geom_jitter(size = 3) +
	geom_hline(yintercept = 7.889) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	
ggsave("figures/light-rstar-inoculation-rfus.pdf", width = 10, height = 6)

all_rfus4 <- all_rfus3 %>% 
	mutate(light_level = str_replace(expt_light_level, "L", "")) %>% 
	mutate(light_level = as.numeric(light_level))

write_csv(all_rfus4, "data-processed/light-rstar-rfus-time.csv")

all_rfus4 %>% 
	ggplot(aes(x = days, y = RFU, color = factor(light_level), group = well_plate)) + geom_point() +
	geom_line() + facet_wrap( ~ population, scales = "free_y") + scale_color_viridis_d(name = "Light level")
ggsave("figures/light-rstar-rfus-time.pdf", width = 15, height = 10)

all_rfus4 %>% 
	ggplot(aes(x = days, y = RFU, color = factor(light_level), group = well_plate)) + geom_point() +
	geom_line() + facet_wrap( ~ light_level, scales = "free_y") + scale_color_viridis_d(name = "Light level")
ggsave("figures/light-rstar-rfus-time-light-level.pdf", width = 15, height = 10)


all_rfus4 %>%
	filter(light_level == 10) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(light_level), group = well_plate)) + geom_point() +
	geom_line() + facet_wrap( ~ light_level + population) + scale_color_viridis_d(name = "Light level")
ggsave("figures/light-rstar-rfus-time-populations.pdf", width = 25, height = 20)


order_to_sample <- all_rfus4 %>% 
	group_by(multitron, plate, population, light_level) %>% 
	summarise(mean_RFU = mean(RFU)) %>% 
	arrange(desc(mean_RFU)) %>% 
	ungroup() %>% 
	distinct(multitron, plate) %>%
	ungroup() %>% 
	select(multitron, plate)
	
write_csv(order_to_sample, "data-processed/order_to_sample_sep20_2018.csv")
