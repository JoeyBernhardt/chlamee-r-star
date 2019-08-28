

### get the cell size info from the gen 5 outputs

library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)
library(plotrix)
library(readxl)
library(janitor)

fraction_imaged <- 0.459987*2/32

bad_photos <- read_excel("data-processed/bad_pictures.xlsx")
bad_photos_p <- bad_photos %>% 
	filter(experiment == "phosphate") %>% 
	mutate(well_plate = paste(well, plate, sep = "_"))

bad_photos_n <- bad_photos %>% 
	filter(experiment == "nitrate") %>% 
	mutate(well_plate = paste(well, plate, sep = "_"))

bad_photos_l <- bad_photos %>% 
	filter(experiment == "light") %>% 
	mutate(well_plate = paste(well, plate, sep = "_"))


treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 

plate_key <- read_excel("data-general/ChlamEE-phosphate-plate-key.xlsx") %>% 
	clean_names() %>% 
	mutate(phosphate_concentration = case_when(n_level == "P1" ~ "0.5",
											   n_level == "P2" ~ "1",
											   n_level == "P3" ~ "2",
											   n_level == "P4" ~ "4",
											   n_level == "P5" ~ "6",
											   n_level == "P6" ~ "8",
											   n_level == "P7" ~ "10",
											   n_level == "P8" ~ "20",
											   n_level == "P9" ~ "35",
											   n_level == "P10" ~ "50"
	)) %>% 
	mutate(phosphate_concentration = as.numeric(phosphate_concentration))

plate_info <- left_join(plate_layout, plate_key, by = c("plate_key")) %>% 
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population)) %>% 
	rename(plate = plate_number) %>% 
	mutate(population_new_color = NA) %>% 
	mutate(population_new_color = case_when(population == 22 & plate == 23 ~ "white",
											population == "anc2" & plate == 23 ~ "yellow",
											TRUE ~ "no")) %>%
	# mutate(population_new_color = ifelse(population == "anc2" & plate == 23, "yellow", NA)) %>% 
	mutate(population = case_when(population_new_color == "white" ~ "anc2",
								  population_new_color == "yellow" ~ "22",
								  TRUE ~ population)) 




RFU_files <- c(list.files("imaging/phosphate-gen5-outputs", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")



all_plates <- map_df(RFU_files, read_excel, range = "B66:N74", .id = "file_name") %>%
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


all_rfus3_images <- all_rfus2 %>% 
	group_by(n_level) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) %>% 
	separate(col = n_level, sep = 1, into = c("p", "phosphate_level")) %>% 
	mutate(phosphate_level = as.numeric(phosphate_level))

sizes <- all_rfus3_images %>% 
	filter(population != "COMBO") %>% 
	mutate(size = as.numeric(RFU)) %>% 
	left_join(., treatments) 


	
	size_sum <- sizes %>% 
		filter(size != 0) %>% 
		group_by(ancestor_id, treatment, phosphate_concentration, population) %>% 
		summarise_each(funs(mean, std.error), size) 

	size_sum %>% 
		ggplot(aes(x = phosphate_concentration, y = mean)) + geom_point() +
		geom_errorbar(aes(x = phosphate_concentration, ymin = mean + std.error, ymax = mean + std.error), color = "purple",  width = 0.2) +
		facet_wrap(~ population, scales = "free") + ylab("Mean cell length") + xlab("Phosphate (uM P)")
	ggsave("figures/cell-sizes-phosphate.png", width = 12, height = 9)
	
	
	### ok let's get the cell density from these files as well
	
	
	
	
	all_plates_counts <- map_df(RFU_files, read_excel, range = "B54:N62", .id = "file_name") %>%
		rename(row = X__1) %>% 
		filter(!grepl("dilution", file_name)) %>% 
		mutate(file_name = str_replace(file_name, " ", ""))
	
	all_plates2_counts <- dplyr::left_join(all_plates_counts, all_times, by = "file_name")
	
	all_temp_RFU_counts <- all_plates2_counts %>% 
		gather(key = column, value = RFU, 3:14) %>% 
		unite(row, column, col = "well", remove = FALSE, sep = "") %>%
		mutate(column = formatC(column, width = 2, flag = 0)) %>% 
		mutate(column = str_replace(column, " ", "0")) %>% 
		unite(col = well, row, column, sep = "") %>% 
		filter(!is.na(RFU)) %>% 
		mutate(plate = as.numeric(plate))
	
	
	
	
	
	all_rfus_raw_counts <- left_join(all_temp_RFU_counts, plate_info, by = c("well", "plate")) %>% 
		mutate(plate = as.numeric(plate)) %>% 
		filter(!is.na(plate_key))
	
	
	
	all_rfus2_counts <- all_rfus_raw_counts %>%
		unite(col = date_time, Date, time, sep = " ") %>%
		mutate(date_time = ymd_hms(date_time)) %>% 
		mutate(population = ifelse(population == "cc1629", "COMBO", population))
	
	
	all_rfus3_images_counts <- all_rfus2_counts %>% 
		group_by(n_level) %>% 
		mutate(start_time = min(date_time)) %>% 
		mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
		unite(col = well_plate, well, plate, remove =  FALSE) %>% 
		separate(col = n_level, sep = 1, into = c("p", "phosphate_level")) %>% 
		mutate(phosphate_level = as.numeric(phosphate_level)) %>% 
		rename(cell_count = RFU) %>% 
		mutate(cell_count = as.numeric(cell_count)) %>% 
		filter(cell_count != 0) %>% 
		filter(population != "COMBO")
	
	
	all_rfus3_images_counts %>% 
		ggplot(aes(x = phosphate_concentration, y = cell_count)) + geom_point() +
		# geom_errorbar(aes(x = phosphate_concentration, ymin = mean + std.error, ymax = mean + std.error), color = "purple",  width = 0.2) +
		facet_wrap(~ population, scales = "free") + ylab("Cell count") + xlab("Phosphate (uM P)")
ggsave("figures/cell-counts-phosphate.png", width = 14, height = 8)


phosphate_counts <- all_rfus3_images_counts %>%
	select(well_plate, phosphate_concentration, population, cell_count) %>% 
	left_join(., sizes) %>% 
	mutate(cells_per_ml = cell_count*(1/fraction_imaged)*(1/0.125)) 
	

write_csv(phosphate_counts, "data-processed/phosphate-cell-counts-sizes-lengths-long.csv")


### the image area is 789 * 583um which is equivalent to 0.459987mm^2, and total well surface area is 32mm^2 (or 34 mm^2)

fraction_imaged <- 0.459987*2/32
1/fraction_imaged
1000/0.125

counts %>% 
	filter(cell_count < 5000) %>%
	mutate(cells_per_ml = cell_count*(1/fraction_imaged)*(1/0.125)) %>% 
	ggplot(aes(x = cells_per_ml, y = size, color = log(phosphate_concentration), group = phosphate_concentration)) + geom_point() +
	ylab("Cell length") + xlab("Cells/mL") + scale_color_viridis_c() + geom_smooth(method = "lm") +
	facet_wrap( ~ phosphate_concentration, scales = "free")
ggsave("figures/phosphate-length-vs-count-slopes-facet.png", width = 12, height = 6)

counts %>% 
	filter(cell_count < 5000) %>%
	mutate(cells_per_ml = cell_count*(1/fraction_imaged)*(1/0.125)) %>% 
	ggplot(aes(x = phosphate_concentration, y = cells_per_ml, group = phosphate_concentration)) + geom_point() +
	ylab("Phosphate (uM P)") + ylab("Cells/mL") +
	facet_wrap( ~ population, scales = "free") + 
	scale_y_continuous(breaks = c(10000,100000,150000, 300000, 500000, 1000000, 2000000)) 
ggsave("figures/phosphate-cells-ml-facet.png", width = 17, height = 12)



# nitrate cell sizes and counts -------------------------------------------

plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 

plate_key <- read_excel("data-general/ChlamEE-nitrate-plate-key.xlsx") %>% 
	clean_names()

plate_info <- left_join(plate_layout, plate_key, by = c("plate_key")) %>% 
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population)) %>% 
	rename(plate = plate_number)



RFU_files <- c(list.files("imaging/nitrate-gen5-outputs", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B66:N74", .id = "file_name") %>%
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


all_rfus3_images_nitrate <- all_rfus2 %>% 
	group_by(n_level) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) %>% 
	separate(col = n_level, sep = 1, into = c("n", "nitrate_level")) %>% 
	mutate(nitrate_level = as.numeric(nitrate_level)) %>% 
	mutate(nitrate_concentration = NA) %>% 
	mutate(nitrate_concentration = case_when(nitrate_level == "1" ~ "5",
											 nitrate_level == "2" ~ "10",
											 nitrate_level == "3" ~ "20",
											 nitrate_level == "4" ~ "40",
											 nitrate_level == "5" ~ "60",
											 nitrate_level == "6" ~ "80",
											 nitrate_level == "7" ~ "100",
											 nitrate_level == "8" ~ "400",
											 nitrate_level == "9" ~ "600",
											 nitrate_level == "10" ~ "1000",
											 TRUE ~ "no")) %>%
	mutate(population_corrected = case_when(plate == "19" & population == "11" ~ "anc5",
											plate == "19" & population == "anc5" ~ "11",
											plate == "24" & population == "5" ~ "30",
											plate == "24" & population == "30" ~ "5",
											TRUE ~ population)) %>% 
	select(-population) %>% 
	rename(population = population_corrected) %>% 
	rename(cell_size = RFU) %>% 
	filter(!is.na(cell_size), cell_size != 0) %>% 
	mutate(cell_size = as.numeric(cell_size)) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	filter(population != "COMBO")

all_rfus3_images_nitrate %>% 
	ggplot(aes(x = nitrate_concentration, y = cell_size)) + geom_point() +
	# geom_errorbar(aes(x = phosphate_concentration, ymin = mean + std.error, ymax = mean + std.error), color = "purple",  width = 0.2) +
	facet_wrap(~ population, scales = "free") + ylab("Cell size") + xlab("Nitrate (uM N)")
ggsave("figures/cell_sizes_nitrate.png", width = 14, height = 10)

### ok let's get the cell density from these files as well


all_plates_counts <- map_df(RFU_files, read_excel, range = "B54:N62", .id = "file_name") %>%
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	mutate(file_name = str_replace(file_name, " ", ""))

all_plates2_counts <- dplyr::left_join(all_plates_counts, all_times, by = "file_name")

all_temp_RFU_counts <- all_plates2_counts %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))





all_rfus_raw_counts <- left_join(all_temp_RFU_counts, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate)) %>% 
	filter(!is.na(plate_key))



all_rfus2_counts <- all_rfus_raw_counts %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time)) %>% 
	mutate(population = ifelse(population == "cc1629", "COMBO", population))


all_rfus3_images_counts_nitrate <- all_rfus2_counts %>% 
	group_by(n_level) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) %>% 
	separate(col = n_level, sep = 1, into = c("n", "nitrate_level")) %>% 
	mutate(nitrate_level = as.numeric(nitrate_level)) %>% 
	rename(cell_count = RFU) %>% 
	mutate(cell_count = as.numeric(cell_count)) %>% 
	filter(cell_count != 0) %>% 
	filter(population != "COMBO") %>% 
	mutate(nitrate_concentration = NA) %>% 
	mutate(nitrate_concentration = case_when(nitrate_level == "1" ~ "5",
											 nitrate_level == "2" ~ "10",
											 nitrate_level == "3" ~ "20",
											 nitrate_level == "4" ~ "40",
											 nitrate_level == "5" ~ "60",
											 nitrate_level == "6" ~ "80",
											 nitrate_level == "7" ~ "100",
											 nitrate_level == "8" ~ "400",
											 nitrate_level == "9" ~ "600",
											 nitrate_level == "10" ~ "1000",
											 TRUE ~ "no")) %>%
	mutate(population_corrected = case_when(plate == "19" & population == "11" ~ "anc5",
											plate == "19" & population == "anc5" ~ "11",
											plate == "24" & population == "5" ~ "30",
											plate == "24" & population == "30" ~ "5",
											TRUE ~ population)) %>% 
	select(-population) %>% 
	rename(population = population_corrected) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	filter(population != "COMBO")



nitrate_cells_ml <- all_rfus3_images_counts_nitrate %>% 
	# filter(cell_count < 5000) %>%
	mutate(cells_per_ml = cell_count*(1/fraction_imaged)*(1/0.125)) 

nitrate_cells_ml %>% 
	ggplot(aes(x = nitrate_concentration, y = cells_per_ml)) + geom_point() +
	# geom_errorbar(aes(x = phosphate_concentration, ymin = mean + std.error, ymax = mean + std.error), color = "purple",  width = 0.2) +
	facet_wrap(~ population, scales = "free") + ylab("Cells/mL") + xlab("Nitrate (uM N)")
ggsave("figures/cell-counts-nitrate.png", width = 14, height = 8)


nitrate_counts <- nitrate_cells_ml %>%
	select(well_plate, nitrate_concentration, population, cells_per_ml) %>% 
	left_join(., all_rfus3_images_nitrate)

nitrate_counts %>% 
	ggplot(aes(x = cells_per_ml, y = cell_size, color = log(nitrate_concentration), group = nitrate_concentration)) + geom_point() +
	scale_color_viridis_c() + geom_smooth(method = "lm") + facet_wrap( ~ nitrate_concentration)
write_csv(nitrate_counts, "data-processed/nitrate-cell-counts-sizes-lengths.csv")


# now the light experiment ------------------------------------------------


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



RFU_files <- c(list.files("imaging/light-gen5-outputs", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


all_plates <- map_df(RFU_files, read_excel, range = "B66:N74", .id = "file_name") %>%
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

all_rfus4 <- all_rfus3 %>% 
	mutate(light_level = str_replace(expt_light_level, "L", "")) %>% 
	mutate(light_level = as.numeric(light_level)) %>% 
	rename(cell_size = RFU) %>% 
	mutate(cell_size = as.numeric(cell_size)) %>% 
	filter(cell_size != 0) %>% 
	filter(population != "COMBO") %>% 
	mutate(percentage = ifelse(percentage == "0.5-0.7", "0.6", percentage)) %>% 
	mutate(percentage = as.numeric(percentage)) %>% 
	mutate(percentage = percentage/100) %>% 
	mutate(light = percentage*250) 

all_rfus4 %>% 
	ggplot(aes(x = light, y = cell_size)) + geom_point() +
	facet_wrap( ~ population, scales = "free") + ylab("Cell size") + xlab(expression("Light" ~ (mu * mol ~ m^{-2} * s^{-1})))
ggsave("figures/light-cell-sizes.png", width = 14, height = 12)

### now bring in the light data counts

all_plates_counts <- map_df(RFU_files, read_excel, range = "B54:N62", .id = "file_name") %>%
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	mutate(file_name = str_replace(file_name, " ", ""))

all_plates2_counts <- dplyr::left_join(all_plates_counts, all_times, by = "file_name")

all_temp_RFU_counts <- all_plates2_counts %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))





all_rfus_raw_counts <- left_join(all_temp_RFU_counts, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate)) %>% 
	filter(!is.na(plate_key))

fraction_imaged <- 0.459987*2/32
1/fraction_imaged
1000/0.125

all_rfus2_counts <- all_rfus_raw_counts %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time)) %>% 
	mutate(population = ifelse(population == "cc1629", "COMBO", population)) %>% 
	rename(cell_count = RFU) %>% 
	filter(population != "COMBO") %>% 
	mutate(percentage = ifelse(percentage == "0.5-0.7", "0.6", percentage)) %>% 
	mutate(percentage = as.numeric(percentage)) %>% 
	mutate(percentage = percentage/100) %>% 
	mutate(light = percentage*250) %>% 
	select(light, population, well, plate, cell_count) %>% 
	filter(cell_count != 0) %>% 
	mutate(cell_count = as.numeric(cell_count)) %>% 
	mutate(cells_per_ml = cell_count*(1/fraction_imaged)*(1/0.125)) 

light_cell_counts_sizes <- left_join(all_rfus4, all_rfus2_counts)
write_csv(light_cell_counts_sizes, "data-processed/light-cell-counts-sizes-lengths.csv")


light_cell_counts_sizes %>% 
	ggplot(aes(x = cells_per_ml, y = cell_size, color = light, group = light)) + geom_point() +
	scale_color_viridis_c() + ylab("Cell size") + xlab("Cells/mL") + geom_smooth(method = "lm")

light_cell_counts_sizes %>% 
	ggplot(aes(x = light, y = cells_per_ml)) + geom_point() + 
	facet_wrap( ~ population, scales = "free")
ggsave("figures/light-cell-count.png", width = 17, height = 12)


# now all cell sizes together --------------------------------------------

### well B3 on plate 14 might have counted some small debris as cells.
## remove well E6 from plate 14, light; C4 plate 13

nitrate_sizes <- read_csv("data-processed/nitrate-cell-counts-sizes.csv") %>% 
	mutate(experiment = "nitrate") %>% 
	filter(!well_plate %in% c(bad_photos_n$well_plate)) %>% 
	select(population, nitrate_concentration, cell_size, cells_per_ml, experiment) %>% 
	filter(nitrate_concentration %in% c(5)) %>% 
	rename(nitrate_area = cell_size)

nitrate_lengths <- read_csv("data-processed/nitrate-cell-counts-sizes-lengths.csv") %>% 
	mutate(experiment = "nitrate") %>% 
	filter(!well_plate %in% c(bad_photos_n$well_plate)) %>% 
	select(population, nitrate_concentration, cell_size, cells_per_ml, experiment) %>% 
	filter(nitrate_concentration %in% c(5)) %>% 
	rename(nitrate_length = cell_size)

phosphate_sizes <- read_csv("data-processed/phosphate-cell-counts-sizes.csv") %>% 
	mutate(experiment = "phosphate") %>% 
	filter(!well_plate %in% c(bad_photos_p$well_plate)) %>% 
	select(population, phosphate_concentration, size, cells_per_ml, experiment) %>% 
	rename(phosphate_area = size) %>% 
	filter(phosphate_concentration %in% c(0.5))

phosphate_lengths <- read_csv("data-processed/phosphate-cell-counts-sizes-lengths.csv") %>% 
	mutate(experiment = "phosphate") %>% 
	filter(!well_plate %in% c(bad_photos_p$well_plate)) %>% 
	select(population, phosphate_concentration, size, cells_per_ml, experiment, treatment) %>% 
	rename(phosphate_length = size) %>% 
	mutate(p_keep = ifelse(phosphate_concentration == 0.5, "yes", "no")) %>% 
	mutate(p_keep = ifelse(phosphate_concentration == 2 & population == "22", "yes", p_keep)) %>% 
	# 
	# mutate(phophate_keep = case_when(phosphate_concentration == "0.5" ~ "yes",
	# 								 phosphate_concentration == "1", treatment == "control" ~ "yes",
	# 								 TRUE ~ "no")) %>% 
	filter(p_keep == "yes") %>% 
	select(-treatment)

light_sizes <- read_csv("data-processed/light-cell-counts-sizes.csv") %>% 
	mutate(experiment = "light") %>% 
	filter(!well_plate %in% c(bad_photos_l$well_plate)) %>% 
	filter(light == 27.50) %>% 
	mutate(weirdos = ifelse(plate == "15" & cells_per_ml > 500000, "weirdo", "fine")) %>% 
	filter(weirdos == "fine") %>% 
	select(population, light, cell_size, cells_per_ml, experiment, well, plate) %>% 
	rename(light_area = cell_size)

light_lengths <- read_csv("data-processed/light-cell-counts-sizes-lengths.csv") %>% 
	mutate(experiment = "light") %>% 
	filter(!well_plate %in% c(bad_photos_l$well_plate)) %>% 
	filter(light == 27.50) %>% 
	mutate(weirdos = ifelse(plate == "15" & cells_per_ml > 500000, "weirdo", "fine")) %>% 
	filter(weirdos == "fine") %>% 
	select(population, light, cell_size, cells_per_ml, experiment, well, plate) %>% 
	rename(light_length = cell_size)

lightb <- left_join(light_sizes, light_lengths) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), light_length, light_area)
phoshpateb <- left_join(phosphate_lengths, phosphate_sizes) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), phosphate_length, phosphate_area)
nitrateb <- left_join(nitrate_sizes, nitrate_lengths) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), nitrate_length, nitrate_area)

all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	mutate(np = 1/cons_vector_new) 
all_sizes3 <- left_join(lightb, phoshpateb) %>% 
	left_join(., nitrateb) %>% 
	left_join(., all_rstars) %>% 
	mutate(nitrate_biovolume = ((nitrate_length/2)^3)*(4/3)*pi) %>% 
	mutate(phosphate_biovolume = ((phosphate_length/2)^3)*(4/3)*pi) %>% 
	mutate(light_biovolume = ((light_length/2)^3)*(4/3)*pi) %>% 
	select(-term, -std.error, -p.value, -statistic, -name) %>% 
	select(ancestor_id, treatment, population, contains("star"), contains("umax"), contains("ks"), contains("biovolume"), everything()) %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS")))
	

write_csv(all_sizes3, "data-processed/all_cell_sizes_low_resources.csv")


all_sizes3 <- read_csv("data-processed/all_cell_sizes_low_resources.csv")
all_sizes3 %>% 
	ggplot(aes(x = nitrate_biovolume, y = n_star)) + geom_point() + 
	geom_smooth(method = "lm", color = "black") + ylab("N* (um N)") + 
	xlab(expression("Cell biovolume" ~ (mu * m^{3}))) 
ggsave("figures/n-star-biovolume.png", width = 8, height =6)

all_sizes3 %>% 
	ggplot(aes(x = phosphate_biovolume, y = p_star)) + geom_point() + 
	geom_smooth(method = "lm", color = "black") + ylab("P* (um P)") + xlab(expression("Cell biovolume" ~ (mu * m^{3}))) 
ggsave("figures/p-star-biovolume.png", width = 8, height =6)

all_sizes3 %>% 
	ggplot(aes(x = light_biovolume, y = i_star)) + geom_point() + 
	geom_smooth(method = "lm", color = "black") + 
	xlab(expression("Cell biovolume" ~ (mu * m^{3}))) +
	ylab(expression("I*" ~ (mu * mol ~ m^{-2} * s^{-1})))
ggsave("figures/i-star-biovolume.png", width = 8, height = 6)
	

library(vegan)
library(rgl)
library(geometry)
library(scatterplot3d)
library(ptinpoly)
library(ggbiplot)


all_sizes3 %>% 
	ggplot(aes(x = phosphate_biovolume, y = nitrate_biovolume)) + geom_point() +
	geom_smooth(method = "lm")


scaled <- all_sizes3 %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	mutate(avg_biovolume = (nitrate_biovolume + phosphate_biovolume)/2) %>%
	mutate(avg_biovolume = scale(avg_biovolume)) %>%
	mutate(nitrate_biovolume = scale(nitrate_biovolume)) %>%
	mutate(n_star = scale(1/n_star)) %>%
	mutate(p_star = scale(1/p_star)) %>%
	mutate(i_star = scale(1/i_star)) %>% 
	mutate(p_umax = scale(p_umax)) %>% 
	mutate(n_umax = scale(n_umax)) %>% 
	mutate(light_umax = scale(light_umax)) %>% 
	mutate(nc = scale(nc)) %>%
	mutate(pc = scale(pc)) %>%
	mutate(phosphate_biovolume = scale(phosphate_biovolume)) %>%
	mutate(light_biovolume = scale(light_biovolume)) %>%
	select(i_star, n_star, p_star, phosphate_biovolume) %>% 
	dplyr::rename(CI = i_star,
		   CN = n_star,
		   CP = p_star,
		   size = phosphate_biovolume)
	
	
	


unscaled <- all_sizes3 %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	select(i_star, n_star, p_star, nitrate_biovolume, treatment, light_biovolume, phosphate_biovolume, population, ancestor_id, treatment) 

pca_size <- prcomp(scaled, scale. = TRUE, center = TRUE)
pca_size2 <- rda(scaled, scale. = TRUE, center = TRUE)

summary(pca_size2)


pcas <- as.data.frame(scores(pca_size2, choices = 1:2)$sites)
pcas1 <- as.data.frame(scores(pca_size2, choices = 1:2)$species) %>% 
	mutate(trait = rownames(.)) %>% 
	mutate(PC1b = PC1 + 0.001) %>% 
	mutate(PC2b = PC1 + 0.001)
pcas$treatment <- unscaled$treatment

pca_plot1 <- pcas %>% 
	ggplot(aes(x = PC1, y = PC2, color = treatment)) + geom_point(size = 3) +
	geom_point(size = 3, color = "black", shape = 1) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, text =  trait), data = pcas1, color = "black",
				 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
	scale_colour_manual(values=cbbPalette, name = "Selection environment") +
	# geom_text(aes(x = PC1, y = PC2, label = trait), data = pcas1, color = "black", size = 4.5, hjust = 0, nudge_x = -0.2, nudge_y = -0.1) +
	ylab("PC2 (35.76% var)") + xlab("PC1 (38.51% var)") 
ggsave("figures/pca-plot2.png", width = 8, height = 6)



library(plyr)
library(dplyr)

summary(pca_size)

pca_plot <- ggbiplot::ggbiplot(pca_size, group = unscaled$treatment,
				   ellipse = FALSE, size = 10, circle = FALSE, scale = 1, varname.size = 5, labels.size=10, color = "black") + 
	scale_color_manual(values = cols_anc, name = "Selection environment") +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) + 
	# xlim(-1.7, 1.2) +
	geom_point(aes(colour=unscaled$treatment), size = 2.5) +
	geom_point(shape = 1, color = "black", size = 2.5) +
	# theme(legend.position = "none") + 
	ylab("PC2 (38.51% var explained)") + xlab("PC1 (35.76% var explained)") 
ggsave("figures/nstar-pstar-istar-biovolume-stoich-pca.png", width = 8, height = 6)
ggsave("figures/nstar-pstar-istar-pbiovolume-pca.png", width = 8, height = 6)
ggsave("figures/nstar-pstar-biovolume-pca.png", width = 8, height = 6)

pca_plot <- ggbiplot::ggbiplot(pca_size, group = unscaled$treatment,
				   ellipse = FALSE, size = 6, circle = TRUE, scale = 1) + 
	scale_color_manual(values = cols_anc, name = "Selection treatment") +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) + xlim(-3, 3) + ylim(-3, 3) + ylab("PC2") + xlab("PC1") + 
	theme(legend.position = "none")


## ok since Mridul doesn't think that interpreting the trade-off from the PCA is a good idea, let's try with mixed effects model


library(lme4)
library(nlme)

unscaled2 <- unscaled %>% 
	mutate(CN = 1/n_star) %>% 
	mutate(CP = 1/p_star) %>% 
	mutate(CI = 1/i_star) %>% 
	dplyr::rename(biovolume = phosphate_biovolume)

scaled2 <- scaled %>% 
	mutate(CN = 1/n_star) %>% 
	mutate(CP = 1/p_star) %>% 
	mutate(CI = 1/i_star) 
light_monod <- read_csv("data-processed/light-monod-params-direct.csv") %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	dplyr::rename(i_ks = ks) %>% 
	dplyr::rename(i_umax = umax)


sizes_sub <- all_sizes3 %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	filter(!is.na(phosphate_biovolume)) %>%
	mutate(CN = 1/n_star) %>% 
	mutate(CP = 1/p_star) %>% 
	mutate(CI = 1/i_star) %>% 
	left_join(., light_monod) %>% 
	select(4:36) 

sizes_sub_meta <- all_sizes3 %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	filter(!is.na(phosphate_biovolume)) %>%
	mutate(CN = 1/n_star) %>% 
	mutate(CP = 1/p_star) %>% 
	mutate(CI = 1/i_star) %>% 
	select(1:3)

size_scaled <- scale(sizes_sub, center = TRUE, scale = TRUE) %>% 
	as.data.frame() %>% 
	bind_cols(., sizes_sub_meta)


library(visreg)
size_scaled_mat <- size_scaled %>% 
	select(CN, CP, CI, nitrate_biovolume, phosphate_biovolume, light_biovolume)

cov_mat <- cov(size_scaled_mat)
pca_res <- prcomp(cov_mat)
summary(pca_res)
ggbiplot::ggbiplot(pca_res)

# mod1 <- lm(CN ~  ancestor_id + n_umax, data = size_scaled)
mod1 <- lm(CN ~  CP + CI + nitrate_biovolume + ancestor_id + n_umax, data = size_scaled)
# mod1 <- lm(CN ~  CP + CI + nitrate_biovolume + ancestor_id, data = size_scaled)
# mod1 <- lm(n_ks ~  p_ks + i_ks + nitrate_biovolume + ancestor_id + n_umax, data = size_scaled)
# mod1 <- lm(n_star ~  p_star + i_star + nitrate_biovolume + ancestor_id + n_umax, data = size_scaled)
summary(mod1)


mod1 <- lm(CN ~  CP + CI + nitrate_biovolume + n_umax + ancestor_id, data = size_scaled)
mod2 <- lm(CP ~  CN + CI + phosphate_biovolume + p_umax + ancestor_id, data = size_scaled)
mod3 <- lm(CI ~  CP + CN + light_biovolume + i_umax + ancestor_id, data = size_scaled)

stargazer(mod1, title="", type = "html",
		  align=TRUE, dep.var.labels= "Competitive ability for nitrogen, CN (1/N*)",
		  covariate.labels=c("CP","CI",
		  				   "Biovolume","umax", "Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "", out="tables/CN-models.htm")

summary(mod2)
stargazer(mod2, title="", type = "html",
		  align=TRUE, dep.var.labels= "Competitive ability for phosphorus, CP (1/P*)",
		  covariate.labels=c("CN","CI",
		  				   "Biovolume","umax", "Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "", out="tables/CP-models.htm")

summary(mod3)
stargazer(mod3, title="", type = "html",
		  align=TRUE, dep.var.labels= "Competitive ability for light, CI (1/I*)",
		  covariate.labels=c("CP","CN",
		  				   "Biovolume","umax", "Anc 3","Anc 4", "Anc 5", "cc1690", "Constant"),
		  omit.stat=c("LL","ser","f"), ci=TRUE, ci.level=0.95, single.row=TRUE, digits = 2, dep.var.caption = "", out="tables/CI-models.htm")




# Partial regression plots ------------------------------------------------

mod1 <- lm(CN ~  CP + CI + nitrate_biovolume + n_umax + ancestor_id, data = size_scaled)
mod2 <- lm(CP ~  CN + CI + phosphate_biovolume + p_umax + ancestor_id, data = size_scaled)
mod3 <- lm(CI ~  CP + CN + light_biovolume + i_umax + ancestor_id, data = size_scaled)




plotn1 <- visreg(mod1, "CP", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CN") + xlab("CP")
plotn2 <- visreg(mod1, "CI", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CN") + xlab("CI")
plotn3 <- visreg(mod1, "n_umax", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CN") + xlab("umax")
plotn4 <- visreg(mod1, "nitrate_biovolume", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CN") + xlab("Biovolume when N-limited")
plotn5 <- visreg(mod1, "ancestor_id", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CN") + xlab("Ancestor")
nplots <- plot_grid(plotn1, plotn2, plotn3, plotn4, plotn5, align = "h", nrow = 3, ncol = 2,
					labels = c("A", "B", "C", "D", "E"))
save_plot("figures/CN-plots.png", nplots,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 0.5, base_width = 3.5, base_height = 2.7
)
mod2 <- lm(CP ~  CN + CI + phosphate_biovolume + p_umax + ancestor_id, data = size_scaled)

plotp1 <- visreg(mod2, "CN", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CP") + xlab("CN")
plotp2 <- visreg(mod2, "CI", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CP") + xlab("CI")
plotp3 <- visreg(mod2, "p_umax", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CP") + xlab("umax")
plotp4 <- visreg(mod2, "phosphate_biovolume", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CP") + xlab("Biovolume when P-limited")
plotp5 <- visreg(mod2, "ancestor_id", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CP") + xlab("Ancestor")
pplots <- plot_grid(plotp1, plotp2, plotp3, plotp4, plotp5, align = "h", nrow = 3, ncol = 2,
					labels = c("A", "B", "C", "D", "E"))
save_plot("figures/CP-plots.png", pplots,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 0.5, base_width = 3.5, base_height = 2.7
)


mod3 <- lm(CI ~  CP + CN + light_biovolume + i_umax + ancestor_id, data = size_scaled)
ploti1 <- visreg(mod3, "CP", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CI") + xlab("CP")
ploti2 <- visreg(mod3, "CN", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CI") + xlab("CN")
ploti3 <- visreg(mod3, "i_umax", gg = TRUE) +
	# geom_point(size = 2) +
	ylab("CI") + xlab("umax")
ploti4 <- visreg(mod3, "light_biovolume", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CI") + xlab("Biovolume when light-limited")
ploti5 <- visreg(mod3, "ancestor_id", gg = TRUE) +
	# geom_point(size = 2) + 
	ylab("CI") + xlab("Ancestor")
iplots <- plot_grid(ploti1, ploti2, ploti3, ploti4, ploti5, align = "h", nrow = 3, ncol = 2,
					labels = c("A", "B", "C", "D", "E"))
save_plot("figures/CI-plots.png", iplots,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 3, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 0.5, base_width = 3.5, base_height = 2.7
)



mod1 <- lm(CN ~  CP + CI + nitrate_biovolume + ancestor_id, data = size_scaled)
summary(mod1)
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
	geom_point(aes(x = change_nstar_mean, y = visregRes, color = treatment), data = resids_in, size = 3) +
	geom_point(aes(x = change_nstar_mean, y = visregRes), data = resids_in, shape = 1, color = "black", size = 3) +
	scale_color_manual(values = c(cols_no_anc), name = "") + 
	theme(legend.position = "none",
		  axis.title = element_text(size = 20),
		  axis.text.x = element_text(size = 20),
		  axis.text.y = element_text(size = 20))



visreg(mod1, "n_umax", gg = TRUE) +
	geom_point(size = 2) +ylab("CN") + xlab("umax")
ggsave("figures/CN-n_umax.png", width = 4, height = 3)

mod2 <- lm(CP ~  CN + CI + biovolume + ancestor_id, data = unscaled2)
mod2 <- lm(CP ~  CN + CI + phosphate_biovolume + ancestor_id, data = size_scaled)
mod2 <- lm(p_ks ~  n_ks + i_ks + phosphate_biovolume + ancestor_id + p_umax, data = size_scaled)
# mod2 <- lm(CP ~  CN + CI + phosphate_biovolume  + p_umax, data = size_scaled)
# mod2 <- lm(p_star ~  n_star + i_star + phosphate_biovolume + ancestor_id + p_umax, data = size_scaled)
summary(mod2)
visreg(mod2)

visreg(mod2, "phosphate_biovolume", gg = TRUE) 
ggsave("figures/CN-biovolume.png", width = 6, height = 4)


mod3 <- lm(CI ~  CP + CN + biovolume + ancestor_id, data = unscaled2)
mod3 <- lm(CI ~  CN + CP + light_biovolume + ancestor_id, data = size_scaled)
mod3 <- lm(i_ks ~  n_ks + p_ks + light_biovolume + ancestor_id + light_umax, data = size_scaled)
summary(mod3)
visreg(mod3, "light_umax", gg=TRUE,
	   points=list(size=2, pch=1))
ggsave("figures/visreg_plot.png")

class(plot1)
ggsave(plot1, "figures/visreg.png", width = 8, height = 6)

cor(log(unscaled2$CN), log(unscaled2$CP))
cor(unscaled2$CI, unscaled2$CN)

library(rayshader)
library(ggplot2)

gg = ggplot(unscaled2, aes(CN, CP, color = biovolume)) +
	geom_point() +
	scale_color_viridis_c(option = "A")
plot_gg(gg,multicore=TRUE,width=5,height=5,scale=250)


plot3d(unscaled2$biovolume, unscaled2$CN,  unscaled2$CP, xlab = "Biovolume when P-limited", ylab = "1/N*", zlab = "1/P*")

mod1 <- lme(log(CN) ~  log(CP) + log(CI) + phosphate_biovolume,  random = ~1|ancestor_id, data = unscaled2)

library(nlme)
mod1 <- lme(CN ~  CP + CI + biovolume,  random = ~1|ancestor_id, data = unscaled2)
summary(mod1)
anova(mod1)
ranef(mod1)

mod2 <- lme(CP ~ CN + CI + biovolume,  random = ~1|ancestor_id, data = unscaled2)
summary(mod2)
anova(mod2)
ranef(mod2)

mod3 <- lme(CI ~ CP + CN + biovolume,  random = ~ 1|ancestor_id, data = unscaled2)
summary(mod3)
anova(mod3)
ranef(mod3)


cols_anc  <- c("black", "#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")

rda_res <- rda(scaled, unscaled$treatment, scale = FALSE)
my.rda <- rda(scaled)
biplot(my.rda)


biplot(my.rda,
	   display = c("sites", 
	   			"species"),
	   type = c("text",
	   		 "points"))

#Add "hulls"
ordihull(my.rda,
		 group = unscaled$treatment)


summary(rda_res)
coef(rda_res)
r2_adj <- RsquareAdj(rda_res)$adj.r.squared
r2 <- RsquareAdj(rda_res)$r.squared

plot(cca_res)
ordihull(cca_res, groups = unscaled$treatment, col = cols_anc, draw = "polygon")


loadings <- as.data.frame(scores(cca_res, choices = c(1,2))$sites) %>% 
	mutate(population = rownames(.)) %>% 
	mutate(population = str_replace(population, "sit", ""))

loadings <- as.data.frame(scores(pca_size, choices = c(1,2))$sites) %>% 
	mutate(population = rownames(.)) %>% 
	mutate(population = str_replace(population, "sit", ""))

loadings$treatment <- unscaled$treatment
loadings$population <- unscaled$population

loadings %>%
	ggplot(aes(x = RDA1, y = RDA2, color = treatment)) + geom_point(size = 3) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = cols_anc, name = "Selection treatment")
	
loadings16 <- scores(pca_res,choices=c(1,2))
summary(pca_size)

ggbiplot::ggbiplot(pca_size, group = unscaled$treatment, ellipse = FALSE, size = 4) + 
	scale_color_manual(values = cols_anc, name = "Selection treatment") +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0)
ggsave("figures/nstar-pstar-biovolume-stoich-pca.png", width = 8, height = 6)
ggsave("figures/nstar-pstar-biovolume-pca.png", width = 8, height = 6)

mod1 <- lm(scaled$n_star ~ scaled$p_star*scaled$nitrate_biovolume)
mod1 <- lm(scaled$n_star ~ scaled$p_star*scaled$phosphate_biovolume)
mod1 <- lm(scaled$n_star ~ scaled$nitrate_biovolume)
summary(mod1)
tidy(mod1, conf.int = TRUE)

mod3 <- lm(scaled$i_star ~ scaled$n_star + scaled$light_biovolume + scaled$p_star)
summary(mod3)

mod2 <- lm(scaled$p_star ~ scaled$n_star*scaled$phosphate_biovolume)
# mod2 <- lm(scaled$p_star ~ scaled$phosphate_biovolume)
summary(mod2)
tidy(mod1, conf.int = TRUE)

unscaled <- all_sizes3 %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	select(n_star, p_star, nitrate_biovolume) 


plot3d(all_sizes3$phosphate_biovolume, all_sizes3$n_star, all_sizes3$p_star, xlab = "Biovolume when P-limited", ylab = "N*", zlab = "P*")
points3d(all_sizes3$nitrate_biovolume, all_sizes3$n_star, all_sizes3$p_star)



all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	# mutate(mean_biovolume_p = mean(phosphate_biovolume)) %>% 
	# mutate(se_biovolume_p = std.error(phosphate_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = nitrate_length, shape = diversity, color = treatment)) + geom_jitter(width = 0.2) +
	# geom_pointrange(alpha = 0.7, aes(color = treatment, x = treatment, y = mean_biovolume_p, ymin = mean_biovolume_p - se_biovolume_p, ymax = mean_biovolume_p + se_biovolume_p),
	# 				position=position_jitterdodge(jitter.width = 1.2, jitter.height = 0,
	# 											  dodge.width = 0.5, seed = 1), size = 0.4) +
	# geom_point(aes(x = treatment, y = mean_biovolume_p), size = 4, shape = 16) +
	# geom_errorbar(aes(x = treatment, ymin = mean_biovolume_p - se_biovolume_p, ymax = mean_biovolume_p + se_biovolume_p), width = 0.1) + 
	scale_color_manual(values = cols_anc, name = "Selection treatment") + ylab("Cell biovolume when P-limited") + xlab("") +
	theme(legend.position = "none")
ggsave("figures/biovolume-p-limited.png", width = 8, height = 6)

all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_n = mean(nitrate_biovolume)) %>% 
	mutate(se_biovolume_n = std.error(nitrate_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = nitrate_biovolume, shape = diversity, color = treatment)) + geom_jitter(width = 0.2) +
	geom_point(aes(x = treatment, y = mean_biovolume_n), size = 4, shape = 16) +
	geom_errorbar(aes(x = treatment, ymin = mean_biovolume_n - se_biovolume_n, ymax = mean_biovolume_n + se_biovolume_n), width = 0.1) + 
	scale_color_manual(values = cols_anc) + ylab("Cell biovolume when N-limited") + xlab("")
ggsave("figures/biovolume-n-limited.png", width = 8, height = 6)


### ok let's get changes in cell size associated with resource limited envts

nsize_anc <- all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_n = mean(nitrate_biovolume)) %>% 
	mutate(se_biovolume_n = std.error(nitrate_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	filter(treatment == "Ancestors") %>% 
	select(ancestor_id, nitrate_biovolume) %>% 
	rename(ancestor_nitrate_biovolume = nitrate_biovolume)

psize_anc <- all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_p = mean(phosphate_biovolume)) %>% 
	mutate(se_biovolume_p = std.error(phosphate_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	filter(treatment == "Ancestors") %>% 
	select(ancestor_id, phosphate_biovolume) %>% 
	rename(ancestor_phosphate_biovolume = phosphate_biovolume)
lsize_anc <- all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_i = mean(light_biovolume)) %>% 
	mutate(se_biovolume_i = std.error(light_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	filter(treatment == "Ancestors") %>% 
	select(ancestor_id, light_biovolume) %>% 
	rename(ancestor_light_biovolume = light_biovolume)


sizes4 <- all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_n = mean(nitrate_biovolume)) %>% 
	mutate(se_biovolume_n = std.error(nitrate_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	left_join(nsize_anc) %>% 
	left_join(psize_anc) %>% 
	left_join(lsize_anc)


sizes5 <- sizes4 %>% 
	mutate(change_nitrate_size = nitrate_biovolume - ancestor_nitrate_biovolume) %>% 
	mutate(change_phosphate_size = phosphate_biovolume - ancestor_phosphate_biovolume) %>% 
	mutate(change_light_size = light_biovolume - ancestor_light_biovolume) 
	

write_csv(sizes5, "data-processed/all-sizes-ancestors-desc.csv")


sizes4 %>% 
	mutate(change_nitrate_size = nitrate_biovolume - ancestor_nitrate_biovolume) %>% 
	group_by(treatment) %>% 
	mutate(mean_change_nitrate_size = mean(change_nitrate_size)) %>% 
	mutate(se_change_nitrate_size = std.error(change_nitrate_size)) %>% 
	ggplot(aes(x = treatment, y = change_nitrate_size, color = treatment)) + geom_point() +
	geom_hline(yintercept = 0) + scale_color_manual(values = cols_anc) + ylab("Change in cell biovolume when N-limited") + xlab("") +
	geom_pointrange(aes(x = treatment, y = mean_change_nitrate_size, ymin = mean_change_nitrate_size - se_change_nitrate_size,
						ymax = mean_change_nitrate_size + se_change_nitrate_size, color = treatment))
ggsave("figures/change-cell-size-nitrate-limited.png", width = 8, height = 6)



sizes4 %>% 
	mutate(change_phosphate_size = phosphate_biovolume - ancestor_phosphate_biovolume) %>% 
	group_by(treatment) %>% 
	mutate(mean_change_phosphate_size = mean(change_phosphate_size)) %>% 
	mutate(se_change_phosphate_size = std.error(change_phosphate_size)) %>% 
	ggplot(aes(x = treatment, y = change_phosphate_size, color = treatment)) + geom_point() +
	geom_hline(yintercept = 0) + scale_color_manual(values = cols_anc) + ylab("Change in cell biovolume when P-limited") + xlab("") +
	geom_pointrange(aes(x = treatment, y = mean_change_phosphate_size, ymin = mean_change_phosphate_size - se_change_phosphate_size,
						ymax = mean_change_phosphate_size + se_change_phosphate_size, color = treatment))
ggsave("figures/change-cell-size-phosphate-limited.png", width = 8, height = 6)




library(tidyverse)
library(plotrix)
library(cowplot)
all_sizes3 %>% 
	mutate(diversity = ifelse(ancestor_id == "cc1690", "Genotypically diverse", "Isoclonal")) %>% 
	# filter(!is.na(phosphate_biovolume)) %>% 
	dplyr::group_by(treatment) %>% 
	mutate(mean_biovolume_l = mean(light_biovolume)) %>% 
	mutate(se_biovolume_l = std.error(light_biovolume)) %>% 
	ungroup() %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>% 
	ggplot(aes(x = treatment, y = light_biovolume, shape = diversity, color = treatment, group = ancestor_id)) + geom_jitter(width = 0.2) +
	geom_point(aes(x = treatment, y = mean_biovolume_l), size = 4, shape = 16) +
	geom_errorbar(aes(x = treatment, ymin = mean_biovolume_l - se_biovolume_l, ymax = mean_biovolume_l + se_biovolume_l), width = 0.1) + 
	scale_color_manual(values = cols_anc) + ylab("Cell biovolume when light-limited") + xlab("")
ggsave("figures/biovolume-light-limited.png", width = 8, height = 6)
