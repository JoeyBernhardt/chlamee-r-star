
library(tidyverse)
library(cowplot)

RFU_files <- c(list.files("data-raw/imaging/plate21", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


all_plates <- map_df(RFU_files, read_excel, .id = "file_name") %>%
	# rename(row = X__1) %>% 
	# filter(!grepl("dilution", file_name)) %>% 
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	clean_names() %>% 
	separate(file_name, into= c("path", "well"), sep = "plate21/", remove = FALSE) %>% 
	# mutate(well = case_when(well == "B2" ~ "B02",
	# 						well == "B3" ~ "B03",
	# 						well == "B4" ~ "B04",
	# 						well == "B5" ~ "B05",
	# 						well == "B6" ~ "B06",
	# 						well == "B7" ~ "B07",
	# 						well == "B8" ~ "B08",
	# 						well == "B9" ~ "B09",
	# 						well == "B10" ~ "B10")) %>% 
	mutate(well_rename = well) %>% 
	mutate(well_rename = str_replace(well_rename, "B", "B0")) %>% 
	mutate(well_rename = str_replace(well_rename, "010", "10")) %>% 
	separate(well_rename, into = c("well_orig", "side"), sep = 3, remove = FALSE) %>% 
	mutate(side = ifelse(!side %in% c("R", "L"), "C", side))

unique(all_plates$well_rename)

all_plates3 <- all_plates %>% 
	filter(area < 120)

plate21_auto <- read_excel("data-raw/imaging/chlamee_phosphate_imaging_20190425_plate21_190429_110636.xls",
						   range = "B101:N109") %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	filter(grepl("B", well)) %>% 
	filter(well != "B11") %>% 
	rename(well_orig = well)
	

all_plates2 <- left_join(all_plates3, plate21_auto, by = "well_orig") %>% 
	left_join(., summary)


unique(plate21_auto$well)
unique(all_plates2$well)

intersect(all_plates2$well, plate21_auto$well)

all_plates2 %>% 
	ggplot(aes(x = mean_chlorophyll_a_445_685, fill = side)) + geom_histogram() +
	geom_vline(aes(xintercept = RFU), color = "black") +
	geom_vline(aes(xintercept = mean), color = "blue") +
	geom_vline(aes(xintercept = median), color = "purple") +
	facet_grid( side ~ well_orig)

library(plotrix)

summary <- all_plates3 %>% 
	group_by(well) %>% 
	summarise_each(funs(mean, std.error, median), mean_chlorophyll_a_445_685)


all_plates2 %>% 
	ggplot(aes(x = RFU, y = mean)) + geom_point() +
	geom_smooth(method = "lm") +
	geom_abline(slope = 1, intercept = 0) +xlim(0, 25000) + ylim(0, 25000) +
	ylab("Mean fluorescence calculated from ind cells") + xlab("Mean fluorescence from bulk export")


phosphate_plate_key <- read_csv("data-processed/chlamee-phosphate-rstar-plate-info.csv")


plate21 <- phosphate_plate_key %>% 
	filter(plate == 21)


plate21_3 <- left_join(plate21_auto, plate21)


plate21_3 %>% 
	group_by(population) %>% 
	summarise_each(funs(mean, std.error), RFU) %>% 
	filter(!is.na(population)) %>% 
	ggplot(aes(x = population, y = mean)) +geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1)



plate21_auto_area <- read_excel("data-raw/imaging/chlamee_phosphate_imaging_20190425_plate21_190429_110636.xls",
						   range = "B77:N85") %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	rename(area = RFU)

plate21_4 <- left_join(plate21_auto_area, plate21)

plate21_4 %>% 
	group_by(population) %>% 
	summarise_each(funs(mean, std.error), area) %>% 
	filter(!is.na(population)) %>% 
	ggplot(aes(x = population, y = mean)) +geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1)


all_traits <- left_join(plate21_4, plate21_3)


all_traits %>% 
	filter(!is.na(plate_key)) %>% 
	ggplot(aes(x = area, y = RFU, color = population)) + geom_point() +
	ylab("Mean per cell fluorescence") + xlab("Mean cell area")




# April 30 exploration ----------------------------------------------------

RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/with-protocol/auto-export", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")

protocol_auto_export <- map_df(RFU_files, read_excel, .id = "file_name", range = "B87:N95") %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 3:14) %>% View
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(RFU != "?????") %>% 
	rename(area = RFU) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll")) %>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))

RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/with-protocol/cell-treated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


protocol_cell_treated <- map_df(RFU_files, read_excel, .id = "file_name") %>%
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	clean_names() %>% 
	separate(file_name, into= c("path", "well", "crap"), sep = c(-6, -3), remove = FALSE) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll"))%>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))


RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/manually-exported", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


manual_cell_treated <- map_df(RFU_files, read_excel, .id = "file_name") %>%
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	clean_names() %>% 
	separate(file_name, into= c("path", "focus_type", "well", "crap"), sep = c(-8, -6, -3), remove = FALSE) %>% 
	mutate(well = str_replace(well, "_", "")) %>% 
	mutate(well = str_replace(well, "B", "B0")) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll")) %>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))


manual_cell_treated2 <- manual_cell_treated %>% 
	select(focus_type, well, area, image_type, magnification) %>% 
	mutate(export_type = "imager") 

protocol_cell_treated2 <- protocol_cell_treated %>% 
	select(well, area, image_type, magnification) %>% 
	mutate(export_type = "protocol") %>% 
	mutate(focus_type = "AF") 
	
protocol_auto_export2 <- protocol_auto_export %>% 
	select(well, area, image_type, magnification) %>% 
	mutate(export_type = "protocol") %>% 
	mutate(focus_type = "AF") %>% 
	mutate(area = as.numeric(area))
	

all_exported <- bind_rows(manual_cell_treated2, protocol_cell_treated2, protocol_auto_export2)
	
	
all2 <- all_exported %>% 
	group_by(export_type, focus_type, well, image_type, magnification) %>% 
	summarise(mean_size = mean(area))
	
	
all2 %>% 
	# filter(!well %in% c("B05", "B07")) %>% 
	ggplot(aes(x = well, y = mean_size, color = focus_type, shape = image_type)) + geom_point() +
	facet_grid(export_type ~ magnification, scales = "free")
ggsave("figures/protocol-comparisons-cell-size.png", width = 8, height = 6)
	

# now look at chorophyll values -------------------------------------------


RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/with-protocol/auto-export", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")

chla_protocol_auto_export <- map_df(RFU_files, read_excel, .id = "file_name", range = "B99:N107") %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(RFU != "?????") %>% 
	rename(chla = RFU) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll")) %>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))

RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/with-protocol/cell-treated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


chla_protocol_cell_treated <- map_df(RFU_files, read_excel, .id = "file_name") %>%
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	clean_names() %>% 
	separate(file_name, into= c("path", "well", "crap"), sep = c(-6, -3), remove = FALSE) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll"))%>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))


RFU_files <- c(list.files("data-raw/imaging/carolina-protocol-tests/excel-outputs/manually-exported", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "") %>% 
	gsub(pattern = ".xls$", replacement = "")


chla_manual_cell_treated <- map_df(RFU_files, read_excel, .id = "file_name") %>%
	mutate(file_name = str_replace(file_name, " ", "")) %>% 
	clean_names() %>% 
	separate(file_name, into= c("path", "focus_type", "well", "crap"), sep = c(-8, -6, -3), remove = FALSE) %>% 
	mutate(well = str_replace(well, "_", "")) %>% 
	mutate(well = str_replace(well, "B", "B0")) %>% 
	mutate(image_type = ifelse(grepl("BF", file_name), "brightfield", "chlorophyll")) %>% 
	mutate(magnification = ifelse(grepl("x10", file_name), "10x", "4x"))


chla_manual_cell_treated2 <- chla_manual_cell_treated %>% 
	rename(fluor = mean_chlorophyll_a_445_685) %>% 
	select(focus_type, well, fluor, image_type, magnification) %>% 
	mutate(export_type = "manual") 

chla_protocol_cell_treated2 <- chla_protocol_cell_treated %>% 
	rename(fluor = mean_chlorophyll_a_445_685) %>% 
	select(well, fluor, image_type, magnification) %>% 
	mutate(export_type = "auto") %>% 
	mutate(focus_type = "AF") 

chla_protocol_auto_export2 <- chla_protocol_auto_export %>% 
	rename(fluor = chla) %>% 
	select(well, fluor, image_type, magnification) %>% 
	mutate(export_type = "auto") %>% 
	mutate(focus_type = "AF") %>% 
	mutate(fluor = as.numeric(fluor))


all_exported_chla <- bind_rows(chla_manual_cell_treated2, chla_protocol_cell_treated2, chla_protocol_auto_export2)


all3 <- all_exported_chla %>% 
	group_by(export_type, focus_type, well, image_type, magnification) %>% 
	summarise(mean_fluor = mean(fluor))


all3 %>% 
	ggplot(aes(x = well, y = mean_fluor, color = focus_type, shape = image_type)) + geom_point() +
	facet_grid(export_type ~ magnification, scales = "free")
ggsave("figures/protocol-comparisons-cell-fluorescence.png", width = 8, height = 6)



