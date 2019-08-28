

### reading in the fluorescence images
library(tidyverse)
library(purrr)
library(readxl)

plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 
plate_info <- read_csv("data-processed/chlamee-phosphate-rstar-plate-info.csv")

RFU_files <- c(list.files("imaging/Fluo_per_cell/results", full.names = TRUE))
RFU_files <- RFU_files[grepl(".csv", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".csv$", replacement = "") %>% 
	gsub(pattern = ".csv$", replacement = "")


all_fl <- map_df(RFU_files, read_csv, .id = "file_name") 


all_fl$file_name[[1]]

all_fl2 <- all_fl %>% 
	separate(col = file_name, into = c("file_path", "well"),
			 sep = c("imaging/Fluo_per_cell/results/")) %>%
	separate(col = well, into = c("well_id", "other"), sep = "_", remove = FALSE) %>% 
	separate(well_id, into = c("letter", "number"), sep = 1, remove = FALSE) %>% 
	mutate(number = as.numeric(number)) %>% 
	mutate(number = formatC(number, width = 2, flag = 0)) %>% 
	unite(col = well, letter, number, sep = "") %>% 
	separate(Label, into = c("stuff","plate"), sep = "plate") %>%
	mutate(plate = str_replace(plate, ".jpg", "")) %>% 
	mutate(plate = as.numeric(plate))


all_fl3 <- left_join(all_fl2, plate_info, by = c("well", "plate"))


all_fl3 %>% 
	ggplot(aes(x = ))