## packages
library(tidyverse)
library(ggsci)

font_import()
## load fonts
extrafont::loadfonts()
?extrafont

library(tidyverse)

df_students <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-05-07/student_teacher_ratio.csv")

df_world_tile <- readr::read_csv("https://gist.githubusercontent.com/maartenzam/787498bbc07ae06b637447dbd430ea0a/raw/9a9dafafb44d8990f85243a9c7ca349acd3a0d07/worldtilegrid.csv") %>% 
	mutate(
		## Namibias two-digit country code is handled as `NA` - let us fix that
		alpha.2 = if_else(name == "Namibia", "NA", alpha.2),
		## We are going to split "Americas" into "North America" and "Sout America"
		region = if_else(region == "Americas", sub.region, region),
		region = if_else(region %in% c("Northern America", "Central America", "Caribbean"), 
						 "North America", region),
		region = if_else(region == "Southern America", "South America", region),
		## to join both data sets, we need a id column
		country_code = alpha.3
	)

df_ratios <- df_students %>% 
	## Let's keep only the most recent data per country
	group_by(country, indicator) %>% 
	filter(year == max(year)) %>% 
	ungroup() %>% 
	# Create `NA`s for countries which do not have any data 2012-2018
	complete(indicator, nesting(country, country_code)) %>% 
	## Let's focus on primary education and keep only countries (coded by letters)
	filter(
		indicator == "Primary Education",
		str_detect(country_code, "[A-Z]")
	) %>% 
	## merge with world tile map data
	full_join(df_world_tile) %>%
	filter(
		!is.na(region),
		!is.na(indicator)
	) %>% 
	group_by(region) %>% 
	mutate(student_ratio_region = median(student_ratio, na.rm = T)) %>% 
	ungroup()

## get data
devtools::source_gist("https://gist.github.com/Z3tt/301bb0c7e3565111770121af2bd60c11")

## tile map as legend
map_regions <- df_ratios %>%
	mutate(region = fct_reorder(region, -student_ratio_region)) %>%
	ggplot(aes(x = x, y = y, fill = region, color = region)) +
	geom_tile(color = "white") +
	scale_y_reverse() +
	scale_fill_uchicago(guide = F) +
	coord_equal() +
	theme_light() +
	theme(line = element_blank(),
		  panel.background = element_rect(fill = "transparent"),
		  plot.background = element_rect(fill = "transparent", color = "transparent"),
		  panel.border = element_rect(color = "transparent"),
		  strip.background = element_rect(color = "grey20"),
		  axis.text = element_blank(),
		  plot.margin = margin(0, 0, 0, 0)) +
	labs(x = NULL, y = NULL)


## calculate worldwide average
world_avg <- df_ratios %>%
	summarize(avg = mean(student_ratio, na.rm = T)) %>%
	pull(avg)

## coordinates for arrows
arrows <- tibble(
	x1 = c(6, 3.65, 1.8, 1.8, 1.8),
	x2 = c(5.6, 4, 2.07, 2.78, 1.1),
	y1 = c(world_avg + 6, 10.5, 9, 9, 76),
	y2 = c(world_avg + 0.1, 18.32, 14.4, 11.85, 83.41195)
)

## final plot
## set seed to fix position of jittered points
set.seed(123)
?theme_light

## final plot
df_ratios %>% 
	mutate(region = fct_reorder(region, -student_ratio_region)) %>%
	ggplot(aes(region, student_ratio, color = region)) +
	geom_segment(aes(x = region, xend = region,
					 y = world_avg, yend = student_ratio_region),
				 size = 0.8) +
	geom_hline(aes(yintercept = world_avg), color = "grey70", size = 0.6) +
	geom_point(aes(region, student_ratio_region), size = 5) +
	geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
	coord_flip() +
	annotate("text", x = 6.3, y = 35, family = "Poppins", size = 2.7, color = "grey20",
			 label = glue::glue("Worldwide average:\n{round(world_avg, 1)} students per teacher")) +
	annotate("text", x = 3.5, y = 10, family = "Poppins", size = 2.7, color = "grey20",
			 label = "Continental average") +
	annotate("text", x = 1.7, y = 11, family = "Poppins", size = 2.7, color = "grey20",
			 label = "Countries per continent") +
	annotate("text", x = 1.9, y = 64, family = "Poppins", size = 2.7, color = "grey20",
			 label = "The Central African Republic has by far\nthe most students per teacher") +
	geom_curve(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
			   arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
			   color = "grey20", curvature = -0.3) +
	annotation_custom(ggplotGrob(map_regions), xmin = 2.5, xmax = 7.5, ymin = 55, ymax = 85) +
	scale_y_continuous(limits = c(0, 90), expand = c(0.005, 0.005), breaks = c(1, seq(20, 80, by = 20))) +
	scale_color_uchicago() +
	labs(x = NULL, y = "Student to teacher ratio",
		 caption = 'Data: UNESCO Institute for Statistics') +
	theme_light(base_size = 15) +
	theme_light(base_size = 15, base_family = "Lao Sangam MN") +
	theme(legend.position = "none",
		  axis.title = element_text(size = 12),
		  axis.text.x = element_text(family = "Lao Sangam MN", size = 10),
		  plot.caption = element_text(size = 9, color = "grey50"),
		  panel.grid = element_blank())

fonts()
