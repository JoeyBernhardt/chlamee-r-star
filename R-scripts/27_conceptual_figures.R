### conceptual figures
library(tidyverse)
library(cowplot)
library(readxl)
library(here)
library(janitor)
theme_set(theme_cowplot())




predict_r <- function(nutrient_concentration){
	r <- umax*(nutrient_concentration/ (ks+ nutrient_concentration))
}



treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")

all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	# read_csv("data-processed/all-rstars.csv") %>% 
	mutate(cons_vector = pc/nc) ## ## pc is the P:C molar ratio, nc is the N:C molar ratio

source(here("R-scripts", "predict_monod.R"))

all_preds_p <- all_rstars %>% 
	dplyr::rename(ks = p_ks,
		   umax = p_umax) %>% 
	split(.$population) %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population") %>% 
	rename(phosphate_concentration = nitrate_concentration.x) %>% 
	filter(phosphate_concentration < 50)

all_preds_n <- all_rstars %>% 
	rename(ks = n_ks,
		   umax = n_umax) %>% 
	split(.$population) %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population") %>% 
	rename(nitrate_concentration = nitrate_concentration.x) 

View(all_rstars)

p_sub <- all_preds_p %>% 
	filter(population %in% c(22, 35))

ggplot() +
	geom_line(aes(x = phosphate_concentration, y = growth_rate, group = population), data = p_sub) 


growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))

monod_function_fast <- function(x, umax = 1.5, ks = 10, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 

monod_function_med <- function(x, umax = 1.2, ks = 5, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
}

monod_function_slow <- function(x, umax = 1, ks = 0.8, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 


rstar_concept <- data.frame(type = c("fast", "med", "slow"), ks = c(10, 5, 0.8), umax = c(1.5, 1.2, 1)) %>% 
	mutate(r_star = ks*0.5/(umax-0.5)) 
	
View(rstar_concept)

rstar_concept %>% 
	ggplot(aes(x = r_star, y = umax, color = type)) + geom_point(size = 6) + labs(y = expression ("Max growth rate"~day^-1)) +
	geom_point(size = 6, shape = 1, color = "black") +
	xlab("P* (uM P)") + scale_color_manual(values = c("deepskyblue4", "orange", "purple"), name = "") +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") 
ggsave("figures/conceptual-monod-trade-off.pdf", width = 6.2, height = 4.6)


gt <- rstar_concept %>% 
	ggplot(aes(x = r_star, y = umax, color = type)) + geom_point(size = 4) + labs(y = expression ("Max growth rate"~day^-1)) +
	geom_point(size = 4, shape = 1, color = "black") +
	xlab("P* (uM P)") + scale_color_manual(values = c("deepskyblue4", "orange", "lightblue"),
																	 name = "") +
 geom_smooth(method = "lm", color = "black", se = FALSE) +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") + 
	panel_border(colour = "black", size = 1) 
	

## "Khmer Sangam MN", Arial


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
gp <- p + stat_function(fun = monod_function_slow, color = "lightblue", size = 2) +
	stat_function(fun = monod_function_med, color = "orange", size = 2) +
	stat_function(fun = monod_function_fast, color = "deepskyblue4", size = 2) +
	xlim(0, 100) + ylim(-0.5, 1) + coord_cartesian() + geom_hline(yintercept = 0) +
	geom_point(aes(x = r_star, y = 0, color = type), data = rstar_concept, size = 5) +
	geom_point(aes(x = r_star, y = 0), data = rstar_concept, size = 5, shape = 1, color = "black") +
	labs(y = expression ("Population growth rate"~day^-1)) +
	xlab("uM P") +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	scale_color_manual(values = c("deepskyblue4", "orange", "lightblue"), name = "") + 
	panel_border(colour = "black", size = 1) 
ggsave("figures/conceptual-monod.pdf", width = 6, height = 4.6)


### gleaner opportunist

monod_function_fast <- function(x, umax = 1.5, ks = 10, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 

monod_function_slow <- function(x, umax = 1, ks = 0.8, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 


	p + stat_function(fun = monod_function_slow, color = "darkgreen", size = 2) +
		stat_function(fun = monod_function_fast, color = "deepskyblue4", size = 2) +
		coord_cartesian() + geom_hline(yintercept = 0) +
		geom_point(aes(x = r_star, y = 0, color = type), data = filter(rstar_concept, type != "med"), size = 7) +
		geom_point(aes(x = r_star, y = 0), data = filter(rstar_concept, type != "med"), size = 7, shape = 1, color = "black") +
		labs(y = expression ("Population growth rate"~day^-1)) +
		xlab("Resource concentration") +
		theme( 
			axis.text = element_text(size=22, color = "black"),
			axis.title=element_text(size=22, color = "black"),
			legend.position = "none") +
		scale_color_manual(values = c("deepskyblue4","darkgreen"), name = "") + 
		geom_hline(yintercept = 1, color = "deepskyblue4", linetype = "dashed") + xlim(0, 100) +
		geom_hline(yintercept = 0.5, color = "darkgreen", linetype = "dashed") 
ggsave("figures/gleaner-opportunist-diagram.pdf", width = 6.6, height = 4.6)
ggsave("figures/gleaner-opportunist-diagram.png", width = 6.6, height = 5)


r_star <-  (0.8*0.5)/(1-0.5)
r_star2 <-  (10*0.5)/(1.5-0.5)

p + geom_point(aes(x = 0.8, y = 0.5), color = "darkgreen", size = 3)  +
	geom_point(aes(x = 5, y = 1), size = 3, color = "deepskyblue4") +ylab("μmax") +
	xlab("R*")
ggsave("figures/gleaner-opportunist-points.pdf", width = 3, height = 2)
ggsave("figures/gleaner-opportunist-points.png", width = 3, height = 2)


## big house fast car
bmonod_function_fast <- function(x, umax = 1.5, ks = 0.8, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 

bmonod_function_med <- function(x, umax = 1.2, ks = 5, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
}

bmonod_function_slow <- function(x, umax = 1, ks = 10, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 


brstar_concept <- data.frame(type = c("fast", "med", "slow"), ks = c(10, 5, 0.8), umax = c(1, 1.2, 1.5)) %>% 
	mutate(r_star = ks*0.5/(umax-0.5)) 



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
bp <-  p + stat_function(fun = bmonod_function_slow, color = "lightblue", size = 2) +
	stat_function(fun = bmonod_function_med, color = "orange", size = 2) +
	stat_function(fun = bmonod_function_fast, color = "deepskyblue4", size = 2) +
	xlim(0, 100) + ylim(-0.5, 1.2) + coord_cartesian() + geom_hline(yintercept = 0) +
	geom_point(aes(x = r_star, y = 0, color = type), data = brstar_concept, size = 5) +
	geom_point(aes(x = r_star, y = 0), data = brstar_concept, size = 5, shape = 1, color = "black") +
	labs(y = expression ("Population growth rate"~day^-1)) +
	xlab("uM P") +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"), name = "") +
	panel_border(colour = "black", size = 1) 

bt <- brstar_concept %>% 
	ggplot(aes(x = r_star, y = umax, color = type)) + geom_point(size = 4) + labs(y = expression ("Max growth rate"~day^-1)) +
	geom_point(size = 4, shape = 1, color = "black") +
	xlab("P* (uM P)") + scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"),
																	 name = "") +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", size = 1) 

## no umax variation

umonod_function_fast <- function(x, umax = 0.7, ks = 0.8, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 

umonod_function_med <- function(x, umax = 0.7, ks = 2, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
}

umonod_function_slow <- function(x, umax = 0.7, ks = 3, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 


urstar_concept <- data.frame(type = c("fast", "med", "slow"), ks = c(3, 2, 0.8), umax = c(0.7, 0.7, 0.7)) %>% 
	mutate(r_star = ks*0.5/(umax-0.5)) 



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
up <- p + stat_function(fun = umonod_function_slow, color = "lightblue", size = 2) +
	stat_function(fun = umonod_function_med, color = "orange", size = 2) +
	stat_function(fun = umonod_function_fast, color = "deepskyblue4", size = 2) +
	xlim(0, 100) + ylim(-0.2, 0.3) + coord_cartesian() + geom_hline(yintercept = 0) +
	geom_point(aes(x = r_star, y = 0, color = type), data = urstar_concept, size = 5) +
	geom_point(aes(x = r_star, y = 0), data = urstar_concept, size = 5, shape = 1, color = "black") +
	labs(y = expression ("Population growth rate"~day^-1)) +
	xlab("uM P") +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"), name = "") +
	panel_border(colour = "black", size = 1) 

ut <- urstar_concept %>% 
	ggplot(aes(x = r_star, y = umax, color = type)) + geom_point(size = 4) + labs(y = expression ("Max growth rate"~day^-1)) +
	geom_point(size = 4, shape = 1, color = "black") +
	xlab("P* (uM P)") + scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"),
																	 name = "") +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", size = 1) 

## no rstar variation
ksa <- 0.8
umaxa <- 1
((0.5*ksa)/(umaxa - 0.5))

ks1 <- 1.2
umax1 <- ((0.5*1.2) + (0.8*0.5))/0.8
((0.5*ks1)/(umax1 - 0.5))

ks2 <- 1.5
umax2 <- ((0.5*1.5) + (0.8*0.5))/0.8
((0.5*ks2)/(umax2 - 0.5))



rmonod_function_fast <- function(x, umax = umaxa, ks = ksa, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 

rmonod_function_med <- function(x, umax = umax1, ks = ks1, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
}

rmonod_function_slow <- function(x, umax = umax2, ks = ks2, m = 0.5){
	growth_rate <- (umax * x / (ks+x)) - m
	growth_rate
} 


rrstar_concept <- data.frame(type = c("fast", "med", "slow"), ks = c(ksa, ks1, ks2), umax = c(umaxa, umax1, umax2)) %>% 
	mutate(r_star = ks*0.5/(umax-0.5)) 



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
rp <- p + stat_function(fun = rmonod_function_slow, color = "lightblue", size = 2) +
	stat_function(fun = rmonod_function_med, color = "orange", size = 2) +
	stat_function(fun = rmonod_function_fast, color = "deepskyblue4", size = 2) +
	xlim(0, 100) + ylim(-0.4,1) + coord_cartesian() + geom_hline(yintercept = 0) +
	geom_point(aes(x = r_star, y = 0, color = type), data = rrstar_concept, size = 5) +
	geom_point(aes(x = r_star, y = 0), data = rrstar_concept, size = 5, shape = 1, color = "black") +
	labs(y = expression ("Population growth rate"~day^-1)) +
	xlab("uM P") +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"), name = "") +
	panel_border(colour = "black", size = 1) 

rt <- rrstar_concept %>% 
	ggplot(aes(x = r_star, y = umax, color = type)) + geom_point(size = 4) + labs(y = expression ("Max growth rate"~day^-1)) +
	geom_point(size = 4, shape = 1, color = "black") +
	xlab("P* (uM P)") + scale_color_manual(values = c("lightblue", "orange", "deepskyblue4"),
																	 name = "") +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=18, family = "Arial", color = "black"),
		axis.title=element_text(size=18, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", size = 1)  


multi_plot <- plot_grid(gp, bp, rp, up, gt, bt, rt, ut, labels = c("A", "B", "C","D", "E", "F", "G", "H"), align = "v", nrow = 2, label_x = 0.1)
save_plot("figures/figure1-conceptual.png", multi_plot,
		  ncol = 4, # we're saving a grid plot of 2 columns
		  nrow = 2, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.1
)

### now make the equivalent R* plot

rstar_concept2 <- data.frame(type = c("fast", "slow"), ks = c(10, 0.8), umax = c(1.5, 1)) %>% 
	mutate(p_star = ks*0.5/(umax-0.5)) 

View(rstar_concept2)
rstar_concept2$n_star <- c(1, 4)
rstar_concept2$nc <- c(0.2057244, 0.6057244)
rstar_concept2$pc <- c(1, 1.3)



rstar_concept2 %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color = type)) + geom_point() +
	geom_segment(aes(color = type, x = max_n_star, y = max_p_star, xend = max_n_star + nc*2, yend = max_p_star + pc*2),
				 size = .5, linetype = 2) +
	geom_segment(aes(color = type, x = max_n_star, y = max_p_star, xend = max_n_star - nc*2, yend = max_p_star - pc*2),
				 size = 0.5, linetype = 2,
				 arrow = arrow(type = "closed", length = unit(0.09, "inches"))) + 
	geom_segment(aes(color = type, x = n_star, y = p_star, xend = n_star, yend = 8),
				 size = 1, linetype = 1) +
	geom_segment(aes(color = type, x = n_star, y = p_star, xend = 6, yend = p_star),
				 size = 1, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	scale_color_manual(values = c("deepskyblue4", "orange"), name = "") +
	# geom_point(aes(x = 1.76, y = 0.097), color = "black") + 
	ylim(0.0, 8) + xlim(0, 6) +
	coord_cartesian() +
	theme( 
		# plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=20, family = "Arial", color = "black"),
		axis.title=element_text(size=20, family = "Arial", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none")+
	panel_border(colour = "black") 
ggsave("figures/conceptual-monod-zngis.pdf", width = 5, height = 4)


