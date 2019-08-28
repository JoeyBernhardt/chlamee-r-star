
combos <- poster_combos[[1]]


find_alphas_combos <- function(combos, SN1, SP1){
	snippet <- all_rstars %>% 
		filter(population %in% c(combos$comp1[[1]], combos$comp2[[1]])) 
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	c1P <- snippet$pc[snippet$population == pop1]
	c2P <- snippet$pc[snippet$population == pop2]
	c1N <- snippet$nc[snippet$population == pop1]
	c2N <- snippet$nc[snippet$population == pop2]
	R1P <- snippet$p_star[snippet$population == pop1]
	R2N <- snippet$n_star[snippet$population == pop2]
	R2P <- snippet$p_star[snippet$population == pop2]
	R1N <- snippet$n_star[snippet$population == pop1]
	D <- 0.5
	
	## how can we define the consumption vector lines?
	
	SN <- SN1
	SP <- SP1
	
	
	r1 <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
	r2 <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])
	
	cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
	cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1], nitrogen_supply = SN1,
						phosphorus_supply = SP1, zone = zone, combination = combos$combination[1])
	alphas_calc2 <- bind_cols(alphas_calc, r_s, comps)
	return(alphas_calc2)
}

all_combos_poster <- read_csv( "data-processed/all_combos_poster.csv")
all_rstars <- read_csv("data-processed/all_rstars_new_stoich.csv") %>% 
	# read_csv("data-processed/all-rstars.csv") %>% 
	mutate(cons_vector = pc/nc) 

poster_combos <- combn(all_combos_poster$population, 2) %>% 
	t() %>% 
	as.data.frame() %>% 
	rename(comp1 = V1,
		   comp2 = V2) %>% 
	mutate(combination = rownames(.)) %>% 
	split(.$combination)


all_combos_allop <- read_csv("data-processed/all_combos_allopatry.csv") %>% 
	filter(combination %in% c(39, 29, 32)) %>% 
	split(.$combination)

ns <- c(1.75)
ns <- seq(1.50, 1.8, by = 0.01)
ps <- c(0.096)
ps <- seq(0.08, 1, by = 0.05)

ns <- seq(min(snip_pop$n_star)-0.01, max(snip_pop$n_star) + 0.01, length.out = 10)
ps <- seq(min(snip_pop$p_star) - 0.001, max(snip_pop$p_star) + 0.001, length.out = 10)

poster_combos[[1]]

snip_pop <- all_rstars %>% 
	filter(population %in% c("anc3", 26, 4))

results <- data.frame()
for(i in ns){
	for(j in ps){
		hold <- map_df(all_combos_allop, find_alphas_combos, SN1 = i, SP1 = j)
		hold$n_supply <- i
		hold$p_supply <- j
		results <- bind_rows(results, hold)
	}}

results3 <- results %>% 
	mutate(rho = sqrt((a12*a21)/(a11*a22))) %>% 
	mutate(fit_ratio = (sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D))) %>% 
	mutate(coexist = rho < fit_ratio &  fit_ratio < 1/rho) %>% 
	mutate(stabil_potential = 1 - rho) %>% 
	mutate(location = "evolved in allopatry") %>% 
	mutate(supply_ratio = n_supply/p_supply) %>% 
	mutate(combination = as.factor(combination)) %>% 
	mutate(competitors = case_when(combination == "29" ~ "low light vs. high salt",
							combination == "32" ~ "ancestor vs. high salt",
							combination == "39" ~ "ancestor vs. low light",
							TRUE ~ "NA")) %>% 
	mutate(unique_combo = rownames(.)) %>% 
	filter(unique_combo %in% c(1, 3))


allo_long %>% 
	# filter(treatment.y %in% c("P", "L", "N")) %>% 
	filter(unique_combination == 29) %>% 
	filter(competition_outcome == "stable coexistence") %>% 
	group_by(combination) %>% 
	mutate(max_n_star = max(n_star)) %>% 
	mutate(max_p_star = max(p_star)) %>% 
	ggplot(aes(x = n_star, y = p_star, color = treatment.y, label = competition_outcome)) + geom_point() +
	geom_segment(aes(color = treatment.y, x = max_n_star, y = max_p_star, xend = max_n_star + nc, yend = max_p_star + pc),
				 size = 1.5, linetype = 2) +
	geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 0.15),
				 size = 1.5, linetype = 1) +
	geom_segment(aes(color = treatment.y, x = n_star, y = p_star, xend = 1.8, yend = p_star),
				 size = 1.5, linetype = 1) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star + nc, yend = p_star + pc),
				 size = 1.5, linetype = 2) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = n_star, yend = 0.15),
				 size = 1.5, linetype = 1) +
	geom_segment(data = symp_snip, aes(color = treatment.y, x = n_star, y = p_star, xend = 1.8, yend = p_star),
				 size = 1.5, linetype = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	geom_point(aes(x = 1.75, y = 0.096), color = "black") + 
	# facet_wrap( ~ unique_combination) + 
	scale_color_manual(values = c(cols[5], "grey", cols[6])) +
	ylim(0.06, 0.15) + xlim(1.2, 1.8) +
	coord_cartesian() +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1.5)  

mylims.x.Chesson <-  1
mylims.y.Chesson  <- 4
x.1 = seq(1, 0, by=-0.001)
x.2 = seq(0, -1, by=-0.001)
y1.1 = 1/((1-x.1))
y2.1 = 1-x.1
y1.2 = 1/((1-x.2))
y2.2 = 1-x.2

cf.df.1 = data.frame(my.x = c(x.1, x.1), 
					 my.y = c(y1.1, y2.1), 
					 whichy = c(rep("y1", times = length(seq(1, 0, by=-0.001))), 
					 		   rep("y2", times = length(seq(1, 0, by=-0.001)))))
rib.dims.1 = data.frame(min.dim = y1.1, max.dim = y2.1, x.dim = x.1)

cf.df.2 = data.frame(my.x = c(x.2, x.2), 
					 my.y = c(y1.2, y2.2), 
					 whichy = c(rep("y1", times = length(seq(0, -1, by=-0.001))), 
					 		   rep("y2", times = length(seq(0, -1, by=-0.001)))))
rib.dims.2 = data.frame(min.dim = y1.2, max.dim = y2.2, x.dim = x.2)


ggplot() + 
	geom_ribbon(data = rib.dims.1, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "grey", 
				alpha = 0.3) +
	geom_line(data = cf.df.1, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") +     
	geom_ribbon(data = rib.dims.2, 
				aes(x = x.dim, 
					ymin = min.dim, 
					ymax = max.dim), 
				fill = "darkgrey", 
				alpha = 0.6) +
	geom_line(data = cf.df.2, 
			  aes(x = my.x, 
			  	y = my.y,
			  	# linetype = whichy,
			  	group = whichy), 
			  col = "black") + 
	geom_abline(intercept = 0, 
				slope = 0, 
				lty = 2, 
				col = "darkgrey") + 
	geom_vline(xintercept = 0, 
			   lty = 2, 
			   col = "darkgrey") +
	panel_border(colour = "black") +
	geom_point(aes(x = stabil_potential, y = fit_ratio, color = competitors), size = 3, data = results3, alpha = 0.75) +
	geom_point(aes(x = stabil_potential, y = fit_ratio), size = 3, data = results3, shape = 1) +
	coord_cartesian(expand = c(0, 0), 
					xlim = c(-1, 1), 
					ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
	scale_y_log10() +
	# Add text for coexistence/priority effect region:
	geom_text(aes(x = 0.0,
				  y = exp(log(4)/2)),
			  label = "Population 2 wins",
			  size = 6.0, 
			  fontface = "bold",
			  family = "Arial Rounded MT Bold") +  
	geom_text(aes(x = 0.0,
				  y = exp(-log(4)/2)),
			  label = "Population 1 wins",
			  size = 6.0, 
			  fontface = "bold",
			  family = "Arial Rounded MT Bold") +  
	geom_text(aes(x = 0.5,
				  y = exp(log(1.7)/18)),
			  label = "Coexistence",
			  size = 6.0, 
			  fontface = "bold",
			  family = "Arial Rounded MT Bold") +  
	geom_text(aes(x = -0.5,
				  y = exp(log(1.7)/18)),
			  label = "Priority effect",
			  size = 6.0, 
			  fontface = "bold", 
			  family = "Arial Rounded MT Bold") +
	scale_color_viridis_d(name = "competitors") +
	
	scale_linetype_manual(values = c("y2" = "solid", 
									 "y1" = "dotted")) +
	xlab(expression(paste("Niche differences ( 1 - ", rho, " )", sep = ""))) + 
	ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
	theme(
		axis.title.y = element_text(angle = 90), 
		axis.title = element_text(size = 14)) +
	theme( 
		plot.margin = unit(c(1,1,1,1), "lines"),
		axis.text = element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		axis.title=element_text(size=20, family = "Arial Rounded MT Bold", color = "black"),
		rect = element_rect(fill = "transparent"),
		legend.position = "none") +
	panel_border(colour = "black", linetype = "solid", size = 1.5)
ggsave("figures/chesson-plot-poster.pdf", width = 6.2, height = 5.2)
