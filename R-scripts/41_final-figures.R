
library(PNWColors)

col_vir <- (viridis_pal(option = "magma", begin = 0.8, end = 0.2)(7))
col_jake <- pnw_palette("Bay",7,type="continuous")
col_vir <- col_jake[c(1,6,4,7,3,2,5)]
cols_no_anc <- col_vir
str(col_jake)

q4 <- qualitative_hcl(7, palette = "Dark2")
cols_no_anc  <- c(q4)
cols_no_anc <- c("#78b98f", "#1b511d", "#7bde3f", "#9a23b1", "#c6c0fe", "#1255d3", "#35c8ef")
library(colorspace)

salt %>% 
	ggplot(aes(x = treatment, y = change_salt_tol_mean, color = treatment, fill = treatment)) + 
	geom_pointrange(alpha = 0.7, aes(shape = diversity, size = diversity2, x = unique_treatment,
									 y = change_salt_tol_mean, ymin = change_salt_tol_lower, ymax = change_salt_tol_upper),
					position=position_dodge2(width = 0.01), data = salt) +
	geom_hline(yintercept = 0) +
	ylab(expression("Change in salt tolerance" ~ (g ~ L^{-1}))) +
	xlab("") + 
	scale_color_manual(values = c("black", cols_no_anc)) +
	scale_fill_manual(values = c("black", cols_no_anc)) +
	# scale_fill_discrete_qualitative(palette = "Dark 2") +
	# scale_color_discrete_qualitative(palette = "Dark 2") +
	theme(legend.position="none",
		  axis.text = element_text(size=19),
		  axis.title=element_text(size=19)) +
	theme(
		axis.text = element_text(size=19),
		axis.title=element_text(size=19)) +
	scale_size(range = c(0.7, 1.2)) +
	scale_shape_manual(values = c(25, 19, 19))  +
	theme(axis.title.x=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	scale_y_continuous(labels = scales::number_format(accuracy = 1)) 
