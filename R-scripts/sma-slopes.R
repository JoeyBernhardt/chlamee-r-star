function (data, x, y, method = "OLS", mod.line = F, group = NA, 
		  same.panel = TRUE, xpos = 0.3, ypos = 0.8, pt.col = "blue", 
		  pt.shape = 1, ...) 
{
	if (is.na(group)) {
		data <- data[complete.cases(data[, c(x, y)]), ]
		xpos <-  0.3
		ypos <- 0.8
		pt.col <- "blue" 
		pt.shape <- 1
		x <-  "change_salt_tol_mean"
		y <- "change_pstar_mean"
		method <-  "SMA"
		data <- all_changes
		fit <- lmodel2(as.formula(paste(y, x, sep = "~")), data)
		n <-  which(method %in% c("OLS", "MA", "SMA"))
		n <- 3
		plot <- ggplot(data, aes_string(x = x, y = y)) + 
			geom_point(colour = "blue") + 
			geom_abline(intercept = fit$regression.results$Intercept[n], slope = fit$regression.results$Slope[n], col = "black") + 
			geom_abline(intercept = fit$confidence.intervals[n, 2], slope = fit$confidence.intervals[n, 5], col = "grey") + 
			geom_abline(intercept = fit$confidence.intervals[n, 3], slope = fit$confidence.intervals[n, 4], col = "grey") + 
			geom_abline(intercept = fit$regression.results$Intercept[1], slope = fit$regression.results$Slope[1], col = "pink") +
			geom_abline(intercept = fit$confidence.intervals[1, 2], slope = fit$confidence.intervals[1, 5], col = "purple") +
			geom_abline(intercept = fit$confidence.intervals[1, 3], slope = fit$confidence.intervals[1, 4], col = "purple") +
			geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
			theme_classic()
		ggsave("figures/RMA-OLS-pstar-nstar.png", width = 5, height = 4)
	}
	if (!is.na(group)) {
		data <- data[complete.cases(data[, c(x, y, group)]), 
					 ]
		regDF <- plyr::ddply(data, plyr::as.quoted(group), lmStat, 
							 form = as.formula(paste(y, x, sep = "~")))
		if (same.panel) {
			eqnDF <- data.frame(text_eqn = plyr::daply(data, 
													   plyr::as.quoted(group), eqn, form = as.formula(paste(y, 
													   													 x, sep = "~")), method = method), text_x = xpos * 
									max(data[[x]], na.rm = T), text_y = spread_num(ypos, 
																				   0.05, length(unique(data[[group]]))) * max(data[[y]], 
																				   										   na.rm = T))
			eqnDF$text_group <- rownames(eqnDF)
			plot <- ggplot(data, aes_string(x = x, y = y, colour = group)) + 
				geom_point() + geom_abline(data = subset(regDF, 
														 Method == method), aes_string(intercept = "Intercept", 
														 							  slope = "Slope", colour = group)) + geom_text(aes(x = text_x, 
														 							  												  y = text_y, colour = text_group, label = text_eqn), 
														 							  											  data = eqnDF, parse = T, show.legend = FALSE) + 
				theme_classic()
		}
		if (!same.panel) {
			eqnDF <- data.frame(text_eqn = plyr::daply(data, 
													   plyr::as.quoted(group), eqn, form = as.formula(paste(y, 
													   													 x, sep = "~")), method = method), text_x = xpos * 
									max(data[[x]], na.rm = T), text_y = ypos * max(data[[y]], 
																				   na.rm = T))
			eqnDF[[group]] <- rownames(eqnDF)
			plot <- ggplot(data, aes_string(x = x, y = y)) + 
				geom_point(aes_string(...)) + geom_abline(data = subset(regDF, 
																		Method == method), aes_string(intercept = "Intercept", 
																									  slope = "Slope"), col = "red") + geom_text(data = eqnDF, 
																									  										   aes(text_x, text_y, label = text_eqn), parse = TRUE, 
																									  										   show.legend = FALSE) + facet_wrap(as.formula(paste("~", 
																									  										   												   group, sep = ""))) + theme_classic()
		}
	}
	if (mod.line) {
		plot = plot + geom_abline(intercept = 0, slope = 1) + 
			geom_abline(intercept = 0, slope = 2, linetype = 2) + 
			geom_abline(intercept = 0, slope = 0.5, linetype = 2)
	}
	return(plot)
}