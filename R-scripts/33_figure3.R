
all_changes <- read_csv("data-processed/all_changes.csv")

all_changes_anc_sum <- all_changes %>% 
	group_by(ancestor_id) %>% 
	summarise_each(funs(mean, std.error), change_istar_mean, change_nstar_mean, change_pstar_mean, change_salt_tol_mean)



plot1a <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_pstar_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_pstar_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_pstar_mean_mean + change_pstar_mean_std.error,
					  ymin = change_pstar_mean_mean - change_pstar_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_pstar_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in P* (um P)") + 
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')')) 

plot1b <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_nstar_mean, y = change_pstar_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_nstar_mean_mean, y = change_pstar_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_nstar_mean_mean, ymax = change_pstar_mean_mean + change_pstar_mean_std.error,
					  ymin = change_pstar_mean_mean - change_pstar_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_nstar_mean_mean - change_nstar_mean_std.error, y = change_pstar_mean_mean,
					   xmax = change_nstar_mean_mean + change_nstar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in P* (um P)") + 
	xlab("Change in N* (um N)") 

plot1c <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=-2, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=2),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_nstar_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_nstar_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_nstar_mean_mean + change_nstar_mean_std.error,
					  ymin = change_nstar_mean_mean - change_nstar_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_nstar_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in N* (um N)") + 
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')'))

plot1d <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_istar_mean, y = change_salt_tol_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_istar_mean_mean, y = change_salt_tol_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_istar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_istar_mean_mean - change_istar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_istar_mean_mean + change_istar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab(bquote('Change in I* ('*mu~'mol'~ m^-2~s^-1*')'))

plot1e <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_pstar_mean, y = change_salt_tol_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_pstar_mean_mean, y = change_salt_tol_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_pstar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_pstar_mean_mean - change_pstar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_pstar_mean_mean + change_pstar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab("Change in P* (uM P)")

plot1f <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_nstar_mean, y = change_salt_tol_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_nstar_mean_mean, y = change_salt_tol_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_nstar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_nstar_mean_mean - change_nstar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_nstar_mean_mean + change_nstar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "none", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab("Change in N* (uM N)")

plot_legend <- ggplot() +
	geom_rect(aes(xmin=0, xmax=2.3, ymin=4, ymax=0),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_rect(aes(xmin=-2.3, xmax=0, ymin=0, ymax=-4),
			  color="transparent", alpha=0.5, fill = "grey") +
	geom_point(aes(x = change_nstar_mean, y = change_salt_tol_mean, color = ancestor_id), data = all_changes, size = 2) +
	
	geom_point(size = 4, aes(x = change_nstar_mean_mean, y = change_salt_tol_mean_mean, color = ancestor_id), data = all_changes_anc_sum) +
	geom_errorbar(aes(x = change_nstar_mean_mean, ymax = change_salt_tol_mean_mean + change_salt_tol_mean_std.error,
					  ymin = change_salt_tol_mean_mean - change_salt_tol_mean_std.error, color = ancestor_id),
				  data = all_changes_anc_sum, width = 0) +
	geom_errorbarh(aes(xmin = change_nstar_mean_mean - change_nstar_mean_std.error, y = change_salt_tol_mean_mean,
					   xmax = change_nstar_mean_mean + change_nstar_mean_std.error, color = ancestor_id),
				   data = all_changes_anc_sum) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	scale_color_manual(values = c(cols_no_anc), name = "") +
	theme(legend.position = "bottom", legend.direction = "horizontal") +
	ylab("Change in salt tolerance (ug/L)") + 
	xlab("Change in N* (uM N)")

legend <- get_legend(plot_legend)

multi_plot_trade_offs <- plot_grid(plot1a, plot1b, plot1c, plot1d, plot1e, plot1f, legend, labels = c("A", "B", "C", "D", "E", "F"), align = "h", ncol = 2, nrow = 4)
save_plot("figures/trait_trade_offs2_anc.png", multi_plot_trade_offs,
		  ncol = 2, # we're saving a grid plot of 2 columns
		  nrow = 4, # and 2 rows
		  # each individual subplot should have an aspect ratio of 1.3
		  base_aspect_ratio = 1.1
)
