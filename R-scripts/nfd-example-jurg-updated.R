# @author: J.W.Spaak
# Example how to compute the ND and FD for a given differential equation setting

# The script is originally written in python, we run python code from within R
# make sure that python is installed on the system used
# This is not part of the install.packages("reticulate") command!
library(reticulate)
library(tidyverse)
library(rootSolve) # to numerically solve a function


# loads the relevant python code
source_python("python/numerical_NFD.py")

## read in the r star and stoichiometry data

rstars <- read_csv("data-processed/all-rstars.csv")
stoich <- read_excel("data-raw/stoichiometry/ChlamEE_140218_stoich_growthrates.xlsx") %>% 
	clean_names() %>% 
	mutate(cN = b_n/b_c) %>% 
	mutate(cP = b_p/b_c) %>% 
	mutate(cNb = b_n/biomass) %>% 
	mutate(cPb = b_p/biomass) %>% 
	mutate(nc = 1/cto_n_molar) %>% ## this is the consumption vector for n
	mutate(pc = 1/cto_p_molar) %>% ## this is the consumption vector for p
	mutate(c_diff = abs(cN - cP)) %>% 
	mutate(ancestor = ifelse(ancestor == "CC 1690", "cc1690", ancestor)) %>% 
	select(cN, cP, c_diff, everything()) %>% 
	mutate(selection_treatment = ifelse(is.na(selection_treatment), "none", selection_treatment)) %>% 
	mutate(ancestor = ifelse(is.na(ancestor), strain_id, ancestor)) %>% 
	select(ancestor, selection_treatment, nc, pc)

stoich2 <- stoich %>% 
	rename(ancestor_id = ancestor, 
		   treatment = selection_treatment,
		   c_N = nc,
		   c_P = pc) %>% 
	mutate(ancestor_id = case_when(ancestor_id == "Anc 2" ~ "anc2",
								   ancestor_id == "Anc 3" ~ "anc3",
								   ancestor_id == "Anc 4" ~ "anc4",
								   ancestor_id == "Anc 5" ~ "anc5",
								   ancestor_id == "cc1690" ~ "cc1690"))
rstars2 <- rstars %>% 
	select(ancestor_id, treatment, p_umax, n_umax, n_ks, p_ks) %>% 
	rename(r_N = n_umax,
		   r_P = p_umax,
		   H_N = n_ks,
		   H_P = p_ks)

all_resource_traits <- left_join(rstars2, stoich2)

snippet <- all_resource_traits %>% 
	filter(treatment %in% c("N", "P"), ancestor_id %in% c("anc4"))

fun <- function(S_P, S_N, df){
	n_spec <- 2 # number of species in the system, must be an integer
	set.seed(0) # set random seed for reproduce ability

	spec <- list()
# df <- snippet
# S_P <- 50
# S_N <- 1000
# parameters for nitrogen competition
spec$c_N <- c(df$c_N[1], df$c_N[2])# consumption of nitrogen
spec$H_N <- c(df$H_N[1], df$H_N[2]) # halfsaturation constant of nitrogen
spec$r_N <- c(df$r_N[1], df$r_N[2]) # maximal growth rate if limited by nitrogen

# parameters for phosphate competition
spec$c_P <- c(df$c_P[1], df$c_P[2]) # identical to nitrogen definitions
spec$H_P <- c(df$H_P[1], df$H_P[2])
spec$r_P <- c(df$r_P[1], df$r_P[2]) 

#parameters for light competition
spec$H_light <- c(0.1,0.1) # halfsaturation constant for light
spec$k <- c(0.008,0.004) # light absorption coefficient
# spec$r_light <- c(1,5) # maximal growth rate for light -- JB note: why is this so much higher?
spec$r_light <- c(100,100) 
# mortality rate of the species
spec$m <- 0.50# assumed to be the washout rate, i.e similar for both

# environmental parameters
res <- list()
res$S_N <- S_N # supply rate of nitrogen
res$S_P <- S_P # supply rate of phosphate

res$d <- spec$m # assumed to be the washout rate, similar to species death rate
res$I_in <- 10000000 # incoming light intensity
res$z <- 1 # mixing depth

# compute the growth rate of the species
growth_function <- function(X, P, N, spec, res, growth_light){
	growth_P <- spec$r_P*P/(P+spec$H_P) # how much they would grow if P limited
	growth_N <- spec$r_N*N/(N+spec$H_N)
	
	# take the minima of all growth rates
	growth <- c(min(growth_P[1], growth_N[1], growth_light[1]),
				min(growth_P[2], growth_N[2], growth_light[2]))
	
	# compute whether N and P are in equilibrium
	change_N <- res$d * (res$S_N-N) - sum(X * spec$c_N * growth)
	change_P <- res$d * (res$S_P-P) - sum(X * spec$c_P * growth)
	return (c(change_N, change_P, growth))
}

# compute the densities of N and P in the experiment
find_N_P <- function(R,parms){
	# extract parameters
	X <- parms$X
	spec <- parms$spec
	res <- parms$res
	growth_light <- parms$growth_light
	P <- R[1]
	N <- R[2]
	return(growth_function(X,P,N,spec,res,growth_light)[1:2])
}

test_f <- function(X){
	# compute growth via light
	if (sum(spec$k*X) > 0){ # species are actually absorbing
		I_out <- res$I_in*exp(-sum(res$z*spec$k*X))
		
		# total light absorbed and converted into energy
		growth_light_tot <- spec$r_light/(spec$k*res$z)*log((res$I_in + spec$H_light)/(I_out + spec$H_light))
		growth_light <- growth_light_tot*spec$k/sum(spec$k*X)
	} else {
		growth_light <- spec$r_light*res$z*res$I_in/(spec$H_light+res$I_in)
	}
	
	# compute, using numerical solving, the densities of N and P
	P_N <- multiroot(find_N_P, c(res$S_P,res$S_N),
					 parms = list(X = X, spec = spec, res = res, growth_light = growth_light))
	P_N <- P_N$root
	
	# use N and P levels to obtain growth rates
	growth <- growth_function(X,P_N[1], P_N[2],spec,res,growth_light)[3:4]
	return (growth - spec$m)
}

# to give an impression whether function works well
test_f(c(1,1))
test_f(c(0,0)) # initial growth rate, should be positive

# compute relevant parameters with python
# the parameter `from_R = TRUE` changes data types from R to python
# the parameters args is to pass additional parameters to the function in python
# in case this seems to be causing problems, just delete args and the spec and res entries in test_f
pars <- NFD_model(test_f, n_spec, from_R = TRUE)
ND <- pars$ND
NO <- pars$NO
FD <- pars$FD
c <- pars$c
results <- data.frame(ND = ND, FD = FD)
return(results)
}


fun(20.10, 50.13, snippet)
fun(20, 20)

results <- data.frame()
for (S_P in seq(0,50, by = 1.13)){
	for (S_N in seq(10,350, by = 5.13)){
		hold <- data.frame(tryCatch(fun(S_P = S_P, S_N = S_N, df = snippet), 
						 error = function(cond){
						 	return(NA)
						 }))
		hold$S_P <- S_P
		hold$S_N <- S_N
		results <- bind_rows(results, hold)
	}
}


library(cowplot)

names(results)
results %>% 
	select(ND, FD, S_P, S_N) %>% 
	filter(!is.na(ND)) %>% View 
	ggplot(aes(x = S_P, y = S_N, color = ND)) + geom_point(size = 4) +
	scale_color_viridis_c()
