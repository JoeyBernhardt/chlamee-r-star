

#### Huisman and Weissing 1994 light limited growth model
library(simecol)


dWRdt <-  function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dW1 = min((r1*R*W1/(M1 + R)), ((k1*W1/(k1*W1 + k2*W2)) * (pmax1/k1) * (log((H1 + Iin)/(H1 + Iout))))) - l1*W1 
			dW2 = min((r2*R*W2/(M2 + R)), ((k2*W2/(k1*W1 + k2*W2)) * (pmax2/k2) * (log((H2 + Iin)/(H2 + Iout))))) - l2*W2 #algae
			dR = l1*(Rin-R) - c1*W1 - c2*W2	#resource
			dIout = Iin*exp(-(k1*W1 + k2*W2))
			return(list(c(dW1, dW2, dR, dIout)))			
		}
	)
}

dWRdt <-  function(t, state, parameters) {
	with(
		as.list(c(state,parameters)),{
			dW2 = (r2*R*W2/(M2 + R)) - l2*W2 
			dW1 = ((k1*W1/(k1*W1 + k2*W2)) * (pmax1/k1) * (log((H1 + Iin)/(H1 + Iout)))) - l1*W1 #algae
			dR = l1*(Rin-R) - c1*W1 - c2*W2	#resource
			dIout = Iin*exp(-sum(k1*W1,k2*W2))
			return(list(c(dW1, dW2, dR, dIout)))			
		}
	)
}

# dR = D * (S - R) - Q * (dN + m * N + D * N)

state <-  c(W1 = 1, W2 = 1, R = 1000, Iout = 1000) ## initial conditions
time <-  seq(from = 0, to = 100, by = 0.001) 
parameters <-  c(Rin = 20, Iin = 20, pmax1 = 8, pmax2 = 10, H1 = 100, H2 = 100, r1 = 2, r2 = 1, M1 = 10, M2 = 10, l1 = 0.5, l2 = 0.5, c1 = 0.02, c2 = 0.04, k1 = 0.008, k2 = 0.004)

hold <- ode(y = state, times = time, func = dWRdt, parms = parameters)
hold2 <- as.data.frame(hold) %>%
	gather(key = type, value = value, 2:5)

hold3 <- as.data.frame(hold) 


hold2 %>% 
	ggplot(aes(x = time, y = value, color = type)) + geom_line() +
	facet_wrap( ~ type, scales = "free")
 