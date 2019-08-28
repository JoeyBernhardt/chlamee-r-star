

predict_monod <- function(df) {
	
	monodcurve<-function(x){
		growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 1000, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(nitrate_concentration.x = x, 
			   growth_rate = preds)
}