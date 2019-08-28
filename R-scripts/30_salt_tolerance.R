


### sigmoid curve



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p +
	stat_function(fun = sigmoid, color = "black") +
	geom_line(aes(x = xes, y = yes))


xes <- seq(-10, 10, length.out = 100)

yes <- sapply(xes, sigmoid)
yes


sigmoid <- function(x){
	1/(1+exp(x*t))
}

t <- 3
x <- seq(-5, 5, 0.01)
plot(x, sigmoid(x), col='blue')

sigmoid2 = function(params, x) {
	params[1] / (1 + exp(-params[2] * (x - params[3])))
}

x = 1:53
y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.18,0.18,0.18,0.33,0.33,0.33,0.33,0.41,
	  0.41,0.41,0.41,0.41,0.41,0.5,0.5,0.5,0.5,0.68,0.58,0.58,0.68,0.83,0.83,0.83,
	  0.74,0.74,0.74,0.83,0.83,0.9,0.9,0.9,1,1,1,1,1,1,1)
ya <- 1/(y + 1)

# fitting code
fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), start=list(a=1,b=.5,c=25))

# visualization code
# get the coefficients using the coef function
params=coef(fitmodel)

y2 <- sigmoid3(params,x)
plot(y2,type="l")
points(y)
points(ya)


sigmoid3 <-  function(params, x) {
	params[1] / (1 + exp(params[2] * (x - params[3])))
}

