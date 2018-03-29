
gen_ccgrid <- function(maxc1sqr,maxc2sqr,n_c1,n_c2){
	c1 <- seq(1,sqrt(maxc1sqr),length = n_c1)
	c2 <- seq(1,sqrt(maxc2sqr),length = n_c2)
	as.data.frame(t(expand.grid(c1 = c1,c2 = c2))) 

}