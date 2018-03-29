
setPar <- function(alpha, Rsq)
{
	sg <- 1  
	
	sh <- sg    
	bxg <- sqrt(alpha)
	byx.g <- sqrt(alpha)
	byg.x <- 1 - alpha
	bys <- alpha
	byh <- alpha
	tau <- (1 - Rsq)/Rsq
	sxg <- alpha * tau   
	ss <- sxg + sg * alpha    
	
	syxgsh <- tau + (alpha ^ 2) * tau * (1 + tau) * (1 + alpha) 
	list(alpha=alpha, Rsq=Rsq, sg=sg, ss=ss, sh=sh, 
			bxg=bxg, 	byx.g=byx.g, byg.x=byg.x, bys=bys, byh=byh, 
			tau=tau, sxg=sxg, syxgsh=syxgsh)
}