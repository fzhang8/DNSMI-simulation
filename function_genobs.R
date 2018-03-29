source('./function_calsigma.R')
source("./function_genlat.R")
source("./function_setpar.R")

genIndicators <- function(alpha, Rsq, nsample, 
													rind, nx, nxPrime, ns, ntranscriptNoise, 
													ng, ngPrime, nh, ngeneNoise, noiseSD){
	
	set.seed(NULL)
	
	set.seed(sample(x = seq(1,20000),size = 1))
	tauInd <- (1 - rind^2)/rind^2
	dfg <- NULL
	dfy <- NULL
	dfx <- NULL
	for (i in 1:nsample){
		
		sigma <- calSigma(setPar(alpha, Rsq))
		latent <- genLatent(sigma)
		
		y <- rnorm(1, latent["Yprime"], sqrt(tauInd * sigma["Yprime", "Yprime"]))
		
		
		x <- c(rnorm(nx, latent["X"], sqrt(tauInd * sigma["X", "X"])), 
						rnorm(nxPrime, latent["Xprime"], sqrt(tauInd * sigma["Xprime", "Xprime"])), 
						rnorm(ns, latent["S"], sqrt(tauInd * sigma["S", "S"])), 
						rnorm(ntranscriptNoise, 0, noiseSD * sqrt(tauInd * sigma["Yprime", "Yprime"])))
		
		
		g <- c(rnorm(ng, latent["G"], sqrt(tauInd * sigma["G", "G"])), 
						rnorm(ngPrime, latent["Gprime"], sqrt(tauInd * sigma["Gprime", "Gprime"])), 
						rnorm(nh, latent["H"], sqrt(tauInd * sigma["H", "H"] )), 
						rnorm(ngeneNoise, 0, noiseSD * sqrt(tauInd * sigma["Yprime", "Yprime"])))
		
		
		dfy <- rbind(dfy, y)
		dfx <- rbind(dfx, x)
		dfg <- rbind(dfg, g)
	}
	truegenes <- c(rep(1, ng), rep(0, ngPrime + nh) )
	truetranscripts <- c(rep(1, nx), rep(0, nxPrime + ns) )
	df <- list(y = dfy - mean(y),
						transcripts = scale(dfx, scale=FALSE), # only centered
						genes = scale(dfg, scale=FALSE), 			# only centered
						truetranscripts = truetranscripts, 
						truegenes = truegenes)
						
	dimnames(df$transcripts) <- NULL
	dimnames(df$genes) <- NULL
	dimnames(df$y) <- NULL
	set.seed(NULL)
	return(df)
}