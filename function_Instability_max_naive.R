
source("./function_gen_ccgrid.R")
source("./function_genW_v5.R")
source("./function_gen_disagree_separately_sum.R")
library(PMA)


InstabilityMax_naive <- function(obs,beta = 0.05,gridc1 = 50,gridc2 = 75, repeatsN = 100){

		w <- W_svd_v5(df = obs,SVDdecomp = FALSE)$w
		######### setting #################
		n <- nrow(w)
		p <- ncol(w)
		c1grid <- gridc1
		c2grid <- gridc2
		subsampling <- floor(0.5 * length(obs$y))
		repeats <- repeatsN
		c1c2grid <- gen_ccgrid(maxc1sqr = n, maxc2sqr = p, n_c1 = c1grid, n_c2 = c2grid)
		
		
		#########################################
		cat("Searchgrid:",length(unique(as.numeric(c1c2grid[1,]))),"x",length(unique(as.numeric(c1c2grid[2,]))),
				"b:",subsampling,"repeats:",repeats,"\n")
		cat("w:",dim(w)[1],"x",dim(w)[2],"\n")
		cat("Generating c1c2grid...")
		cat("Done\n")
		cat("computing Dhat list...\n")
		
		###############  ####################
		cat("generating sample W list...")
		totaltime11 <- Sys.time()
		sampleW_list <- lapply(X = seq(1,repeats), FUN = function(i){
											sapl <- sample(x = seq(1,length(obs$y)),size = subsampling,replace = FALSE)
											sampleobs <- list(y = obs$y[sapl],
																	transcripts = obs$transcripts[sapl,],
																	genes = obs$genes[sapl,])
											w <- W_svd_v5(df = sampleobs,SVDdecomp = FALSE)$w
								})

		totaltime22 <- Sys.time()

		cat("Done ")
		cat("|time: ",
		round(difftime(totaltime22,totaltime11,"auto"),2),
		units(difftime(totaltime22,totaltime11,"auto")),"\n")
		#######################################

		cat("\n1. computing D...\n")
		counts <- 0
		Dlist <- apply(X = c1c2grid, MARGIN = 2, FUN = function(cc){
														totaltime1 <- Sys.time()
														counts <<- counts + 1
														cat("point:",counts) 
														Duvhat <- gen_disagree_UandV_separate_sum(c1 = cc[1], c2 = cc[2], W_list = sampleW_list)
														totaltime2 <- Sys.time()
														cat("|time: ",
																round(difftime(totaltime2,totaltime1,"auto"),3),
																units(difftime(totaltime2,totaltime1,"auto")),"\n")
																
														return(Duvhat)
							})
							
		
		Duvec <- c()
		Dvvec <- c()
		Up1meanvec <- c()
		Vp1meanvec <- c()
		for(i in 1:length(Dlist)){
    	Duvec <- c(Duvec,Dlist[[i]][["DU"]])
    	Dvvec <- c(Dvvec,Dlist[[i]][["DV"]])
    	Up1meanvec <- c(Up1meanvec,Dlist[[i]][["Up1mean"]])
    	Vp1meanvec <- c(Vp1meanvec,Dlist[[i]][["Vp1mean"]])
		}
		Dmaxvec <- pmax(Duvec, Dvvec)
		
		cat(" Done\n")
		###############################################
		
		
		################  ########################
		cat("\n6. Dbar...")
		c1c2posmx <- matrix(seq(1,length(c1c2grid)),c1grid,c2grid)
		c1c2posgrid <- as.data.frame(t(expand.grid(c1pos = seq(1,c1grid),c2pos = seq(1,c2grid))))

							
		Dbarmax <- sapply(X = c1c2posgrid, FUN = function(pp){
								max(Dmaxvec[as.vector(c1c2posmx[1:pp[1],1:pp[2]])])
							})
		

		cat(" Done\n")
		#################################################
		
		
#############  ######################

		c1c2betacandidates <- c1c2grid[,Dbarmax <= beta,drop = FALSE]
		c1c2betaselected <- c1c2betacandidates[,length(c1c2betacandidates)]

##############################################################
		
		PMDout <- PMD(w,type = "standard",sumabs = NULL,
								sumabsu = ifelse(round(as.numeric(c1c2betaselected[1]),5) == round(sqrt(dim(w)[1]),5),
																	as.numeric(c1c2betaselected[1]) - 0.000001,as.numeric(c1c2betaselected[1])),
								sumabsv = ifelse(round(as.numeric(c1c2betaselected[2]),5) == round(sqrt(dim(w)[2]),5),
																	as.numeric(c1c2betaselected[2]) - 0.000001,as.numeric(c1c2betaselected[2])),
								trace = FALSE)
	
	return(list(U_spas = PMDout$u, V_spas = PMDout$v, d_spas = PMDout$d,
							opt_c1 = as.numeric(c1c2betaselected[1]),
							opt_c2 = as.numeric(c1c2betaselected[2]),
							w = w,
							sampleW_list = sampleW_list))

}

