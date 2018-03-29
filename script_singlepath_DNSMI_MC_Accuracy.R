
load("./selectW_modified.Rdata")
rm(list = setdiff(ls(),"pars"))
mpars <- pars
rm(pars)

currentfolder <- "nw10" # choose between nw1, nw10, nw15, nw1010, nw1015, nw10s15, nw10s25
dir.create(path = paste("./",currentfolder,sep = ""))

###############
if(currentfolder == "nw1"){mpars <- mpars[1]
}else if(currentfolder == "nw15"){mpars <- mpars[15]
}else{mpars <- mpars[10]}

if(currentfolder == "nw1010"){ mpars$nw10$nsample <- 1000
}else if(currentfolder == "nw1015"){ mpars$nw10$nsample <- 1500
}



if(currentfolder == "nw10s15"){
mpars$nw10$ng <- 30
mpars$nw10$ngPrime <- 56
mpars$nw10$nh <- 56
mpars$nw10$ngeneNoise <- 58
mpars$nw10$nx <- 45
mpars$nw10$nxPrime <- 85
mpars$nw10$ns <- 85
mpars$nw10$ntranscriptNoise <- 85
}else if(currentfolder == "nw10s25"){
mpars$nw10$ng <- 50
mpars$nw10$ngPrime <- 49
mpars$nw10$nh <- 49
mpars$nw10$ngeneNoise <- 52
mpars$nw10$nx <- 75
mpars$nw10$nxPrime <- 75
mpars$nw10$ns <- 75
mpars$nw10$ntranscriptNoise <- 75
}
#################

library(parallel)

#### set multi core #################################
		cat("Setting parallel cores...")
		setDefaultCluster(NULL)
		setDefaultCluster(cl <- makePSOCKcluster(detectCores(logical = TRUE) - 0)) 
		cat(length(cl))
		clusterEvalQ(expr = {
  		source("./function_Instability_max_naive.R")
  		source("./function_genW_v5.R")
			source("./function_gen_eff_obs_V2.R")
  		source("./function_kappafrPMD.R")

  		library(PMA)
  		library(stringr)

		})
		cat(" Done\n")
		##############################################


####### record time ################
totaltime1 <- Sys.time()
####################################

############## generte W settings ##################

commands <- list()
for(i in 1:length(mpars)){
				commands[[names(mpars[i])]] <- unlist(mpars[i])
				names(commands[[names(mpars[i])]]) <- names(mpars[[i]])
}

rm(mpars)

######################################################
beta <- 0.05	# This is the instability threshold
MC <- 3				# This is Monte Carlo number, can be set to any positive number
k <- 1				# Do not change this number, keep it as 1
rdataname <- "Accuracy single path"
mx <- currentfolder
set.seed(NULL)
vsn <- paste(MC, " Monte Carlo", " k",k," v",
						sample(x = seq(1,999999),size = 10,replace = FALSE)[sample(seq(1,10),1)],
						" beta",beta,sep = "")
						
############################################
cat("MC:",MC,", k:",k,", vsn:",vsn,"\n")

############### Monte carlo ###########

clusterExport(varlist=c("commands","k","vsn","mx","currentfolder","beta"))

montecarloout <- parSapply(X = seq(1,MC),FUN = function(jj){

	for(i in 1:length(commands)){
		cat("\n",names(commands[i]),"\n")
		
		################# generate obs ###################
		wvsn <- names(commands[i])
		outpt <- gen_eff_obs_v2(
				alpha = as.numeric(commands[[i]]["alpha"]),
				Rsq = as.numeric(commands[[i]]["Rsq"]),
				nsample = as.numeric(commands[[i]]["nsample"]),
				rind = as.numeric(commands[[i]]["rind"]),
				nx = as.numeric(commands[[i]]["nx"]),
				nxPrime = as.numeric(commands[[i]]["nxPrime"]),
				ns = as.numeric(commands[[i]]["ns"]),
				ntranscriptNoise = as.numeric(commands[[i]]["ntranscriptNoise"]),
				ng = as.numeric(commands[[i]]["ng"]),
				ngPrime = as.numeric(commands[[i]]["ngPrime"]),
				nh = as.numeric(commands[[i]]["nh"]),
				ngeneNoise = as.numeric(commands[[i]]["ngeneNoise"]),
				noiseSD = as.numeric(commands[[i]]["noiseSD"])
				)
		
		####### observed #########
		obs <- outpt$Obs
}
	
	cat("permute",jj,"\n")

	checkfile <- paste("./",currentfolder,"/",mx,"_",vsn," MC check.txt",sep="")
	#######################
		


for(kk in 1:k){
	cat("k:",kk,"\n")
	write(x = paste("MC:",jj,"k:",kk,"w dim:",ncol(obs$genes),"x",ncol(obs$transcripts),"\n"),
				file = checkfile,append = TRUE)

	if(kk == 1){
		sparseout <- InstabilityMax_naive(obs = obs,beta = beta,gridc1 = 50,gridc2 = 75, repeatsN = 100)
	
		
		# gridc1 is the number of search gird on c1 marginal dimension
		# gridc2 is the number of search gird on c2 marginal dimension
		# repeatsN is the number of subsamples for DNSMI
		
		
		
		accuracyU <- kappafrPMD(PMDoutvec = sparseout$U_spas,trueposinvec = seq(1,outpt$Input$actvGcolmn))
		Ukappa1 <- accuracyU$Kappa
		Usensi <- accuracyU$Sensi
		Uspeci <- accuracyU$Speci
		
		
		accuracyV <- kappafrPMD(PMDoutvec = sparseout$V_spas,trueposinvec = seq(1,outpt$Input$actvXcolmn))
		Vkappa1 <- accuracyV$Kappa
		Vsensi <- accuracyV$Sensi
		Vspeci <- accuracyV$Speci
		
		
		nonzeroU_N <- sum(sparseout$U_spas != 0)
		nonzeroV_N <- sum(sparseout$V_spas != 0)
	}


}

write(x = paste("MC:",jj,"completed",kk,"\n"), file = checkfile,append = TRUE)


return(c(Ukappa1 = Ukappa1,Vkappa1 = Vkappa1,UVkappaSum1 = Ukappa1 + Vkappa1,
				Usensi = Usensi, Uspeci = Uspeci, Vsensi = Vsensi, Vspeci = Vspeci,
					SelectU_N = nonzeroU_N, SelectV_N = nonzeroV_N))
})





save.image(file = paste("./",currentfolder,"/",
									rdataname,"_",mx,"_",vsn,".Rdata",sep=""))


########## time used ###############
totaltime2 <- Sys.time()
cat("\nTotal time used: ",
		round(difftime(totaltime2,totaltime1,"auto"),2),
		units(difftime(totaltime2,totaltime1,"auto")),"\n")
###################################


