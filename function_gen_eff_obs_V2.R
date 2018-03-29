source("./function_genobs.R")

gen_eff_obs_v2 <- function(alpha = 0.8,
												Rsq = 0.7,
												nsample = 1000,
												rind = 0.8,
												nx = 60,
												nxPrime = 50,
												ns = 50,
												ntranscriptNoise = 50,
												ng = 50,
												ngPrime = 50,
												nh = 50,
												ngeneNoise = 50,
												noiseSD = 1){
	cat("generating obs...")


	input <- list(alpha = alpha,
								Rsq = Rsq,
								nsample = nsample,
								rind = rind,
								nx = nx,
								nxPrime = nxPrime,
								ns = ns,
								ntranscriptNoise = ntranscriptNoise,
								ng = ng,
								ngPrime = ngPrime,
								nh = nh,
								ngeneNoise = ngeneNoise,
								noiseSD = noiseSD)
							
		
	input <- c(input, Gcolmn = input$ng + input$ngPrime + input$nh + input$ngeneNoise,
									actvGcolmn = input$ng,
									Xcolmn = input$nx + input$nxPrime + input$ns + input$ntranscriptNoise,
									actvXcolmn = input$nx)

	
	obs <- genIndicators(alpha = input$alpha,
															Rsq = input$Rsq,
															nsample = input$nsample,
															rind = input$rind,
															nx = input$nx,
															nxPrime = input$nxPrime,
															ns = input$ns,
															ntranscriptNoise = input$ntranscriptNoise,
															ng = input$ng,
															ngPrime = input$ngPrime,
															nh = input$nh,
															ngeneNoise = input$ngeneNoise,
															noiseSD = input$noiseSD)

	cat("Done\n")
	return(list(Obs = obs, Input = input))
}