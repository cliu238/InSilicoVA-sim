# load packages
library(openVA)
library(MCMCpack)
source("functions.R")

# read PHMRC data from local CSV file
PHMRC <- read.csv("../PHMRC/IHME_PHMRC_VA_DATA_ADULT_Y2013M09D11_0.csv")
causes <- as.character(unique(PHMRC$gs_text34))
N <- dim(PHMRC)[1] # 7841
N.train <- round(N * 0.75) # 5881
N.test <- N - N.train # 1960

# number of causes and list of causes
C <- 34 
causes <- as.character(unique(PHMRC$gs_text34))

# set seed and repeat 500 times with and without HCE each
metrics <- NULL
for(rep in 1:1000){
	if(rep > 500){
		rep <- rep - 500
		has.HCE <- FALSE
	}else{
		has.HCE <- TRUE
	}
	set.seed(123 + as.integer(rep))

	is.train <- sample(1:N, N.train)
	test <- PHMRC[-is.train, ]
	train <- PHMRC[is.train, ]

	## resample within test data
	#  obtain a new CSMF
	csmf.resample <- rdirichlet(1, rep(1, C))
	names(csmf.resample) <- causes

	#  remove causes if not exist in test data
	noexist.causes <- which(causes %in% test$gs_text34 == FALSE)
	if(length(noexist.causes) > 0){
		csmf.resample[noexist.causes] <- 0
		csmf.resample <- csmf.resample / sum(csmf.resample)
	}
	#  find the new list of causes
	causes.test <-  sample(size = dim(test)[1], x = 1:length(causes), replace = TRUE, prob = csmf.resample)
	#  create a pool of observations for each cause
	testlist <- NULL
	for(i in 1:C){
	  testlist[[i]] <- which(test$gs_text34 == causes[i])
	} 
	#  create list of observations ~ new CSMF
	test.index <- rep(NA, dim(test)[1])
	for(i in 1:dim(test)[1]){
	  test.index[i] <- sample(testlist[[causes.test[i]]], 1)
	}
	test <- test[test.index, ]
	csmf.true <- table(test$gs_text34)
	csmf.true <- csmf.true[causes] / sum(csmf.true)

	if(!has.HCE){
		binary <- ConvertData.phmrc(train, test, phmrc.type = "adult", cause = "gs_text34")
		reduced <- paste0("a1_01_", 1:14)
		binary$output <- binary$output[, -which(colnames(binary$output) %in% reduced)]
		binary$output.test <- binary$output.test[, -which(colnames(binary$output.test) %in% reduced)]
		## run InSilicoVA
		insilico <- codeVA(data = binary$output.test, 
					data.type = "customize", 
					model = "InSilicoVA",
	                data.train = binary$output, 
	                causes.train = "Cause", 
	                jump.scale = 0.05, 
	                convert.type = "fixed",
	                Nsim=10000, auto.length = FALSE)

	}else{
		 insilico <- codeVA(data = test,
                    data.type = "PHMRC",
                    model = "InSilicoVA",
                    data.train = train,
                    causes.train = "gs_text34",
                    phmrc.type = "adult",
                    jump.scale = 0.05,
                    convert.type = "fixed",
                    Nsim=10000, auto.length = FALSE)
	}

	
	## get CSMF and individual assignment
	csmf <- getCSMF(insilico)[, "Mean"]
	csmf <- csmf[causes]
	topcause <- getTopCOD(insilico)
	## get CSMF accuracy
	csmfacc <- getCSMF_accuracy(csmf= csmf, truth = csmf.true)
	cccsmfacc <- (csmfacc - (1 - exp(-1))) / (1 - (1 - exp(-1)))
	## get CCC
	ccc <- getCCC(est = topcause[,2], truth = as.character(test$gs_text34), causes = causes)
	names(ccc) <- causes
	metrics <- list(ccc = ccc, csmfacc = csmfacc, cccsmfacc = cccsmfacc)
	if(has.HCE){
		save(metrics, file = paste0("rda/full-HCE-rep", rep, ".rda"))
	}else{
		save(metrics, file = paste0("rda/full-noHCE-rep", rep, ".rda"))
	}
}