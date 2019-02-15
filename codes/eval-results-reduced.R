# Combine simulation results
mets <- vector("list", 500)
mets2 <- vector("list", 500)
N <- 1960
C <- 34
for(i in 1:500){
	if(file.exists(paste0("rda-r/reduced-HCE-rep", i, ".rda"))){
	load(paste0("rda-r/reduced-HCE-rep", i, ".rda"))
	mets[[i]] <- metrics
	}else{
		print(i)
	}
	if(file.exists(paste0("rda-r/reduced-noHCE-rep", i, ".rda"))){
	load(paste0("rda-r/reduced-noHCE-rep", i, ".rda"))
	mets2[[i]] <- metrics
	}else{
		print(i+500)
	}
}

csmfacc <- cccsmfacc <- csmfacc2 <- cccsmfacc2 <- rep(NA, 500)
ccc <- ccc2 <- matrix(NA, 500, 34)

causes <- names(mets[[1]]$ccc)
for(i in 1:500){
	if(length(mets[[i]]) > 0){
		# with HCE
		csmfacc[i] <- mets[[i]]$csmfacc
		cccsmfacc[i] <- mets[[i]]$cccsmfacc
		ccc[i, ] <- mets[[i]]$ccc[causes]
	}
	if(length(mets2[[i]]) > 0){
		# without HCE
		csmfacc2[i] <- mets2[[i]]$csmfacc
		cccsmfacc2[i] <- mets2[[i]]$cccsmfacc
		ccc2[i, ] <- mets2[[i]]$ccc[causes]
	}
}
# CSMF accuracy (HCE): 0.66, 0.76, 0.83
# CSMF accuracy (no HCE): 0.62, 0.72, 0.81
quantile(csmfacc, c(0.025, 0.5, 0.975),  na.rm=TRUE)
quantile(csmfacc2, c(0.025, 0.5, 0.975),  na.rm=TRUE)

# InSilicoVA results in Flaxman et al, 2018: 
#   no HCE: 0.5, 2.1, 3.9
#      HCE: 12.6, 13.9, 15.5
# Our results:
#   no HCE: -4.3, 23.2, 47.5
# some HCE: 7.8, 33.6, 53.0
# Tariff 2.0:
#   no HCE: 21.6, 23.1, 24.3
#      HCE: 36.5, 37.6, 38.9    
quantile(cccsmfacc, c(0.025, 0.5, 0.975),  na.rm=TRUE)
quantile(cccsmfacc2, c(0.025, 0.5, 0.975),  na.rm=TRUE)
# InSilicoVA results in Flaxman et al, 2018: 
#   no HCE: 28.3, 28.5, 28.7
#      HCE: 33.9, 34.1, 34.5
# Our results:
#   no HCE: 28.5, 33.3, 37.5
# some HCE: 33.7, 38.8, 43.7           
# Tariff 2.0:
#   no HCE: 37.6, 37.8, 37.9
#      HCE: 50.2, 50.5, 50.7    
quantile(apply(ccc, 1, mean), c(0.025, 0.5, 0.975),  na.rm=TRUE)
quantile(apply(ccc2, 1, mean), c(0.025, 0.5, 0.975),  na.rm=TRUE)

ccctab <- t(apply(ccc, 2, quantile, c(0.025, 0.5, 0.975),  na.rm=TRUE))
ccctab2 <- t(apply(ccc2, 2, quantile, c(0.025, 0.5, 0.975),  na.rm=TRUE))
ccctab <- data.frame(ccctab)
ccctab2 <- data.frame(ccctab2)
rownames(ccctab) <- causes
ccctab$X50.noHCE <- ccctab2$X50.


compare <- read.csv("../Flexman_2018/Table3.csv")
ccctab$cause <- tolower(causes)
ccctab$cause[ccctab$cause == "tb"] <- "tuberculosis"
compare[,1] <- tolower(compare[,1])
sum(ccctab$cause %in% compare[,1])

# Plot: results v.s. published InSilicoVA results
ccctab <- merge(ccctab, compare, by.x = "cause", by.y = "PHMRC.cause")
plot(ccctab[, "InSilico.NoHCE.Median"], ccctab[, "X50.noHCE"] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="InSilicoVA (Tariff 2.0) by Flexman et al. (2018)", 
	ylab = "InSilicoVA", main = "CCC by cause (without HCE)")
abline(c(0, 1), col = "red", lty=2)

# Plot: results v.s. published Tariff 2.0 results 
plot(ccctab[, "Tariff.NoHCE.Median"], ccctab[, "X50.noHCE"] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="Tariff 2.0 from Siena et al. (2015)", 
	ylab = "InSilicoVA", main = "CCC by cause (without HCE)")
abline(c(0, 1), col = "red", lty=2)


# Plot: results v.s. published InSilicoVA results
plot(ccctab[, "InSilico.HCE.Median"], ccctab[, "X50."] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="InSilicoVA (Tariff 2.0) by Flexman et al. (2018)", 
	ylab = "InSilicoVA", main = "CCC by cause (with HCE)")
abline(c(0, 1), col = "red", lty=2)

# Plot: results v.s. published Tariff 2.0 results 
plot(ccctab[, "Tariff.HCE.Median"], ccctab[, "X50."] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="Tariff 2.0 from Siena et al. (2015)", 
	ylab = "InSilicoVA", main = "CCC by cause (with HCE)")
abline(c(0, 1), col = "red", lty=2)
