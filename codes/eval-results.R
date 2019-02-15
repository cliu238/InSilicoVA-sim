# Combine simulation results
mets <- vector("list", 500)
mets2 <- vector("list", 500)
N <- 1960
C <- 34
for(i in 1:500){
	if(file.exists(paste0("rda/full-HCE-rep", i, ".rda"))){
	load(paste0("rda/full-HCE-rep", i, ".rda"))
	mets[[i]] <- metrics
	}else{
		print(i)
	}
	if(file.exists(paste0("rda/full-noHCE-rep", i, ".rda"))){
	load(paste0("rda/full-noHCE-rep", i, ".rda"))
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
# CSMF accuracy (no HCE): 0.63, 0.72, 0.81
# CSMF accuracy (HCE): 0.66, 0.75, 0.83
quantile(csmfacc, c(0.025, 0.5, 0.975),  na.rm=TRUE)
quantile(csmfacc2, c(0.025, 0.5, 0.975),  na.rm=TRUE)

# InSilicoVA results in Flaxman et al, 2018: 
#   no HCE: 0.5, 2.1, 3.9
#      HCE: 12.6, 13.9, 15.5
# Our results:
#   no HCE: -0.6, 23.8, 48.1
# some HCE: 6.5, 31.4, 52.5
# Tariff 2.0:
#   no HCE: 21.6, 23.1, 24.3
#      HCE: 36.5, 37.6, 38.9    
quantile(cccsmfacc, c(0.025, 0.5, 0.975),  na.rm=TRUE)
quantile(cccsmfacc2, c(0.025, 0.5, 0.975),  na.rm=TRUE)
# InSilicoVA results in Flaxman et al, 2018: 
#   no HCE: 28.3, 28.5, 28.7
#      HCE: 33.9, 34.1, 34.5
# Our results:
#   no HCE: 27.4, 32.6, 37.8
# some HCE: 34.1, 39.7, 44.9          
# Tariff 2.0:
#   no HCE: 37.6, 37.8, 37.9
#      HCE: 50.2, 50.5, 50.7    
quantile(apply(ccc, 1, mean, na.rm=TRUE), c(0.025, 0.5, 0.975))
quantile(apply(ccc2, 1, mean,na.rm=TRUE), c(0.025, 0.5, 0.975))


## Create the latex table
cccm <- apply(ccc, 1, mean, na.rm=TRUE)
cccm2 <- apply(ccc2, 1, mean, na.rm=TRUE)
tab1 <- data.frame(matrix(NA, 4, 6))
tab1[1, 1] <- round(median(cccsmfacc2) * 100, 1)
tab1[1, 2] <- paste0("(", paste(round(quantile(cccsmfacc2, c(0.025, 0.975))*100, 1), collapse = ", "), ")")
tab1[2, 1] <- round(median(cccsmfacc) * 100, 1)
tab1[2, 2] <- paste0("(", paste(round(quantile(cccsmfacc, c(0.025, 0.975))*100, 1), collapse = ", "), ")")
tab1[3, 1] <- round(median(cccm2) * 100, 1)
tab1[3, 2] <- paste0("(", paste(round(quantile(cccm2, c(0.025, 0.975))*100, 1), collapse = ", "), ")")
tab1[4, 1] <- round(median(cccm) * 100, 1)
tab1[4, 2] <- paste0("(", paste(round(quantile(cccm, c(0.025, 0.975))*100, 1), collapse = ", "), ")")
tab1[, 3] <- c(2.1, 13.9, 28.5, 34.1)
tab1[, 4] <- c("(0.5, 3.9)", "(12.6, 15.5)", "(28.3, 28.7)", "(33.9, 34.5)")
tab1[, 5] <- c(23.1, 37.6, 37.8, 50.5)
tab1[, 6] <- c("(21.6, 24.3)", "(36.5, 38.9)", "(37.6, 37.9)", "(50.2, 50.7)")
rownames(tab1) <- c("CCCSMF: no HCE", "CCCSMF:    HCE", "CCC: no HCE", "CCC:    HCE")
colnames(tab1) <- rep(c("median", "95% UI"), 3)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- paste0(paste0('& \\multicolumn{2}{c}{', c("InSilicoVA", "InSilicoVA (reported)", "Tariff 2.0"), '}', collapse=''), '\\\\')

library(xtable)
xtab1 <- xtable(tab1, align="rrrrrrr")
digits(xtab1) <- 1
print(xtab1, add.to.row = addtorow)


# Some more comparisons
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


# Plot: both results v.s. published Tariff 2.0 results 
plot(ccctab[, "Tariff.NoHCE.Median"], ccctab[, "X50.noHCE"] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="Tariff 2.0 from Siena et al. (2015)", 
	ylab = "InSilicoVA", main = "CCC by cause (without HCE)", pch=19, col="blue")
points(ccctab[, "Tariff.NoHCE.Median"], ccctab[, "InSilico.NoHCE.Median"], pch=19, col="orange")
segments(x0=ccctab[, "Tariff.NoHCE.Median"], x1=ccctab[, "Tariff.NoHCE.Median"], y0=ccctab[, "X50.noHCE"] * 100, y1=ccctab[, "InSilico.NoHCE.Median"])
abline(c(0, 1), col = "red", lty=2)


# Plot: results v.s. published InSilicoVA results
plot(ccctab[, "InSilico.HCE.Median"], ccctab[, "X50."] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="InSilicoVA (Tariff 2.0) by Flexman et al. (2018)", 
	ylab = "InSilicoVA", main = "CCC by cause (with HCE)")
abline(c(0, 1), col = "red", lty=2)

# Plot: results v.s. published Tariff 2.0 results 
plot(ccctab[, "Tariff.HCE.Median"], ccctab[, "X50."] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="Tariff 2.0 from Siena et al. (2015)", 
	ylab = "InSilicoVA", main = "CCC by cause (with HCE)")
abline(c(0, 1), col = "red", lty=2)


# Plot: both results v.s. published Tariff 2.0 results 
plot(ccctab[, "Tariff.HCE.Median"], ccctab[, "X50."] * 100, xlim = c(0, 100), ylim = c(0, 100), xlab="Tariff 2.0 from Siena et al. (2015)", 
	ylab = "InSilicoVA", main = "CCC by cause (without HCE)", pch=19, col="blue")
points(ccctab[, "Tariff.HCE.Median"], ccctab[, "InSilico.HCE.Median"], pch=19, col="orange")
segments(x0=ccctab[, "Tariff.HCE.Median"], x1=ccctab[, "Tariff.HCE.Median"], y0=ccctab[, "X50."] * 100, y1=ccctab[, "InSilico.HCE.Median"])
abline(c(0, 1), col = "red", lty=2)

