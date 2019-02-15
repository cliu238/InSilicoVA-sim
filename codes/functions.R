getCCC <- function(est, truth, causes){
	C <- length(causes)
	N <- length(truth)
	ccc <- rep(NA, C)
	for(i in 1:C){
		TP <- length(intersect(which(est == causes[i]), which(truth == causes[i])))
		FN <- length(intersect(which(est != causes[i]), which(truth == causes[i])))
    	ccc[i] = (TP / (TP + FN) - 1/C) / (1 - 1/C)
	}
	return(ccc)
}

# just in case, this function returns the removed symptoms in the PHMRC shortened questionnaire.
# However, for this analysis, the full questionnaire was used in Tariff 2.0, so this is not needed
getRemoved <- function(){
c("a2_01", "a2_03", "a2_06", "a2_08", "a2_12", "a2_16", "a2_17", "a2_18", "a2_19", "a2_20", "a2_23", "a2_24", "a2_28", "a2_33", "a2_37", "a2_38", "a2_38_s1", "a2_40", "a2_41", "a2_42", "a2_45", "a2_46", "a2_46a", "a2_46a_s1", "a2_46b", "a2_48", "a2_49", "a2_54", "a2_69", "a2_70", "a2_71", "a2_76", "a2_78", "a2_79", "a2_80", "a2_81", "a2_86", "a3_20", "a4_03")
}

