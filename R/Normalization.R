
NormalizeLog2 <- function(mat, mean.center = T) {
	log2Scaled <- log2(mat)
	for (i in 1:nrow(log2Scaled)) {
		ind <- which(is.finite(log2Scaled[i,]) == FALSE)
		log2Scaled[i,ind] <- 0

		if (mean.center) {
			log2Scaled[i,] <- log2Scaled[i,]-mean(log2Scaled[i,])
		}
	}

	return(as.matrix(log2Scaled))
}