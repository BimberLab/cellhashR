# Reverse the vector x and return the value at the Nth index. If N is larger
# than the length of the vector, return the last value in the reversed vector.
#
# @param x vector of interest
# @param N index in reversed vector
#
# @return returns element at given index
#
MaxN <- function(x, N = 2){
	len <- length(x)
	if (N > len) {
		warning('N greater than length(x).  Setting N=length(x)')
		N <- length(x)
	}
	sort(x, partial = len - N + 1)[len - N + 1]
}

# Local maxima estimator
#
# Finding local maxima given a numeric vector
#
# @param x A continuous vector
#
# @return Returns a (named) vector showing positions of local maximas
#
# @author Tommy
# @references \url{https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima}
#
# @examples
# x <- c(1, 2, 9, 9, 2, 1, 1, 5, 5, 1)
# LocalMaxima(x = x)
#
LocalMaxima <- function(x) {
	# Use -Inf instead if x is numeric (non-integer)
	y <- diff(x = c(-.Machine$integer.max, x)) > 0L
	y <- cumsum(x = rle(x = y)$lengths)
	y <- y[seq.int(from = 1L, to = length(x = y), by = 2L)]
	if (x[[1]] == x[[2]]) {
		y <- y[-1]
	}
	return(y)
}

GetPlotColors <- function(total, palette = 'Set1') {
 	return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, total))), "Set1")(total))
}

GetTotalPlotPages <- function(totalValues, perPage = 4) {
	return(ceiling(totalValues / perPage))
}

SimplifyHtoNames <- function(v) {
	return(sapply(v, function(x){
		x <- unlist(strsplit(x, '-'))
		if (length(x) > 1) {
			x <- x[-(length(x))]
		}

		paste0(x, collapse = "-")
	}))
}
