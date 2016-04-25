#' Selecting the best 2-Band combinations for Normalized Difference Vegetation Index (NDVI)
#'
#' This function developes a optimization algorithm based on correlation analysis between spectral parameter 'X' and y, which
#' determines the best spectral bands of the full spectrum that are most predictive for 'y'.
#'
#' @param x data set of spectral matrix data
#' @return bands(b1 b2) for NDVI = (b1 - b2)/(b1 + b2)
#' @details
#' This function runs, throughly, a calculation of NDVI = (Xi-Xj)/(Xi+Xj) using all the possible pairs/combinations of any two bands (Xi,Xj)
#' within the full spectrum range. A correlation analysis is then performed between the y and all possible NDVIs (Xi,Xj), and it calculates
#' the correlation coefficients (r) which indicates the predictive performance of each NDVI and its corresponding two-band combination. The
#' output is the wavelength (nm) indicating the best two bands that produce the highest value of r.
#' @seealso \code{cor}
#' @examples
#' ndvi()
#' @export

ndvi <- function(){

  library(ggplot2)
  library(Matrix)
  library(reshape2)
  library(colorRamps)
  library(RColorBrewer)
  # library(scales)

  #  ------------------------------------------------------------------------

  data <- read.csv('exampleData.csv')[1:20,1:952]
  wavelength <- as.numeric(gsub("X","",names(data)[(1+1):952]))

  data <- data.matrix(data) # devtools::use_data(data)

  # data <- matrix(as.numeric(unlist(data)),nrow=nrow(data))

  aa <- matrix(data[, 1]) # Variable of interest, e.g., Chl, N, LAI
  ss <- data[, 2:ncol(data)]; # Reflectance spectra

  n <- dim(ss)[2] # Returns the Number of wavebands

  # (Rj-Ri)/(Rj+Ri)

  R2 <- Matrix(0,n,n,sparse = TRUE)  # Zero sparse matrix to storage Z coordinates

  Rj <- ss

  ones <- matrix(1,1,n)

  for (cI in 1:n){
    Ri <- ss[,cI]
    Ri <- Ri %*% ones  # Turn Ri values to dimensional matrix
    # VI formular
    V <- (Rj-Ri)/(Rj+Ri)

    # Squared values (R2) of the Pearson Correlation coefficients
    Rcorr <- (cor(aa, V))^2
    # To store the value of R2 (Z)
    # R2=R2+sparse([1:n],i,Rcorr,n,n);
    spR <- sparseMatrix(i = c(1:n),j = rep(cI,n), x = as.numeric(Rcorr), dims = c(n,n))
    R2 <- R2 + spR
  }

  max(R2, na.rm = TRUE)
  ZZ <- as.matrix(R2)
  str(ZZ)
  max(ZZ, na.rm = TRUE)


  #--------------------------------------------------------------------------
  # identify the max R2 and its corresponding bands
  #--------------------------------------------------------------------------

  R2max <- max(ZZ, na.rm = TRUE)
  print(paste('The maximun value of R^2 is', as.character(round(R2max,4))))

  ind_max <- which(ZZ == R2max, arr.ind = TRUE)
  ind_max

  bestBands = wavelength[ind_max[1,]]


  #----------------------------------
  # plot ndvi
  #----------------------------------

  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  zp1 <- ggplot(melt(ZZ), aes(Var1+350-1,Var2+350-1, fill=value))
  zp1 <- zp1 + geom_tile()
  zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(100))
  # zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
  # zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
  zp1 <- zp1 + coord_equal()
  zp1 <- zp1 + theme_bw()
  print(zp1)

}


