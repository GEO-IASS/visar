# devtools::use_data_raw()
# * Add data creation scripts in data-raw
# * Includes data-raw/ on .Rbuildignore

exampleData <- read.csv('data-raw/exampleData.csv', header=FALSE)[1:20,1:952]
# wavelength <- as.numeric(gsub("X","",names(exampleData)[(1+1):952]))

exampleData <- data.matrix(exampleData) # Convert a Data Frame to a Numeric Matrix
exampleData <- matrix(as.numeric(unlist(exampleData)),nrow=nrow(exampleData))

visa.spectra <- exampleData

devtools::use_data(exampleData, overwrite = TRUE)
devtools::use_data(visa.spectra, overwrite = TRUE)

# Internal data
# Sometimes functions need pre-computed data tables
# devtools::use_data(exampleData, internal = TRUE)


load("data/exampleData.rda")
aa <- matrix(exampleData[, 1]) # Variable of interest, e.g., Chl, N, LAI
ss <- exampleData[, 2:ncol(exampleData)] # Reflectance spectra
