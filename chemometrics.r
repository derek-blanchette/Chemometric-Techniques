# Page numbers refer to
# "Chemometric Techniques for Quantitative Analysis" by Richard Kramer
# http://www.amazon.com/Chemometric-Techniques-Quantitative-Analysis-Richard/dp/0824701984/

#Overview of process: 
# -Read in the spectra data.
# -Read in the concentration data.
# -Estimate/fit the calibration.
# -Validate the calibration against the unknowns.

#This library is not needed for the estimation.
#Input was provided in Excel format.
library(xlsx)

#Read in the standard series (known concentrations).

#Read in each spectra and cbind them together
spectra <- NULL
for(i in seq(1 , 26)){
  sample <- read.xlsx("STD Series.xlsx", i)
  keep <- sample[,2]
  spectra <- cbind(spectra, keep)
}

spec2 <- spectra[1:191,]

concs <- read.csv("Concentrations.csv",
                  sep = ",", header = TRUE)

concs <- as.matrix(t(concs))

FirstPart <- spec2 %*% t(concs)
SecondPart <- chol2inv( chol((concs %*% t(concs))) )

###############################################
# Page 52 - This is our calibration
##############################################

EstimatedPureSpectra <- FirstPart %*% SecondPart

K <- EstimatedPureSpectra

#Benzocaine
plot(seq(190,380), EstimatedPureSpectra[,1] , type = "l", 
     main = "Estimated Pure Benzocaine",
     ylab = "Intensity",
     xlab = "Wavelength")

# PVP
plot(seq(190, 380), EstimatedPureSpectra[,2], type = "l",
     main = "Estimated Pure PVP",
     ylab = "Intensity",
     xlab = "Wavelength")

#Concentration Space
plot(c(0,t(concs)[,1]), c(0,t(concs)[,2]),
     main = "Concentration Combinations for Standards",
     xlab = "Benzocaine",
     ylab = "PVP",
     pch = 15,
     col = "blue")


#######################################################
#   Test the Unknowns
#######################################################

#Read in the unknown data

#Read in each spectra and cbind them together
#I think I can just do the concentration determination all at once

unknowns <- NULL

for(i in seq(1 , 8)){
  sample <- read.xlsx("SP Series.xlsx", i)
  keep <- sample[,2]
  unknowns <- cbind(unknowns, keep)
}

#190 rows analysis
unknowns2 <- unknowns[1:190,]
samp1 <- unknowns2[,1]
K <- K[-191,]

#Cholesky factorization is used to find inverse 
#to increase numerical stability!

FirstPart <- chol2inv( chol(t(K) %*% K) )
#Here is where we plug in the spectra, opposite side of K.
SecondPart <- t(K) %*% unknowns2
Concentration <- FirstPart %*% SecondPart

##########################################
#170 rows analysis
unknowns2 <- unknowns[1:170,]
samp1 <- unknowns2[,1]
K <- K[1:170,]

FirstPart <- chol2inv( chol(t(K) %*% K) )
#Here is where we plug in the spectra, opposite side of K.
SecondPart <- t(K) %*% unknowns2
Concentration <- FirstPart %*% SecondPart

UnknownConcs <- t(Concentration)
colnames(UnknownConcs) <- c("Benzocaine", "PVP")

ValidationConcs <- read.csv("Unknown Concentrations.csv", sep = ",", header = TRUE)

#Differences
#Absolute Error / Absolute Value = Relative Error
Differences <- ValidationConcs - UnknownConcs
Scaled <- abs(Differences / ValidationConcs)

paste("Benzocaine")
mean(Scaled[,1])
sd(Scaled[,1])

paste("PVP")
mean(Scaled[,2])
sd(Scaled[,2])
