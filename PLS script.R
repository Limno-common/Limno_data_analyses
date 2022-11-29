####Karla's PLS script
##I use my own data in this script. I have already standardized and centered
##everything in Excel, but you can also do this in R, if your dataset is big.

#This analysis requires the package "pls"
library(pls)

###loading my data from the work directory; set your own directory with setwd()
data <- read.csv("pls_data_std.csv")

###Define the y-variable for the PLS-analysis. In this case, I chose c13. 
y <- data$c13
y[is.na(y)] <- 0		#missing values have to be set to 0 for the PLS

###Define the x-variables for the PLS analysis. In this case, it's all variables
###except for c13 and om.autoch
###Take note: the x-variables have to be converted into a matrix and dataframe
###so that the structure works in the PLS! 
x <- data		#define all x-variables
x$om.autoch <- NULL	#remove variable
x$c13 <- NULL		#remove y-variable
x <- as.matrix(x)		#transform the x-variables to a matrix
frame <- data.frame(y, x=I(x))	#transform the matrix into a dataframe
str(frame)			#check the structure

###The dataframe for x-variables should look like this now: 
#variables:
 #$ y: num  0.31 0.619 -0.112 -1.075 -0.625 ...
 #$ x: 'AsIs' num [1:24, 1:14] -0.412 -1.042 -1.138 -0.898 -1.117 ...
 # ..- attr(*, "dimnames")=List of 2
 # .. ..$ : NULL
 # .. ..$ : chr [1:14] "gonyo.density" "chl" "secchi" "sechi.epi.ratio" ...


###Do the PLS using the plsr() command
#y~x <- your y-variable and x-variable-matrix
#data <- the dataframe that contains the x-variable-matrix
#if you want to use only part of the x-matrix, you can specify the number of components with ncomp=n. The ncomp-default is the number of columns in your x-variable-matrix
#validation <- internal validation for the PLS-model. You can choose "CV" (uses cross-validation segments) or "LOO" (leave-one-out, cross-validation leaving out one segment).
#method <- method used to fit the PLS-model. Here I used orthogonal scores ("oscorespls"). It's a classical PLS-method, but check out the package-info for more on this.

gonyo.pls <- plsr(y~x, data=frame, validation="LOO", method = "oscorespls")
summary(gonyo.pls)	#summary of PLS (includes validation and training of PLS-model)

###Let's make a list of the x-variables and their VIP-scores
library(plsVarSel)	#load this library
wtf <- which.min(gonyo.pls$validation$PRESS)	#compiles the PLS-results for the different components
vip <- VIP(gonyo.pls, wtf)	#calculates VIP-scores from the PLS-results
sort(vip, decreasing=TRUE)	#sorts the VIP scores from highest to lowest. 

###VIP scores: every variable with a VIP>1 is considered to significantly contribute to the PLS-model
 