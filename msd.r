################################################################################
#  Modified Stahel-Donoho Estimators for Multivariate Outlier Detection  
#    Ver.1.6 2009/07/14
#    Ver.1.7 2018/10/19   Modify gso function to stop warning messages
################################################################################
#   by WADA, Kazumi (National Statistics Center of Japan)
#   Published Ver. 1.6 at http://www.stat.go.jp/training/2kenkyu/pdf/ihou/67/wada1.pdf 
################################################################################

msd <- function(inp, nb=0, sd=0, tm="EUR") {

inp_d <- ncol(inp)            # number of variables
inp_n <- nrow(inp)            # number of observations

##############################
# create orthogonal bases
##############################

if (sd != 0) set.seed(sd)

## number of orthogonal bases
if (nb == 0) bb_n <- trunc(exp(2.1328+0.8023*inp_d) / inp_d)
   else bb_n <- nb
rn <- bb_n * inp_d^2           # number of random variables for the bases
basis <- array(runif(rn), c(inp_d, inp_d, bb_n))

## orthonormalization by function "gso"
basis <- apply(basis, 3, gso)
basis <- array(basis, c(inp_d, inp_d, bb_n))

##############################
# projection and residual computation
##############################

prj <- array(0, c(inp_n, inp_d, bb_n))  # for projection
res <- array(0, c(inp_n, inp_d, bb_n))  # resudal
wt <- array(0, c(inp_n, inp_d, bb_n))   # weight by observations x variable x bases
wts <- array(0, c(inp_n, bb_n))         # weight by observations x bases
bwt <- rep(0, inp_n)                    # the smallest weight by observations
kijun <- qchisq(0.95, inp_d)            # reference for trimming

Fprj <- function(pj) t(pj %*% t(inp))   # projection
prj <- apply(basis, 3, Fprj)
prj <- array(prj, c(inp_n, inp_d, bb_n)) 

medi <- apply(prj, c(2, 3), median)     # median
madx <- apply(prj, c(2, 3), mad)        # median absolute deviation (MAD) / 0.674 

for (i in 1:bb_n) {                     # robust standardization of residuals
    res[,,i] <- t(abs(t(prj[,,i]) - medi[,i]) / madx[,i])
}

### trimming weight 
if (tm == "CAN") {
   k0 <- which(res <= 1.75)
   k1 <- which(res > 1.75 & res <= 3.5)
   k2 <- which(res > 3.5)
   wt[k0] <- 1
   wt[k1] <- 1.75 / res[k1]
   wt[k2] <- 0
}
else {                                  # trimming by a Huber-like weight function
   k0 <- which(res <= sqrt(kijun))
   k1 <- which(res > sqrt(kijun))
   wt[k0] <- 1
   wt[k1] <- kijun / (res[k1]^2)
}
wts <- apply(wt, c(1,3), prod)          
bwt <- apply(wts, 1, min)               # selecting the smallest weight 

### initial robust covariance matrix
u1 <- apply(inp * bwt, 2, sum) / sum(bwt)
V1 <- t(t(t(inp) - u1) * bwt) %*% (t(t(inp) - u1) * bwt) / sum(bwt^2)

### avoiding NaN error
u1 <- ifelse(is.nan(u1), 0, u1)
V1 <- ifelse(is.nan(V1), 0, V1)

### robust PCA (LAPACK)
eg <- eigen(V1, symmetric=TRUE)         
ctb <- eg$value / sum(eg$value)         # contribution ratio

##############################
# projection pursuit (PP)
##############################
res2 <- array(0, c(inp_n, inp_d))       # residuals
wt2  <- array(0, c(inp_n, inp_d))       # weight by observations x variables
wts2 <- array(0, inp_n)                 # final weight by observations

prj2 <- t(eg$vector %*% (t(inp) - u1))  # projection
medi2 <- apply(prj2, 2, median)         # median and 
madx2 <- apply(prj2, 2, mad)            # MAD for standardization
res2 <- t(abs(t(prj2) - medi2) / madx2) # standardized residuals

# trimming
if (tm == "CAN") {
   k0 <- which(res2 <= 1.75)
   k1 <- which(res2 > 1.75 & res2 <= 3.5)
   k2 <- which(res2 > 3.5)
   wt2[k0] <- 1
   wt2[k1] <- 1.75 / res2[k1]
   wt2[k2] <- 0
}
else {                                  # trimming by the Huber-like weight function
   k0 <- which(res2 <= sqrt(kijun))
   k1 <- which(res2 > sqrt(kijun))
   wt2[k0] <- 1
   wt2[k1] <- kijun / (res2[k1]^2)
}

wts2 <- apply(wt2, 1, prod)  

if (tm == "EUR") wts2 <- pmin(wts2, bwt)

##############################
# final mean vector and covariance matrix
##############################

u2 <- apply(inp * wts2, 2, sum) / sum(wts2)
V2 <- t(t(t(inp) - u2) * wts2) %*% (t(t(inp) - u2) * wts2) / sum(wts2^2)

return(list(u1=u1, V1=V1, bwt=bwt, u2=u2, V2=V2, wts2=wts2, eg=eg, ctb=ctb))
}

###########################################################
# gso: Gram-Schmidt Orthonormalization for function "msd"
###########################################################

gso <- function(basis) {
    bd <- ncol(basis) 	
    bn <- nrow(basis)	
    # basis[1,] <- basis[1,] / sqrt(t(basis[1,]) %*% basis[1,])
    basis[1,] <- as.vector(basis[1,]) / sqrt(as.vector(t(basis[1,]) %*% basis[1,]))
    for (i in 2 : bd ) {
    #    wk1 <- basis[i,]
        wk1 <- as.vector(basis[i,])
        for (j in 1:(i-1)) {
#            wk2 <- basis[j,]
            wk2 <- as.vector(basis[j,])
#            basis[i,] <- basis[i,] - (t(wk1) %*% wk2) * wk2
            basis[i,] <- as.vector(basis[i,]) - (as.vector(t(wk1) %*% wk2) * wk2)
        }
#        basis[i,] <- basis[i,] / sqrt(t(basis[i,]) %*% basis[i,])
        basis[i,] <- as.vector(basis[i,]) / sqrt(as.vector(t(basis[i,]) %*% basis[i,]))
    }
    return(basis)
} 

################################################################################

