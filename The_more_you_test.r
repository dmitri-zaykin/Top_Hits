# "The more you test, the more you find:
# The smallest P-values become increasingly
# enriched with real findings as more tests are conducted"
#
# (Vsevolozhskaya et al; submitted)
#
#
# Scripts in this directory are to be continuously updated
# Stay tuned
#
# This script computes values such as given by Table 3
# The "simulation" method can be very slow!
# If you trust the formulas, set the simulation number (oo)
# to some small value to get a rough check
############################
# Command line example: 
# R CMD BATCH --vanilla --slave --no-timing '--args j1=1 j2=25 j3=50 alpha=0.25 gamma1=0.05 oo=500' Table_3.r /dev/null &
#
#
# default values
# note: these can be overwritten by values given in the command line arguments (example at the top of this file)
alpha <- 0.1 # alpha=0 means "no truncation"
gamma1 <- 0 # point null?
# in Table 3, these are the three values of "i"
j1 <- 1
j2 <- 25
j3 <- 50
# in Table 3, these are the three values of "k"
Ki = c(500, 1000, 2000)

# number of simulations
oo <- 1500


args=(commandArgs(TRUE))
if(length(args) > 0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

jj = c(j1, j2, j3)

# Output file names
# append alpha, gammas[1] to the file name
# and do everything for 1 value of alpha
File1 <- paste("Table_3.Full", "_", jj[1], ".", jj[2], ".", jj[3], "_a", alpha, ".txt", sep="")
File2 <- paste("Table_3", "_", jj[1], ".", jj[2], ".", jj[3], "_a", alpha, ".txt", sep="")


Df <- 1
Hgam0 <- function(x, w, gamX) {
    sum(w*pchisq(x, df=Df, ncp=gamX))
}
Hgam <- Vectorize(Hgam0, vectorize.args=c("x"))

Hgam.ord0 <- function(x, j, n, w, gamX) {
    1 - pbinom(j-1, n, Hgam(x, w, gamX))
}
Hgam.ord <- Vectorize(Hgam.ord0, vectorize.args=c("x"))
Hgam.inv0 <- function(p, j, n, ww, gm, br=c(-10000, 10000)) {
    G = function(x) Hgam.ord(x, j, n, ww, gm) - p
    return( uniroot(G,br)$root )
}
Hgam.inv <- Vectorize(Hgam.inv0, vectorize.args=c("p"))
CHgam0 <- function(x, w, gam) {
    1 - Hgam(x,w,gam)
}
CHgam <- Vectorize(CHgam0, vectorize.args=c("x"))
CHgam.ord0 <- function(x, j, n, w, gam) {
    1 - pbinom(j-1, n, CHgam(x, w, gam))
}
CHgam.ord <- Vectorize(CHgam.ord0, vectorize.args=c("x"))
#
fy <- function(x2, nc, df=1) { dchisq(x2, df=df, ncp=nc) }
Post.tab0 <- function(z, nc, b) {
  nn = length(b)
  dn = fy(z, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy(z, nc[i])*b[i] / dn }
  pst
}
Post.tab <- Vectorize(Post.tab0, vectorize.args=c("z"))

Nc <- 10
Shape = 1
Scale = Nc/Shape
Trunc <- 100
Nb <- 600
step <- Trunc/Nb

Compute.mixture = TRUE # if FALSE, read from the file "Weights.RData"
#
if(Compute.mixture) {
    require(actuar)
    w = discretize(pgamma(x, shape=Shape, scale=Scale), method = "upper", from=step, to=Trunc+step, step=step)
    gammas <- seq(from=step, to=Trunc+step,  length.out=Nb)
    w <- c(NA, w)
    Chop <- Trunc
    w <- w[1:Chop]
    w[1] = 0.9 # force high probability for the first bin
    v = w[-1]
    v = (1-w[1]) * v/sum(v)
    w[-1] = v
    w <- w/sum(w)
    (w <- w/sum(w))
    (gammas <- gammas[1:Chop])
    gammas <- 10*gammas/max(gammas)
    #save(file="Weights.RData", w, gammas)
} else {
    load(file="Weights.RData")
}

gammas[1] = gamma1 # Point Null?

postE0 <- postE1 <- postE2 <- rep(0, length(Ki))
expE0 <- expE1 <- expE2 <- rep(0, length(Ki))
i=1
for(k in Ki) { 
    Fb <- qchisq(1-alpha/k, df=1)
    pFb <- Hgam(Fb, w, gammas)
    Hgam0.Trun <- function(x, w, gamX) { sum(w*pchisq(x, df=Df, ncp=gamX)) / pFb }
    Hgam.Trun <- Vectorize(Hgam0.Trun, vectorize.args=c("x"))
    Hgam.Trun.ord0 <- function(x, j, n, w, gamX) { 1 - pbinom(j-1, n, Hgam.Trun(x, w, gamX)) }
    Hgam.Trun.ord <- Vectorize(Hgam.Trun.ord0, vectorize.args=c("x"))
    fu0 <- function(x, j, n, w, gamX) { 1-Hgam.Trun.ord(x, j, n, w, gamX) }
    fu.Trunc <- Vectorize(fu0, vectorize.args=c("x"))
    x1 = integrate(fu.Trunc, 0, Fb, k-jj[1]+1, k, w, gammas, rel.tol = .Machine$double.eps^0.5, stop.on.error = FALSE)
    expE0[i] = x1$value
    postE0[i] <- sum(gammas*Post.tab(x1$value, gammas, w))
    x1 = integrate(fu.Trunc, 0, Fb, k-jj[2]+1, k, w, gammas, rel.tol = .Machine$double.eps^0.5, stop.on.error = FALSE)
    expE1[i] = x1$value
    postE1[i] <- sum(gammas*Post.tab(x1$value, gammas, w))
    x1 = integrate(fu.Trunc, 0, Fb, k-jj[3]+1, k, w, gammas, rel.tol = .Machine$double.eps^0.5, stop.on.error = FALSE)
    expE2[i] = x1$value
    postE2[i] <- sum(gammas*Post.tab(x1$value, gammas, w))
    cat(k, Fb, postE0[i], postE1[i], postE2[i], expE0[i], expE1[i], expE2[i], "\n")
    i <- i+1
}



if(Df == 1) Sgammas = sqrt(gammas)

Num.Bins <- length(w)

cat(file=File1, "k & Av(post) & Av(gamma) & Av(X2) & E(gamma) & E(X2)", "\\\\", "\n", append=TRUE)
cat(file=File1, "jj =", jj[1], jj[2], jj[3], "\n", append=TRUE)

cat(file=File2, "jj =", jj[1], jj[2], jj[3], "\n", append=TRUE)
cat(file=File2, "k & Av(gamma) & Av(X2) & E(gamma) & E(X2)", "\\\\", "\n", append=TRUE)


ii=1
for(k in Ki) {
    Fb = qchisq(1-alpha/k, df=1)
    postE.emp0 <- postE.emp1 <- postE.emp2 <- rep(0, oo)
    expE.emp0 <- expE.emp1 <- expE.emp2 <- rep(0, oo)
    expMx0 <- expMx1 <- expMx2 <- rep(0, oo)
    for(i in 1:oo) {
        repeat {
            idx <- sample(1 : Num.Bins, prob=w, size=k, replace=TRUE)
            if(Df==1) {
                x2 <- (rnorm(k, mean = Sgammas[idx]))^2 # for df=1 much faster than rchisq
            } else {
                x2 <- rchisq(k, df=Df, ncp=gammas[idx])
            }
            if(max(x2) < Fb) break;
        }
        mxi <- order(x2)[k - jj[1] + 1] # index of jj-th largest X2
        expMx0[i] <- (gammas[idx])[mxi]
        expE.emp0[i] <- x2[mxi]
        postE.emp0[i] <- sum(gammas*Post.tab(x2[mxi], gammas, w))
        mxi <- order(x2)[k - jj[2] + 1] # index of jj-th largest X2
        expMx1[i] <- (gammas[idx])[mxi]
        expE.emp1[i] <- x2[mxi]
        postE.emp1[i] <- sum(gammas*Post.tab(x2[mxi], gammas, w))
        mxi <- order(x2)[k - jj[3] + 1] # index of jj-th largest X2
        expMx2[i] <- (gammas[idx])[mxi]
        expE.emp2[i] <- x2[mxi]
        postE.emp2[i] <- sum(gammas*Post.tab(x2[mxi], gammas, w))
        if(i %% 100 == 0) cat(".", i, ".", sep=""); if(i %% 1000 == 0) cat("\r")
    }
    cat("\n")
    cat(file=File1, k, "&",
        round(mean(postE.emp0),2), "&",
        round(mean(postE.emp1),2), "&",
        round(mean(postE.emp2),2), "&",
        round(mean(expMx0),2), "&",
        round(mean(expMx1),2), "&",
        round(mean(expMx2),2), "&",
        round(mean(expE.emp0),2), "&",
        round(mean(expE.emp1),2), "&",
        round(mean(expE.emp2),2), "&",
        round(postE0[ii],2), "&",
        round(postE1[ii],2), "&",
        round(postE2[ii],2), "&",
        round(expE0[ii],2), "&",
        round(expE1[ii],2), "&",
        round(expE2[ii],2), "&",
        "\\\\", "\n", append=TRUE)
    cat(file=File2, k, "&",
        round(mean(expMx0),2), "&",
        round(mean(expMx1),2), "&",
        round(mean(expMx2),2), "&",
        round(postE0[ii],2), "&",
        round(postE1[ii],2), "&",
        round(postE2[ii],2), "&",
        round(expE0[ii],2), "&",
        round(expE1[ii],2), "&",
        round(expE2[ii],2), "&",
        "\\\\", "\n", append=TRUE)
    ii <- ii+1
}
