make_dists <- function(cruz, segs){

  # calculate the density surface
  # use true relationships from paper
  f_elev <- function(x) x
  f_chap <- function(x) x + x^2

  # truth is the sum of these functions per grid cell
  true_surface <- function(elevation, chaparral){
    f_elev(elevation) + f_chap(chaparral)
  }

  # rescale covariates
  cruz$elevation <- scale(cruz$elevation)[,1]
  cruz$elevation <- cruz$elevation - min(cruz$elevation)
  cruz$chaparral <- scale(cruz$chaparral)[,1]
  cruz$chaparral <- cruz$chaparral - min(cruz$chaparral)

  # calculate proportion of population at each grid cell
  cruz$p <- true_surface(cruz$elevation, cruz$chaparral)
  cruz$p <- cruz$p/sum(cruz$p)

  # spread out true N over the density
  trueN <- 2500
  cruz$N <- trueN * cruz$p


  # now generate the locations of birds
  point_pattern <- c()

  # number of birds per cell is poisson process
  cruz$birds <- rpois(rep(1,nrow(cruz)), lambda=cruz$N)
  Nbirds <- sum(cruz$birds)

  truncation <- 300
  # now create the bird distribution jittering each observation from
  # the point centre
  # we assume uniform distribution within a grid cell (reasonable for
  # distance sampling assumptions!)
  for(i in 1:nrow(cruz)){
    if(cruz$birds[i] > 0){
      for(j in 1:cruz$birds[i]){
        point_pattern <- rbind(point_pattern,
                               cbind(cruz$x[i]+rnorm(1, 0, truncation),
                                     cruz$y[i]+rnorm(1, 0, truncation)))
      }
    }
  }
  point_pattern <- as.data.frame(point_pattern)
  names(point_pattern) <- c("x" , "y")

  # need to calculate distances between points and birds
  # https://stat.ethz.ch/pipermail/r-devel/2017-June/074488.html
  # fast Euclidean distances
  euc_dist <- function(m) {
    mtm <- Matrix::tcrossprod(m)
    sq <- rowSums(m*m)
    sqrt(outer(sq,sq,"+") - 2*mtm)
  }
  # concatenate points and birds locations
  pts <- as.matrix(rbind(segs[,c("x","y")], point_pattern[,c("x","y")]))
  # calculate all idstances
  dists <- euc_dist(pts)
  # first segs rows and pts columns gives lookup we need
  dists <- dists[1:nrow(segs), (nrow(segs)+1):nrow(dists)]


  return(list(dists=dists,
              Nbirds=Nbirds))
}


# get the observations by sampling at the segments and
# applying the detection function
make_sample <- function(dists, detfn, betas, segs){
  samp <- data.frame(distance=NA,
                     Sample.Label=NA,
                     chaparral=NA)
  # loop over segments
  for(i in 1:nrow(dists)){
    # get distances for this segment
    these_dists <- dists[i, ]
    # only include those within truncation
    candidate_distances <- which(these_dists <= 300)
    # what's the probability of that bird being observed
    p <- detfn(these_dists[candidate_distances], betas,
               segs$chaparral[i])
    # randomly sample it
    ind <- which(runif(length(p)) < p)
    if(length(ind) > 0){
      # create data
      samp <- rbind(samp,
                    data.frame(distance=these_dists[candidate_distances[ind]],
                               Sample.Label=i,
                               chaparral=segs$chaparral[i]))
    }
  }

  samp$size <- 1
  samp$object <- 1:nrow(samp)
  return(samp[-1,])
}

# define half-normal detection function with one covariate
detfn <- function(x, betas, cov){
  sigma <- exp(betas[1] + betas[2]*cov)
  exp(-x^2/(2*sigma^2))
}

# what quantile of the estimate does truth lie in?
get_N_quantile <- function(N, Nhat, cvN){

  meantrans <- function(m, cv) log(m)-0.5*log(cv^2+1)
  setrans <- function(cv) sqrt(log(cv^2+1))

  plnorm(N, meantrans(Nhat, cvN), setrans(cvN))
}

# as above but using the empirical CDF, optionally weighted
get_e_quantile <- function(N, Ns, wts=NULL){
  if(is.null(wts)){
    ecdf(Ns)(N)
  }else{
    ecdfvals <- Hmisc::wtd.Ecdf(Ns, wts)
    rval <- approxfun(ecdfvals$x, ecdfvals$ecdf,
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    rval(N)
  }
}

# importance sampling based on code from mgcv
gam.imp <- function(b, ns=10000){

  beta <- coef(b)
  Vb <- vcov(b)
  X <- model.matrix(b)

  # beta proposals
  bs <- rmvn(ns, beta, Vb)

  # log proposal density
  lfp <- dmvn(t(bs), beta, Vb)
  # log density for penalized MLE
  lfp1 <- dmvn(matrix(beta, ncol=1), beta, Vb)

  # storage
  wts <- rep(0, ns)

  # loglik for penalized MLE
  lpl0 <- lpl(b, beta, X)

  for (i in 1:ns) {
    # loglik of proposal...
    lpl1 <- lpl(b, bs[i,], X)
    # store weight
    wts[i] <- exp(lfp1-lfp[i]+lpl1-lpl0)
  }
  list(bs=bs, wts=wts)
}


# taken from mgcv source code
## Simple post fit mcmc for mgcv.
## (c) Simon N. Wood (2020)
## some functions to extract important components of joint density from
## fitted gam...
## evaluate penalty for fitted gam, possibly with new beta
# patched to include parapen
bSb <- function(b,beta=coef(b)) {
  bSb <- k <-  0
  sp <- if (is.null(b$full.sp)) b$sp else b$full.sp ## handling linked sp's


  # the parapen bits
  # need to do something clever with L at some point
  if(!is.null(b$paraPen)){
    for (i in 1:length(b$paraPen$S)) {
      k <- k + 1
      # get indices
      ii <- b$paraPen$off[i]
      ii <- ii:(ii+ncol(b$paraPen$S[[i]])-1)
      # add to penalty
      bSb <- bSb + sp[b$paraPen$full.sp.names[i]]*
                    (t(beta[ii])%*%b$paraPen$S[[i]]%*%beta[ii])
    }
  }


  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + sp[k]*(t(beta[ii])%*%b$smooth[[i]]$S[[j]]%*%beta[ii])
      }
    }
  }


  bSb
} ## bSb
devg <- function(b,beta=coef(b),X=model.matrix(b)) {
## evaluate the deviance of a fitted gam given possibly new coefs, beta
## for general families this is simply -2*log.lik
  if (inherits(b$family,"general.family")) {
    -2*b$family$ll(b$y,X,beta,b$prior.weights,b$family,offset=b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y,b$family$linkinv(X%*%beta+b$offset),b$prior.weights))
  }
} ## devg
lpl <- function(b,beta=coef(b),X=model.matrix(b)) {
## log joint density for beta, to within uninteresting constants
  -(devg(b,beta,X)/b$sig2+bSb(b,beta)/b$sig2)/2
}

