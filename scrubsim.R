# main simulation code
# this runs all simulations, in parallel where possible

library(dsm)
library(mgcv)
library(Distance)
library(parallel)


# load detection function parameters
load("df_pars.RData")

simit <- function(these_betas, truncation, betaset, id){

  load("issj.RData")
  # generate a new distribution of birds, calculate distances
  # between birds and samplers and calculate true population size
  dists <- make_dists(cruz, segs)
  Nbirds <- dists$Nbirds
  dists <- dists$dists

  # which of those birds were observed?
  obs <- make_sample(dists, detfn, these_betas, segs)

  # fit detection function
  df_fit <- try(suppressMessages(ds(obs, transect="point", formula=~chaparral,
                                    truncation=truncation)))

  # handle the case where the detection function didn't fit
  if("try-error" %in% class(df_fit)){
    return(data.frame(model  = NA,
                      Nbirds = NA,
                      Nhat   = NA,
                      var    = NA,
                      n_obs  = NA))
  }

  # trim the observations according to the truncation
  obs <- obs[obs$distance <= df_fit$ddf$meta.data$width, ]

  # fit spatial model
  dsm_b <- dsm(count~s(x,y, k=20, bs="tp"),
               observation.data=obs, segment.data=segs, ddf.obj=df_fit,
               transect="point", family=poisson())

  # make a prediction/variance estimate (varprop)
  vp <- dsm_varprop(dsm_b, cruz)
  # same but assuming independence
  vg <- dsm.var.gam(dsm_b, cruz, cruz$off.set)

  # do some gam posterior sampling for the varprop model
  br <- gam.imp(ns=1000, vp$refit)
  cruz$XX <- matrix(0, ncol=2, nrow=nrow(cruz))
  Xp <- predict(vp$refit, cruz, type="lpmatrix")
  vp_Nhat_post <- colSums(cruz$off.set * exp(Xp%*%t(br$bs)))

  vp_Nhat_post <- cbind(vp_Nhat_post, br$wts)


  nimfit <- try(nimble_dsm(dsm_b, obs, segs, cruz, vp))

  if("try-error" %in% class(nimfit)){
    return(data.frame(model  = NA,
                      Nbirds = NA,
                      Nhat   = NA,
                      var    = NA,
                      n_obs  = NA,
                      id     = id))
  }

  # write out the last 1000 samples for fully Bayesian...
  write.csv(file=paste0("posteriors/fb/out-", id, ".csv"),
            nimfit,
            #[(nrow(nimfit)-1000):nrow(nimfit), ],
            row.names=FALSE)
  # and varprop
  write.csv(file=paste0("posteriors/vp/out-", id, ".csv"), vp_Nhat_post,
            row.names=FALSE)

  # median and variance of the Nhats
  nimfit <- list(Nhat     = mean(nimfit[,1]),
                 Nhat_var = var(nimfit[,1]))


  # return summary statistics
  data.frame(model   = c("varprop", "delta", "onestage"),
             Nbirds  = Nbirds,
             Nhat    = c(sum(vp$pred), vg$pred[[1]], nimfit$Nhat),
             var     = c(vp$var[1,1], vg$pred.var, nimfit$Nhat_var),
             n_obs   = nrow(obs),
             betaset = betaset,
             id      = id)
}




# function to run simulations in parallel for a given scenario

# need to make a PSOCK parallel process for the compilation to be safe
# but this means the (parent) environment of par_func is empty when we start
# sourcing 1_nimble_setup.R gives us the objects we need
# https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html
par_func <- function(X){
  library(dsm)
  library(mgcv)
  library(Distance)
  library(nimble)

  # get fitting function
  source("nimble_dsm.R")
  source("support_functions.R")

  res <- simit(these_betas, truncation="10%", betaset=names(all_betas)[b],
               id=paste(b, X, sep="-"))

  write.table(res, file=paste0("results_csvs/results-", names(all_betas)[b], ".csv"),
             append=TRUE, col.names=FALSE, sep=", ")

  gc()
  return(res)
}

# how many simulations per scenario

nsim <- 200

# run the simulation
for(b in seq_along(all_betas)){
  these_betas <- all_betas[[b]]

  # setup a cluster
  this_cluster <- makeCluster(5, outfile="cluster_out.txt")

  # export some objects to the cluster environment (others are loaded inside)
  clusterExport(cl=this_cluster, c("simit", "these_betas", "all_betas", "b"))

  # run the simulation
  results <- parLapply(cl = this_cluster, X = 1:nsim,
                       fun = par_func)

  # stop the cluster
  stopCluster(this_cluster)
}


