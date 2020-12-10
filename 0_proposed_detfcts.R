# detection functions that will be used in the simulation
# parameters here are saved for use later

load("issj.RData")
source("support_functions.R")

# plot candidate detection functions
plot_pdf <- function(betas, title="", truncation){
  xx <- seq(0, truncation, length.out=250)
  # plot levels from paper
  chaplevs <- c(24, 47.5, 71)
  # fiddle to rescale
  chaplevs_scaled <- (chaplevs/100-s_chaparral_detail[1])/s_chaparral_detail[2]

  nu <- sig <- rep(0,3)
  for(i in seq_along(chaplevs)){
    sig[i] <- exp(betas[1] + betas[2]*chaplevs_scaled[i])
    nu[i] <- 2 * pi * sig[i]^2 *
              (1 - exp(-(truncation^2)/(2*sig[i]^2)))
  }

  plot(xx, xx*detfn(xx, betas, chaplevs_scaled[2])/nu[2], type="l",
       main=paste0(paste(round(betas, 2), collapse=", "), " ", title),
       xlab="Distance (m)", ylab="Probability of detection")
  lines(xx, xx*detfn(xx, betas, chaplevs_scaled[1])/nu[1], lty=2)
  lines(xx, xx*detfn(xx, betas, chaplevs_scaled[3])/nu[3], lty=3)
}


# parameter combinations with very silly names
good_betas <- c(4.68, -0.2)
bad_betas <- c(3.5, -0.2)
middle_betas <- c(4, -0.2)
middle2_betas <- c(4.3, -0.2)
middle3_betas <- c(3.75, -0.2)
middle4_betas <- c(3.875, -0.2)
vgood_betas <- c(5, -0.2)

all_betas <- list(good_betas, bad_betas, middle_betas, middle2_betas,
                  middle3_betas, middle4_betas, vgood_betas)
names(all_betas) <- c("good", "bad", "middle", "middle2", "middle3", "middle4", "vgood")

# save the parameters
save(all_betas, file="df_pars.RData")

# rough plots, see supplementary material for better ones
#par(mfrow=c(2,3))
#lapply(all_betas, plot_pdf, title="", truncation=300)


