# Simulate and thin high-resolution Langevin tracks (used by simulation_study/).

# simulate langevin track (High-res)
sumLTrack <- function(dt, gamma2, covlist, beta, loc0, nobs) {
  x = mvnfast::rmvn(nobs, rep(0,2), dt * gamma2 * diag(1,2,2))
  x[1, ] <- loc0 
  for (i in 2:nrow(x)) {
    grad = bilinearGradVec(matrix(x[i-1, 1:2], nrow=1), covlist)
    x[i, ]  = x[i, ] + x[i-1, ] + (dt * gamma2 / 2)*beta %*% grad[,1,]
  }
  return(x)
}

# thin high-resolution track
thinTrack <- function(x, thin) {
  n = nrow(x)
  x = x[seq(1, n, thin), ]
  return(x)
}

# simulate thinned Langevin movement model
simLMM <- function(dt, gamma2, covlist, beta, loc0, nobs){ 
  x <- sumLTrack(dt, gamma2, covlist, beta, loc0, nobs)  # Simulate high-resolution track
  #x <- thinTrack(x, delta) # Thin the track
  return(x)
}