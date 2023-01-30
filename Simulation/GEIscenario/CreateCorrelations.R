# Simulate a correlation matrix between environments
#============================================================
nEnvPerYear <- nEnvs*nCycles
rho <- 0.5 # baseline correlation in k-th group
epsilon <- 0.99 - max(rho)
Cmat <- simCmat(groups = 1, 
                size = nEnvPerYear,
                rho = rho,
                delta = 0, 
                epsilon = epsilon,
                eidim = 6, 
                skew = 0.6)
min <- summary(Cmat[upper.tri(Cmat)])["Min."]
max <- summary(Cmat[upper.tri(Cmat)])["Max."]
while ((min > 0 & min < 0.1 & max > 0.9) == FALSE) {
  Cmat <- simCmat(groups = nCycles, 
                  size = nEnvPerYear,
                  rho = rho,
                  delta = 0, 
                  epsilon = epsilon,
                  eidim = 6, 
                  skew = 0.6)
  min <- summary(Cmat[upper.tri(Cmat)])["Min."]
  max <- summary(Cmat[upper.tri(Cmat)])["Max."]
}
print(summary(Cmat[upper.tri(Cmat)]))
hist(Cmat[upper.tri(Cmat)])
# plotCMAT(cmat = Cmat, order = T)
