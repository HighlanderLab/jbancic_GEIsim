# Simulate genetic correlation matrix with moderate GEI
Ce <- simCmat(n_envs   = n_envs,
              mean_cor = 0.2,
              epsilon  = 0.95-0.2,
              rank = 6,
              skew = 0.5)
Ce[Ce > 1] <- 1
# plotCmat(Ce, den_order = T)$heatmap
# plotCmat(Ce, den_order = T)$hist
# Ge <- sqrt(diag(De)) %*% Ce %*% sqrt(diag(De))
# measureGEI(Ge, prop = T,disentangle = F)
# measureGEI(Ge, prop = T,disentangle = T)