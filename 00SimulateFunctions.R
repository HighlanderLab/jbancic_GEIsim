# Script name: Functions to simulate GxE interaction
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# This script contains all functions required to simulate GxE interaction


################################################

# This function simulates a between-environment genetic correlation matrix

simCmat = function (n_groups = 1,   # no. environment n_groups
                    n_envs = 20,    # no. environments within each group
                    mean_cor = 0.5, # base correlation within n_groups
                    delta = 0,      # base correlation between n_groups
                    epsilon,        # range in correlations within and across groups
                    rank = 6,       # rank of the noise (no. eigenvectors)
                    skew = 0){      # amount of skewness, no skewness = 1 
  
  if(length(n_envs) == 1 & n_groups > 1){n_envs <- rep(n_envs, n_groups)}
  if(length(mean_cor) == 1 & n_groups > 1){mean_cor <- rep(mean_cor, n_groups)}
  
  if(skew < -1 | skew > 1){
    stop("The skew value must be between -1 and 1")
  }
  skew.abs <- abs(skew)
  tot_envs <- sum(n_envs)
  base_cor <- matrix(rep(delta, tot_envs * tot_envs), ncol = tot_envs)
  for(i in 1:n_groups){
    cor <- matrix(rep(mean_cor[i], n_envs[i] * n_envs[i]), ncol = n_envs[i])
    if (i == 1){base_cor[1:n_envs[1], 1:n_envs[1]] <- cor}
    if (i != 1){base_cor[(sum(n_envs[1:(i - 1)]) + 1):sum(n_envs[1:i]), 
                         (sum(n_envs[1:(i - 1)]) + 1):sum(n_envs[1:i])] <- cor}
  }
  diag(base_cor) <- 1 - epsilon
  eivect <- c()
  for (i in 1:tot_envs){
    ei <- runif(rank, -1, 1*(1-skew.abs))
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  Cmat <- base_cor + t(eivect) %*% eivect
  if (skew < 0){
    Cmat <- -Cmat
  }
  return(Cmat)
  # bend non-positive definite matrices
}


################################################

# This function quantifies the variance explained by non-crossover 
# and crossover interaction in a between-environment genetic correlation matrix

measureGEI <- function(cov_mat,
                       prop = TRUE,
                       disentangle = FALSE,
                       groups = NULL,
                       old_measure = FALSE){
  if(!isSymmetric(cov_mat)) {stop("cov_mat is not symmetric")}
  if(!mean(cov_mat) > 0) {stop("cov_mat is not positive (semi) definite")}
  tot_var1 <- mean(diag(cov_mat))
  tot_var2 <- 1
  if(prop == F){
    tot_var1 <- 1
    tot_var2 <- mean(diag(cov_mat))}
  if (disentangle == FALSE) {
    Gvar   <- sum(cov_mat)/ncol(cov_mat)^2/tot_var1 
    # equivalent to: mean(cov_mat)/tot_var1
    GEIvar <- tot_var2 - Gvar
    var <- list(Gvar = Gvar, GEIvar = GEIvar)
  } 
  if (disentangle == T) {
    if (old_measure == F) {
      Gvar   <- sum(cov_mat)/ncol(cov_mat)^2/tot_var1
      Hetvar <- sum((sqrt(diag(cov_mat)) - mean(sqrt(diag(cov_mat))))^2)/ncol(cov_mat)/tot_var1 # het of variance
      nGEI <- Gvar + Hetvar
      cGEI <- tot_var2 - nGEI
    } else {
      # if(length(unique(sign(ss$u[,1]))) > 1){stop("First principal component captures cGEI, not just nGEI")}
      nGEI <- svd(cov_mat,nu = 1)$d[1]/ncol(cov_mat)/tot_var1
      cGEI <- tot_var2 - nGEI
    }
    var <- list(nGEI = nGEI, cGEI = cGEI)
  }
  return(var)
  # groups is a list with items in each group, can be names or indexes
  # initial check to see if all indexes are in cov_mat
  # if(!is.null(groups)){
  #   # within
  #
  #   # between
  #   bGvar <- sum(cov_mat)/ncol(cov_mat)/tot_var
  #   var <- list(#within = list(Gvar = wGvar, GEIvar = wGEIvar),
  #               #between = list(Gvar = bGvar, GEIvar = bGEIvar),
  #               overall = list(Gvar = Gvar, GEIvar = GEIvar))
  # }
}
# measureGEI(cov_mat = Ge, prop = T, disentangle = F)
# measureGEI(cov_mat = Ge, prop = T, disentangle = T)
# measureGEI(cov_mat = Ge, prop = T, disentangle = T, old_measure = T)

# to add: disentangle GxY, GxE,
# within and across, overall



################################################

# This functions plots and orders between-environment genetic correlation matrix

plotCmat <- function(cor_mat, 
                     den_order = TRUE, 
                     groups = NULL,
                     axis_title = "Environment"){
  require(cluster)
  require(tidyr)
  require(ggplot2)
  
  tot_envs <- ncol(cor_mat)
  env_order1 <- 1:tot_envs
  if(is.null(colnames(cor_mat)) | is.null(rownames(cor_mat))){colnames(cor_mat) <- rownames(cor_mat) <- env_order1}
  cor_mat <- as.matrix(cor_mat)
  if(isSymmetric(cor_mat) != TRUE){stop("Matrix is not symmetric")}
  Vmat = FALSE
  if (length(table(diag(cor_mat))) > 1) {
    Vmat = TRUE; print("Matrix does not have 1s on the diagonal. Assuming this is a variance matrix")
  }
  env_order2 <- colnames(cor_mat)
  if (!is.null(groups)){
    if(length(unlist(groups)) != tot_envs){stop("Groups specified do not conform with 'cor_mat'")}
    if(!is.list(groups)){stop("'groups' must be a list")}
    n_groups <- length(groups)
    n_envs <- unlist(lapply(groups, function(x) length(x)))
  }
  if(den_order){
    if (!is.null(groups)){
      dis_mat_list <- lapply(groups, function(x) 1 - cor_mat[x,x])
      env_order1 <- unlist(lapply(dis_mat_list, function(x) colnames(x)[agnes(x = x, diss = TRUE, method = "average")$order]))
      env_order2 <- colnames(cor_mat)[as.numeric(env_order1)]
    }
    if(is.null(groups)){
      dis_mat <- 1 - cor_mat
      env_order1 <- agnes(x = dis_mat, diss = TRUE, method = "average")$order
      env_order2 <- colnames(cor_mat)[as.numeric(env_order1)]
    } 
  }
  cor.df <- gather(as.data.frame(cor_mat[env_order2, env_order2]))
  names(cor.df) <- c('Env1', "Cor")
  cor.df$Env1 <- factor(cor.df$Env1, levels = rev(env_order2))
  cor.df$Env2 <- factor(rep(env_order2, times = tot_envs), levels = env_order2)
  cor.df <- cor.df[order(cor.df$Env1, cor.df$Env2),]
  rownames(cor.df) <- NULL
  cor.df <- cor.df[,c("Env1","Env2","Cor")]
  if (Vmat == FALSE) {cor.df$Cor[cor.df$Env1 == cor.df$Env2] <- NA}
  
  # hh <- rev(rainbow(256, start = 0, end = 2/3))
  hh <- rev(rainbow(256, start = 0, end = 2/3))[c(5:100,150:250)]
  pp <- ggplot(data = cor.df, aes(y= Env1, x = Env2, fill=Cor)) + geom_tile() +
    labs(x = axis_title, y = axis_title) +
    scale_fill_gradientn(colours = hh, na.value = "white", limits= c(-1.1,1.1)) +
    theme(axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -1))
  if (Vmat == TRUE) {
    # hh <- hsv(seq(0.5,0.65,length.out = 100))
    hh <- hsv(seq(0.5,0.68,length.out = 200), 1, 0.9)
    pp <- pp + scale_fill_gradientn(colours = hh, na.value = "white") +
      labs(fill = "Var")
  }
  
  cor.df2 <- data.frame(Source = "All correlations",
                        Cor = cor.df[order(cor.df$Env1, rev(cor.df$Env2)),][lower.tri(cor_mat, diag= F),]$Cor)
  cor.df2$Mean <- mean(cor.df2$Cor)
  qq <- ggplot(data = cor.df2, aes(Cor)) + geom_histogram(bins = 20, color = "black", fill = "grey") +
    labs(x = paste0("Pair-wise correlations between ", tolower(axis_title),"s")) + 
    geom_vline(aes(xintercept = Mean), linetype = 2, colour = "steelblue", linewidth = 0.6) + 
    theme(axis.title.y = element_text(vjust = 2),
          axis.title.x = element_text(vjust = -1))
  
  if (!is.null(groups)){
    frames <- data.frame(Env1 = cumsum(c(0.5, n_envs))[1:n_groups],
                         Env2 =  cumsum(c(tot_envs + 0.5, -n_envs))[1:n_groups])
    pp <- pp + geom_rect(data = frames, linewidth = 0.5, fill = NA, colour = "black",
                         aes(xmin = Env1, xmax = Env1+n_envs,
                             ymin = Env2, ymax = Env2-n_envs))
    
    cor.df3 <- data.frame(Source = "Correlations within groups",
                          Cor = unlist(lapply(groups, function(x) cor_mat[x,x][upper.tri(cor_mat[x,x])])))
    cor.df3$Mean <- mean(cor.df3$Cor)
    cor_mat2 <- cor_mat
    for(i in 1:n_groups){cor_mat2[groups[[i]],groups[[i]]] <- NA}
    cor.df4 <- data.frame(Source = "Correlations between groups",
                          Cor = cor_mat2[upper.tri(cor_mat2)])
    cor.df4 <- cor.df4[!is.na(cor.df4$Cor),]
    cor.df4$Mean <- mean(cor.df4$Cor)
    cor.df5 <- rbind(cor.df2, cor.df3, cor.df4)
    cor.df5$Source <- factor(cor.df5$Source, levels = c("Correlations within groups", "Correlations between groups", "All correlations"))
    qq <- ggplot(data = cor.df5, aes(Cor)) + geom_histogram(bins = 20, color = "black", fill = "grey") +
      labs(x = paste0("Pair-wise correlations between ", tolower(axis_title),"s")) + 
      geom_vline(aes(xintercept = Mean), linetype = 2, colour = "steelblue", linewidth = 0.6) + 
      theme(axis.title.y = element_text(vjust = 2),
            axis.title.x = element_text(vjust = -1)) + facet_wrap(~Source)
  }
  temp <- list(heatmap = pp, hist = qq, data = cor.df, order = env_order1)
  return(temp)
}

################################################

# This function a between-environment genetic correlation matrix

simPheno <- function(n_genos,            # number of genotypes
                     Cmat,               # between-environment genetic correlation matrix
                     De,                 # environment genetic variance matrix
                     n_reps = 2,         # number of replications per environment
                     mu = 1,             # phenotypic mean
                     varE = 1,           # error variance
                     h2 = 0.5,           # narrow-sense heritabilty 
                     subset_envs = NULL, 
                     rank = NULL, 
                     prop_spatial = 0,   # proportion of variance explained by spatial variation
                     n_cols = NULL, 
                     n_rows = NULL) 
{
  require(FieldSimR)
  # Sim parameters -------------
  Ce <- as.matrix(Cmat)
  if (isSymmetric(Ce) == F) {
    stop("Cmat not symmetric")
  }
  n_envs <- dim(Ce)[1]
  
  # Simulate environmental main effects -------------
  tau <- scale(rnorm(n_envs))*sqrt(varE)  # environment (E) main effects
  hist(tau,20, main = "Environment means")
  X <- kronecker(diag(n_envs), rep(1,n_genos*n_reps))  # design matrix
  
  # Simulate genetic effects -------------
  # add to call if simulating variances inside the function meanG = 1, varG = 0.2,
  # De <- diag(abs(rnorm(n_envs))) # diagonal genetic variance matrix between environments
  # De <- diag(c(scale(rgamma(n_envs, shape = 2, scale = 0.3))*sqrt(varG) + meanG))
  # De <- diag(rep(1,n_envs)) # homozygous
  if (length(De) != n_envs) {stop("The length of De is not equal to dimension of Cmat")}
  hist(De,20, main = "Genetic variances")
  De <- diag(De)
  Ge <- sqrt(De) %*% Ce %*% sqrt(De)  # genetic variance matrix
  # j=1
  # while(length(unique(sign(svd(Ge)$u[,1]))) > 1){
  #   cat("Resampling genetic variances:",j,"\n")
  #   De <- diag(c(scale(rgamma(n_envs, shape = 2, scale = 0.3))*sqrt(varG) + meanG))
  #   hist(diag(De),20, main = "Genetic variances")
  #   Ge <- sqrt(De) %*% Ce %*% sqrt(De)  # genetic variance matrix
  #   j=j+1
  # }
  # cumsum(svd(Ge)$d)/sum(svd(Ge)$d)
  S <- svd(Ge)$u
  D <- diag(svd(Ge)$d)
  slopes <- scale(matrix(rnorm(n_genos*n_envs), ncol = n_envs))  # standardise the slopes for each term
  slopes <- slopes %*% sqrt(D)  # assign the appropriate variances to each term
  u <- c(slopes %*% t(S)) # genotype by environment (GxE) interaction effects
  
  # Data frame of genetic effects -------------
  ge_df <- data.frame(env = factor(rep(1:n_envs, each = n_genos*n_reps)),
                      rep = factor(1:n_reps),
                      id  = factor(rep(1:n_genos, each = n_reps)),
                      ge.Trait.1 = rep(u, each = n_reps))
  ge_df <- ge_df[order(ge_df$env, ge_df$rep),]
  # simple randomisation of id's within reps in envs
  ge_df$ord <- unlist(lapply(1:(n_envs*n_reps), function(x) sample(1:n_genos)))  
  ge_df <- ge_df[order(ge_df$env, ge_df$rep, ge_df$ord),]  # reorder based on randomisation
  rownames(ge_df) <- NULL
  # Z <- model.matrix(ge.Trait.1 ~ id:env - 1, data  = ge_df)  # design matrix
  
  # Simulate spatial plot errors using FieldSimR
  if ((is.null(n_cols) & is.null(n_rows)) == T) {stop("Specify n_cols and n_rows")}
  hist(h2 <- abs(rnorm(n_envs, h2, 0.1)), main = "Heritabilities")
  # h2 <- rep(0.4, n_envs) # HOMOGENEOUS
  var_R <- diag(De)/h2 - diag(De)
  # hist(diag(De)/(diag(De) + var_R))         # plot level heritability's for each environment
  # hist(diag(De)/(diag(De) + var_R/n_reps))  # line mean heritability's for each environment
  # Line mean overall
  error_df <- field_trial_error(n_envs   = n_envs,
                                n_traits = 1,
                                n_blocks = n_reps,
                                n_cols   = n_cols,
                                n_rows   = n_rows,
                                var_R    = var_R,
                                prop_spatial = prop_spatial,
                                ext_dir  = "row")
  # plot_effects(error_df[error_df$env == 1,], effect = "e.Trait.1") # quick plot of the errors
  e <- error_df$e.Trait.1
  
  # Create phenotypes -------------
  # y <- mu + X %*% tau + Z %*% u + e
  Xtau <- rep(tau, each = n_genos*n_reps)
  y <- mu + Xtau + ge_df$ge.Trait.1 + e
    
  # Create MET dataset -------------
  met.df <- data.frame(error_df[, c("env", "block", "col", "row")], 
                       id = ge_df$id,
                       u = ge_df$ge.Trait.1,
                       e,
                       y)
  # Return object
  return(df <- list(Ce      = Ce,
                    Ge      = Ge,
                    trueG   = rowMeans(matrix(u, ncol = n_envs, byrow = F)), # main effs
                    trueGE = matrix(u, ncol = n_envs, byrow = F), # genotype by environment (GxE) interaction effects
                    varG    = diag(De),
                    varE    = var_R,
                    met.df  = met.df))
}


# test <- simPheno(Ce, n_genos = 100, mu = 1, varE = 1, n_reps = 2, h2 = 0.3, prop_spatial = 0, n_cols = 20, n_rows = 10)
# str(test)
# cor(test$y,test$u)