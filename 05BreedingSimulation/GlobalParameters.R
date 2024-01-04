# Global parameters script contains key parameters for simulation of
# a breeding program  and founder parents

#=======================================================================
##-- Stage-specific parameters
#=======================================================================
# Crossing block
nParents  = 30  # no. of parents (and founders)
nCrosses  = 100 # no. of crosses per year

# HDRW
nDH       = 100 # no. of DH lines produced per cross
nEnvHDRW  = 1   # h2=0.07
nRepHDRW  = 1
varE_HDRW = 8 * 1.474541 # 1.474541 is the mean of simulated variance

# PYT
famMax   = 10 # max no. of DH lines per cross
nPYT     = 500
selPYT   = "estMean" # Selection strategy
nEnvPYT  = 2       # h2=0.2
nRepPYT  = 1
varE_PYT = 4 * 1.474541 # 1.474541 is the mean of simulated variance

# AYT
nAYT     = 50
selAYT   = "estMean" # Selection strategy
nEnvAYT  = 5       # h2=0.66
nRepAYT  = 2
varE_AYT = 4 * 1.474541 # 1.474541 is the mean of simulated variance

# EYT
nEYT     = 10
selEYT   = "estMean" # Selection strategy
nEnvEYT  = 20      # h2=0.90
nRepEYT  = 2
varE_EYT = 2 * 1.474541 # 1.474541 is the mean of simulated variance

##-- Number of cycles and environments in simulation
nCycles = 30  
nEnvs   = nEnvEYT # max no. of envs in a year
nYearsTP = 5 # number of years of training population

##-- Founder parents mean and variance
# Decompose between environment variance matrix into k terms 
# kTerms should be lower than the maximum number of environments at any
# given time in the simulation
kTerms = nEnvs*nYearsTP

# Mean and variance of additive effects
initMeanG = rep(0, times = kTerms) # kTerms for scores
initVarG  = rep(1, times = kTerms)
corA      = diag(kTerms)

# Number of QTL per chromosome
nQtl = 300

# Number of SNP per chromosome
nSnp = 600

#=======================================================================
#-- Output variables
#=======================================================================
#-- Initialize a placeholder for tracked variables
trackParams = data.frame(matrix(NA,nrow = 0,ncol = 0))
