# Global parameters script contains key parameters for simulation of
# a breeding program  and founder parents

##-- Stage-specific parameters
# Crossing block
nParents  = 30  # no. of parents (and founders)
nCrosses  = 100 # no. of crosses per year

# HDRW
nDH       = 100 # no. of DH lines produced per cross
nEnvHDRW  = 1   # h2=0.07
nRepHDRW  = 1
varE_HDRW = 8

# PYT
famMax   = 10 # max no. of DH lines per cross
nPYT     = 500
nEnvPYT  = 2  # h2=0.2
nRepPYT  = 1
varE_PYT = 8

# AYT
nAYT     = 50
nEnvAYT  = 5  # h2=0.66
nRepAYT  = 2
varE_AYT = 8

# EYT
nEYT     = 10
nEnvEYT  = 20  # h2=0.90
nRepEYT  = 2
varE_EYT = 8

##-- Number of cycles and environments in simulation
nCycles = 20  
nEnvs   = nEnvEYT # max no. of envs in a year

##-- Founder parents mean and variance
# Additive effects
initMeanG  = rep(0, times = nEnvs)
initVarG   = rep(1, times = nEnvs)
corA       = diag(nEnvs)

# Dominance effects (assumed 0)
initMeanD  = rep(0, times = nEnvs)
initVarD   = rep(0, times = nEnvs)
corD       = diag(nEnvs)

# Number of QTL per chromosome
nQtl = 300

# Number of SNP per chromosome
nSnp = 600