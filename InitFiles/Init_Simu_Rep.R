library(HDMT)
library(MASS)
library(qvalue)
library(DACT)
library(locfdr)
library(fdrtool)
library(doMC)
library(ks)
library(kernelboot)
library(ssa)
library(adaFilter)
library(fdrtool)
library(repfdr)
#library(knockoff)
library(doParallel)
#library(RMediation)
# Load methods
# source("Methods/Rect.R")
# source("Methods/Rect2.R")
# source("Methods/Oracle_FDR.R")
# source("Methods/cutoff_gen.R")
# source("Methods/Symmetric_fdr.R")
# source("Methods/Empirical_fdr2.R")
# source("Methods/Mirror_fdr_est.R")
# source("Methods/TwoDMirror.R")

# Load necessary Tools
source("Tools/Basic.R")
source("Tools/Fit_Med.R")
source("Tools/mainFun.R")
source("Tools/Data_Generating_Process.R")

# Load Methods
source("Methods/DACT_z_stat.R")
source("Methods/Sobel_Test.R")
source("Methods/Pmax_Test.R")
source("Methods/MT_Comp.R")
source("Methods/conditionPC.R")
Rcpp::sourceCpp("Methods/JointMirror.cpp")
Rcpp::sourceCpp("Methods/modules.cpp")
source("Methods/loadJointMirror.R")
source("Methods/JointMirrorR.R")

## Read the information in bash
args <- commandArgs(trailingOnly = TRUE)
DGP_name = args[1] 
H1_type <- args[2]
alter_frac <- args[3]
alter_frac <- as.numeric(alter_frac)
is.oneside = (H1_type=="OneSide")

# Basic Setting
num_cores <- 20
nsim <- 100
beta0 <- 0.3
# simu_n <- 250
# simu_p <- 5000

set.seed(222222)
seeds <- sample(20000:99999999,nsim,replace=FALSE)

trgt.fdr.level <- 0.2

# Parameters for Symmetric Algorithm
num_Split <- 50
prop_Split <- 1/3


# Parameters for generating M X and Y
X_prob <- 0.2
X_type <- "binom"
sigma_X <- 5
sigma.M <- 1
sigma.Y <- 1
Sigma.MY <- diag(2)

# Parameters for Replicate K experiments 
simuM <- 10000

## Create Result Folder
dir.fold <- paste0("Res/",paste0(DGP_name,"_",H1_type,"_",alter_frac))
if(!dir.exists(dir.fold)){
  dir.create(dir.fold)
}

dir.subfold <- paste0(dir.fold,"/","M_",simuM)
if(!dir.exists(dir.subfold)){
  dir.create(dir.subfold)
}
