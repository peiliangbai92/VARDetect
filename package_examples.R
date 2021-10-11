# !! Code contains commands to clear workspace objects. Make sure to save any changes before execution !!

rm(list=ls())
gc()

######## Loading Packages #######################
# library("mvtnorm")
# library("MTS")
# library("lattice")
# library("igraph")
# library("pracma")
# library("sparsevar")
# library("Rcpp")
# library("RcppArmadillo")
# 
# ######## Call Functions from the source code #########################
# source("VARDetect/R/Functions_BSS_SGL.R")
# source("VARDetect/R/Functions_LpS.R")
# sourceCpp("VARDetect/src/Functions_BSS_SGL.cpp")


######## load package from R CRAN
# Please uncomment the following code to load "VARDetect" package
# install.packages("VARDetect")
# library("VARDetect")
require(VARDetect)
stopifnot(packageVersion("VARDetect") >= "0.1.5")


###################
### Section 3.4 ###
###################

# section 3.4 and Figure 4
nob <- 4000; p <- 15
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m <- length(brk); q.t <- 1
sp_density <- rep(0.05, m * q.t)
signals <- c(-0.6, 0.6, -0.6)
try <- simu_var(method = "sparse", nob = nob, k = p, lags = q.t,
                brk = brk, sigma = diag(p), signals = signals,
                sp_density = sp_density, sp_pattern = "random", seed = 1)
data <- as.matrix(try$series)
fit <- tbss(data, method = "sparse", q = 1, refit = TRUE)

plot(fit, display = "cp")
plot(fit, display = "param")
plot(fit, display = "granger", threshold = 0.2, layout = "nicely")
plot(fit, display = "density", threshold = 0.1)

print(fit)
summary(fit)


###################
### section 3.5 ###
###################
nob <- 4000; p <- 15
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m <- length(brk); q.t <- 1
signals <- c(-0.6, 0.8, -0.6)
try_simu <- simu_tbss(nreps = 5, simu_method = "sparse", nob = nob, k = p, lags = q.t,
          brk = brk, sigma = diag(p), signals = signals,
          sp_pattern = "off-diagonal", est_method = "sparse", q = q.t, refit = TRUE)
summary(try_simu, critical = 5)

###################
### Section 4.1 ###
###################

######## Random sparse
nob <- 4000; p <- 15
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1);
m <- length(brk)
q.t <- 2
signals <- c(-0.6, -0.4, 0.6, 0.4, -0.6, -0.4)

sp_density <- rep(0.05, m * q.t)
try <- simu_var(method = "sparse", nob = nob, k = p, lags=q.t,
                brk = brk, sigma = diag(p), signals = signals,
                sp_density = sp_density, sp_pattern = "random", seed = 1)
print(plot_matrix(do.call("cbind", try$model_param), m * q.t))



######## 1-off diagonal:
nob <- 4000; p <- 15
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1);
m <- length(brk)
q.t <- 2
signals <- c(-0.6, -0.4, 0.6, 0.4, -0.6, -0.4)
try <- simu_var(method = "sparse", nob = nob, k = p, lags = q.t,
                brk = brk, sigma = diag(p), signals = signals,
                sp_pattern =  "off-diagonal", seed = 1)
# display the true transtion matrices
print(plot_matrix(do.call("cbind", try$model_param), m * q.t))
data <- try$series
data <- as.matrix(data)
MTS::MTSplot(data)

#run the tbss method
fit <- tbss(data, method = "sparse", q = q.t)

#display the estimated break points
print(fit)
summary(fit)

# display the estiamted transtion matrices
plot(fit, display = "cp")
plot(fit, display = "param")


####### Sparse with different time lags:
nob <- (10^3); p <- 15
brk <- c(floor(nob / 2), nob + 1)
m0 <- length(brk) -1
q.t <- 2
m <- m0 + 1
signals <- c(-0.8, 0.6, 0.4)
try <- simu_var('sparse', nob = nob, k = p, brk = brk, lags_vector = c(1,2),
              sp_pattern = "off-diagonal", signals = signals)
print(plot_matrix(do.call("cbind", try$model_param), m * q.t))
data <- try$series
data <- as.matrix(data)
MTS::MTSplot(data)

# lags selection (unknown lags)
library(sparsevar)
d_full <- c(1, 2, 3, 4)
BIC_full <- rep(0, length(d_full))
for(i in 1:length(d_full)){
  d <- d_full[i]
  #run the tbss method
  fit <- tbss(data, method = "sparse", q = d, refit = TRUE)
  sparse_mats <- fit$sparse_mats
  cp_est <- fit$cp   
  cp_full <- c(1, cp_est, nob+1)
  BIC <- 0
  for(j in  1:(length(cp_est) + 1)){
    data_temp <- as.matrix(data[(cp_full[j]): (cp_full[j+1]-1), ])
    n_temp <- dim(data_temp )[1]
    sparse_mats_temp <- sparse_mats[[j]]
    residual <- c()
    for(t in ( (d+1) :n_temp) ){
      y_pred <- 0
      for(dd in 1:d){
        phi <- sparse_mats_temp[, ( (dd-1)*p +1) :( (dd)*p ) ]
        y_pred <- y_pred +  phi%*% (data_temp[t-dd,])
      }
      residual <- cbind(residual, data_temp[t,] - y_pred)
    }
    sigma.hat <- 0*diag(p);
    for(t in 1: (n_temp -d)){sigma.hat <- sigma.hat +  residual[, t]%*%t(residual[, t]);  }
    sigma.hat <- (1/(n_temp - d))*sigma.hat;
    log.det <- log(det(sigma.hat));
    count <- sum(sparse_mats_temp !=0)
    BIC <- BIC + log.det + log((n_temp - d))*count/(n_temp - d)
  }
  BIC_full[i] <- BIC
  
}
BIC_full
#choose the one with the smallest BIC
d_full[which.min(BIC_full)]

#run the tbss method
fit <- tbss(data, method = "sparse", q = q.t, refit = TRUE)
#detected change points:
print(fit)

plot(fit, display = "cp")
plot(fit, display = "param")

###################
### Section 4.2 ###
###################


######column-wise separate across all lags
nob <- 4000
p <- 20
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m <- length(brk)
q.t <- 2


signals <- c(-0.8, -0.4, 0.6, -0.4, -0.8, -0.4)
num_group <- 3
group_index <- vector('list', num_group)
group_index[[1]] <- c(1, 5)
group_index[[2]] <- c(31)
try <- simu_var(method = "group sparse", nob = nob, k = p, lags = q.t,
                brk = brk, sigma = diag(p), signals = signals,
                group_index = group_index, group_type = "columnwise")
print(plot_matrix(do.call("cbind", try$model_param), m * q.t ))
data <- try$series
data <- as.matrix(data)
MTS::MTSplot(data)

# run the tbss method
fit <- tbss(data, method = "group sparse", q = q.t,
            group.case = "columnwise", group.index = as.list(c(0: (p * q.t - 1))))

#display the estimated break points
print(fit)

# display the estiamted transtion matrices
plot(fit, display = "cp")
plot(fit, display = "param")



######## row-wise simultaneously across all lags
nob <- 4000
p <- 20
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m <- length(brk)
q.t <- 2


signals <- c(-0.8, 0.4, 0.6, -0.3, -0.8, 0.4)
num_group <- q.t + 1
group_index <- vector('list', num_group)
group_index[[1]] <- c(1, 3)
group_index[[2]] <- c(1, 3) + p
try <- simu_var(method = "group sparse", nob = nob, k = p, lags = q.t,
                sigma = diag(p), brk = brk, signals = signals,
                group_index = group_index, group_type = "rowwise")
print(plot_matrix(do.call("cbind",try$model_param), m * q.t ))
data <- try$series
data <- as.matrix(data)
MTS::MTSplot(data)

#run the tbss method
group.index <- vector("list", p)
for(i in 1:p){
  group.index[[i]] <- rep(i - 1, q.t) + seq(0, p * (q.t - 1), p)
}
fit <- tbss(data, method = "group sparse", q = q.t,
            group.case = "rowwise", group.index = group.index)

#display the estimated break points
print(fit)

# display the estiamted transtion matrices
plot(fit, display = "cp")
plot(fit, display = "param")


###### Hierarchical row-wise  across all lags
nob <- 4000
p <- 10
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m <- length(brk)
q.t <- 2
signals <- c(-0.4, -0.4, 0.4, -0.4, -0.4, -0.4)

num_group <- q.t + 1
group_index <- vector('list', num_group)
group_index[[1]] <- c(1, 3, 10)
group_index[[2]] <- c(3 + p)
try <- simu_var(method = "group sparse", nob = nob, k=p, lags = q.t, sigma = diag(p),
                brk = brk, signals = signals, seed = 1,
                group_index = group_index, group_type = "rowwise")
print(plot_matrix(do.call("cbind",try$model_param), m * q.t ))
data <- try$series
data <- as.matrix(data)
ts.plot(data)


#run the tbss method
group.index <- vector("list", p * q.t)
for(i in 1:p){
  for(j in 1:q.t){
    if(j == 1){
      group.index[[(j - 1) * p + i]] <- c((q.t - 1) * p + i) - 1

    }else{
      group.index[[(j - 1) * p + i]] <-rep(i - 1, q.t) + seq(0, p * (q.t - 1), p)
    }
  }
}
#run the tbss method
fit <- tbss(data, method = "group sparse", q = q.t,
            group.case = "rowwise", group.index = group.index)

#display the estimated break points
print(fit)

# display the estiamted transtion matrices
plot(fit, display = "cp")
plot(fit, display = "param")


###################
### Section 4.3 ###
###################
### fixed low rank plus time-varying sparse structure
### new example in Section 4.5
nob <- 300; p <- 15
brk <- c(floor(nob/3), floor(2*nob/3), nob+1)
m <- length(brk)
signals <- c(-0.7, 0.85, -0.7)
rank <- rep(2, m)
singular_vals <- c(1, 0.75)
info_ratio <- rep(0.35, 3)
q.t <- 1
try <- simu_var(method = "fLS", 
                nob = nob, 
                k = p, 
                lags = q.t, 
                brk = brk,
                sigma = as.matrix(diag(p)), 
                signals = signals,
                seed = 1,
                rank = rank,
                singular_vals = singular_vals, 
                info_ratio = info_ratio,
                sp_pattern = "off-diagonal", 
                spectral_radius = 0.9)
data <- as.matrix(try$series)
MTS::MTSplot(data)
print(plot_matrix(do.call("cbind", try$model_param), m * q.t))
fit <- tbss(data, method = "fLS", mu = 150)
print(fit)
plot(fit, display = "cp")
plot(fit, display = "param")


###################
### Section 4.4 ###
###################
### Low rank plus sparse structure
nob <- 300
p <- 20
brk <- c(floor(nob / 3), floor(2 * nob / 3), nob + 1)
m0 <- length(brk)-1
q.t <- 1
m <- m0 + 1
signals <- c(-0.7, 0.8, -0.7)
rank <- c(1, 3, 1)
singular_vals <- c(1, 0.75, 0.5)
info_ratio <- rep(0.35, 3)

# generating data process
try <- simu_var(method = "LS", nob = nob, k = 20, lags = 1, brk = brk, sigma = as.matrix(diag(p)), signals = signals,
                rank = rank, singular_vals = singular_vals, info_ratio = info_ratio, sp_pattern = "off-diagonal",
                spectral_radius = 0.9)
data <- try$series
MTS::MTSplot(data)
print(plot_matrix(do.call("cbind", try$model_param), m*q.t))

lambda1 = lambda2 = lambda3 <- c(2.5, 2.5)
mu1 = mu2 = mu3 <- c(15, 15)

# fit without cv selection
fit <- lstsp(data, lambda.1 = lambda1, mu.1 = mu1,
             lambda.2 = lambda2, mu.2 = mu2,
             lambda.3 = lambda3, mu.3 = mu3, alpha_L = 0.25,
             step.size = 5, niter = 20, skip = 5,
             cv = FALSE, verbose = FALSE)
print(fit)
plot(fit, display = "cp")
plot(fit, display = "param")


#############################################
### Section 4.5: real dataset application ###
#############################################
####### Stock-return dataset
### loading data
data(weekly)  # if VARDetect package is already installed, please use data(...) function to load data
#load("VARDetect/data/weekly.RData")

# use sparse to fit
set.seed(100)
lambda.1.max <- 1e-2
nlam <- 20
lambda.1.min <-  lambda.1.max * 1e-3
delata.lam <- (log(lambda.1.max) - log(lambda.1.min)) / (nlam - 1)
lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min * exp(delata.lam * (nlam - jjj)))
fit <- tbss(weekly, method = "sparse",
            lambda.1.cv = lambda.1.cv,
            lambda.2.cv = 0.05,
            block.size = 8, an.grid = c(10, 15, 18))
print(fit)
plot(fit, display = "cp")

# use group sparse to fit
lambda.1.max <- 0.01
nlam <- 20
lambda.1.min <-  lambda.1.max * 1e-4
delata.lam <- (log(lambda.1.max) - log(lambda.1.min)) / (nlam - 1)
lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min * exp(delata.lam * (nlam - jjj)))
fit_group <- tbss(weekly, method = "group sparse",
                  group.case = "columnwise", max.iteration = 50,
                  lambda.1.cv = lambda.1.cv, lambda.2.cv = 0.1,
                  block.size = 8, an.grid = c(10, 15, 18))
print(fit_group)
plot(fit_group, display = "cp")

### plot Granger causal networks for partitioned segments
# for sparse structure, we only present segments 2, 3, and 6
plot(fit, display = "granger", threshold = 0.2)

# for group sparse structure, we present segments 2, 3, and 7
plot(fit_group, display = "granger", threshold = 0.2)

### plot sparsity levels for segments
plot(fit, display = "density", threshold = 0.2)
plot(fit_group, display = "density", threshold = 0.2)

# use LSTSP to fit
### use LSTSP to investigate the change points for the dataset
data <- weekly
n <- dim(data)[1]
k <- dim(data)[2]

lambda.1 <- c(0.015, 0.015)
mu.1 <- c(0.05, 0.05)
lambda.2 <- c(0.01, 0.01)
mu.2 <- c(0.05, 0.05)
lambda.3 <- c(0.008, 0.008)
mu.3 <- c(0.0335, 0.0335)
N <- n - 1
omega <- (0.0015) * (((log(N))^1)*log(k))

fit <- lstsp(data, lambda.1 = lambda.1, mu.1 = mu.1, 
             lambda.2 = lambda.2, mu.2 = mu.2, 
             lambda.3 = lambda.3, mu.3 = mu.3, 
             h = 80, step.size = 40, omega = omega,
             niter = 20, skip = 5, verbose = TRUE)

print(fit)
plot(fit, display = 'cp')
plot(fit, display = 'density', threshold = 0.05)

### comparison with SBS algorithm
library("hdbinseg")  ### load SBS package

time.stamp <- data[,1]
fit <- sbs.alg(t(weekly), cp.type = 2)
idx <- fit$ecp
cp_times <- time.stamp[idx]
MTS::MTSplot(weekly)
abline(v = idx, lwd = 2, col = 'red')



############################################# EEG dataset
### loading data
#load("VARDetect/data/eeg.RData")
load(url("https://github.com/peiliangbai92/VARDetect/blob/main/data/eeg.RData?raw=true"))

### BIC select time lag by using TBSS
bic_scores <- rep(0, 5)
lambda.1.cv <- c(0.1)
lambda.2.cv <- c(0.001)
ts = data
for(l in 1:5){
  fit <- tbss(as.matrix(data), method = "sparse", q = l, 
              lambda.1.cv = lambda.1.cv, 
              lambda.2.cv = lambda.2.cv, 
              max.iteration = 20, 
              block.size = floor(0.8*sqrt(n)), 
              an.grid = c(120, 150, 200))
  selected_cps <- fit$cp
  ### refit the functional connectivity network
  segments <- c(1, selected_cps, n)
  est_mats <- NULL
  res <- 0
  for(j in 1:(length(segments)-1)){
    s <- segments[j]+50
    e <- segments[j+1]-50
    
    seg_data <- ts[s:e,]
    refit <- fitVAR(seg_data, l, nfolds=5)
    res <- res + log(det(refit$sigma))
    df <- 0
    for(ll in 1:l){
      est_mat <- refit$A[[ll]]
      for(row in 1:p){
        for(col in 1:p){
          if(abs(est_mat[row, col]) > 0.01){
            df <- df + 1
          }
        }
      }
    }
    res <- res + log(e-s)/(e-s) * (df^2)
  }
  bic_scores[l] <- res
  print('===============================================================')
}
q.t <- which.min(bic_scores)

### use TBSS
lambda.1.cv <- c(0.1)
# lambda.2.cv <- c(1)*sqrt(log(p)/n)
lambda.2.cv <- c(0.001)
fit <- tbss(as.matrix(data), method = "sparse", q = 1, 
            lambda.1.cv = lambda.1.cv, 
            lambda.2.cv = lambda.2.cv, 
            block.size = floor(0.8*sqrt(n)), 
            an.grid = c(120, 150, 200), refit = FALSE)
print(fit)
plot(fit, display = 'cp')
plot(fit, display = 'density', threshold = 0.25)
plot(fit, display = 'granger', threshold = 0.75)


### use LSTSP
lambda.1 = lambda.2 <- c(0.5, 0.5)
mu.1 = mu.2 <- c(5, 5)
lambda.3 <- c(0.35, 0.35)
mu.3 <- c(200, 200)
N <- n-1
lambda.2 <- c((1/1)*(log(N)*log(p))/N, (1/1)*(log(N)*log(p))/N)
mu.2 <- c(1, 1)
omega <- (250)*(((log(N))^1)*log(p))
h <- 8*floor(sqrt(n))+1
steps <- floor(0.45*h)

fit <- lstsp(as.matrix(data), lambda.1 = lambda.1, mu.1 = mu.1, 
             lambda.2 = lambda.2, mu.2 = mu.2, 
             lambda.3 = lambda.3, mu.3 = mu.3, 
             omega = omega, h = h, step.size = steps, skip = 125, 
             verbose = TRUE)
ranks <- rep(0, length(fit$cp)+1)
for(i in 1:(length(fit$cp)+1)){
  ranks[i] <- qr(fit$lowrank_mats[[i]])$rank
  print(ranks[i])
}
print(fit)
plot(fit, display = 'cp')
plot(fit, display = 'density', threshold = 0.25)
plot(ranks, ylab = 'rank', type = 'o')

