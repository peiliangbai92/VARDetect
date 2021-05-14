## usethis namespace: start
#' @useDynLib VARDetect, .registration = TRUE
## usethis namespace: end
NULL


##########################################################################################
################ Main function: the main function for LpS detection algo #################
##########################################################################################
#' Main function for the low rank plus sparse structure VAR model
#' @description Main function for the low-rank plus sparse structure VAR model
#' @param data A n by p dataset matrix
#' @param lambda.1 tuning parameter for sparse component for the first step
#' @param mu.1 tuning parameter for low rank component for the first step
#' @param lambda.1.seq a sequence of lambda to the left segment for cross-validation, it's not mandatory to provide
#' @param mu.1.seq a sequence of mu to the left segment, low rank component tuning parameter
#' @param lambda.2 tuning parameter for sparse for the second step
#' @param mu.2 tuning parameter for low rank for the second step
#' @param lambda.3 tuning parameter for estimating sparse components
#' @param mu.3 tuning parameter for estimating low rank components
#' @param alpha_L a positive numeric value, indicating the restricted space of low rank component, default is 0.25
#' @param omega tuning parameter for information criterion, the larger of omega, the fewer final selected change points
#' @param h window size of the first rolling window step
#' @param step.size rolling step
#' @param tol tolerance for the convergence in the second screening step, indicates when to stop
#' @param niter the number of iterations required for FISTA algorithm
#' @param backtracking A boolean argument to indicate use backtrack to FISTA model
#' @param skip The number of observations need to skip near the boundaries
#' @param cv A boolean argument, indicates whether the user will apply cross validation to select tuning parameter, default is FALSE
#' @param nfold An positive integer, the number of folds for cross validation
#' @param verbose If is TRUE, then it will print all information about current step.
#' @return A list object including 
#' \describe{
#'     \item{cp}{Final estimated change points}
#'     \item{sparse_mats}{Final estimated sparse components}
#'     \item{lowrank_mats}{Final estimated low rank components}
#'     \item{est_phi}{Final estimated model parameter, equals to sum of low rank and sparse components}
#' }
#' @export
#' @examples
#' \donttest{
#' nob <- 100
#' p <- 15
#' brk <- c(50, nob+1)
#' rank <- c(1, 3)
#' signals <- c(-0.7, 0.8)
#' singular_vals <- c(1, 0.75, 0.5)
#' info_ratio <- rep(0.35, 2)
#' try <- simu_var(method = "LS", nob = nob, k = p, lags = 1, brk = brk, 
#'                 sigma = as.matrix(diag(p)), signals = signals, 
#'                 rank = rank, singular_vals = singular_vals, info_ratio = info_ratio, 
#'                 sp_pattern = "off-diagonal", spectral_radius = 0.9)
#' data <- try$series
#' 
#' lambda1 = lambda2 = lambda3 <- c(2.5, 2.5)
#' mu1 = mu2 = mu3 <- c(15, 15)
#' fit <- lstsp(data, lambda.1 = lambda1, mu.1 = mu1, 
#'              lambda.2 = lambda2, mu.2 = mu2, 
#'              lambda.3 = lambda3, mu.3 = mu3, alpha_L = 0.25, 
#'              step.size = 5, niter = 20, skip = 5, 
#'              cv = FALSE, verbose = FALSE)
#' print(fit$cp)
#' }
lstsp <- function(data, lambda.1 = NULL, mu.1 = NULL, lambda.1.seq = NULL, mu.1.seq = NULL, 
                  lambda.2, mu.2, lambda.3, mu.3, alpha_L = 0.25, omega = NULL, h = NULL, 
                  step.size = NULL, tol = 1e-4, niter = 100, backtracking = TRUE, skip = 5, 
                  cv = FALSE, nfold = NULL, verbose = FALSE){
    ###### model dimensions
    n <- dim(data)[1]
    p <- dim(data)[2]
    if(is.null(h)){
        h <- 2*floor(sqrt(n))+1
    }
    
    ###### Step 1: rolling-window find candidate change points ######
    candi_pts <- first.step.detect(data = data, h = h, step.size = step.size, 
                                   lambda = lambda.1, mu = mu.1, alpha_L = alpha_L, 
                                   skip = skip, lambda.1.seq = lambda.1.seq, mu.1.seq = mu.1.seq,
                                   cv = cv, nfold = nfold, verbose = verbose)
    
    ###### Step 2: screening out redundant candidate and obtain the final change points ######
    N <- n - 1
    if(is.null(omega)){
        omega <- (1/0.05) * (((log(N))^1)*log(p))
    }
    screen_select <- second.step.detect(data, candi_pts, omega, lambda.2, mu.2, 
                                        alpha_L = alpha_L, verbose = verbose)
    final.pts <- screen_select$pts
    
    ###### Step 3: estimate model parameters by using final change points ######
    pts <- c(1, final.pts, n)
    est_phi = sp_mats = lr_mats <- vector('list', length(pts)-1)
    for(i in 1:(length(pts)-1)){
        s <- pts[i]
        e <- pts[i+1]
        interval_data <- data[s:e, ]
        ns <- dim(interval_data)[1]
        
        X <- interval_data[1:(ns-1), ]
        y <- interval_data[2:ns, ]
        fit <- fista.LpS(X, y, lambda = lambda.3, mu = mu.3, alpha_L = alpha_L, 
                         niter = niter, backtracking, x.true = diag(p))
        sp_mats[[i]] <- t(fit$sparse.comp)
        lr_mats[[i]] <- t(fit$lr.comp)
        est_phi[[i]] <- sp_mats[[i]] + lr_mats[[i]]
    }
    # S3 class
    final_result <- list(cp = final.pts, sparse_mats = sp_mats, lowrank_mats = lr_mats, est_phi = est_phi)
    class(final_result) <- append(class(final_result), "LpSDetected")
    return(final_result)
}




#' Single change point detection for low-rank plus sparse model structure
#' @description Single change point detection for low-rank plus sparse model using minimizing SSE method.
#' @param data A n by p dataset, n is the number of observations, p is the number of variables
#' @param lambda A 2-d vector recording tuning parameters for sparse components left/right
#' @param mu A 2-d vector recording tuning parameters for low rank components left/right
#' @param alpha_L A numeric value, a positive number indicating the constraint space for low rank matrix 
#' @param skip The number of observations at the boundaries we should skip 
#' @return A list object, which includes the followings:
#' \describe{
#'   \item{cp}{A vector recording all estimated change points}
#'   \item{S_hat1}{Estimated sparse component for the left-handed side}
#'   \item{S_hat2}{Estimated sparse component for the right-handed side}
#'   \item{L_hat1}{Estimated low rank component for the left-handed side}
#'   \item{L_hat2}{Estimated low rank component for the right-handed side}
#'   \item{sse}{The sum of squared errors for all passed time points}
#' }
detect.LpS <- function(data, lambda, mu, alpha_L, skip = 50){
    sse <- c()
    n <- dim(data)[1]
    k <- dim(data)[2]
    
    lr_est <- vector("list", n - 2*skip)
    sp_est <- vector('list', n - 2*skip)
    for(t in (skip+1):(n-skip)){
        ###### segmentation of dataset ######
        seg1 <- data[1:t,]
        seg2 <- data[((t+1):n),]
        
        ### fitting the left-side dataset
        X1 <- seg1[1:(t-1),]
        Y1 <- seg1[2:t,]
        fit1 <- fista.LpS(X1, Y1, lambda = lambda[1], mu = mu[1], alpha_L = alpha_L, 
                          niter = 100, backtracking = TRUE, x.true = diag(k))
        
        ### fitting the right-side dataset
        X2 <- seg2[1:(n-t-1),]
        Y2 <- seg2[2:(n-t),]
        fit2 <- fista.LpS(X2, Y2, lambda = lambda[2], mu = mu[2], alpha_L = alpha_L, 
                          niter = 100, backtracking = TRUE, x.true = diag(k))
        
        ###### find the residual and full SSE ######
        phi.hat1 <- t(fit1$lr.comp + fit1$sparse.comp)
        phi.hat2 <- t(fit2$lr.comp + fit2$sparse.comp)
        
        lr_est[[t-skip]] <- cbind(t(fit1$lr.comp), t(fit2$lr.comp))
        sp_est[[t-skip]] <- cbind(t(fit1$sparse.comp), t(fit2$sparse.comp))
        
        pred.err1 <- (1/t) * norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        pred.err2 <- (1/(n-t)) * norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        sse <- c(sse, pred.err1 + pred.err2)
    }
    
    idx <- which.min(sse) + skip
    
    L_hat1 <- lr_est[[idx-skip]][,c(1:k)]; L_hat2 <- lr_est[[idx-skip]][,-(c(1:k))]
    S_hat1 <- sp_est[[idx-skip]][,c(1:k)]; S_hat2 <- sp_est[[idx-skip]][,-(c(1:k))]
    return(list(cp = idx, S_hat1 = S_hat1, S_hat2 = S_hat2, L_hat1 = L_hat1, L_hat2 = L_hat2, sse = sse))
}



#' Single change point detection for low-rank plus sparse model with cross-validation
#' @description single change point detection with cross-validation
#' @param data a n by p datatset matrix
#' @param lambda.1.seq a numeric vector, the sequence of tuning parameter for sparse component for the left segment
#' @param mu.1.seq a numeric vector, the sequence of tuning parameter for low rank component for the left segment
#' @param alpha_L a numeric value, indicates the constraint space of low rank component
#' @param skip a positive integer, indicates the number of time points should be skipped close to the boundary
#' @param nfold a positive integer, the number of folds for cross validation
#' @return A list object, including
#' \describe{
#'   \item{cp}{Estimated change point}
#'   \item{S_hat1}{Estimated sparse component for the left-handed side}
#'   \item{S_hat2}{Estimated sparse component for the right-handed side}
#'   \item{L_hat1}{Estimated low rank component for the left-handed side}
#'   \item{L_hat2}{Estimated low rank component for the right-handed side}
#'   \item{sse}{The sum of squared errors for all passed time points}
#' }
cv.detect.LpS <- function(data, lambda.1.seq = NULL, mu.1.seq = NULL, 
                          alpha_L = 0.25, skip = 50, nfold = 5){
    sse <- c()
    n <- dim(data)[1]
    k <- dim(data)[2]

    if(is.null(lambda.1.seq)){
        lambda.1.seq <- sqrt(log(k)/n) * exp(seq(-2, 5, length.out = nfold))
    }
    if(is.null(mu.1.seq)){
        mu.1.seq <- sqrt((k)/n) * exp(seq(-2, 5, length.out = nfold))
    }
    
    lr_est <- vector("list", n-2*skip)
    sp_est <- vector('list', n-2*skip)
    for(t in (skip+1):(n-skip)){
        ###### segmentation of dataset ######
        seg1 <- data[1:t,]
        seg2 <- data[((t+1):n),]
        
        ###### tuning parameter selection for the left-side segment ######
        ret <- cv.tuning.selection(seg1, lambda.1.seq, mu.1.seq, 
                                   alpha_L = alpha_L, nfold = nfold)
        lambda1 <- ret$lambda
        mu1 <- ret$mu
        
        ### fitting the left-side dataset
        X1 <- seg1[1:(t-1),]
        Y1 <- seg1[2:t,]
        fit1 <- fista.LpS(X1, Y1, lambda = lambda1, mu = mu1, alpha_L = alpha_L, 
                          niter = 100, backtracking = TRUE, x.true = diag(k))
        
        ###### tuning parameter selection for the right-side segment ######
        ret <- cv.tuning.selection(seg2, lambda.1.seq, mu.1.seq, 
                                   alpha_L = alpha_L, nfold = nfold)
        lambda2 <- ret$lambda
        mu2 <- ret$mu
        
        # fitting the right-side dataset
        X2 <- seg2[1:(n-t-1),]
        Y2 <- seg2[2:(n-t),]
        fit2 <- fista.LpS(X2, Y2, lambda = lambda2, mu = mu2, alpha_L = alpha_L, 
                          niter = 100, backtracking = TRUE, x.true = diag(k))
        
        ###### find the residual and full SSE ######
        phi.hat1 <- t(fit1$lr.comp + fit1$sparse.comp)
        phi.hat2 <- t(fit2$lr.comp + fit2$sparse.comp)
        
        lr_est[[t-skip]] <- cbind(t(fit1$lr.comp), t(fit2$lr.comp))
        sp_est[[t-skip]] <- cbind(t(fit1$sparse.comp), t(fit2$sparse.comp))
        
        pred.err1 <- (1/n) * norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        pred.err2 <- (1/n) * norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        sse <- c(sse, pred.err1 + pred.err2)
    }
    
    idx <- which.min(sse) + skip
    L_hat1 <- lr_est[[idx-skip]][,c(1:k)]; L_hat2 <- lr_est[[idx-skip]][,-(c(1:k))]
    S_hat1 <- sp_est[[idx-skip]][,c(1:k)]; S_hat2 <- sp_est[[idx-skip]][,-(c(1:k))]
    return(list(cp = idx, S_hat1 = S_hat1, S_hat2 = S_hat2, L_hat1 = L_hat1, L_hat2 = L_hat2, sse = sse))
}


########################################################################################
################## Cross validation related functions ##################################
########################################################################################
#' cross-validation index function, separate train and test sets
#' @description a function to setup training and testing datasets, select data point in test every nfold time points
#' @param period an input dataset with size of n by p
#' @param nfold an positive integer, indicates the number of folds in cv
#' @return a list object, includes
#' \describe{
#'     \item{train}{training dataset}
#'     \item{test}{test dataset}
#' }
cv.separate <- function(period, nfold = 5){
    n <- dim(period)[1]
    if(n < 10){
        # print("Warning: too few data points!")
        # break
        stop("Warning: too few data points!")
    }else{
        cv.index <- c()
        for(i in 1:n){
            if(i %% nfold == (nfold-1)){
                cv.index <- c(cv.index, i)
            }
        }
    }
    data_train <- period[-cv.index,]
    data_test <- period[cv.index,]
    return(list(train = data_train, test = data_test))
}


#' a function to apply cross-validation to select tuning parameter by minimizing SSE
#' @param data a n by p dataset matrix
#' @param lambda.seq a numeric vector, indicates the sequence of tuning parameters of sparse components
#' @param mu.seq a numeric vector, the sequence of tuning parameters of low rank components
#' @param alpha_L a positive numeric value, indicating the constraint space of low rank components
#' @param nfold a positive integer, the number of folds for cv
#' @return a list of object, including
#' \describe{
#'     \item{grid}{the grid of lamdbas and mus}
#'     \item{lambda}{final selected tuning parameter for sparse}
#'     \item{mu}{final selected tuning parameter for low rank}
#' }
cv.tuning.selection <- function(data, lambda.seq, mu.seq, alpha_L = 0.25, nfold = 5){
    ### separate data time series into training and testing set
    separate <- cv.separate(data, nfold = nfold)
    data_train <- separate$train
    data_test <- separate$test
    
    ntrain <- dim(data_train)[1]
    ntest <- dim(data_test)[1]
    
    X_train <- data_train[1:(ntrain-1),]
    Y_train <- data_train[2:ntrain,]
    
    X_test <- data_test[1:(ntest-1),]
    Y_test <- data_test[2:ntest,]
    
    ### grid search part: find the tuning parameter that minimizes the prediction error on test set
    n <- dim(data)[1]
    k <- dim(data)[2]
    grid <- matrix(0, nrow = length(lambda.seq), ncol = length(mu.seq))
    for(r in 1:length(lambda.seq)){
        for(c in 1:length(mu.seq)){
            lambda <- lambda.seq[r]
            mu <- mu.seq[c]
            fit <- fista.LpS(X_train, Y_train, lambda = lambda, mu = mu, alpha_L = alpha_L, 
                             niter = 100, backtracking = TRUE, diag(k))
            if(qr(fit$lr.comp)$rank == 0 | qr(fit$lr.comp)$rank > k/2){
                grid[r,c] <- Inf
            }else{
                residual <- Y_test - X_test %*% (fit$sparse.comp + fit$lr.comp)
                grid[r,c] <- norm(residual, "F")^2
            }
        }
    }
    idx <- which(grid == min(grid), arr.ind = TRUE)
    final.lambda <- lambda.seq[idx[1]]
    final.mu <- mu.seq[idx[2]]
    return(list(grid = grid, lambda = final.lambda, mu = final.mu))
}



#' Auxiliary function to calculate loss at the estimated change points
#' @description Function to calculate loss function at the estimated change points
#' @param data A n by p dataset, n is the number of observations, p is the number of variables
#' @param pts a vector of change points
#' @param lambda tuning parameter for sparse component
#' @param mu tuning parameter for low-rank component
#' @param alpha_L a positive number, indicating the value of constraint space of low rank component, default is 0.25
#' @return A list object, including
#' \describe{
#'  \item{L.n}{Value of objective function}
#' }
break.var.lps <- function(data, pts, lambda, mu, alpha_L = 0.25){
    n <- dim(data)[1]
    k <- dim(data)[2]
    m <- length(pts)
    L.n <- rep(0, m+1)
    if(m == 0){
        pts.temp <- c(1, n+1)
    }else{
        pts.temp <- c(1, pts, n+1)
    }
    for(mm in 1:(m+1)){
        data.temp <- data[(pts.temp[mm]):(pts.temp[mm+1]-1),]
        n.temp <- dim(data.temp)[1]
        x.temp <- data.temp[1:(n.temp-1),]
        y.temp <- data.temp[2:n.temp,]
        
        ### estimate the coefficient matrices and calculate the corresponding prediction error
        try <- fista.LpS(x.temp, y.temp, lambda = lambda, mu = mu, alpha_L = alpha_L, 
                         niter = 20, backtracking = TRUE, diag(k))
        est.coef <- t(try$lr.comp) + t(try$sparse.comp)
        pred.error <- y.temp - x.temp %*% t(est.coef)
        L.n[mm] <- sum(pred.error^2)
    }
    return(list(L.n = sum(L.n)))
}


#' Backward selection function for the second screening step
#' @description Backward elimination algorithm
#' @param data A n by p dataset matrix
#' @param pts A numeric vector, which includes all candidate change points obtained from the first step
#' @param lambda Tuning parameter for sparse component
#' @param mu Tuning parameter for low-rank component
#' @param alpha_L Constraint space for low rank component, default is 0.25
#' @return A list object, containing
#' \describe{
#'     \item{L.n}{Value of objective function}
#'     \item{L.n.current}{Current value of objective function}
#' }
backward.selection <- function(data, pts, lambda, mu, alpha_L = 0.25){
    n <- dim(data)[1]
    k <- dim(data)[2]
    m <- length(pts)
    L.n <- rep(0, m)
    L.n.current <- rep(0, 1)
    
    try <- break.var.lps(data, pts, lambda, mu)
    L.n.current <- try$L.n
    for(mm in 1:m){
        pts.temp <- pts[-mm]
        try <- break.var.lps(data, pts.temp, lambda, mu, alpha_L)
        L.n[mm] <- try$L.n
    }
    return(list(L.n = L.n, L.n.current = L.n.current))
}




##########################################################################################
################# First step: rolling window method to select candidates #################
##########################################################################################
#' First step rolling window function
#' @description Rolling window scheme function for the first step
#' @param data the whole data matrix
#' @param h window size
#' @param step.size rolling step size, default is NULL. If Null, the step size is 1/4 of the window size
#' @param lambda a 2-d vector of tuning parameters for sparse components, available when cv is FALSE
#' @param mu a 2-d vector of tuning parameters for low rank components, available when cv is FALSE
#' @param alpha_L a numeric value, indicates the size of constraint space of low rank component
#' @param skip the number of observations we should skip near the boundaries, default is 3
#' @param lambda.1.seq the sequence of sparse tuning parameter to the left segment, only available when cv is TRUE
#' @param mu.1.seq the sequence of low rank tuning, only available for cv is TRUE
#' @param cv a boolean argument, indicates whether use cross validation or not
#' @param nfold a positive integer, indicates the number of folds of cross validation
#' @param verbose if TRUE, then all information for current stage are printed
#' @return A vector which includes all candidate change points selected by rolling window
first.step.detect <- function(data, h, step.size = NULL, lambda, mu, alpha_L = 0.25, 
                              skip = 3, lambda.1.seq = NULL, mu.1.seq = NULL,
                              cv = FALSE, nfold = NULL, verbose = FALSE){
    n <- dim(data)[1]
    k <- dim(data)[2]
    
    if(is.null(step.size)){
        step.size <- floor(0.25 * h)
    }
    s <- 1
    e <- s + h
    candi_cp <- c()
    while(e <= n-1){
        interval_data <- data[s:e, ]
        
        ### use single change point detection method to find the candidates
        if(!cv){
            fit <- detect.LpS(interval_data, lambda = lambda, mu = mu, alpha_L = alpha_L, skip = skip)
            current_cp <- fit$cp
            candi_cp <- c(candi_cp, current_cp + s)
        }else{
            fit <- cv.detect.LpS(interval_data, lambda.1.seq = lambda.1.seq, mu.1.seq = mu.1.seq, 
                                 alpha_L = alpha_L, skip = skip, nfold = nfold)
            current_cp <- fit$cp
        }
        s <- s + step.size
        if(s + h <= n){
            e <- s + h
        }else{
            e <- n
        }
        if(verbose){
            print(paste("Finished interval:", s, e))
        }
    }
    first.step.cp <- candi_cp[-1]
    first.step.cp <- first.step.cp[order(first.step.cp)]
    
    ### remove repeat candidate change points
    m <- length(first.step.cp)
    pts <- c()
    for(mm in 1:(m-1)){
        if(abs(first.step.cp[mm] - first.step.cp[mm+1]) > 2){
            pts <- c(pts, first.step.cp[mm])
        }
    }
    pts <- c(pts, first.step.cp[m])  ## final points after the first rolling window selection.
    return(pts)
}

##########################################################################################
################ Second step: screening with backward selection algorithm ################
##########################################################################################
#' Backward elimination algorithm for screening in the second step
#' @description Backward elimination algorithm function for screening
#' @param data a n by p dataset matrix
#' @param pts a vector includes all candidate change points obtained by the first step
#' @param omega tuning parameter for the information criterion function
#' @param lambda tuning parameter for sparse component estimation
#' @param mu tuning parameter for low rank component estimation
#' @param alpha_L a numeric value, indicates the size of constraint space of low rank component
#' @param verbose if TRUE, then it provides all information for current stage
#' @return A list object includes
#' \describe{
#'     \item{pts}{Final selected change points}
#'     \item{ic}{Values of information criterion}
#' }
second.step.detect <- function(data, pts, omega, lambda, mu, alpha_L = 0.25, verbose = FALSE){
    m <- length(pts)
    if(m == 0){
        stop('no change point!')
    }
    mm <- 0; ic <- 0
    while(mm < m){
        mm <- mm + 1
        try <- backward.selection(data, pts, lambda, mu, alpha_L)
        L.n <- try$L.n; L.n.curr <- try$L.n.current
        if(min(L.n) + (m-1)*omega >= L.n.curr + m*omega){
            ic <- L.n.curr + (m - mm + 1) * omega
            break
        }else{
            pts <- pts[-which(L.n == min(L.n))]
            if(verbose){
                print(pts)
            }
        }
    }
    return(list(pts = pts, ic = ic))
}


########################################################################################
######################## FISTA estimation functions ####################################
########################################################################################
#' A function to solve low rank plus sparse model estimation using FISTA algorithm
#' @description A function to solve low rank plus sparse model estimation
#' @param A A design matrix with size of n by p
#' @param b A matrix, (or vector) with size of n by p (or n by 1)
#' @param lambda A positive numeric value, indicating the tuning parameter for sparse component
#' @param mu A positive numeric value, indicating the tuning parameter for low rank component
#' @param alpha_L The constraint coefficient of low rank component, default is 0.25
#' @param niter The maximum number of iterations required for FISTA 
#' @param backtracking A boolean argument, indicating that use backtracking in the FISTA
#' @param x.true A p by p matrix, the true model parameter. Only available for simulation.
#' @return A list object, including the followings
#' \describe{
#'   \item{sparse.comp}{Estimated sparse component}
#'   \item{lr.comp}{Estimated low-rank component}
#'   \item{obj.val}{Values of objective function}
#'   \item{rel.err}{Relative errors compared with the true model parameters if available}
#' }
#'
fista.LpS <- function(A, b, lambda, mu, alpha_L = 0.25, 
                      niter = 100, backtracking = TRUE, x.true){
    tnew = t <- 1
    p <- dim(A)[2]
    x1 <- matrix(0, nrow = p, ncol = p)
    xnew1 = xnew2 <- x1
    y1 = y2 <- x1
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b
    
    if(backtracking == TRUE){
        L <- norm(A, "F")^2 / 5
        gamma <- 2
    }else{
        L <- norm(A, "F")^2
    }
    
    obj.val = rel.err <- c()
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            found <- FALSE
            while(found == FALSE){
                y <- y1 + y2
                prox1 <- prox.sparse.func(y1, y, A, b, 2*L.bar, lambda, AtA, Atb)
                prox2 <- prox.nuclear.func(y2, y, A, b, 2*L.bar, mu, AtA, Atb)
                
                ### Restricted solution space
                for(j in 1:p){
                    for(k in 1:p){
                        if(abs(prox2[j,k]) > alpha_L){
                            prox2[j,k] <- alpha_L * sign(prox2[j,k])
                        }
                    }
                }
                
                prox <- prox1 + prox2
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    found <- TRUE
                }else{
                    L.bar <- L.bar * gamma
                }
            }
            L <- L.bar
        }
        x1 <- xnew1 
        x2 <- xnew2
        xnew1 <- prox1
        xnew2 <- prox2
        t = tnew
        tnew <- (1 + sqrt(1 + 4*t^2))/2
        y1 <- xnew1 + (t - 1) / tnew * (xnew1 - x1)
        y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
        xnew <- xnew1 + xnew2
        
        obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
        rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
    }
    return(list(sparse.comp = xnew1, lr.comp = xnew2, obj.val = obj.val, rel.err = rel.err))
}



#' Shrinkage function for sparse soft-thresholding
#' @description Shrinkage function for sparse soft-thresholding
#' @param y A matrix, or a vector for thresholding
#' @param tau A positive number, threshold
#' @return A thresholded matrix, or vector
shrinkage <- function(y, tau){
    z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for(i in 1:nrow(y)){
        for(j in 1:ncol(y)){
            z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
        }
    }
    return(z)
}

#' Shrinkage function for low-rank soft-thresholding
#' @description Shrinkage function for low-rank soft-thresholding
#' @param y A matrix, or a vector for thresholding
#' @param tau A positive number, threshold
#' @return A thresholded matrix, or vector
shrinkage.lr <- function(y, tau){
    z <- rep(0, length(y))
    for(i in 1:length(y)){
        z[i] <- sign(y[i]) * max(0, abs(y[i]) - tau)
    }
    return(z)
}


#' Gradient function of quardratic loss
#' @description Gradient function of quardratic loss
#' @param x A vector, or matrix, indicating the model parameter
#' @param AtA A p by p Gram matrix for corresponding design matrix A
#' @param Atb An inner product for design matrix A and corresponding matrix (vector) b
#' @return Value of gradients
gradf.func <- function(x, AtA, Atb){
    return(AtA %*% x - Atb)
}


#' Nuclear norm penalty for low-rank component
#' @description Nuclear norm penalty for low-rank component
#' @param x Model parameter
#' @param lambda Tuning parameter
#' @return Value of nuclear norm penalty term
nuclear.pen <- function(x, lambda){
    d <- svd(x)$d
    return(lambda * sum(d))
}


#' L1-norm penalty for sparse component
#' @description L1-norm penalty for sparse component 
#' @param x Model parameter
#' @param lambda Tuning parameter
#' @return Value of l1-norm penalty term
sparse.pen <- function(x, lambda){
    return(lambda*sum(x))
}


#' Main loss function for quardratic loss
#' @description Main loss function
#' @param x Model parameters
#' @param A Design matrix with size of n by p
#' @param b Correspond vector or matrix
#' @return Value of objective function
f.func <- function(x, A, b){
    return(0.5 * norm(A %*% x - b, "F")^2)
}


#' An auxiliary function in FISTA algorithm
#' @description Auxiliary function for FISTA implementation
#' @param x Model parameter for previous update
#' @param y Model parameter for updating
#' @param A An n by p design matrix
#' @param b A correspond vector, or matrix with size of n by 1 or n by p
#' @param L Learning rate
#' @param AtA Gram matrix for design matrix A
#' @param Atb Inner product for design matrix A and correspond vector b
#' @return Value of function Q
Q.func <- function(x, y, A, b, L, AtA, Atb){
    return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + 0.5 * L * norm(x - y, "F")^2)
}


#' Proximal function with nuclear norm penalty updating
#' @description Proximal function with nuclear norm
#' @param w1 previously updated model parameter
#' @param y updated model parameter
#' @param A design matrix
#' @param b correspond vector, or matrix
#' @param L learning rate
#' @param lambda tuning parameter for low-rank component
#' @param AtA Gram matrix of design matrix A
#' @param Atb inner product of design matrix A and correspond vector b
#' @return Value of proximal function with nuclear norm penalty
prox.nuclear.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    d <- shrinkage.lr(svd(Y)$d, 2*lambda / L)
    return(svd(Y)$u %*% diag(d) %*% t(svd(Y)$v))
}


#' Proximal function with l1-norm penalty updating
#' @description Proximal function with l1-norm
#' @param w1 previously updated model parameter
#' @param y updated model parameter
#' @param A design matrix
#' @param b correspond vector, or matrix
#' @param L learning rate
#' @param lambda tuning parameter for sparse component
#' @param AtA Gram matrix of design matrix A
#' @param Atb inner product of design matrix A and correspond vector b
#' @return Value of proximal function with l1-norm penalty
prox.sparse.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    return(shrinkage(Y, 2*lambda / L))
}



#' Objective function
#' @description objective function, main loss function and penalties
#' @param x.lr low-rank component
#' @param x.sparse sparse component
#' @param A design matrix 
#' @param b correspond vector
#' @param lambda a tuning parameter for sparse component
#' @param mu a tuning parameter for low-rank component
#' @return value of objective function
obj.func <- function(x.lr, x.sparse, A, b, lambda, mu){
    ### x.sparse is a list
    m <- length(x.sparse)
    loss <- 0
    for(i in 1:m){
        loss <- loss + f.func((x.lr[[i]] + x.sparse[[i]]), A, b) + sparse.pen(x.sparse[[i]], lambda) + nuclear.pen(x.lr[[i]], mu)
    }
    return(loss)
}