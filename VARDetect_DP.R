#' Function represents the objective function Ln corresponding to the interval [s,e]
#' @param data a numeric matrix, indicating the dataset with size of n by p
#' @param s a positive integer, representing the start time point for the series
#' @param e a positive integer, representing the end time point for the series
#' @param lambda a numeric value, indicating the tuning parameter for the sparse component
#' @param mu a numeric value, indicating the tuning parameter for the low rank component
#' @return A numeric value indicating the value of objective function Ln
L.n <- function(data, s, e, lambda, mu){
    temp_data <- data[s:e, ]
    n.temp <- dim(temp_data)[1]
    p <- dim(temp_data)[2]

    X <- temp_data[1:(n.temp-1),]
    y <- temp_data[2:n.temp, ]
    try <- fista.LpS(X, y, lambda, mu, niter = 50, backtracking = TRUE, diag(p))
    est.coef <- t(try$lr.comp) + t(try$sparse.comp)
    pred.error <- y - X %*% t(est.coef)

    res <- sum(pred.error^2)
    return(res)
}

#' The main function for dynamical programming estimation
#' @param data a numeric matrix, indicating the dataset with size of n by p
#' @param s a numeric value, starting time point for the data time point
#' @param t a numeric value, indicating the time point for the series
#' @param flag a boolean value, for noting the local minima of the dp method
#' @param lambda a tuning parameter for sparse component
#' @param mu a tuning parameter for low rank component
#' @return A list, containing all estimated change points
#' \describe{
#'     \item{final.cps}{finally estimated change points}
#' }
#' @export
#' @examples
#' library("VARDetect")
#' nob <- 300
#' p <- 15
#' brk <- c(floor(nob/3), floor(2*nob/3), nob+1)
#' lambda <- 0.1; mu <- 1
#' try <- simu_var(nob = nob, k = p, lags = 1, sigma = 0.01*diag(p), brk = brk,
#'                 sp_pattern = "off-diagonal")
#' data <- as.matrix(try$series)
#' fit <- dp.detect(data, s = 1, t = 5, flag = 0, lambda = lambda, mu = mu)
#' print(fit$final.cps)
dp.detect <- function(data, s, t, flag, lambda, mu){
    ret <- c()
    n <- dim(data)[1]
    p <- dim(data)[2]

    ### penalty tuning parameter
    gamma <- log(n)*p*0.05

    ### use dynamical programming method
    while(s < n-1){
        s <- s+1
        while(t < n && flag == 0){
            t <- t+1
            L.temp <- c()
            L.total <- ifelse(t-s >= gamma, L.n(data, s, t, lambda, mu), 0)
            for(l in (s+2):(t-2)){
                ### separately consider the left, right segments for [s,l] and [l+1,t], and total [s,t]
                L.left <- ifelse(t-s >= gamma, L.n(data, s, l, lambda, mu), 0)
                L.right <- ifelse(t-s >= gamma, L.n(data, l, t, lambda, mu), 0)
                L.temp <- c(L.temp, L.left + L.right + gamma)
            }
            if(min(L.temp) < L.total){
                local_min <- which.min(L.temp) + s
                ret <- c(ret, local_min)
                s <- local_min
                flag <- 1
            }
            cat(paste("subinterval time point:", t, "\n", sep = " "))
        }
        flag <- 0
        cat(paste("current time point:", s, "\n", sep = " "))

    }
    return(list(final.cps = ret))
}
