# Stochastic Point Method
SI.SPM <- function(
    h,       # density h(x) to be integrated
    from,    # left end point of x
    to,      # right end point of x
    M,       # the upper bound of h(x) in [from,to]
    N        # the number of points to be generated
){
    x <- runif(N, min = from, max = to)
    y <- runif(N, min = 0, max = M)
    hx <- h(x)
    p_hat <- mean(y<=hx)
    I <- p_hat*M*(to-from)
    Var <- I*(M*(to-from)-I)/N
    return(list(I,Var))
}


# Mean Value Method
SI.MVM <- function(
    h,       # density h(x) to be integrated
    from,    # left end point of x
    to,      # right end point of x
    N        # the number of points to be generated
){
    x <- runif(N, min = from, max = to)
    hx <- h(x)
    I <- (to-from)* mean(hx)
    Var <- mean(((to-from)*hx-I)^2)/N
    return(list(I,Var))
}


# Important Sampling Method
SI.ISM <- function(
    h,       # density h(x) to be integrated
    g,       # test density g(x)
    G_inv,   # inverse function of G(x)
    N,       # the number of points to be generated
    min_G = 0,
    max_G = 1
){
    U <- runif(N,min_G,max_G)
    X <- G_inv(U)
    I <- mean(h(X)/g(X)*(max_G-min_G))
    Var <- mean((h(X)/g(X)*(max_G-min_G)-I)^2)/N
    return(list(I,Var))
}

# Stratified Sampling Method
SI.SSM <- function(
    h,       # density h(x) to be integrated
    from,    # left end point of x
    to,      # right end point of x
    level,
    N        # the number of points to be generated
){
    interval <- seq(from,to,length.out=level+1)
    I <- 0
    Var <- 0
    for (i in c(1:level)){
        MVMresult <- SI.MVM(h,interval[i],interval[i+1],as.integer(N/level))
        I <- I + MVMresult[[1]]
        Var <- Var + MVMresult[[2]]
    }
    return(list(I,Var))
}

