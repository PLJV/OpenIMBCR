#' build a series of response plots for a fitted model and input table
partialPredict <- function(m=NULL, var=NULL, type='state', plot=T, nCores=NULL,
                           xlab=NULL, ylab=NULL, main=NULL, xTransform=NULL){
  require(parallel)
  nCores <- ifelse(is.null(nCores), parallel::detectCores()-1, nCores)
      cl <- parallel::makeCluster(nCores)
  # fetch our site-level covariates from the training data
  xlab <- ifelse(is.null(xlab), var, xlab)
  ylab <- ifelse(is.null(ylab), "Density (birds/km2)", ylab)
  t <- m@data@siteCovs
  # run our model for each unique value of x, averaging the predictions across all non-focal variables as we go
  x <- seq((min(t[,var])-sd(t[,var])),(max(t[,var])+sd(t[,var])),length.out=100)
  y <- matrix(NA, nrow=length(x), ncol=4)
  partial <- function(x){
    require(unmarked)
    run_table <- t
      run_table[, var] <- x
    p <- predict(m, run_table, type = type)
      y <- colMeans(p)
  }
  y <- do.call(rbind, parLapply(cl, x, fun=partial))
  colnames(y) <- colnames(p)
  parallel::stopCluster(cl)
  if(!is.null(xTransform)){
    x <- eval(parse(text = xTransform))
  }
  if(plot){
    dev.new()
    plot(y=y[,1], x=x, xlab=xlab, ylab=ylab, main=main, col="white")
    grid(lwd=1.4)
    lines(y=y[,1], x=x, col="black")
      lines(y=y[,4], x=x, col="grey")
        lines(y=y[,3], x=x, col="grey")
  } else {
    return(list(x=x, y=y))
  }
}
