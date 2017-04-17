#' derive a data.frame of all possible combinations of vars=
mCombinations <- function(vars=NULL,verbose=T){
  if(verbose) cat(" -- deriving model combinations:")
  # calculate : combinations w/o repetition (n!/(r!(n-r)!)... or 2^n
  m_len <- vector(); for(i in 1:length(vars)) { m_len <- append(m_len,dim(combn(vars,m=i))[2]) }
    m_len <- sum(m_len)
  # define the model space of all possible combinations of predictors
  models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
  k <- 1 # row of our models data.frame
  for(i in 1:length(vars)){
   combinations <- combn(vars,m=i)
   # build a formula string with lapply comprehension
   f <- function(j){
     paste("~doy~",paste(combinations[,j],collapse="+"),collapse="")
   }
   # lapply over all combinations of our current n covariates
   # (avoiding a slow, nested for-loop)
   models[k:(k+ncol(combinations)-1),1] <-
     unlist(lapply(1:ncol(combinations),FUN=f))
   k <- k+ncol(combinations)
   if(verbose) cat(".");
  };
  if(verbose) cat("\n");
  return(models)
}
#' perform a random walk on an unmarked dataframe with a user-specified unmakred funtion
#' @export'
randomWalk_dAIC <- function(vars=NULL, step=1000, umdf=NULL,
                            umFunction=unmarked::distsamp, nCores=NULL,retAll=FALSE){
  require(parallel)
  if(!require(unmarked)){ stop("function requires the unmarked package is installed") }
  nCores <- ifelse(is.null(nCores), parallel::detectCores()-1, nCores)
      cl <- parallel::makeCluster(nCores)
  models <- OpenIMBCR:::mCombinations(vars)
  # parallelize our runs across nCores processors (defined at top)
  total_runs <- 1:nrow(models)
  # prime the pump
  cat(" -- starting a random walk:\n")
  # assign a null model AIC to beat (below)
  m <- distsamp(~doy~1, umdf,keyfun="hazard",output="density",unitsOut="kmsq") # use AIC to determine an optimal detection function
  minimum <- data.frame(formula="~doy~1",AIC=m@AIC) # begin with our null (intercept) model
  # iterate over total_runs and try and minimize AIC as you go
  while(length(total_runs)>1){
    # randomly sample total_runs that the cluster will consider for this run
    focal_runs <- sample(total_runs,
                         replace=F,
                         size=ifelse(length(total_runs) > step, step, length(total_runs))
                         )
    # build models for this run
    runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
      runs <- parLapply(cl=cl, runs, fun=umFunction, data=umdf,keyfun="hazard", output="density", unitsOut="kmsq") # this should change based on user-specified function
        runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
    # if we beat the running lowest AIC, append it to the random walk table
    if(runs[which(runs == min(runs))[1]] < min(minimum$AIC)){
      minimum <- rbind( minimum,
                        data.frame( formula=models[focal_runs[which(runs == min(runs))[1]],'formula'],
                                        AIC=runs[which(runs == min(runs))[1]] )
                      )
    }
    total_runs <- total_runs[!(total_runs %in% focal_runs)]
    cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
  };
  cat("\n");
  parallel::stopCluster(cl)
  if(retAll){
    return(
        data.frame( formula=models[focal_runs,'formula'],
                    AIC=runs )
          )
  } else {
    return(minimum)
  }
}
#' use a globally-sensitive numerical optimization procedure to select covariates for
#' inclusion in our model. This should be considerably faster than randomWalk_dAIC()
simulatedAnnealing_dAIC <- function(m){
  return(NA)
}
