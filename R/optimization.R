#' stump code that will perform a random walk on an unmarked model object.
#' Needs to be clean-up and validated.
#' @export'
randomWalkUnmarked <- function(m){
  if(!require(unmarked)){ stop("function requires the unmarked package is installed") }
  # calculate : combinations w/o repetition (n!/(r!(n-r)!)... or 2^n
  m_len <- vector(); for(i in 1:length(n)) { m_len <- append(m_len,dim(combn(n,m=i))[2]) }
   m_len <- sum(m_len)
  # define the model space of all possible combinations of predictors
  models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
  k <- 1
  cat(" -- deriving model combinations:")
  for(i in 1:length(n)){
   combinations <- combn(n,m=i)
   for(j in 1:ncol(combinations)){
     models[k,1] <- paste("~doy~",paste(combinations[,j],collapse="+"),collapse="")
     k <- k+1
   }; cat(".");
  }; cat("\n");
  # parallelize our runs across nCores processors (defined at top)
  step <- 1000
  total_runs <- 1:nrow(models)
  # prime the pump
  cat(" -- starting a random walk:\n")
  # assign a null model AIC to beat (below)
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
      runs <- parLapply(cl=cl, runs, fun=distsamp, data=umf,keyfun="hazard", output="density", unitsOut="kmsq")
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
  }; cat("\n");
}
