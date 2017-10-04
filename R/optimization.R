#' downsample records in a continuous 'distances' distribution using a
#' quantile cut-off. Distances can be literal, or the output of a Mahalanobis or Euclidean distance function (see: the 'FNN'
#' package). This function is in testing.
#' @export
downsample_by_quantile <- function(distances=NULL,
                                   p=0.1,
                                   lte=F,
                                   gte=F,
                                   byid=F){

  if(lte)
    ret <- distances < quantile(distances,p=p)
  if(gte)
    ret <- distances > quantile(distances,p=p)
  if(byid){
    return(which(ret))
  }
  return( distances[ ret ] )
}
#' downsample records in a continuous 'distances' distribution using a
#' normal distribution function (can be truncated). Distances can be literal,
#' or the output of a Mahalanobis or Euclidean distance function (see: the 'FNN'
#' package). This function is in testing. It behaves oddly sometimes and needs to be refactored.
#' Use with caution.
#'
#' @param shape is multiplier applied to the SD of the distances
#' vector. Multipliers from 1 -> 0 will restrict the variance (increase the
#' shoulder character of our output).
#' @export
downsample_by_normal_dist <- function(distances=NULL, bins=11,
                                      shape=1,
                                      byid=F,
                                      calc_median_bin=F,
                                      use_mean=F){
  counts <- hist(distances,breaks=bins,plot=F)$counts
  # pad our counts by one class so they agree with the length of counts
  counts <- append(counts,0)
  breaks <- hist(distances,breaks=bins,plot=F)$breaks

  # should we calculate a CT consistent with our bins?
  if(calc_median_bin){
    central_tendency <- hist(distances,
        breaks=10,plot=F)$density
    central_tendency <- breaks[ which(central_tendency ==
        max(central_tendency)) ]
  # or just calculate the mean?
  } else if(use_mean){
    central_tendency <- mean(distances)
  # by default, assume the shoulder is the central tendency
  } else {
    central_tendency <- breaks[1]
  }

  probs <- sapply(breaks, FUN=function(x){
      dnorm(x=x,mean=central_tendency,
        sd=sd(distances)*shape
      )
  })
  # min/max normalize our probabilities to sampling densities --
  # ignore warnings about size of counts. The last count is often
  # junk (0).
  strata_densities <- ceiling(round(
        probs - min(probs),20) / diff(range(probs)
      ) * counts )
  # iterate over our breaks, downsampling accordingly
  keep <- vector()
  for(i in 1:(length(breaks)-1)){
    if(strata_densities[i] != 0){
      sample_bin <- distances[distances < breaks[i+1] &
        distances >= breaks[i]]
      miss <- length(sample_bin)-strata_densities[i]
      if(miss < 0){
        warning(paste(miss,"fewer samples in bin",
          "than requested for this strata -- this",
          "shouldn't happen. Adjusting..."))
        strata_densities[i] = length(sample_bin)
      }
      keep <- append(keep,
        which(distances %in%
          sample(sample_bin, size=strata_densities[i])
        )[1:strata_densities[i]]
      )
    }
  }
  # sanity-check
  if( (1 - length(unique(na.omit(keep)))/sum(counts)) > 0.4){
    warning(paste(
        "we lost",
        100*( 1 - length(unique(na.omit(keep))) / sum(counts) ),
        "% of the records in downsampling"
      ))
  }
  if(byid){
    ret <- unique(na.omit(keep))
  } else {
    ret <- distances[unique(na.omit(keep))]
  }
  return(ret)
}
#' hidden function to derive a data.frame of all possible combinations of vars=
mCombinations <- function(siteCovs=NULL,availCovs=NULL,detCovs=NULL,verbose=T){
  if(verbose) cat(" -- deriving model combinations:")
  # calculate : combinations w/o repetition (n!/(r!(n-r)!)... or 2^n
  m_len <- vector(); for(i in 1:length(siteCovs)) { m_len <- append(m_len,dim(combn(siteCovs,m=i))[2]) }
    m_len <- sum(m_len)
  # define the model space of all possible combinations of predictors
  models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
  k <- 1 # row of our models data.frame
  for(i in 1:length(siteCovs)){
   combinations <- combn(siteCovs,m=i)
   # build a formula string with lapply comprehension
   f <- function(j){
     paste(
       paste("~",paste(combinations[,j],collapse="+"), sep=""),
       "~1",
       paste("~",paste(detCovs,collapse="+"),sep=""),
       collapse=""
     )
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
#' perform a random walk on an unmarked dataframe with a user-specified unmarked funtion
#' @export
randomWalk_dAIC <- function(siteCovs=NULL, availCovs=NULL, detCovs=NULL,
                            step=100, umdf=NULL,
                            umFunction=unmarked::distsamp,
                            nCores=NULL,retAll=FALSE, ...){
  require(parallel)
  if(!require(unmarked)){ stop("function requires the unmarked package is installed") }
  nCores <- ifelse(is.null(nCores), parallel::detectCores()-1, nCores)
      cl <- parallel::makeCluster(nCores)
  models <- OpenIMBCR:::mCombinations(
      siteCovs=siteCovs,
      availCovs=availCovs,
      detCovs=detCovs
    )
  # parallelize our runs across nCores processors (defined at top)
  total_runs <- 1:nrow(models)
  # prime the pump
  cat(" -- starting a random walk:\n")
  # assign a null model AIC to beat (below)
  m <- unmarked::gdistsamp(
      lambdaformula=~1+offset(log(effort)), # abundance
      phiformula=~1,                     # availability
      pformula=formula(paste("~",paste(detCovs, collapse="+"),sep="")), # detection
      data=umdf,
      ...
    )
  # begin with our null (intercept) model
  minimum <- data.frame(formula="~1+offset(log(effort))~1~doy+starttime",AIC=m@AIC)
  # iterate over total_runs and try and minimize AIC as you go
  while ( length(total_runs) > 1 ){
    # randomly sample total_runs that the cluster will consider for this run
    focal_runs <- sample(total_runs,
                         replace=F,
                         size=ifelse(length(total_runs) > step, step, length(total_runs))
                         )
    # build models for this run across our cluster
    if(identical(umFunction, unmarked::gdistsamp)){
      # split the formula comprehension into many arguments
      functionFactory <- function(x,data=NULL,...){
        formulas <- gsub(na.omit(unlist(strsplit(
            Reduce(paste, deparse(x)),
            split="~"))),
            pattern=" |[\n]",
            replacement=""
          )
        lambda <- paste("~",formulas[2],sep="")
        phi <- paste("~",formulas[3],sep="")
        p <- paste("~",formulas[4],sep="")
        return(unmarked::gdistsamp(
            lambdaformula=lambda,
            phiformula=phi,
            pformula=p,
            data=data,
            ...
          ))
      }
    } else if(identical(umFunction, unmarked::distsamp)) {
      functionFactory <- function(x,data=NULL,...){
        formulas <- na.omit(unlist(strsplit(
          Reduce(paste, deparse(x)),
          split="~"))
        )
        formula=as.formula(paste(
            paste("~",formulas[3],sep=""),
            paste("~",formulas[2],sep="")
          ))
        return(umFunction(formula, data=data, ...))
      }
    } else {
      functionFactory <- umFunction
    }
    runs <- lapply(as.list(models[focal_runs,1]),FUN=as.formula)
    runs <- parLapply(
        cl=cl,
        runs,
        fun=functionFactory,
        data=umdf

      )
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
#' inclusion in our model. This should be considerably faster for walking large variable
#' space than randomWalk_dAIC()
simulatedAnnealing_dAIC <- function(m){
  return(NA)
}
#' Frequentist slope intercept test first described by Bartuszevige. Does the
#' confidence interval of a given variable cross the intercept (i.e., x_n=0)?
#' returns pass/fail by default.
#' @param alpha alpha value for our test (default is 0.975)
#' @export
bartuszevige_intercept_test <- function(m=NULL, var=NULL, alpha=0.975){
  se <- OpenIMBCR::SE(m@opt)
    se <- se[1:max(which(grepl(rownames(m@opt$hessian),pattern="ntercept")))-1] # state covs
      se <- se * qnorm(alpha)
        se <- if(!is.null(var)) se[var] else se
  if(is.null(var)){
    crosses_zero <-
      matrix(
        c(
          m@estimates@estimates$state@estimates + se,
          m@estimates@estimates$state@estimates - se
        ),
        ncol=2)
    crosses_zero <- apply(crosses_zero,MARGIN=1,FUN=prod) < 0
      names(crosses_zero) <- rownames(m@opt$hessian)[1:max(which(grepl(rownames(m@opt$hessian),pattern="ntercept")))-1]
  } else {
    crosses_zero <-
    prod(range(m@estimates@estimates$state@estimates[var] + se,
         m@estimates@estimates$state@estimates[var] - se)
        ) < 0
  }
  return(!crosses_zero)
}
