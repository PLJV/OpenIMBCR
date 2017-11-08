# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2016
# Description : various tasks related to model selection and optimization

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' downsample records in a continuous 'distances' distribution using a
#' quantile cut-off. Distances can be literal, or the output of a Mahalanobis 
#' or Euclidean distance function (see: the 'FNN'package). This function is 
#' in testing.
#' @export
downsample_by_quantile <- function(
  distances=NULL,
  p=0.1,
  lte=F,
  gte=F,
  byid=F)
{

  if(lte)
    ret <- distances < quantile(distances,p=p)
  if(gte)
    ret <- distances > quantile(distances,p=p)
  if(byid){
    return(which(ret))
  }
  return( distances[ ret ] )
}
#' testing : downsample records in a continuous 'distances' distribution using a
#' normal distribution function (can be truncated). Distances can be literal,
#' or the output of a Mahalanobis or Euclidean distance function (see: the 'FNN'
#' package). This function is in testing. It behaves oddly sometimes and 
#' needs to be refactored. Use with caution.
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
#' testing : hidden function to derive a data.frame of all possible 
#' combinations of vars=
mCombinations <- function(siteCovs=NULL,availCovs=NULL,detCovs=NULL,
                          offset=NULL, verbose=T){
  if(verbose) cat(" -- deriving model combinations:")
  # calculate : combinations w/o repetition (n!/(r!(n-r)!)... or n^2
  m_len <- sum(unlist(lapply(
      1:length(siteCovs), 
      FUN=function(i) dim(combn(siteCovs,m=i))[2] 
    )))
  # define the model space of all possible combinations of predictors
  models <- data.frame(formula=rep(NA,m_len),AIC=rep(NA,m_len))
  k <- 1 # row of our models data.frame
  for(i in 1:length(siteCovs)){
   combinations <- combn(siteCovs,m=i)
   # build a formula string with lapply comprehension
   f <- function(j){
     paste(
       if(is.null(siteCovs)){
         "~1"
       } else {
         paste(
           paste("~",paste(combinations[,j],collapse="+"), sep=""),
           ifelse(is.null(offset),NULL,paste("+",offset,sep=""))
         )
       },
       if(is.null(availCovs)){
         "~1"
       } else {
         paste("~",paste(availCovs,collapse="+"), sep="")
       },
       if(is.null(detCovs)){
         "~1"
       } else {
         paste("~",paste(detCovs,collapse="+"),sep="")
       },
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
#' testing: perform a random walk on an unmarked dataframe with a
#' user-specified unmarked function and some depth for exploring
#' our model covariates. This allows for exploring models across
#' high-dimensional space without having to fit millions of models.
#' But it's controversial, because Kyle made it up. There are no papers
#' in the statistics literature demonstrating that this is a valid way
#' of exploring high-dimensional datasets. So, let's write one.
#' @export
randomWalk_dAIC <- function(
  siteCovs=NULL,
  availCovs=NULL,
  detCovs=NULL,
  step=100,
  umdf=NULL,
  offset=NULL,
  depth=1,
  umFunction=unmarked::distsamp,
  nCores=NULL, ...){
  # define our workspace and set-up parallel
  if(!require(unmarked)){ stop("function requires the unmarked package is installed") }
  nCores <- ifelse(
      is.null(nCores),
      parallel::detectCores()-1,
      nCores
    )
  cl <- parallel::makeCluster(nCores)
  models <- OpenIMBCR:::mCombinations(
      siteCovs=siteCovs,
      availCovs=availCovs,
      detCovs=detCovs,
      offset=offset
    )
  if(depth<1){
    sample_size <- round(nrow(models) * depth)
    models <- models[sample(1:nrow(models), size=sample_size), ]
  }
  # append a null model to the top of our models data.frame
  null_model <- gsub(paste(
    ifelse(
        !is.null(offset),
        paste(paste("~1+",offset,sep=""),"~1~"),
        "~1~1~"
      ),
      paste(detCovs,collapse="+"),
      sep=""
    ),
    pattern=" ",
    replacement=""
  )
  models <- rbind(data.frame(formula=null_model, AIC=NA), models)
  # parallelize our runs across nCores processors (defined at top)
  total_runs <- 1:nrow(models)
  cat(" -- starting a random walk:\n")
  # begin with a null (intercept) model
  # iterate over total_runs and try and minimize AIC as you go
  while ( length(total_runs) > 1 ){
    # randomly sample total_runs that the cluster will consider for this run
    focal_runs <- sample(
        total_runs,
        replace=F,
        size=ifelse(length(total_runs) > step, step, length(total_runs))
      )
    # use factory comprehension to determine function handling for
    # different unmarked model-fitting calls
    if(identical(umFunction, unmarked::gdistsamp)){
      # split the formula comprehension into many arguments
      functionFactory <- function(x,data=NULL,offset=NULL,...){
        formulas <- gsub(na.omit(unlist(strsplit(
            Reduce(paste, deparse(x)),
            split="~"))),
            pattern=" |[\n]",
            replacement=""
          )
        lambda <- paste("~",formulas[2],sep="")
        phi <- paste("~",formulas[3],sep="")
        p <- paste("~",formulas[4],sep="")
        return(tryCatch(umFunction(
              lambdaformula=lambda,
              phiformula=phi,
              pformula=p,
              data=data,
              ...
            ),
            error = function(e) NA
          ))
      }
    } else if(identical(umFunction, unmarked::distsamp)) {
      functionFactory <- function(x,data=NULL,offset=NULL,...){
        formulas <- na.omit(unlist(strsplit(
          Reduce(paste, deparse(x)),
          split="~"))
        )
        formula=as.formula(paste(
            paste("~",formulas[3],sep=""),
            paste("~",formulas[2],sep="")
          ))
        return(tryCatch(umFunction(
              formula,
              data=data,
              ...
            ),
            error = function(e) NA
          ))
      }
    } else {
      # uncaught haymaker
      functionFactory <- umFunction
    }
    # set-up our model runs by chunking the formulas table
    # across our 'focal_runs'
    runs <- lapply(as.list(as.character(
        models[focal_runs, 1])),
        FUN = as.formula
      )
    runs <- parLapply(
        cl=cl,
        runs,
        fun=functionFactory,
        data=umdf,
        ...
      )
    # drop any models that unmarked failed to fit
    runs <- unlist(lapply(runs, na.omit))
    # fetch our AIC's for those models we retained
    runs <- unlist(lapply(runs,FUN=function(x){x@AIC}))
    # if we beat the running lowest AIC, append it to the random walk table
    if(!exists("minimum")){
      minimum <- data.frame(
        formula=models[focal_runs[which(runs == min(runs))[1]],'formula'],
        AIC=runs[which(runs == min(runs))[1]]
      )
    } else if(runs[which(runs == min(runs))[1]] < min(minimum$AIC)){
      minimum <- rbind(minimum,
                       data.frame(
                         formula=models[focal_runs[which(runs == min(runs))[1]],'formula'],
                         AIC=runs[which(runs == min(runs))[1]]
                       ))
    }
    total_runs <- total_runs[!(total_runs %in% focal_runs)]
    cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
  };
  cat("\n");
  parallel::stopCluster(cl)
  return(minimum)
}
#'
#'
allCombinations_dAIC <- function(
  siteCovs=NULL,
  availCovs=NULL,
  detCovs=NULL,
  step=100,
  umdf=NULL,
  offset=NULL,
  ic=OpenIMBCR:::AICc,
  umFunction=unmarked::distsamp,
  nCores=NULL,
  ...){
  # define our workspace and set-up parallel
  if(!require(unmarked)){ stop("function requires the unmarked package is installed") }
  nCores <- ifelse(
      is.null(nCores),
      parallel::detectCores()-1,
      nCores
    )
  cl <- parallel::makeCluster(nCores)
  models <- OpenIMBCR:::mCombinations(
      siteCovs=siteCovs,
      availCovs=availCovs,
      detCovs=detCovs,
      offset=offset
    )
  # append a null model to the top of our models data.frame
  null_model <- gsub(paste(
    ifelse(
        !is.null(offset),
        paste(paste("~1+",offset,sep=""),"~1~"),
        "~1~1~"
      ),
      paste(detCovs,collapse="+"),
      sep=""
    ),
    pattern=" ",
    replacement=""
  )
  models <- rbind(data.frame(formula=null_model, AIC=NA), models)
  # parallelize our runs across nCores processors (defined at top)
  total_runs <- 1:nrow(models)
  cat(" -- building models across all covariate combinations:\n")
  # begin with a null (intercept) model
  # iterate over total_runs and try and minimize AIC as you go
  model_selection_table <- data.frame(
      formula=NA,
      AIC=NA
    )
  while ( length(total_runs) > 1 ){
    # randomly sample total_runs that the cluster will consider for this run
    focal_runs <- sample(
        total_runs,
        replace=F,
        size=ifelse(length(total_runs) > step, step, length(total_runs))
      )
    # use factory comprehension to determine function handling for
    # different unmarked model-fitting calls
    if(identical(umFunction, unmarked::gdistsamp)){
      # split the formula comprehension into many arguments
      functionFactory <- function(x,data=NULL,offset=NULL,...){
        formulas <- gsub(na.omit(unlist(strsplit(
            Reduce(paste, deparse(x)),
            split="~"))),
            pattern=" |[\n]",
            replacement=""
          )
        lambda <- paste("~",formulas[2],sep="")
        phi <- paste("~",formulas[3],sep="")
        p <- paste("~",formulas[4],sep="")
        return(tryCatch(umFunction(
              lambdaformula=lambda,
              phiformula=phi,
              pformula=p,
              data=data,
              ...
            ),
            error = function(e) NA
          ))
      }
    } else if(identical(umFunction, unmarked::distsamp)) {
      functionFactory <- function(x,data=NULL,offset=NULL,...){
        formulas <- na.omit(unlist(strsplit(
          Reduce(paste, deparse(x)),
          split="~"))
        )
        formula=as.formula(paste(
            paste("~",formulas[3],sep=""),
            paste("~",formulas[2],sep="")
          ))
        return(tryCatch(umFunction(
              formula,
              data=data,
              ...
            ),
            error = function(e) NA
          ))
      }
    } else {
      # uncaught haymaker
      functionFactory <- umFunction
    }
    # set-up our model runs by chunking the formulas table
    # across our 'focal_runs'
    runs <- lapply(as.list(as.character(
        models[focal_runs, 1])),
        FUN = as.formula
      )
    runs <- parLapply(
        cl=cl,
        runs,
        fun=functionFactory,
        data=umdf,
        ...
      )
    # drop any models that unmarked failed to fit
    keep <- !unlist(lapply(runs, is.na))
    runs <- runs[keep]
    # fetch our AIC's for those models we retained
    runs <- unlist(lapply(
        runs,
        FUN=ic
      ))
    # append our AIC's and formulas to the output table
    if(sum(keep)>0){
      model_selection_table <- rbind(
        model_selection_table,
        data.frame(
            formula=models[focal_runs * keep, 'formula'],
            AIC=runs
          )
        )
    }
    total_runs <- total_runs[!(total_runs %in% focal_runs)]
    cat(paste("[jobs remaining:",length(total_runs),"]",sep=""));
  };
  cat("\n");
  parallel::stopCluster(cl)
  return(na.omit(
      model_selection_table[order(model_selection_table[,2], decreasing=F),]
    ))
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
#' estimates a K parameter (upper-bound of integration) as a pre-cursor for
#' various count-based mixture models (e.g., poisson and negative binomial)
calc_k <- function(df=NULL, multiplier=1){
  return(floor(max(
      rowSums(unmarked:::getY(df)))*multiplier
    ))
}
#' testing: there is some ambiguity in selecting a good upper-bounds for integration
#' in gdistamp (the K parameter). This builds an intercept model around a range
#' of K values and selects the smallest K value for which a decrease in AIC is no
#' longer observed. Might also consider checking for stability in the beta parameters
#' for this as well.
gdistsamp_find_optimal_k <- function(
  df=NULL,
  allHabitatCovs=NULL,
  allDetCovs=NULL,
  mixture=NULL,
  keyfun="halfnorm",
  K=NULL)
{
  if(is.null(K)){
    K <- unlist(lapply(
        seq(1, 2.5, by=0.25),
        FUN=function(x) calc_k(df, multiplier=x)
      ))
  }
  all_covs_m <- unlist(lapply(
     K,
     FUN=function(x){
       m <- try(suppressMessages(unmarked::gdistsamp(
          # abundance
          as.formula(paste(
            "~",
            paste(allHabitatCovs, collapse="+"),
            "+offset(log(effort))",
            sep=""
          )),
          ~1, # availability
          as.formula(paste("~",paste(allDetCovs, collapse="+"))), # detection
          data=imbcr_df,
          keyfun="halfnorm",
          mixture=mixture,
          unitsOut="kmsq",
          se=T,
          K=x
         )))
       if(class(m) == "try-error"){
           return(NA)
       } else {
           return(OpenIMBCR:::AIC(m))
       }
     }
  ))
  # bug-fix : did we not find a single valid K parameter?
  if(all(is.na(all_covs_m))){
    stop(
        "couldn't find convergence for any K parameters in our sequence -- this shouldn't happen. Maybe there's",
        " zero-inflation of our input data set?"
      )
  # bug-fix : did we only find one working value? Use it
  } else if( sum(!is.na(all_covs_m)) == 1 ){
    return(
      list(
        K=K[which(!is.na(all_covs_m))],
        AIC=all_covs_m[which(!is.na(all_covs_m))]
      )
    )
  }
  # bug-fix swap out our NA values with noise that
  # our second derivative won't settle on
  if(any(is.na(all_covs_m))){
    all_covs_m[is.na(all_covs_m)] <-
      which(is.na(all_covs_m)) * 2 * max(all_covs_m, na.rm=T)
  }
  # at what array position does the change in AIC start to bottom out?
  tail <- quantile(abs(round(diff(diff(all_covs_m)),3)), p=0.6) # median, but more tail-ey
  min_k <- min(which(abs(round(diff(diff(all_covs_m)),2)) < tail))
  return(
    list(
      K=K[min_k],
      AIC=all_covs_m[min_k]
    )
  )
}
gdistsamp_find_optimal_key_func <- function(x=NULL){
  return(NA)
}
#' accepts a variable name an an input data.frame and 
#' conducts a Pearson's product moment correlation test
#' against all variables in the data.frame. 
#' @param x data.frame or matrix containing covariates to test
#' @param var variable name (as.character)
check_peasons_product_moment <- function(x=NULL, var=NULL){

}
#' testing :
check_correlation_matrix <- function(
  var=NULL,
  x_correlation_matrix=NULL,
  imbcr_df=NULL,
  correlation_threshold=0.5)
{
  correlation <- x_correlation_matrix[var,'PC1']
  # is this positive correlation? It shouldn't be
  if (correlation > 0) {
    warning(
      paste(
        "positive correlation observed between",var,
        "and our fragmentation metric. They should be the",
        "inverse of each other. Consider -1*PCA transformation")
    )
  }
  if (abs(correlation) > correlation_threshold){
      warning(paste("dropping",var,
        "-- it's strongly correlated with our fragmentation PCA")
      )
      imbcr_df@siteCovs <- imbcr_df@siteCovs[ , 
        !grepl(colnames(imbcr_df@siteCovs), pattern=var) ]
  }
  return(imbcr_df)
}
#' testing : drop some zero transects so that zero transects are not greater
#' than some multiple of our non-zero transects (default is 2X).
#' @param multiple multiple specifying how many times greater our zero 
#' transects can be relative to our non-zero transects.
balance_zero_transects <- function(imbcr_df=NULL, multiple=2){
  nonzero_transects <- which(rowSums(imbcr_df@y)!=0)
  zero_transects <- which(rowSums(imbcr_df@y)==0)
  if(length(zero_transects)/length(nonzero_transects) > multiple){
    downsample_to <- floor(length(nonzero_transects)*multiple)
    zero_transects <- sample(zero_transects, size=downsample_to)
    # resample our input table by throwing out some of our zero transects
    imbcr_df@y <- imbcr_df@y[c(nonzero_transects,zero_transects) , ]
    imbcr_df@siteCovs <- imbcr_df@siteCovs[c(nonzero_transects,zero_transects) , ]
    imbcr_df@tlength <- rep(1, nrow(imbcr_df@y))
  }
  return(imbcr_df)
}
#' testing : given a user-specified mixture distribution, try
#' and find an optimal K value from a sequence of possible K
#' values describing the upper-bounds of an integration (using
#' gdistsamp_find_optimal_k()). Often this will fail if we have
#' a lot of zeros in our input dataset. This function will attempt
#' a sequence of different downsampling densities to get around the
#' zero-inflation problem
gdistsamp_find_optimal_k_with_transect_downsampling <- function(
  imbcr_df=NULL,
  allHabitatCovs=NULL,
  allDetCovs=NULL,
  dist=NULL,
  K=NULL)
{
  if(is.null(K)){
    K <- unlist(lapply(
        seq(from=1,to=2.5,by=0.15),
        FUN=function(x) calc_k(imbcr_df, multiplier=x)
      ))
  }

  get_optim_k_aic <- function(
    x=NULL,
    imbcr_df=NULL,
    allHabitatCovs=NULL,
    allDetCovs=NULL,
    dist=NULL)
  {
     if(is.null(x)){
       downsample <- imbcr_df
     } else {
       downsample <- balance_zero_transects(imbcr_df, multiple=x)
     }
     # try to find an optimal K
     mixture_aic <- try(gdistsamp_find_optimal_k(
        df=downsample,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        mixture=dist,
        K=K
      ))
     # if we failed, try and downsample to nonzero*multiple density
     # and see if we can converge
     if (class(mixture_aic) == "try-error"){
       return(NA)
     } else {
       return(mixture_aic)
     }
  }

  # default : will attempt to find an optimal K without
  # downsampling first
  mixture <- get_optim_k_aic(x=NULL, imbcr_df, allHabitatCovs, allDetCovs, dist)

  if(is.na(mixture)){
    # generate a sequence of downsampling
    # densities and take the one that has
    # the lowest AIC -- this can take some
    # time
    multiple <- seq(from=0.9,to=0.3,by=-0.1)
    m <- lapply(
      X=multiple,
      FUN=function(x) get_optim_k_aic(x, imbcr_df, allHabitatCovs, allDetCovs, dist)
    )
    if(all(unlist(lapply(m, is.na)))){
      stop(
        "we failed to find convergence from a number of ",
        "different multiples -- quitting"
        )
    } else {
      multiple <- multiple[!unlist(lapply(m, is.na))]
      m <- m[!unlist(lapply(m, is.na))]
    }
    lowest_aic <- which.min(
        unlist(lapply(m, FUN=function(x) return(x$AIC)))
      )
    return(list(
        fit=m[[lowest_aic]],
        downsample_multiple=multiple[lowest_aic]
      ))
  } else {
    return(list(
      fit=mixture,
      downsample_multiple=NA
    ))
  }
}
