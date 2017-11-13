# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2016
# Description : various tasks related to fitting models through optimization 
# or pre-canned modeling interfaces (as provided by 'unmarked').

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

#' expit link function
expit <- function(x) 1/(1+exp(-x))
#' logit link function
logit <- function(x) log(x/(1-x))
#' calculate AIC from a likelihood function optimization procedure
AIC <- function(m){
  if(inherits(m, "unmarkedFit")){
    return(m@AIC)
  # sometimes a calling function will just pass a single numeric value, 
  # for instance AICc() when passed a vector of standard AIC values
  } else if(is.numeric(m)) {
    return(m) 
  # default: if it isn't an unmarked model object, assume we rolled our own likelhood with nlm()
  } else {
    return((m$minimum*2) + (2*length(m$estimate)))
  }
}
#' calculate AICc from a likelihood function optimization procedure
AICc <- function(m=NULL, p=NULL, n=NULL){
  if(inherits(m, "unmarkedFit")){
    # num abundance parameters minus intercept
    p <- unlist(strsplit(paste(as.character(
          m@formlist[[1]]),
          collapse=""
        ),
        "[+]"
      ))
    # don't include effort
    p <- sum(!grepl(p, pattern="effort"))
    # number of site observations
    n <- nrow(unmarked:::getY(m@data))
  } else {
    p <- length(m$estimate)
    n <- m$size
  }
  return(OpenIMBCR:::AIC(m) + ( (2*p*(p+1)) / (n-p-1) ) )
}
#' calculated BIC from a likelihood function optimization procedure
BIC <- function(m) (m$minimum*2) + (q*log(m$size))
#'
#' calculate SE from a likelihood function optimization procedure
SE <- function(m) sqrt(diag(solve(m$hessian)))
#' calculate Akaike weights from a vector of AIC values
akaike_weights <- function(aic_values=NULL, precision=5){
  weights <- exp( -0.5 * (aic_values-min(aic_values)) )
  return(round(weights/sum(weights), precision))
}
#' recursively calculate all possible permutations of an input table n
permutations <- function(n){
   if(n==1){
       return(matrix(1))
   } else {
       sp <- permutations(n-1)
       p <- nrow(sp)
       A <- matrix(nrow=n*p,ncol=n)
       for(i in 1:n){
           A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
       }
       return(A)
   }
}
#' re-fit a model specified by a formula object
gdistsamp_refit_model <- function(formula=NULL, imbcr_df=NULL, K=NULL, mixture=NULL, keyfun=NULL){
  formula <- unlist(strsplit(formula, split="~"))
  return(unmarked::gdistsamp(
      lambdaformula=as.formula(paste(
        "~",
        formula[2],
        sep=""
      )),
      phiformula=~1,
      pformula=as.formula(paste(
        "~",
        formula[4]
      )),
      data=imbcr_df,
      keyfun=keyfun,
      mixture=mixture,
      unitsOut="kmsq",
      se=T,
      K=K
    ))
}
#' fit the single-season distance model of Royle (2004) to pooled (station-level)
#' IMBCR data
#' @param habitatCovs vector of habitat covariates to fit
#' @param detCovs detection covariates to fit
#' @export
imbcr_single_season_pooled_1km_distance <- function(
    habitatCovs=NULL, 
    detCovs=NULL, 
    imbcr_df=NULL
  )
{
  K <- unlist(lapply(
      X=seq(1, 2, by=0.25),
      FUN=function(x) OpenIMBCR:::calc_k(imbcr_df, multiplier=x)
    ))
  #
  # Testing : select an optimal detection function
  #
  key_function <- c("halfnorm","hazard","exp","uniform")
  # test AIC for a sequence of K values for the Poisson mixture distribution
  pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
      df=imbcr_df,
      allHabitatCovs=habitatCovs,
      allDetCovs=detCovs,
      mixture="P",
      K=K
    ))
  # if there is zero-inflation, K selection may fail. Try to downsample our
  # zeros to something managable and see if we still fail to fit a model
  if (class(pois_aic) == "try-error"){
      imbcr_df <- OpenIMBCR:::balance_zero_transects(imbcr_df, multiple=0.5)
      pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
          df=imbcr_df,
          allHabitatCovs=habitatCovs,
          allDetCovs=detCovs,
          mixture="P",
          K=K
        ))
      # if downsampling didn't help -- quit with an error. We don't
      # want to try to interpret a model with questionable data.
      if (class(pois_aic)=="try-error"){
        stop("poisson : couldn't find an optimal K value,",
             "even after downsampling over-abundant zeros")
      } else {
        K_pois <- pois_aic$K
        pois_aic <- pois_aic$AIC
      }
  } else {
    K_pois <- pois_aic$K
    pois_aic <- pois_aic$AIC
  }
  # now try the same thing for the negative-binomial mixture distribution
  negbin_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
      df=imbcr_df,
      allHabitatCovs=habitatCovs,
      allDetCovs=detCovs,
      mixture="NB",
      K=K
    ))
  # did we fail because of zero-inflation? Try to re-fit.
  if (class(negbin_aic) == "try-error"){
      # this will do nothing if we've already downsampled to 'multiple'
      # for the Poisson mixture
      imbcr_df <- OpenIMBCR:::balance_zero_transects(imbcr_df, multiple=0.75)
      pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
          df=imbcr_df,
          allHabitatCovs=habitatCovs,
          allDetCovs=detCovs,
          mixture="NB",
          K=K
        ))
      # if downsampling didn't help -- quit with an error. We don't
      # want to try to interpret a model with questionable data.
      if (class(negbin_aic)=="try-error"){
        stop("poisson : couldn't find an optimal K value,",
             "even after downsampling over-abundant zeros")
      } else {
        K_negbin <- negbin_aic$K
        pois_aic <- negbin_aic$AIC
      }
  } else {
    K_negbin <- negbin_aic$K
    negbin_aic <- negbin_aic$AIC
  }
  # is there support for the negative binomial being a superior mixture dist
  # to the poisson?
  if(diff(c(negbin_aic, pois_aic)) > 4){
    K <- K_negbin
    mixture_dist <- "NB"
    # if we don't see a large improvement in AIC from using the
    # negative binomial, favor the simpler Poisson mixture
  } else {
    K <- K_pois
    mixture_dist <- "P"
  }
  # do model selection for our input covariates
  model_selection_table <- OpenIMBCR:::allCombinations_dAIC(
    siteCovs=habitatCovs,
    detCovs=detCovs,
    step=500,
    umdf=imbcr_df,
    umFunction=unmarked::gdistsamp,
    mixture=mixture_dist,
    unitsOut="kmsq",
    K=K,
    se=T,
    keyfun="halfnorm",
    offset="offset(log(effort))"
  )

}
#' testing: fit a single-season occupancy model  assumes a constant probability of species
#' detection across transects, which is probably inappropriate for IMBCR data and will
#' lead to inaccurate predictions of occupancy. This is "model m0" from the literature
#' and loosely follows Andy Royle's (2008) model specification. It is designed so that
#' parameter estimates can be derived through an optimization procedure like nlm() and
#' variables can be easily selected/un-selected, such as for calculating AIC.
#'
#' @param table a data.frame of detection histories ('det') and covariates on
#' occupancy ('b0') and probability of detection ('a0')
#' @param det_parameters proposed parameter values for detection, including 1 intercept (a0)
#' @param occ_parameters proposed paramter values for occupancy, including 1 intercept (b0)
#'
#' @export
singleSeasonOccupancy <- function(det_parameters=NULL, occ_parameters=NULL, p=NULL, table=NULL){
  if(is.null(table)){
    stop("table= argument requires an input table containing covariates and a 'det' field defining detection histories for your sites.")
  }
  # assign a vars vector containing a composite of our det_parameters and occ_parameters
  if(is.null(p)){
    vars <- c("a0", colnames(table[,2:(length(det_parameters)+1)]), "b0", colnames(table[,(length(det_parameters)+2):ncol(table)]))
  } else {
    vars <- c("a0",colnames(table[,2:(det_parameters+1)]), "b0", colnames(table[,(det_parameters+2):ncol(table)]))
  }
  # The p-parameter is typically assigned by nlm() and should include proposed values for our two intercept terms (in addition to our covariates)
  if(!is.null(p)){
    det_parameters = p[1:(det_parameters+1)]
    occ_parameters = p[(length(det_parameters)+1):length(p)]
  }
  if(length(det_parameters)+length(occ_parameters) > ncol(table)+2-1 ){ # accounting for intercept terms, but ignoring the first column (detections)
    stop("length of paramters is greater than the number of columns specified by table=")
  }
  # define the number of transects, number of stations per transect, and target variables we are considering for this iteration
            M <- nrow(table) # number of sites (IMBCR transects)
    nStations <- nchar(as.character(table[1,1])) # number of stations in each transect (should be 16)
   covarNames <- vars
  # re-build a consistent table (t) of detection histories, intercepts, and covariate data we can work with
  y <- table[,1] # assume our 'detection' field is always the first column
    y <- suppressWarnings(matrix(as.numeric(matrix(unlist(strsplit(as.character(y),split="")))),nrow=M,ncol=nStations))
  t <- matrix(rep(0,M*length(covarNames)),ncol=length(covarNames)) # build a table of zeros for our covariate data
    colnames(t) <- covarNames
     t[,"a0"] <- rep(1,M) # fill our intercept columns
     t[,"b0"] <- rep(1,M)
  # assign values for our focal table from user-specific source table
  t <- data.frame(t)
  table <- data.frame(table)
  t[,vars[!grepl(vars,pattern=0)]] <- table[,vars[!grepl(vars,pattern=0)]] # assign all covariate data (minus intercept data) from source table
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # assign proposed parameters for this optimization step, from nlm()
  tryCatch(coeffs[vars] <- c(det_parameters,occ_parameters),
           warning=function(w){
             stop("caught a warning assigning detection + occurrence parameters + intercept terms. Did you remember to generate proposed paramter values for the intercept terms?")
           })
  detection_coeffs <- 1:(which(grepl(vars,pattern="b0$"))-1)
  occupancy_coeffs <- which(grepl(vars,pattern="b0$")):length(vars)
  # parameters on detection
  # matrix operation : prob <- expit(a0*t[,'a0'] + a1*t[,"tod"] + a2*t[,"doy"] + a3*t[,"intensity"])
  prob <- sweep(t[,vars[detection_coeffs]], MARGIN=2, coeffs[detection_coeffs],`*`)
    prob <- expit(rowSums(prob))
  # matrix operation : psi <- expit(b0*t[,'b0'] + b1*t[,"perc_ag"] + b2*t[,"perc_grass"] + b3*t[,"perc_shrub"] + b4*t[,"perc_tree"] + b5*t[,"perc_playa"])
  psi <- sweep(t[,vars[occupancy_coeffs]], MARGIN=2, coeffs[occupancy_coeffs],`*`)
     psi <- expit(rowSums(psi))
  # solve for likelihood
  likelihood <- rep(NA,M)
  for(i in 1:M){
    detections <- y[i,] # individual detections for our repeat visits at site i
    na_det <- is.na(detections) # any NA values?
    nd <- sum(detections[!na_det]) # check for zero-detections
    p <- prob[i] # what is the predicted probability of detections for this site, given our parameters?
    # calculate likelihood of occupancy, given our calculated probability of detection for the focal transect
    cp <- (p^detections)*((1-p)^(1-detections))
      cp[na_det] <- 1 # set any NA values to 1
    likelihood[i] <- log(prod(cp)*psi[i] + ifelse(nd==0,1,0)*(1-psi[i])) # joint probability across detections, e.g.: http://stats.stackexchange.com/questions/211848/likelihood-why-multiply
  }
  sum(-1*likelihood)
}
#' testing : fit the single-season, two-state occupancy model of Royle-Nichols (the "RN" model) to IMBCR data. This model can account for abundance-induced heterogeneity of detection and capitalize
#' on replicate observations of abundance within sampling units. Royle first implemented the algorithm for use with BBS data, but we can theoretically extend the model
#' within-season repeat observations made at IMBCR stations. Requires the user specify an upper_bound parameter for density within each 1-km transect. Like the single-season occupancy model (model m0), the RN model
#' does not allow for heterogeneity in detection within sites, which may result in biased estimates of detection in instances where an observer detects a species early in the sampling process.
#' Can be accomodated with a removal design.
#'
#' @export
singleSeasonRN <- function(det_parameters=NULL, occ_parameters=NULL, p=NULL, upper_bound=NULL, quasi_binom=FALSE, table=NULL){
  if(is.null(table)){
    stop("table= argument requires an input table containing covariates and a 'det' field defining detection histories for your sites.")
  }
  if(is.null(upper_bound)){
    stop("upper_bound= parameter (e.g., a rational prediction of an absolute maximum bird density per IMBCR transect) must be defined")
  }
  # assign a vars vector containing a composite of our det_parameters and occ_parameters
  if(is.null(p)){
    vars <- c("a0", colnames(table[,2:(length(det_parameters)+1)]), "b0", colnames(table[,(length(det_parameters)+2):ncol(table)]))
  } else {
    vars <- c("a0",colnames(table[,2:(det_parameters+1)]), "b0", colnames(table[,(det_parameters+2):ncol(table)]))
  }
  # The p-parameter is typically assigned by nlm() and should include proposed values for our two intercept terms (in addition to our covariates)
  if(!is.null(p)){
    det_parameters = p[1:(det_parameters+1)]
    occ_parameters = p[(length(det_parameters)+1):length(p)]
  }
  if(length(det_parameters)+length(occ_parameters) > ncol(table)+2-1 ){ # accounting for intercept terms, but ignoring the first column (detections)
    stop("length of paramters is greater than the number of columns specified by table=")
  }
  # define the number of transects, number of stations per transect, and target variables we are considering for this iteration
            M <- nrow(table) # number of sites (IMBCR transects)
    nStations <- nchar(as.character(table[1,1])) # number of stations in each transect (should be 16)
   covarNames <- vars
  # re-build a consistent table (t) of detection histories, intercepts, and covariate data we can work with
  y <- table[,1] # assume our 'detection' field is always the first column
    y <- suppressWarnings(matrix(as.numeric(matrix(unlist(strsplit(as.character(y),split="")))),nrow=M,ncol=nStations))
  t <- matrix(rep(0,M*length(covarNames)),ncol=length(covarNames)) # build a table of zeros for our covariate data
    colnames(t) <- covarNames
     t[,"a0"] <- rep(1,M) # fill our intercept columns
     t[,"b0"] <- rep(1,M)
  # assign values for our focal table from user-specific source table
  t <- data.frame(t)
  table <- data.frame(table)
  t[,vars[!grepl(vars,pattern=0)]] <- table[,vars[!grepl(vars,pattern=0)]] # assign all covariate data (minus intercept data) from source table
  # by default, set our coefficients = 0 for this run
  coeffs <- rep(0,length(covarNames))
    names(coeffs) <- covarNames
  # assign proposed parameters for this optimization step, from nlm()
  tryCatch(coeffs[vars] <- c(det_parameters,occ_parameters),
           warning=function(w){
             stop("caught a warning assigning detection + occurrence parameters + intercept terms. Did you remember to generate proposed paramter values for the intercept terms?")
           })
  detection_coeffs <- 1:(which(grepl(vars,pattern="b0$"))-1)
  occupancy_coeffs <- which(grepl(vars,pattern="b0$")):length(vars)
  # parameters on detection (the rate parameter of our poisson)
  # matrix operation : prob <- expit(a0*t[,'a0'] + a1*t[,"tod"] + a2*t[,"doy"] + a3*t[,"intensity"])
  r <- sweep(t[,vars[detection_coeffs]], MARGIN=2, coeffs[detection_coeffs],`*`)
    r <- expit(rowSums(r))
  # matrix operation : psi <- expit(b0*t[,'b0'] + b1*t[,"perc_ag"] + b2*t[,"perc_grass"] + b3*t[,"perc_shrub"] + b4*t[,"perc_tree"] + b5*t[,"perc_playa"])
  lambda <- sweep(t[,vars[occupancy_coeffs]], MARGIN=2, coeffs[occupancy_coeffs],`*`)
     psi <- exp(rowSums(lambda))
  # solve for likelihood
  likelihood <- rep(NA,M)
  for(i in 1:M){
    # determine a truncated (by upper_bound) probability for our Poisson count,
    gN <- dpois(0:upper_bound, lambda[i])
      gN <- gN/sum(gN)
    # individual detections for our repeat visits at site i (e.g., [1,1,1,3,2,1])
    detections <- y[i,]
    na_det <- is.na(detections) # any NA values?

    pmat <- 1 - outer((1-r[i,]),0:upper_bound,"^")

    focal <- t((pmat^detections)*(1-pmat)^(1-detections))
      focal[,na_det] <- 1
        focal <- apply(focal,1,prod)

    likelihood[i] <- sum(focal*gN)
 }
  sum(-1*likelihood)
}

