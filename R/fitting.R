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
akaike_weights <- function (aic_values = NULL, precision = 5){
  weights <- exp(-0.5 * (aic_values - min(aic_values, na.rm=T)))
  return(round(weights/sum(weights, na.rm=T), precision))
}
