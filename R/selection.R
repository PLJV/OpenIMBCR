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
#' estimates a K parameter (upper-bound of integration) as a pre-cursor for
#' various count-based mixture models (e.g., poisson and negative binomial)
calc_k <- function(df=NULL, multiplier=1){
  return(floor(max(
      rowSums(unmarked:::getY(df)))*multiplier
    ))
}

