# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2017
# Description : various tasks related to model weighting, averaging, and
# prediction for large model sets.

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

#' testing: use the parallel package to predict across a large input
#' table (with unmarked) using chunking.
#' @export
par_unmarked_predict <- function(run_table=NULL, m=NULL){

      steps <- seq(0, nrow(run_table), by=100)
  run_table <- data.frame(run_table@data)

  if(steps[length(steps)] != nrow(run_table)){
    steps <- append(steps, nrow(run_table))
  }

  cl <- parallel::makeCluster(parallel::detectCores()-1)

  parallel::clusterExport(
      cl,
      varlist=c("run_table","m","steps"),
      envir=environment()
    )

  predicted_density <- parallel::parLapply(
    cl=cl,
    X=1:(length(steps)-1),
    fun=function(x){
      unmarked::predict(
        m,
        newdata=run_table[(steps[x]+1):steps[x+1],],
        type="lambda"
      )
    }
  )
  # bind list results and return table to user
  return(do.call(rbind, predicted_density))
}
#' testing : us Akaike weights to predict across a series of models
#' defined in a user-specified model selection table
akaike_predict <- function(
  mod_sel_tab=NULL, 
  train_data=NULL,
  pred_data=NULL, 
  daic_cutoff=2,
  K=NULL, 
  mixture=NULL,
  keyfun=NULL,
  drop_intercept_m=T)
{
  keep <- mod_sel_tab$AIC < min(mod_sel_tab$AIC) + daic_cutoff
  mod_sel_tab <- mod_sel_tab[keep,]
  # sometimes an intercept only model shows up within 2AIC of the top model --
  # by default, let's not use it for model averaging
  if(drop_intercept_m){
    mod_sel_tab <- mod_sel_tab[
        !grepl(as.character(mod_sel_tab[,1]), pattern="~1+offset") , 
      ]
  }
  models <- lapply(
    X=as.character(mod_sel_tab$formula),
    FUN=function(x){
        OpenIMBCR:::gdistsamp_refit_model(
            formula = x, 
            imbcr_df = train_data, 
            K = K, 
            mixture = mixture,
            keyfun=keyfun
        )
      }
  )
  return(lapply(
      X=1:nrow(mod_sel_tab), 
      function(x) list(model=models[[x]], weight=mod_sel_tab$weight[x])
    ))
}
#' testing : accepts a list of akaike predictions (as returned by akaike_predict
#' ) and uses the corresponding $weight list-item for each model to build
#' a data.frame of weighted-average predictions. 
#' @param col column of the prediction table to average. Typically there are th
#' ree columns from predict, specifying prediction, upper-bound, and lower-bound
akaike_weight_predictions <- function(akaike_list=NULL, col=1){
    table <- (do.call(cbind, lapply(
        X=akaike_list,
        FUN=function(x){
            x$prediction[,col]
        }
    )))
    weights <- unlist(lapply(
        X=akaike_list,
        FUN=function(x){
            x$weight
        }
    ))
    return(apply(
        table,
        MARGIN=1,
        FUN=weighted.mean,
        w=weights
    ))
}
