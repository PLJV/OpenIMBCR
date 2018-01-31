# Author: Kyle Taylor <kyle.taylor@pljv.org>
# Year : 2016
# Description : various tasks related to making variable response plots 
# and visualizations from unmarked and other objects.

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

#' build a response plot for a given variable using a fitted unmarked model
#' @export
partialPredict <- function(m=NULL, var=NULL, type='state', plot=T, nCores=NULL,
                           xlab=NULL, ylab=NULL, main=NULL, xlim=NULL,
                           ylim=NULL, xTransform=NULL){
  require(parallel)
  nCores <- ifelse(is.null(nCores), parallel::detectCores()-1, nCores)
      cl <- parallel::makeCluster(nCores)
  # fetch our site-level covariates from the training data
  t <- na.omit(m@data@siteCovs)
  # run our model for each unique value of x, averaging the predictions across all non-focal variables as we go
  x <- seq((min(t[,var])-sd(t[,var])),(max(t[,var])+sd(t[,var])),length.out=300)
  y <- matrix(NA, nrow=length(x), ncol=4)
  partial <- function(x){
    require(unmarked)
    run_table <- t
      run_table[, var] <- x
    p <- predict(m, run_table, type = type)
      y <- colMeans(p)
  }
  y <- do.call(rbind, parLapply(cl, x, fun=partial))
    colnames(y) <- c("predicted", "se", "lower", "upper")
  parallel::stopCluster(cl)
  if(!is.null(xTransform)){
    x <- eval(parse(text = xTransform))
  }
  # define our plot parameters, if asked
  if(plot){
    xlab <- if (is.null(xlab)) var else xlab
    ylab <- if (is.null(ylab)) "Density (birds/km2)" else ylab
    xlim <- if (is.null(xlim)) range(x) else xlim
    ylim <- if (is.null(ylim)) range(y) else ylim
    dev.new()
    plot(y=y[,1], x=x, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
         main=main, col="white")
    grid(lwd=1.4)
    lines(y=y[,1], x=x, col="black")
      lines(y=y[,4], x=x, col="grey")
        lines(y=y[,3], x=x, col="grey")
  } else {
    return(list(x=x, y=y))
  }
}
#' plot a half-normal detection function from an unmarked object
plot_hn_det <- function(x=NULL, breaks=NULL){
  param <- exp(coef(x, type = "det"))
  plot(
    function(x) unmarked:::gxhn(x, param), 0, max(breaks),
    xlab = "Distance (m)", ylab = "Detection probability"
  )
  grid(); grid();
}