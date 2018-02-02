load(commandArgs(trailingOnly = T)[1])
# consider using options(error=traceback)
options(warn = -1, error=traceback)

if (file.exists(gsub(
    r_data_file,
    pattern = "pois_glm",
    replacement = "hds_pois"
  ))){
  cat("previous workspace for this bird found in CWD (skipping)\n")
  q("no")
}

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65

#
# Local accessory functions, some of which may overload what's
# saved in the loaded Rdata file
#

plot_hn_det <- function(x=NULL, breaks=NULL){
  param <- exp(coef(x, type = "det"))
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability"
    )
  grid(); grid();
}

plot_model_pi <- function(tests = NULL, unmarked_m = NULL){
  hinge_m_pred <- density(as.vector(floor(
      predict(tests@objects[[1]], type = "response")))
    )
  hds_m_pred <- density(as.vector(floor(
      unmarked::predict(unmarked_m, type = "state")[, 1]))
    )

  xlim <- range(c(hinge_m_pred$x, hds_m_pred$x))
  ylim <- range(c(hinge_m_pred$y, hds_m_pred$y))

  dev.new();

  plot(
      hinge_m_pred,
      xlim = c(xlim[1] - (diff(xlim) * 0.1), xlim[2] + (diff(xlim) * 0.1)),
      ylim = c(ylim[1], ylim[2] + 0.01),
      col = "red",
      main = ""
    )

  lines(
      hds_m_pred,
      col = "blue"
    )

  grid(); grid();
}

fit_distsamp <- function(formula = NULL, data = NULL, keyfun = "halfnorm"){
  return(unmarked::distsamp(
        formula = as.formula(paste(
          "~1~",
          formula,
          "+offset(log(effort))",
          sep = ""
          )),
        #formula = 1+offset(log(effort))~1+offset(log(effort)),
        data = data,
        keyfun = keyfun,
        unitsOut = "kmsq",
        output = "abund",
        #starts = c(c(as.vector(coef(m))),1,1),
        se = T
  ))
}

quadratics_to_keep <-function(m){

     linear_terms <- grepl(names(unmarked::coef(m)), pattern = ")1")
  quadratic_terms <- grepl(names(unmarked::coef(m)), pattern = ")2")

  # test : are we negative and are we a quadratic term?
  keep <- (unmarked::coef(m) < 0) * quadratic_terms
    keep <- names(unmarked::coef(m))[keep==1]
      keep <- gsub(keep, pattern = ")2", replacement = ")")
  if(length(keep)>0){
    keep <- gsub(keep, pattern = "lam[(]", replacement = "")
      keep <- gsub(keep, pattern = "[)][)]", replacement = ")")
    return(keep)
  } else {
    return(NULL)
  }
}

calc_all_distsamp_combinations <- function(vars = NULL){
  formulas <- gsub(OpenIMBCR:::mCombinations(
          siteCovs = vars,
          availCovs = NULL,
          detCovs = NULL,
          offset = "offset(log(effort))")[,1],
        pattern = " ~1 ~1",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = " ",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = "[+]offset[(]log[(]effort[)][)]",
      replacement = ""
    )
  formulas <- gsub(formulas, pattern = "~", replacement = "")
  return(formulas)
}

calc_emp_dispersion_statistic <- function(x = NULL){
  observed <- sd(x)/round(mean(x))
  predicted <- rpois(n = length(x), round(mean(x)))
    predicted <- sd(predicted)/round(mean(predicted))
  return(observed/predicted)
}

glm_to_distsamp <- function(m=NULL, umdf=NULL){
  umdf@siteCovs <- m$data

  to_use <- names(coefficients(m))
    to_use <- to_use[2:length(to_use)]
  # drop any )1 or )2 poly() postfixes
  to_use <- unique(gsub(to_use, pattern="[)][0-9]", replacement=")"))

  to_use <- paste(
      "~1~",
      paste(to_use, collapse="+"),
      "+offset(log(effort))",
      sep=""
    )

  return(unmarked::distsamp(
      formula=as.formula(to_use),
      data=umdf,
      se=T,
      keyfun="halfnorm",
      unitsOut="kmsq",
      output="abund"
    ))
}

#
# MAIN
#

# top models to re-fit
keep <- tests@crits - min(tests@crits) < AIC_SUBSTANTIAL_THRESHOLD

# get weights
aic_weights <- OpenIMBCR:::akaike_weights(tests@crits)[keep]

# re-fit
unmarked_models <- lapply(
    X=tests@objects[keep],
    FUN=glm_to_distsamp,
    umdf=umdf
  )

# predict across as many cores as we can

units@data$effort <- median(effort)

predict_df <- units@data

require(parallel)
cl <- parallel::makeCluster(parallel::detectCores() - 1)

parallel::clusterExport(cl=cl, varlist = c("predict_df"))

# version 1: don't predict in row chunks -- parallelize across models
predict_df <- parallel::parLapply(
  cl=cl,
  X=unmarked_models,
  fun=function(model){
    return(unmarked::predict(
        model,
        type = 'state',
        newdata = predict_df,
        se = F
      ))
})

parallel::stopCluster(cl); rm(cl); rm(keep);

predict_df <- lapply(
    X=1:length(aic_weights),
    FUN=function(x){ predicted[[x]]$Predicted }
  )

if (length(predict_df) > 1){
  # join our predictions across models into a single
  # matrix that we can apply a weight across
  predict_df <- do.call(cbind, predict_df)
  predict_df <- sapply(
      1:nrow(predict_df),
      FUN=function(i){
        weighted.mean(x=predict_df[i, ], w = aic_weights)
      }
    )
} else {
  predict_df <- unlist(predict_df)
}

predicted@data <- data.frame(pred = predicted);

cat(" -- writing to disk\n")

rgdal::writeOGR(
  predicted,
  ".",
  tolower(paste(argv[2],
      "_imbcr_hds_pois_prediction_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      sep="")
  ),
  driver="ESRI Shapefile",
  overwrite=T
)

r_data_file <- gsub(
    r_data_file,
    pattern="pois_glm",
    replacement="hds_pois"
  )

save(
    compress=T,
    list=ls(pattern="[a-z]"),
    file=r_data_file
  )
