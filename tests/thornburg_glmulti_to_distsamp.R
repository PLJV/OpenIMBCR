
# consider using options(error=traceback)
options(warn = -1)
load(commandArgs(trailingOnly = T)[1])

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
CORRELATION_THRESHOLD     <- 0.55

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

#
# MAIN
#

source <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile("RawData_PLJV_IMBCR_20161201.csv"),
    four_letter_code = toupper(argv[2])
  )

# drop observations that are beyond our cut-off (but keep our NA's)
source <- source[ is.na(source$radialdistance) | (source$radialdistance <= cutoff) , ]

# overload 'breaks' and do some re-binning of our distance breaks to force
# a shoulder (yes, I know -- grumble, grumble)

breaks <- c(
  0,
  as.vector(quantile(
      source$radialdistance,
      probs = seq(0.5, 1, length.out = 8), na.rm = T)
    )
  )


# breaks <- round(seq(
#     from = 0,
#     to = quantile(source$radialdistance, 1, na.rm = T),
#     length.out = 7
#   ))

y <- do.call(rbind, lapply(
  X = unique(source$transectnum),
  FUN = function(x){
    matrix(table(cut(
        source[source$transectnum == x,]$radialdistance,
        breaks = breaks,
        labels = paste("dst_class_",1:(length(breaks)-1),sep = ""),sep = "")
      ), nrow = 1)
  }))

# calculate transect-level centroid (lat/lon) and merge it into our
# site-level covariates table for unmarked
coords <- do.call(rbind, lapply(
  X = unique(source$transectnum),
  FUN = function(x){
    d <- rgeos::gCentroid(sp::spTransform(
        s[s$transect == x ,],
        sp::CRS(raster::projection("+init=epsg:4326"))
      ))@coords
    return(data.frame(lon = d[1],lat = d[2]))
  }))

siteCovs <- do.call(rbind, lapply(
  X = unique(source$transectnum),
  FUN = function(x){
    s[s$transect == x ,]@data
  }))

siteCovs <- cbind(siteCovs, coords)
siteCovs$effort <- calc_transect_effort(source)

unmarked_data_frame <- unmarked::unmarkedFrameDS(
     y = y,
     siteCovs = siteCovs,
     #yearlySiteCovs=NULL,
     survey = "point",
     unitsIn = "m",
     dist.breaks = breaks
     #numPrimary=1
   )

# define the covariates we are going to use for our modeling
vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct","lat","lon")

# drop strongly-correlated variables
x_cor_matrix <- cor(siteCovs[,vars])

cor_threshold <- abs(x_cor_matrix) > CORRELATION_THRESHOLD

passing <- vector()
failing <- vector()

for(i in 1:ncol(cor_threshold)){
  not_correlated <- sum(which(cor_threshold[,i])!=i) == 0
  if(not_correlated){
    passing <- append(passing, colnames(cor_threshold)[i])
  } else {
    failing <- append(failing, colnames(cor_threshold)[i])
  }
}

rm(cor_threshold);

if(length(failing)>0){
  cat("-- these variables dropped due to colinearity problems:", failing,"\n")
  dev.new();
  corrplot::corrplot(x_cor_matrix)
  vars <- unlist(strsplit(
      readline("enter a comma-separated list of vars you want to keep:"),
      split=","
    ))
} else {
  vars <- passing
}

# Find an appropriate key function for our detection sub-model

keyfunction <-
  data.frame(
      key = c("hazard",
        "exp",
        "uniform",
        "halfnorm")
    )

keyfunction$aic <- as.vector(sapply(
  X = as.vector(keyfunction$key),
  FUN = function(x){
    unmarked::distsamp(
      formula = ~1~1+offset(log(effort)),
      keyfun = x,
      unitsOut = "kmsq",
      output = "abund",
      se = F,
      data = unmarked_data_frame
      )@AIC
  }))+AIC_RESCALE_CONST

# test: do we substantially improve on the half-normal detection function
# by using an alternative? If no, stick to half-normal
if( keyfunction[which.min(keyfunction$aic) , 'aic'] <
    keyfunction$aic[keyfunction$key == "halfnorm"]-(AIC_SUBSTANTIAL_THRESHOLD*200)
  ){
  keyfunction <- as.vector(keyfunction[which.min(keyfunction$aic) , 'key'])
} else {
  keyfunction <- "halfnorm"
}

# Fit a full model to our passing variables
unmarked_m_full <- unmarked::distsamp(
      formula = as.formula(paste(
        "~1~",
        paste(paste("poly(",vars,",2)",sep = ""),collapse = "+"),
        "+offset(log(effort))",
        sep = ""
        )),
      data = unmarked_data_frame,
      keyfun = keyfunction,
      unitsOut = "kmsq",
      output = "abund",
      se = T
)

# param <- exp(coef(unmarked_m_full, type="det"))
# plot(function(x) gxhn(x, param), 0, max(breaks),
# 	xlab="Distance (m)", ylab="Detection probability")
# grid(); grid()

# Drop unintuitive quadratic terms
quadratics <- quadratics_to_keep(unmarked_m_full)

if(length(quadratics)>0){
  vars_to_lin_poly <- function(x){
    if(length(!grepl(x, pattern = "1)"))>0){
      already_scaled <- x[grepl(x, pattern="1)")]
      regulars <- paste("poly(",x[!grepl(x, pattern="1)")]," ,1)",sep="")
      return(c(already_scaled,regulars))
    } else {
      return(x)
    }
  }
  vars <- vars[
    !as.vector(sapply(vars, FUN=function(p){
        sum(grepl(x=quadratics, pattern=p))>0
      }))
  ]
  # use AIC to justify our proposed quadratic terms
  for(q in quadratics){
    lin_var <- gsub(
        q,
        pattern=", 2[)]",
        replacement=", 1)"
      )

    m_lin_var <- OpenIMBCR:::AIC(fit_distsamp(
        formula=paste(
            c(ifelse(length(vars)>0, vars_to_lin_poly(vars) ,""),
              lin_var,
              quadratics[!(quadratics %in% q)]
              ),
              collapse="+"
            ),
        data=unmarked_data_frame,
        keyfun=keyfunction
      ))+AIC_RESCALE_CONST

    m_quad_var <- OpenIMBCR:::AIC(fit_distsamp(
        formula=paste(
        c( ifelse(length(vars)>0, vars_to_lin_poly(vars) ,""),
           quadratics
           ),
           collapse="+"
         ),
        data=unmarked_data_frame,
        keyfun=keyfunction
      ))+AIC_RESCALE_CONST
    # if we don't improve our AIC with the quadratic by at-least X aic units
    # (pretty substatial support), keep the linear version
    if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD){
      quadratics <- quadratics[!(quadratics %in% q)]
      vars <- c(vars,lin_var)
    }
  }

  # now convert any remaining linear vars to a scale consistent with poly()
  if(sum(!grepl(vars,pattern="poly"))>0){
    vars <- vars_to_lin_poly(vars)
  }

  vars <- c(vars, quadratics)

  # refit our full model minus un-intuitive quadratics
  unmarked_m_full <- fit_distsamp(
      formula=paste(vars, collapse="+"),
      unmarked_data_frame,
      keyfun=keyfunction
    )
}

# use model selection across our candidate variables

formulas <- calc_all_distsamp_combinations(vars)

unmarked_tests <- lapply(
    formulas,
    FUN = fit_distsamp,
    data = unmarked_data_frame,
    keyfun = keyfunction
  )

unmarked_tests <- unmarked::fitList(unmarked_tests)
  names(unmarked_tests) <- formulas

model_selection_table <- unmarked::modSel(unmarked_tests)

unmarked_top_models <- model_selection_table@Full$formula[
    model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD
  ]

aic_weights <- model_selection_table@Full$AICwt[
    model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD
  ]

formulas <- sapply(
  strsplit(sapply(strsplit(unmarked_top_models, split = "~"),
      FUN=function(x){ x[3] }), split = "[+]"),
      FUN=function(x) { return(paste(x[1:(length(x) - 1)], collapse = "+")) }
    )

# ok -- in order to predict with these models, we need to abandon the scaling
# used by poly during our model selection work so that our output
# is numerically consistent with observations taken across our 1
# km units data set. We can do this by appending raw=T to poly()

# linear terms with scale
formulas <- unlist(lapply(
  X=formulas,
  FUN=function(x){
    x <- strsplit(x, split = "[+]")
    x <- lapply(x[[1]],
    FUN=function(i){
      if(grepl(i, pattern = ", 1")){
          i <- gsub(i, pattern="poly[(]", replacement = "scale(")
          i <- gsub(i, pattern="[,] ", replacement = "^")
          i <- gsub(i, pattern="\\^1", replacement = "")
      }
        return(i)
      })
    paste(x, collapse = "+")
  }
))

unmarked_tests <- lapply(
    formulas,
    FUN=fit_distsamp,
    data=unmarked_data_frame,
    keyfun=keyfunction
  )

# add a fake effort offset
units$effort <- mean(unmarked_data_frame@siteCovs$effort)

cat(" -- calculating lat/lon for grid cell centroids\n")

coords <- rgeos::gCentroid(units, byid = T)

coords <- sp::spTransform(
    coords,
    sp::CRS(raster::projection("+init=epsg:4326"))
  )@coords[, c(1, 2)]

colnames(coords) <- c("lon","lat")

units@data <- cbind(
    units@data[, !grepl(names(units), pattern = "lon|lat")],
    coords
  )

# do some weighted model averaging
cat(
    " -- model averaging across prediction units",
    "table (this could take some time)\n"cd
  )

predict_df <- units@data

#steps <- round(seq(1,nrow(predict_df), length.out=parallel::detectCores()-1))
# add 1 to the last step to accomodate our lapply splitting
#steps[length(steps)] <- steps[length(steps)]+1

require(parallel)
cl <- parallel::makeCluster(parallel::detectCores() - 1)

#parallel::clusterExport(cl=cl, varlist=c("predict_df","steps"))
parallel::clusterExport(cl=cl, varlist = c("predict_df"))

# version 1: don't predict in row chunks -- parallelize across models
predicted <- parallel::parLapply(
  cl=cl,
  X=unmarked_tests,
  fun=function(model){
    return(unmarked::predict(
        model,
        type = 'state',
        newdata = predict_df,
        se = F,
        engine = 'C'
        #newdata=predict_df[seq(steps[i], (steps[i+1]-1)),]
      ))
  })

# version 2: don't predict across models -- parallelize across row chunks
# which is much faster but might distort model results a bit due to scale()
# issues
# predicted <- lapply(
#  X=unmarked_tests,
#  FUN=function(model){
#    predicted <- parallel::parLapply(
#      cl=cl,
#      X=1:(length(steps)-1),
#      fun=function(i, type='state'){
#        return(unmarked::predict(
#            model,
#            type='state',
#            se=F,
#            engine='C',
#            newdata=predict_df[seq(steps[i], (steps[i+1]-1)),]
#          ))
#      })
#    }
#)

# clean-up our 'parallel' call
rm(predict_df);
parallel::stopCluster(cl); rm(cl);
# join the individual rows from within each model run
# and return a single data.frame for each model
#predicted <- lapply(
#    predicted,
#    FUN=function(x) do.call(rbind, x)
#  )
# each model run will have a predicted column and
# a confidence interval -- we are only interested in
# the predicted column right now
predicted <- lapply(
    X=1:length(aic_weights),
    FUN=function(x){ predicted[[x]]$Predicted }
  )
# if we have more than one model in the top models
# table, let's average the results across our models
# using AIC weighting parameter taken from 'unmarked'.
if (length(predicted) > 1){
  # join our predictions across models into a single
  # matrix that we can apply a weight across
  predicted <- do.call(cbind, predicted)
  predicted <- sapply(
      1:nrow(predicted),
      FUN=function(i){
        weighted.mean(x=predicted[i, ], w = aic_weights)
      }
    )
} else {
  predicted <- unlist(predicted)
}
# copy our units shapefile for our predictions
pred_units <- units;
pred_units@data <- data.frame(pred = predicted);

rm(predicted);

# above, scale(units) compresses our explanatory variables slightly relative to
# the data used to fit each model -- below is a hack that uses min-max normalization
# to re-scale our predictions relative to our full model's predicted max(N)

PRED_MAX <- max(unmarked::predict(
  unmarked_m_full, type="state")[, 1], na.rm = T)

pred_units$pred[pred_units$pred < 1] <- 0
pred_units$pred[pred_units$pred > PRED_MAX] <- PRED_MAX

pred_units$pred <-
  ( pred_units$pred - min(pred_units$pred) ) /
  ( max(pred_units$pred, na.rm=T) - min(pred_units$pred, na.rm=T) )

pred_units$pred <- pred_units$pred * PRED_MAX

rgdal::writeOGR(
  pred_units,
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
