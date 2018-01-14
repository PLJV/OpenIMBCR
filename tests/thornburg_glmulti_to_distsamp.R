
options(error=traceback)
load(commandArgs(trailingOnly=T)[1])

if(file.exists(gsub(
    r_data_file, 
    pattern="pois_glm", 
    replacement="hds_pois"
  ))){
  stop("previous workspace for this bird found in CWD (skipping)")
}

#
# Local accessory functions, some of which may overload what's 
# saved in the loaded Rdata file
#

plot_model_pi <- function(tests=NULL, unmarked_m=NULL){
  hinge_m_pred <- density(as.vector(floor(predict(tests@objects[[1]], type="response"))))
  hds_m_pred <- density(as.vector(floor(unmarked::predict(unmarked_m, type="state")[,1])))

  xlim <- range(c(hinge_m_pred$x, hds_m_pred$x))
  ylim <- range(c(hinge_m_pred$y, hds_m_pred$y))

  dev.new();
  
  plot(
      hinge_m_pred, 
      xlim=c(xlim[1]-(diff(xlim)*0.1),xlim[2]+(diff(xlim)*0.1)),
      ylim=c(ylim[1], ylim[2]+0.01),
      col="red",
      main=""
    )
    
  lines(
      hds_m_pred,
      col="blue"
    )

  grid(); grid();
}

fit_distsamp <- function(formula=NULL, data=NULL){
  return(unmarked::distsamp(
        formula=as.formula(paste(
          "~1+offset(log(effort))~",
          formula,
          "+offset(log(effort))",
          sep=""
          )),
        #formula=1+offset(log(effort))~1+offset(log(effort)),
        data=data,
        keyfun="halfnorm",
        unitsOut="kmsq",
        output="abund",
        #starts=c(c(as.vector(coef(m))),1,1),
        se=T
  ))
}

quadratics_to_keep <-function(m){

     linear_terms <- grepl(names(unmarked::coef(m)), pattern=")1")
  quadratic_terms <- grepl(names(unmarked::coef(m)), pattern=")2") 

  # test : are we negative and are we a quadratic term?
  keep <- (unmarked::coef(m) < 0) * quadratic_terms
    keep <- names(unmarked::coef(m))[keep==1]
      keep <- gsub(keep, pattern=")2", replacement=")")
  if(length(keep)>0){
    keep <- gsub(keep, pattern="lam[(]", replacement="")
      keep <- gsub(keep, pattern="[)][)]", replacement=")")
    return(keep)
  } else {
    return(NULL)
  }
}

calc_all_distsamp_combinations <- function(vars=NULL){
  formulas <- gsub(OpenIMBCR:::mCombinations(
          siteCovs=vars, 
          availCovs=NULL, 
          detCovs=NULL, 
          offset="offset(log(effort))")[,1], 
        pattern=" ~1 ~1", 
      replacement=""
    )
  formulas <- gsub(
      formulas, 
      pattern=" ", 
      replacement=""
    )
  formulas <- gsub(
      formulas, 
      pattern="[+]offset[(]log[(]effort[)][)]", 
      replacement=""
    )
  formulas <- gsub(formulas, pattern="~", replacement="")
  return(formulas)
}

calc_emp_dispersion_statistic <- function(x=NULL){
  observed <- sd(x)/round(mean(x))
  predicted <- rpois(n=length(x), round(mean(x)))
    predicted <- sd(predicted)/round(mean(predicted))  
  return(observed/predicted)
}

#
# MAIN
#

source <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile("RawData_PLJV_IMBCR_20161201.csv"),
    four_letter_code=toupper(argv[2])
  )

# drop observations that are beyond our cut-off (but keep our NA's)
source <- source[ is.na(source$radialdistance) | (source$radialdistance <= cutoff) , ]

y <- do.call(rbind, lapply(
  X=unique(source$transectnum),
  FUN=function(x){
    matrix(table(cut(
        source[source$transectnum == x,]$radialdistance, 
        breaks=breaks, 
        labels=paste("dst_class_",1:(length(breaks)-1),sep=""),sep="")
      ), nrow=1)
  }))

siteCovs <- do.call(rbind, lapply(
  X=unique(source$transectnum),
  FUN=function(x){
    s[s$transect == x ,]@data
  }))
  
siteCovs$effort <- calc_transect_effort(source)
  
unmarked_data_frame <- unmarked::unmarkedFrameDS(
     y=y,
     siteCovs=siteCovs,
     #yearlySiteCovs=NULL,
     survey="point",
     unitsIn="m",
     dist.breaks=breaks
     #numPrimary=1
   )

vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct")

# drop strongly-correlated variables
x_cor_matrix <- cor(s@data[,vars])

cor_threshold <- abs(x_cor_matrix)>0.55

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

# Fit a full model to our passing variables
unmarked_m_full <- unmarked::distsamp(
      formula=as.formula(paste(
        "~1+offset(log(effort))~",
        paste(paste("poly(",vars,",2)",sep=""),collapse="+"),
        "+offset(log(effort))",
        sep=""
        )),
      #formula=1+offset(log(effort))~1+offset(log(effort)),
      data=unmarked_data_frame,
      keyfun="halfnorm",
      unitsOut="kmsq",
      output="abund",
      #starts=c(c(as.vector(coef(m))),1,1),
      se=T
)

# Drop unintuitive quadratic terms
quadratics <- quadratics_to_keep(unmarked_m_full)

if(length(quadratics)>0){
  vars <- vars[
    !as.vector(sapply(vars, FUN=function(p){ 
        sum(grepl(x=quadratics, pattern=p))>0  
      }))
  ]
  # use AIC to justify our proposed quadratic terms
  for(q in quadratics){
    lin_var <- gsub(
        gsub(q,pattern="poly[(]", replacement=""),
        pattern=", 2[)]",
        replacement=""
      )
    
    m_lin_var <- OpenIMBCR:::AIC(fit_distsamp(
        formula=paste(
            c(vars,lin_var,quadratics[!(quadratics %in% q)]), collapse="+"
            ), 
        unmarked_data_frame
      ))+100000
      
    m_quad_var <- OpenIMBCR:::AIC(fit_distsamp(
        formula=paste(c(vars, quadratics), collapse="+"), 
        unmarked_data_frame
      ))+100000
    # if we don't improve our AIC with the quadratic by at-least 6 aic units
    # (pretty substatial support), keep the linear version
    if( m_lin_var-m_quad_var < 7){
      quadratics <- quadratics[!(quadratics %in% q)]
      vars <- c(vars,lin_var)
    }
  }
    
  vars <- c(vars, quadratics)

  # refit our full model minus un-intuitive quadratics
  unmarked_m_full <- fit_distsamp(
      formula=paste(vars, collapse="+"),
      unmarked_data_frame
    )
}

# use model selection across our candidate variables

formulas <- calc_all_distsamp_combinations(vars)

unmarked_tests <- lapply(
    formulas, 
    FUN=fit_distsamp, data=unmarked_data_frame
  )

names(unmarked_tests) <- formulas
  unmarked_tests <- unmarked::fitList(unmarked_tests)

model_selection_table <- unmarked::modSel(unmarked_tests)

unmarked_top_models <- model_selection_table@Full$formula[
    model_selection_table@Full$delta < 6
  ]

aic_weights <- model_selection_table@Full$AICwt[
    model_selection_table@Full$delta < 6
  ]

formulas <- sapply(
  strsplit(sapply(strsplit(unmarked_top_models, split="~"), 
      FUN=function(x){ x[3] }), split="[+]"), 
      FUN=function(x) { return(paste(x[1:(length(x)-1)], collapse="+")) }
    )

unmarked_tests <- lapply(
    formulas, 
    FUN=fit_distsamp, data=unmarked_data_frame
  )

units$effort <- mean(unmarked_data_frame@siteCovs$effort)

# do some weighted model averaging

cat(
    " -- model averaging across prediction units",
    "table (this could take some time)\n"
  )

predict_df <- units@data

steps <- round(seq(1,nrow(predict_df), length.out=100))
# add 1 to the last step to accomodate our lapply splitting
steps[length(steps)] <- steps[length(steps)]+1

require(parallel)
cl <- parallel::makeCluster(parallel::detectCores()-1)

parallel::clusterExport(cl=cl, varlist=c("predict_df","steps"))

predicted <- lapply(
  X=unmarked_tests, 
  FUN=function(model){
    predicted <- parallel::parLapply(
      cl=cl,
      X=1:(length(steps)-1), 
      fun=function(i, type='state'){ 
        return(unmarked::predict(
            model, 
            type=type, 
            state=state, 
            newdata=predict_df[seq(steps[i], (steps[i+1]-1)),] 
          )) 
      })
    }
)

predicted <- lapply(predicted, FUN=function(x) do.call(rbind, x))

rm(predict_df); 
parallel::stopCluster(cl);
rm(cl);

predicted <- sapply(
    X=1:length(aic_weights), 
    FUN=function(x){ predicted[[x]]$Predicted }
  )
  
predicted <- sapply(
    1:nrow(predicted), 
    FUN=function(x){ weighted.mean(x=predicted[i,], w=aic_weights) }
  )

# copy our units shapefile for our predictions
pred_units <- units;
pred_units@data <- data.frame(pred=predicted); 

rm(predicted);

# censor any predictions greater than K (max)

cat(
    " -- number of sites with prediction greater than max(K):",
    sum(pred_units$pred > max(rowSums(unmarked_data_frame@y))),
    "\n"
  )

pred_units$pred[(pred_units$pred > max(rowSums(unmarked_data_frame@y)))] <- 
  max(rowSums(unmarked_data_frame@y))

pred_units$pred[pred_units$pred<1] <- 0

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

r_data_file <- gsub(r_data_file, pattern="pois_glm", replacement="hds_pois")

save(
    compress=T,
    list=ls(),
    file=r_data_file
  )

#plot_model_pi_densities(tests, unmarked_m)
#round(sd(s$detections))/round(mean(s$detections))
#mean(unmarked::getP(unmarked_m))-mean(pDet)

