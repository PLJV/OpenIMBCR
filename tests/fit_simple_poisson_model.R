options(warn = -1, error=traceback)

require(OpenIMBCR)
require(glmulti)
require(unmarked)

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65

calc_transect_summary_detections <- function(s=NULL, name=NULL, field='est_abund'){
  return(OpenIMBCR:::calc_route_centroids(
      s = s,
      four_letter_code = toupper(name),
      use=field
    ))
}

quadratics_to_keep <- function(m){
  quadratic_terms <- grepl(names(coefficients(m)), pattern = ")2")
  # test : are we negative and are we a quadratic term?
  keep <- (coefficients(m) < 0) * quadratic_terms
  quadratic_terms <- names(coefficients(m))[keep == 1]
    quadratic_terms <- gsub(quadratic_terms, pattern = ")2", replacement = ")")
  # test: are we a positive linear term and a negative quadratic 
  direction_of_coeffs <- coefficients(m)/abs(coefficients(m))
  direction_of_coeffs <- direction_of_coeffs[ gsub(names(direction_of_coeffs), pattern="[)].", replacement=")") %in% quadratic_terms ]
  steps <- seq(1, length(direction_of_coeffs), by = 2)
  keep <- names(direction_of_coeffs[ which ( direction_of_coeffs[steps] + direction_of_coeffs[steps+1]  == 0 ) ])
    quadratic_terms <- gsub(keep, pattern = ")1", replacement = ")")
  # test are both our linear and quadratic terms negative? drop if so 
  if (length(quadratic_terms) > 0){
    return(quadratic_terms)
  } else {
    return(NULL)
  }
}

pred_hn_det_from_distance <- function(x=NULL, dist=NULL){
  param <- exp(coef(x, type = "det"))
  return(as.vector(unmarked:::gxhn(x=dist, param)))
}

#
# MAIN
#

# Read-in our IMBCR transect data

argv <- commandArgs(trailingOnly = T)

cat(" -- fitting a model for :", argv[2], "\n")

r_data_file <- tolower(paste(
      tolower(argv[2]),
      "_imbcr_pois_glm_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
      ".rdata",
      sep = ""
    ))

# Do we have lurking output in the CWD?
if (file.exists(r_data_file)) {
   cat(" -- found existing rdata file:", r_data_file, "(skipping)\n");
   stop();
}

s <- OpenIMBCR:::scrub_imbcr_df(
     OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
     ),
     four_letter_code = toupper(argv[2])
  )

if( sum(s$radialdistance>0, na.rm=T) < 180){
  cat(" -- not enough detections to fit a meaningful model for", argv[2],"(quitting)\n")
  q("no")
}

detections <- OpenIMBCR:::calc_dist_bins(s)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))

# fit an intercept-only detection function in unmarked

umdf <- unmarked::unmarkedFrameDS(
    y=as.matrix(detections$y), 
    siteCovs=data.frame(effort=effort), 
    dist.breaks=detections$breaks, 
    survey="point", 
    unitsIn="m"
  )

intercept_m <- unmarked::distsamp(
    formula = ~1 ~1+offset(log(effort)), 
    data = umdf, 
    se = T, 
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "abund"
  )

per_obs_det_probabilities <- round(sapply(
    s$radialdistance, 
    function(x) pred_hn_det_from_distance(intercept_m, dist=x)), 
    2
  )

# don't allow a 0.0 detection probability
if(min(per_obs_det_probabilities, na.rm=T)==0){
  min <- as.numeric(quantile(per_obs_det_probabilities, na.rm=T, p=seq(0,0.2,by = 0.001)))
    min <- min(min[min > 0])
  per_obs_det_probabilities[per_obs_det_probabilities == 0] <- min
}


s$est_abund <- round(s$cl_count / per_obs_det_probabilities)

s <- calc_transect_summary_detections(
    s=s, 
    name=toupper(argv[2]), 
    field='est_abund'
  )

s$effort <- effort

# read-in habitat covariates
units <- OpenIMBCR:::readOGRfromPath(argv[1])

# spatial join of transect detection data with habitat covariates
# summarized by-unit
s <- OpenIMBCR:::spatial_join(s, units)

# add latitude and longitude
cat(" -- calculating spatial covariates\n")
 
    s <- OpenIMBCR:::calc_lat_lon(s)
units <- OpenIMBCR:::calc_lat_lon(units)

# define the covariates we are going to use in our analysis

vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct","lat","lon")

# drop strongly-correlated variables
x_cor_matrix <- cor(s@data[,vars])

cor_threshold <- abs(x_cor_matrix)>CORRELATION_THRESHOLD

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
vars <- unlist(strsplit(readline("enter a comma-separated list of vars you want to keep:"), split=","))
} else {
  vars <- passing
}



# fit our full model
m <- glm(
    formula=paste("est_abund~", paste(paste("poly(",vars,",2)",sep=""), collapse="+",sep=""),"+offset(log(effort))",sep=""),
    family=poisson,
    data=s@data
  )

# drop un-intuitive quadratics
quadratics <- quadratics_to_keep(m)

if(length(quadratics)>0){
  vars <- vars[!as.vector(sapply(vars, FUN=function(p){ sum(grepl(x=quadratics, pattern=p))>0  }))]
  # use AIC to justify our proposed quadratic terms
  for(q in quadratics){
    lin_var <- gsub(
        q,
        pattern=", 2[)]",
        replacement=", 1)"
      )

    m_lin_var <- AIC(glm(
        formula=paste("est_abund~",
          paste(
            c(paste("poly(", paste(vars[1:2], ", 1)", sep=""), sep=""),lin_var,quadratics[!(quadratics %in% q)]), collapse="+"),
            "+offset(log(effort))",
            sep=""
          ),
        family=poisson,
        data=s@data
      ))+AIC_RESCALE_CONST
    m_quad_var <- AIC(glm(
        formula=paste("est_abund~",
            paste(c(paste("poly(", paste(vars[1:2], ", 1)", sep=""), sep=""), quadratics), collapse="+"),
            "+offset(log(effort))",
            sep=""
          ),
        family=poisson,
        data=s@data
      ))+AIC_RESCALE_CONST
    # if we don't improve our AIC with the quadratic by at-least 8 aic units
    # (pretty substatial support), keep the linear version
    if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD){
      quadratics <- quadratics[!(quadratics %in% q)]
      vars <- c(vars,gsub(gsub(lin_var, pattern="poly[(]", replacement=""), pattern=", 1[)]", replacement=""))
    }
  }

  vars <- c(paste("poly(", paste(vars, ", 1)", sep=""), sep=""), quadratics)

  # refit our full model minus un-intuitive quadratics
  m <- glm(
      formula=paste("est_abund~", paste(vars, collapse="+"), "+offset(log(effort))", sep=""),
      family=poisson,
      data=s@data
    )
# we are not keeping any of our quadratics? still scale with poly()
} else {
  vars <- paste("poly(", paste(vars, ", 1)", sep=""), sep="")
}

# use model selection with interactions across our candidate variables
tests <- glmulti::glmulti(m, intercept=T, family=poisson, level=1)

# predict across our full run dataset
predicted <- units

vals <- as.vector(predict(
      tests,
      select=AIC_SUBSTANTIAL_THRESHOLD,
      newdata=units@data,
      type="response"))

# average across all models within a dAIC threshold of the top model
if(class(vals)!="numeric"){
  predicted@data <- data.frame(
      pred=as.vector(floor(predict(
        tests,
        select=AIC_SUBSTANTIAL_THRESHOLD,
        newdata=units@data,
        type="response")$averages)
      )
    )
# if there was only one top model, averaging won't work
} else {
  predicted@data <- data.frame(
    pred=as.vector(floor(predict(
      tests@objects[[1]],
      newdata=units@data,
      type="response"))
    )
  )
}

rm(vals)

# censor any predictions greater than K (max)
cat(
    " -- number of sites with prediction greater than predicted max(K):",
    sum(predicted$pred > max(round(predict(m, type="response")))),
    "\n"
  )

predicted$pred[( predicted$pred > max(round(predict(m, type="response"))) )] <- max(s$est_abund)
predicted$pred[predicted$pred<1] <- 0

rgdal::writeOGR(
  predicted,
  ".",
  tolower(paste(argv[2],
      "_imbcr_pois_glm_prediction_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      sep="")
  ),
  driver="ESRI Shapefile",
  overwrite=T
)

save(
    compress=T,
    list=ls(),
    file=r_data_file
  )

