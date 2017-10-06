require(unmarked)
require(OpenIMBCR)
require(parallel)

r_data_file <- commandArgs(trailingOnly=T)

stopifnot(file.exists(argv[1]))

#' re-fit a model specified by a formula object
fit_final_model <- function(formula=NULL, imbcr_df=NULL){
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
      keyfun="halfnorm",
      mixture="P",
      se=T,
      K=500,
    ))
}
calc_total_area <- function(table=NULL,
                            total_area_filter=NULL,
                            total_area_suffix="_ar$"){
  if(inherits(table, "Spatial")){
    table <- table@data
  } else if(inherits(table, "unmarked")){
    table <- table@siteCovs
  }
  total_area <- colnames(table)[
      grepl(colnames(table), pattern=total_area_suffix)
    ]
  if(!is.null(total_area_filter)){
    total_area <- total_area[
        !grepl(total_area, pattern=total_area_filter)
      ]
  }
  return(as.vector(apply(
      table[ , total_area],
      MARGIN=1,
      FUN=sum,
      na.rm=T
    )))
}
#' testing: fragmentation metric calculator. Using a new input table, do a
#' partial PCA reconstruction using a previously fit PCA. Then re-calculate a
#' new PCA from the retained components and accept the first axis (PC1) as our
#' "fragmentation" metric.
calc_fragmentation_metric <- function(table=NULL,
                                      m_pca=NULL,
                                      total_area_filter=NULL,
                                      total_area_suffix="_ar$"){
  if(inherits(table, "Spatial")){
    table <- table@data
  } else if(inherits(table, "unmarked")){
    table <- table@siteCovs
  }
  # calculate total area
  table$total_area <- calc_total_area(
      table,
      total_area_filter=total_area_filter,
      total_area_suffix=total_area_suffix
    )
  table <- scale(table)
  # build a score matrix from our initial PCA
  x_hat <- pca_m[[2]]$x[,keep_components] %*% t(pca_m[[2]]$rotation[,keep_components])
  x_hat <- x_hat[,c("mn_p_ar","pat_ct","inp_dst")]
  x_hat <- -1*scale(
    x_hat,
    center = colMeans(table[,c("mn_p_ar","pat_ct","inp_dst")]),
    scale = T
  )
  # Re-calculate a PCA from our partial reconstruction
  pca_2 <- prcomp(x_hat, scale.=T, center=T)
  # subset the scores matrix ($x) for our single retained component
  scores_matrix <- as.matrix(pca$x[,1])
  colnames(scores_matrix) <- "PC1"
  # drop our lurking configuration metrics
  table <- cbind(
      table,
      scores_matrix
    )
  table <- table[ ,
      !grepl(
          colnames(x@siteCovs),
          pattern="mn_p_ar$|pat_ct$|inp_dst$|total_area"
        )
    ]
  return(table)
}
#' testing : attempt to scale a dataset using the mean-variance
#' of a previously recorded dataset
recenter_input_table <- function(table=NULL, summary_table=NULL){
  vars_to_recenter <- colnames(table)
  # bug-fix don't re-center a PC score if we
  # are using fragmentation
  vars_to_recenter <- vars_to_recenter[
      !grepl(vars_to_recenter, pattern="PC")
    ]
  for(var in vars_to_recenter){
    table[,var] <- as.vector(scale(
      table[,var],
      center=summary_table[var,'mean'],
      scale=summary_table[var,'sd']
    ))
  }
  return(table)
}
#' testing: use the parallel package to predict across a large input
#' table (with unmarked) using chunking.
par_unmarked_predict <- function(run_table=NULL, m_final=NULL){

      steps <- seq(0, nrow(run_table), by=100)
  run_table <- data.frame(run_table@data)

  if(steps[length(steps)] != nrow(run_table)){
    steps <- append(steps, nrow(run_table))
  }

  cl <- parallel::makeCluster(parallel::detectCores()-1)

  parallel::clusterExport(
      cl,
      varlist=c("run_table","m_final","steps"),
      envir=environment()
    )

  predicted_density <- parLapply(
    cl=cl,
    X=1:(length(steps)-1),
    fun=function(x){
      unmarked::predict(
        m_final,
        newdata=run_table[(steps[x]+1):steps[x+1],],
        type="lambda"
      )
    }
  )
  # bind list results and return table to user
  return(do.call(rbind, predicted_density))
}

#
# MAIN
#

load(argv[1])

# find the minimum model from the random walk table

final_model_formula <- gsub(
    as.character(
      model_selection_table$formula[which.min(model_selection_table$AIC)]
    ),
    pattern=" ",
    replacement=""
  )

# bug fix : drop spatial variables?
# final_model_formula <- gsub(
#   final_model_formula, pattern="[+]lon|[+]lat", replacement=""
# )

cat(" -- re-fitting our final model (selected by minimum AIC)\n")

m_final <- fit_final_model(
    final_model_formula,
    imbcr_df
  )

cat(" -- reading our input vector data containing covariates for predict()\n")

units <- OpenIMBCR:::readOGRfromPath(
    "/global_workspace/thornburg/vector/units_attributed.shp"
  )

cat(" -- calculating centroids for each USNG unit")
centroids <- rgeos::gCentroid(
    as(
      sp::spTransform(units, sp::CRS(raster::projection("+init=epsg:4326")))
      ,'SpatialPolygons'),
    byid=T
  )@coords

colnames(centroids) <- c("lon", "lat")
units@data <- cbind(units@data, centroids)

# are the ranges of conditions we are predicting into wildly different
# than the conditions we trained our model with?

# scale our input table and then pass-off our configuration metrics
# to our pca model for re-scaling if needed

if(grepl(final_model_formula, pattern="+PC")){
  cat(" -- calculating a fragmentation metric from our input configuration data")
  units@data <- calc_fragmentation_metric(units, total_area_filter="_rd_")
}

variables_to_use <- unlist(lapply(
    colnames(units@data),
    function(x) grepl(final_model_formula, pattern=x))
  )

cat(
    " -- mean-variance centering our input dataset so it's consistent with",
    "model training data\n"
  )

# re-centering using stored mean/sd values doesn't work very well for some
# reason, but the initial centering was done using the same USNG units across
# the pljv region. Let's just center using what we've got in our units dataset
# and see if the models fit a little better
# units@data <- recenter_input_table(
#     units@data[,variables_to_use],
#     habitat_vars_summary_statistics
#   )


units@data <- units@data[,unlist(lapply(units@data[1,], is.numeric))]
units@data <- scale(units@data)

# append a fake effort offset that assumes each unit was sampled
# across all stations
units@data$effort <- mean(m_final@data@siteCovs$effort)

cat(" -- predicting against optimal model:")
predicted_density <- par_unmarked_predict(units, m_final)

if(mean(predicted_density[,1])>500){ # outliers dragging our PI past K?
  warning(
    "mean prediction is larger than K -- censoring, but this ",
    "probably shouldn't happen"
  )
  nonsense_filter <- quantile(predicted_density[,1], probs=0.99)
  predicted_density[,1][predicted_density[,1]>nonsense_filter] <- NA
}

# write our prediction to the attribute table
spp_name <- strsplit(r_data_file[1], split="_")[[1]][1]

units@data[, spp_name] <-
  as.vector(predicted_density[,1])

units@data <- as.data.frame(
    units@data[,spp_name ]
  )

colnames(units@data) <- spp_name

cat(" -- writing to disk\n")

rgdal::writeOGR(
  units,
  dsn=".",
  layer=paste(spp_name,"_pred_density_1km", sep=""),
  driver="ESRI Shapefile",
  overwrite=T
)

save(
    compress=T,
    list=c("argv","habitat_vars_summary_statistics",
           "model_selection_table",
           "imbcr_df_original",
           "imbcr_df","allHabitatCovs","intercept_m","pca_m",
           "kitchen_sink_m", "predicted_density"),
    file=paste(
      tolower(r_data_file[1]),
      sep="")
    )
