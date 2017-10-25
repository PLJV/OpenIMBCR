require(unmarked)
require(OpenIMBCR)
require(parallel)

r_data_file <- commandArgs(trailingOnly=T)

stopifnot(file.exists(r_data_file))

#' re-fit a model specified by a formula object
refit_model <- function(formula=NULL, imbcr_df=NULL, K=NULL, mixture=NULL){
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
      mixture=mixture,
      unitsOut="kmsq",
      se=T,
      K=K
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
#'
#'
which_component_max_area <- function(m=NULL, area_metric="total_area"){
  # find the eigenvector rotation maximums for each component
  return(as.vector(which.max(abs(m$rotation[area_metric,]))))
}
#' testing: fragmentation metric calculator. Using a new input table, do a
#' partial PCA reconstruction using a previously fit PCA. Then re-calculate a
#' new PCA from the retained components and accept the first axis (PC1) as our
#' "fragmentation" metric.
calc_fragmentation_metric <- function(table=NULL,
                                      m_pca=NULL,
                                      total_area_filter=NULL,
                                      total_area_suffix="_ar$",
                                      keep_total_area=T){
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

  table <- table[ , as.vector(unlist(lapply(table[1,], FUN=is.numeric)))]
  table <- data.frame(scale(table), stringsAsFactors=F)
  scores_matrix <- predict(pca_m[[2]], newdata=table)

  # build a score matrix from our initial PCA
  total_area_component <- which_component_max_area(m=pca_m[[2]])
  keep_components <- !( 1:ncol(pca_m[[2]]$x) %in% total_area_component )

  x_hat <- scores_matrix[,keep_components] %*% t(pca_m[[2]]$rotation[,keep_components])
  x_hat <- x_hat[,c("mn_p_ar","pat_ct","inp_dst")]

  x_hat <- scale(
    x_hat,
    center = colMeans(table[,c("mn_p_ar","pat_ct","inp_dst")]),
    scale = T
  )

  # Re-calculate a PCA from our partial reconstruction
  pca <- prcomp(x_hat, scale.=T, center=T)
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
          colnames(table),
          pattern=ifelse(keep_total_area,"mn_p_ar$|pat_ct$|inp_dst$","mn_p_ar$|pat_ct$|inp_dst$|total_area")
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

akaike_predict <- function(
  mod_sel_tab=NULL, 
  train_data=NULL,
  pred_data=NULL, 
  daic_cutoff=4,
  K=NULL, 
  mixture=NULL){
  keep <- mod_sel_tab$AIC < min(mod_sel_tab$AIC) + daic_cutoff
  mod_sel_tab <- mod_sel_tab[keep,]
  models <- lapply(
    X=as.character(mod_sel_tab$formula),
    FUN=function(x){
        refit_model(
            formula = x, 
            imbcr_df = train_data, 
            K = K, 
            mixture = mixture
        )
      }
  )
  return(lapply(
      X=1:nrow(mod_sel_tab), 
      function(x) list(model=models[[x]], weight=mod_sel_tab$weight[x])
    ))
}

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

#
# MAIN
#

load(r_data_file)

# find the minimum model from the random walk table

top_model_formula <- gsub(
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

cat(" -- re-fitting our top model (selected by minimum AIC) and our series of Akaike weighted models\n")

top_model_m <- refit_model(
    top_model_formula,
    intercept_m@data,
    K=K,
    mixture=mixture_dist
  )

akaike_models_m <- akaike_predict(
    model_selection_table, 
    train_data = intercept_m@data, 
    K=K, 
    mixture = mixture_dist
  )

cat(" -- reading our input vector data containing covariates for predict()\n")

units <- OpenIMBCR:::readOGRfromPath(
    "/global_workspace/thornburg/vector/units_attributed.shp"
  )

cat(" -- calculating centroids for each USNG unit\n")
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

cat(
    " -- mean-variance centering our input dataset so it's consistent with",
    "model training data\n")


using_fragmentation_metric <- grepl(
    paste(unlist(lapply(akaike_models_m,
    FUN=function(x) as.character(x$model@formula))), collapse=" "),
    pattern="PC"
  )

if(using_fragmentation_metric){
  cat(" -- calculating a fragmentation metric from our input configuration data\n")
  units@data <- calc_fragmentation_metric(units, total_area_filter="_rd_")
  units@data <- units@data[,unlist(lapply(units@data[1,], is.numeric))]
  # to be consistent with training data, don't scale our PC
  not_pc <- !grepl(colnames(units@data),pattern="PC1")
  units@data[,not_pc] <- data.frame(
      scale(units@data[,not_pc]),
      stringsAsFactors=F
    )
} else{
  units@data <- units@data[,unlist(lapply(units@data[1,], is.numeric))]
  units@data <- data.frame(scale(units@data), stringsAsFactors=F)
}

# append a fake effort offset that assumes each unit was sampled
# across all stations
units@data$effort <- mean(top_model_m@data@siteCovs$effort)

cat(" -- predicting against top model\n")
predicted_density_top_model <- par_unmarked_predict(units, top_model_m)

if(length(akaike_models_m)>1){
    cat(" -- model averaging against akaike-weighted selection of models\n")
    predicted_density_akaike_models <- lapply(
        X=akaike_models_m,
        FUN=function(x) {
        gc() # forces us to drop an old cluster if it's lurking
        prediction <- par_unmarked_predict(units, x$model)
        return(list(prediction=prediction, weight=x$weight))
    })
    # merge our tables
    predicted_density <- akaike_weight_predictions(
      predicted_density_akaike_models,
      col=1 # we're not averaging across stderr here
    )
} else {
    predicted_density <- as.vector(predicted_density_top_model[,1])
}

# do some sanity checks and report weird predictions

if(mean(predicted_density, na.rm=T)>top_model_m@K){ # outliers dragging our PI past K?
  warning(
    "mean prediction is larger than K -- censoring, but this ",
    "probably shouldn't happen"
  )
  nonsense_filter <- median(predicted_density)*3
  predicted_density[predicted_density>nonsense_filter] <- NA
}

cat(paste(
    " -- ", 
    sum(predicted_density > (2*top_model_m@K), na.rm=T)/length(predicted_density),
    "% of our predictions were more than 2X the K upper-bounds of our integration", 
    sep=""
  ), "\n")

cat(
    " -- total number of predicted birds across the JV:", 
    sum(predicted_density, na.rm = T),
    "\n"
  )

# write our prediction to the attribute table
spp_name <- strsplit(r_data_file[1], split="_")[[1]][1]

units@data[, spp_name] <-
  as.vector(predicted_density)

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
    list=c("all_covs_m",
           "argv",
           "habitat_vars_summary_statistics",
           "imbcr_df_original",
           "intercept_m",
           "intercept_m_negbin_aic",
           "intercept_m_pois_aic",
           "K",
           "top_model_m",
           "akaike_models_m",
           "mixture_dist",
           "model_selection_table",
           "pca_m",
           "predicted_density",
           "predicted_density"),
    file=paste(
      tolower(r_data_file[1]),
      sep="")
    )
