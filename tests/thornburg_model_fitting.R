#
# accepts two arguments at runtime -- (1) a full path to the attributed
# USNG units training dataset and (2) the four-letter bird code for the
# species we are fitting our model to.
#

require(raster)
require(rgdal)
require(rgeos)

stopifnot(grepl(
    tolower(.Platform$OS.type), pattern = "unix"
  ))

system("clear")

#
# Define some useful local functions for manipulating IMBCR data
#
#' hidden function that greps for four-letter-codes
birdcode_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^bird")])
}
#' hidden function that greps for four-letter-codes
commonname_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^c.*.[.]n.*.")])
}
#' hidden function that greps for the distance field name
distance_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^rad")])
}
#' hidden function that greps for the transect field name
transect_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^tran")])
}
#' hidden function that greps for the timeperiod field name
timeperiod_fieldname <- function(df=NULL){
  return(names(df)[grepl(tolower(names(df)),pattern="^tim")])
}
#' kludging to back-fill any transect stations in an imbcr data.frame
#' that were sampled, but where a focal species wasn't observed, with
#' NA values
#' @export
scrub_imbcr_df <- function(df,
                           allow_duplicate_timeperiods=F,
                           four_letter_code=NULL,
                           back_fill_all_na=F,
                           drop_na="none"){
  # throw-out any lurking 88 values, count before start values, and
  # -1 distance observations
  df <- df[!df@data[, timeperiod_fieldname(df)] == 88, ]
  df <- df[!df@data[, timeperiod_fieldname(df)] == -1, ]
  df <- df[!df@data[, distance_fieldname(df)]   == -1, ]
  # build a dataframe for our detections
  detected <- toupper(df@data[, birdcode_fieldname(df)]) ==
    toupper(four_letter_code)
  # define a pool of potential non-detections
  not_detected <- df[!detected, ]
  not_detected@data[, distance_fieldname(df)] <- NA
  not_detected@data[, birdcode_fieldname(df)] <- toupper(four_letter_code)
  not_detected@data[, commonname_fieldname(df)] <-
    as.character(df@data[which(detected == T)[1],commonname_fieldname(df)])
  not_detected@data[, 'cl_count']             <- 0 # not used, but stay honest
  # allow a single NA value for each station, but only keep the NA values
  # if we didn't observe the bird at that point -- by default, don't allow
  # duplicate NA's within a time-period
  not_detected <- not_detected[!duplicated(not_detected@data[,
                    c(transect_fieldname(not_detected), "year", "point",
                      if (allow_duplicate_timeperiods)
                        timeperiod_fieldname(not_detected)
                      else NULL
                    )
                  ]), ]
  # allow multiple detections at stations
  detected <- df[detected, ]
  # take the merge of detections and non-duplicated, non-detections as
  # our new data.frame
  transect_heuristic <- function(x=NULL){
    x <- x@data[, c(transect_fieldname(df), 'year', 'point')]
    return(round(sqrt(as.numeric(x[,1])) + sqrt(x[,2]) + sqrt(x[,3]),5))
  }
  not_detected <- not_detected[
      !transect_heuristic(not_detected) %in%
      transect_heuristic(detected),
    ]
  df <- rbind(y=detected, x=not_detected)
  # zero-inflation fix (1) : don't drop any transects
  if(grepl(tolower(drop_na), pattern="none")){
      df <- sp:::rbind.SpatialPointsDataFrame(detected, not_detected)
  }
  # zero-inflation fix (2) : drop transects without at least one detection
  else if(grepl(tolower(drop_na),pattern="some")){
    valid_transects <- unique(detected@data[,transect_fieldname(df)])
    df <- df[df@data[,transect_fieldname(df)] %in% valid_transects,]
  }
  # zero-inflation fix (3) : drop all NA values
  else if(grep(tolower(drop_na), pattern="all")){
    df <- detected
  }
  df[order(sqrt(as.numeric(df$transectnum))+sqrt(df$year)+sqrt(df$point)),]
}
#' hidden function used to clean-up an unmarked data.frame (umdf) by dropping
#' any NA columns attributed by scrub_imbcr_df(), mean-center (scale) our site
#' covariates (but not sampling effort!), and do some optional quantile filtering
#' that drops covariates with low variance, which is a useful 'significance
#' pruning' precursor for principal components analysis. Prefer dropping the NA
#' bin here (rather than in scrub_imbcr_df), so that we still have an accurate
#' account of total sampling effort to attribute in scrub_unmarked_dataframe().
scrub_unmarked_dataframe <- function(x=NULL, normalize=T, prune_cutoff=NULL){
  row.names(x@y) <- NULL
  row.names(x@siteCovs) <- NULL
  x@y <- x@y[,!grepl(colnames(x@y), pattern="_NA")]
  x@obsToY <- matrix(x@obsToY[,1:ncol(x@y)],nrow=1)
  # do some quantile pruning of our input data, selectively dropping
  # an arbitrary number of variables based on a user-specified
  # low-variance threshold
  if(!is.null(prune_cutoff)){
    # e.g., what is the total variance for each cov across all sites?
    # drop those standardized variables with < prune_cutoff=0.05 variance
    effort_field <- ifelse(
        sum(grepl(colnames(x@siteCovs), pattern="effort")),
        "effort",
        NULL
      )
    vars_to_scale <- colnames(x@siteCovs)[
        !grepl(tolower(colnames(x@siteCovs)), pattern="effort")
      ]
    # bug-fix : only try to prune numeric variables
    is_numeric <- apply(
        x@siteCovs[1,vars_to_scale],
        MARGIN=2,
        FUN=function(x) !is.na(suppressWarnings(as.numeric(x)))
      )
    if(length(vars_to_scale)!=sum(is_numeric)){
      warning(paste(
          "the following input variables are not numeric and cannot be",
          "filtered by quantile and will not be pruned:",
          paste(
              vars_to_scale[!is_numeric],
              collapse=", "
            )
        ))
      vars_to_scale <- vars_to_scale[is_numeric]
    }
    # calculate relative variance across all sites for each variable (column)
    variance <- apply(
      x@siteCovs[,vars_to_scale],
      MARGIN=2,
      FUN=function(x) ( (x - min(x)) / (max(x)-min(x)) ) # quick min-max normalize
    )
    # min-max will return NA on no variance (e.g., divide by zero)
    variance[is.na(variance)] <- 0
    variance <- apply(
        variance,
        MARGIN=2,
        FUN=var
      )
    # drop variables that don't meet our a priori variance threshold
    dropped <- as.vector(variance < quantile(variance, p=prune_cutoff))
    if(sum(dropped)>0){
      warning(paste(
        "prune_cutoff dropped these variables due to very small variance: ",
        paste(colnames(x@siteCovs[,vars_to_scale])[dropped], collapse=", "),
        sep=""
      ))
      keep <- unique(c(
        names(is_numeric[!is_numeric]),
        effort_field,
        vars_to_scale[!dropped]
      ))
      x@siteCovs <- x@siteCovs[, keep]
    }
  }
  # normalize our site covariates?
  if(normalize){
    # don't try to normalize non-numeric values -- drop these as site covs
    x@siteCovs <-
      x@siteCovs[ , as.vector(unlist(lapply(x@siteCovs[1,], FUN=is.numeric)))]
    # don't normalize the "effort" field
    vars_to_scale <- colnames(x@siteCovs)[
        !grepl(tolower(colnames(x@siteCovs)), pattern="effort")
      ]
    # scaling call
    x@siteCovs[,vars_to_scale] <- as.data.frame(
        scale(x@siteCovs[,vars_to_scale])
      )
    # sanity check : do some variable pruning based on variance
    # from our normalization step -- drop variables with low variance
    # from consideration and report dropped variables to user
    dropped <- as.vector(unlist(lapply(x@siteCovs[1,], FUN=is.na)))
    if(sum(dropped)>0){
      warning(paste(
        "scale() dropped these variables due to very small variance: ",
        paste(colnames(x@siteCovs)[dropped], collapse=", "),
        sep=""
      ))
      x@siteCovs <- x@siteCovs[,!dropped]
    }
  }
  return(x)
}
#'
#' @export
calc_transect_effort <- function(df=NULL){
  transects <- unique(as.character(df@data[,transect_fieldname(df)]))
  for(t in transects){
    for(y in unique(df@data[df@data[, transect_fieldname(df)] == t, "year"])){
      focal <- df@data[, transect_fieldname(df)] == t & df@data[, "year"] == y
      df@data[focal, 'effort'] <- length(unique(df@data[focal,'point']))
    }
  }
  return(df)
}
#'
#' @export
calc_time_of_day <- function(df=NULL){
  df$tod <- as.numeric(df$starttime)
}
#'
#' @export
calc_day_of_year <- function(df=NULL){
  if(inherits(df,"Spatial")){
     s <- df
    df <- df@data
  }
  df$doy <- as.numeric(strftime(as.POSIXct(as.Date(as.character(
    df$date), "%m/%d/%Y")),format="%j"))
  if(exists("s")){
    s@data <- df
    return(s)
  } else {
    return(df)
  }
}
#' accepts a Spatial*DataFrame object with a distance observation
#' field name. Will attempt to arbitrarily bin the distances into
#' a number of distance bin intervals (specified by breaks). Will
#' reclass raw distances to bin identifier (e.g., distance class 3).
#' @export
calc_dist_bins <- function(df=NULL, p=0.90, breaks=10){
  if(inherits(df,"Spatial")){
     s <- df
    df <- df@data
  }
  # define our bin intervals from breaks
  if(length(breaks) == 1){
    bin_intervals <- seq(
        from=0,
        to=quantile(df[,distance_fieldname(df)], p=p, na.rm=T),
        length.out=breaks+1
      )
  } else {
    bin_intervals <- breaks
  }
  # build a distance class using the our calculated breaks
  df[,'dist_class'] <- 0
  for (j in length(bin_intervals):2){
    match <- which(df[, distance_fieldname(df)] <= bin_intervals[j])
    df[match, 'dist_class'] <- as.character(j-1)
  }
  # if we haven't matched but a radial distance was recorded, it
  # belongs in the furthest distance bin
  match <- df[,'dist_class'] == 0 & !is.na(df[, distance_fieldname(df)])
    df[match,'dist_class'] <- as.character(length(breaks)-1)
  # assume all remaining unmatched values are non-detections
  df[df[,'dist_class'] == 0, 'dist_class'] <- NA
  # return the breaks and the processed data.frame
  # back to user for inspection
  if(exists("s")){
    s@data <- df
    return(list(distance_breaks=bin_intervals,processed_data=s))
  } else {
    return(list(distance_breaks=bin_intervals,processed_data=df))
  }
}
#' hidden function that summarizes imbcr transect covariate data and metadata
#' by year (with list comprehension). This allows you to calculate covariates
#' at the IMBCR station level and then pool (summarize) the observations by
#' transect and year
pool_by_transect_year <- function(x=NULL, df=NULL, breaks=NULL, covs=NULL,
                                  summary_fun=median){
  breaks <- length(breaks)
  transect_year_summaries <- data.frame()
  # summarize focal_transect_year by breaking into counts within
  # distance classes and binding effort, year, and covs calculated
  # at the transect scale
  years <- sort(unique(df[df[,transect_fieldname(df)] == x, "year"]))
  for(year in years){
    focal_transect_year <- df[
      df[,transect_fieldname(df)] == x & df$year == year, ]
    # pre-allocate zeros for all bins
    distances <- rep(0,(breaks-1))
      names(distances) <- 1:(breaks-1)
    # build a pivot table of observed bins
    dist_classes <- sort(focal_transect_year$dist_class) # drop NA's
      dist_classes <- table(as.numeric(dist_classes))
    # merge pivot with pre-allocate table and add an NA bin
    distances[names(distances) %in% names(dist_classes)] <- dist_classes
      distances <- append(distances,
                          sum(is.na(focal_transect_year$dist_class)))
    distances <- as.data.frame(matrix(distances,nrow=1))
      names(distances) <- paste("distance_",c(1:(breaks-1),"NA"),sep="")
    # summarize each of the covs across the transect-year
    summary_covs <- matrix(rep(NA,length(covs)),nrow=1)
      colnames(summary_covs) <- covs
    for(cov in covs){
      # some covs are year-specific; filter accordingly
      cov_year <- names(focal_transect_year)[
          grepl(names(focal_transect_year),pattern=cov)
        ]
      if(length(cov_year)>1){
          cov_year <- cov_year[grepl(cov_year,pattern=as.character(year))]
        }
      summary_covs[,cov_year] <- summary_fun(
          focal_transect_year[,cov_year],
          na.rm=T
        )
    }
    # post-process pooled transect-year
    # keep most of our vars intact, but drop those that lack meaning at
    # the transect scale or that we have summarized above
    meta_vars <- colnames(df)[!colnames(df) %in%
              c(transect_fieldname(df), "year", "dist_class",
                distance_fieldname(df), "timeperiod", "point", "how",
                  "FID", "visual", "migrant", "cl_count", "cl_id",
                    "ptvisitzone", "ptvisiteasting", "ptvisitnorthing",
                      "rank", covs)]
    # build our summary transect-year data.frame
    focal_transect_year <- cbind(
        focal_transect_year[1, meta_vars],
        data.frame(transectnum=x, year=year),
        distances,
        summary_covs
      )
    # merge into our annual summary table
    transect_year_summaries <-
      rbind(transect_year_summaries,focal_transect_year)
  }
  return(transect_year_summaries)
}
#' accepts a formatted IMBCR SpatialPointsDataFrame and builds an
#' unmarkedFrameGDS data.frame that we can use for modeling with
#' the unmarked package. Will optionally add latitude and longitude
#' attributes (WGS84).
#' @export
build_unmarked_gds <- function(df=NULL,
                               numPrimary=1,
                               distance_breaks=NULL,
                               covs=NULL,
                               unitsIn="m",
                               summary_fun=median,
                               drop_na_values=T
                               ){
  if(inherits(df, "Spatial")){
    df <- df@data
  }
  # determine distance breaks / classes, if needed
  if(is.null(distance_breaks)){
    distance_breaks  = df$distance_breaks
    distance_classes = append(sort(as.numeric(unique(
                            df$dist_class))),
                            NA
                          )
  } else {
    distance_classes = append(1:length(distance_breaks)-1, NA)
  }
  # parse our imbcr data.frame into transect-level summaries
  # with unmarked::gdistsamp comprehension
  transects <- unique(df[,transect_fieldname(df)])
  # pool our transect-level observations
  transects <- do.call(rbind,
      lapply(
          transects,
          FUN=pool_by_transect_year,
          df=df, breaks=distance_breaks,
          covs=covs
        )
    )
  # bug fix : drop entries with NA values before attempting PCA or quantile pruning
  if(drop_na_values){
    transects <- transects[ !as.vector(rowSums(is.na(transects)) > 0) , ]
    transects <- transects[ , !grepl(colnames(transects), pattern="_NA")]
  }
  # build our unmarked frame and return to user
  return(unmarked::unmarkedFrameGDS(
      # distance bins
      y=transects[,grepl(names(transects),pattern="distance_")],
      # covariates that vary at the site (transect) level
      siteCovs=transects[,!grepl(colnames(transects),pattern="distance_")],
      # not used (covariates at the site-year level)
      yearlySiteCovs=NULL,
      survey="point",
      unitsIn=unitsIn,
      dist.breaks=distance_breaks,
      numPrimary=numPrimary # should be kept at 1 (no within-season visits)
    ))
}
#' pavlacky fragmentation pca
pca_reconstruction <- function(x,
                               frag_covs=NULL,
                               total_area_filter=NULL,
                               total_area_suffix="_ar$",
                               drop_total_area=T,
                               scale=T,
                               center=T,
                               test=1)
{
  which_component_max_area <- function(m=NULL, area_metric="total_area"){
    # find the eigenvector rotation maximums for each component
    return(as.vector(which.max(abs(m$rotation[area_metric,]))))
  }
  if(sum(grepl(colnames(x@siteCovs),pattern=total_area_suffix))==0){
    stop("couldn't find an area suffix in our input table")
  }
  if(is.null(frag_covs)){
    fragmentation_metrics <- x@siteCovs[ ,
        grepl(colnames(x@siteCovs), pattern="mn_p_ar$|pat_ct$|inp_dst$")
      ]
  } else {
    fragmentation_metrics <- x@siteCovs[, frag_covs]
  }
  if(sum(grepl(colnames(x@siteCovs), pattern="total_area"))==0){
    warning("no total_area metric found in the input table -- calculating ",
    "one internally. This might behave poorly if you've mean-centered your ",
    "data.")
    s@siteCovs[,'total_area'] <- calc_total_area(
        x@siteCovs,
        total_area_filter=total_area_filter,
        total_area_suffix=total_area_suffix
      )
  }
  # bug check : do we need to mean-center our total area metric?
  if(abs(mean(x@siteCovs$total_area)/sd(x@siteCovs$total_area))>1){
    warning("there was a lot of variance in the total_area metric, suggesting ",
    "that it wasn't mean-variance centered. We are going to scale it before using ",
    "it in our PCA.")
    x@siteCovs[,'total_area'] <- scale(x@siteCovs[,'total_area'])
  }

  # bind our fragmentation metrics and our new "total area"
  # metric
  fragmentation_metrics <- cbind(
      fragmentation_metrics,
      total_area=x@siteCovs[ , 'total_area']
    )
  pca <- prcomp(fragmentation_metrics, scale.=T, center=T)
  # which component captures the greatest variance for "total_area"
  total_area_component <- which_component_max_area(pca)
  keep_components <- !( 1:ncol(pca$x) %in% total_area_component )
  # drop our total area component and collapse our remaining components
  # into a single variable reconstruction representing fragmentation

  # test (1; dim reduction) : take the component that captures the greatest
  # remaining variance after dropping the 'total_area' component. This is
  # often PC2 (if total_area was represented best by PC1), but not always.
  # So we have to dig to make sure.
  if(test==1){
    # update our keep components to whatever has the next highest variance
    keep_components <- which(
        round(pca$sdev, 6) ==
        round(
          pca$sdev[keep_components][which.max(pca$sdev[keep_components])],
          6
        )
      )
    # subset the scores matrix ($x) for our single retained component
    scores_matrix <- as.matrix(pca$x[,keep_components])
    colnames(scores_matrix) <- paste("PC",keep_components,sep="")
    # drop our lurking configuration metrics
    x@siteCovs <- cbind(
        x@siteCovs,
        scores_matrix
      )
    x@siteCovs <- x@siteCovs[ ,
        !grepl(
            colnames(x@siteCovs),
            pattern="mn_p_ar$|pat_ct$|inp_dst$|total_area"
          )
      ]
    # return unmarked df and pca model to user
    return(list(x, pca))
  }
  # test (2; dim reduction) : take the cross product of our 3 retained
  # components after dropping 'total_area'; this is essentially an interaction
  # term for our three retained components
  else if(test==2){
    # subset our scores matrix for our "keeper" components
    scores_matrix <- pca$x[,keep_components]
    cross_product <- matrix(apply(
        scores_matrix,
        #pca$x[, keep_components] %*% t(pca$rotation[,keep_components]),
        MARGIN=1,
        FUN=prod
      ))
    colnames(cross_product) <- "fragmentation"
    # drop our lurking configuration metrics
    x@siteCovs <- cbind(
        x@siteCovs,
        cross_product
      )
    x@siteCovs <- x@siteCovs[ ,
        !grepl(
            colnames(x@siteCovs),
            pattern="mn_p_ar$|pat_ct$|inp_dst$|total_area"
          )
      ]
    # return unmarked df and pca model to user
    return(list(x, pca))
  }
  # test (3; reconstruction) : reconstruct total area from our three retained
  # components (after dropping the total_area component). What does this
  # represent?
  else if(test==3){
    return(NULL)
  }
  # test (4; reconstruction) : reconstruct our three fragmentation metrics
  # (after dropping the total area component) amd then re-calculate a PCA,
  # accepting the first component as representative of "fragmentation"
  else if(test==4){
    # we are going to use two PCA objects, here
    # keep a copy of both for later validation and prediction
    pca_1 <- pca
    x_hat <- pca$x[,keep_components] %*% t(pca$rotation[,keep_components])
    x_hat <- x_hat[,c("mn_p_ar","pat_ct","inp_dst")]
    x_hat <- -1*scale(
      x_hat,
      center = colMeans(x@siteCovs[,c("mn_p_ar","pat_ct","inp_dst")]),
      scale = T
    )
    if(!drop_total_area){
        total_area <- pca$x[,total_area_component] %*% t(pca$rotation[,total_area_component])
        x@siteCovs$total_area <- total_area[,'total_area']
    }
    # Re-calculate a PCA from our partial reconstruction
    pca_2 <- prcomp(x_hat, scale.=T, center=T)
    # subset the scores matrix ($x) for our single retained component
    scores_matrix <- as.matrix(pca$x[,1])
    colnames(scores_matrix) <- "PC1"
    # drop our lurking configuration metrics
    x@siteCovs <- cbind(
        x@siteCovs,
        scores_matrix
      )
    x@siteCovs <- x@siteCovs[ ,
        !grepl(
          colnames(
            x@siteCovs),
            pattern=ifelse(drop_total_area,"mn_p_ar$|pat_ct$|inp_dst$|total_area","mn_p_ar$|pat_ct$|inp_dst$")
        )
      ]
    # return unmarked df and pca model to user
    return(list(x,pca_1,pca_2))
  }
}
#' A vanilla implementation of PCA that accepts a user-specified variance
#' threshold for proportion of variance explained and drops all trailing
#' components after a threshold is met. Meant to be used on all covariates
#' considered by a model as an alternative to model selection via IC.
pca_dim_reduction <- function(x,
                              covs=NULL,
                              scale=T,
                              center=T,
                              var_threshold=0.9,
                              force=F)
{
  # find the minimum number of components needed to explain our
  # variance threshold
  find_min_variance_explained <- function(x){
    var_explained <- round(diffinv(x$sdev/sum(x$sdev)), 2)
      var_explained <- var_explained[2:length(var_explained)]
    return(min(which(var_explained >= var_threshold)))
  }
  if (is.null(covs)){
    stop("covs= argument must specify input covariates names for our PCA")
  }
  # bug-fix drop columns that have zero variance, otherwise prcomp() will
  # fail when you attempt to mean-center
  x@siteCovs[is.na(x@siteCovs)] <- 0
  no_variance <- floor(apply(
      x@siteCovs,
      MARGIN=2,
      FUN=var
    )) == 0
  x@siteCovs <- x@siteCovs[ , !no_variance]
  # by default, assume that the user has not mean-centered siteCovs
  # with scrub_unmarked_dataframe(). Duplicate re-scaling here shouldn't
  # change anything.
  if ( sum(!covs %in% colnames(x@siteCovs))>0 ){
    warnings(paste(
        "the following were not found in the source data.frame: covs=",
        paste(covs[!covs %in% colnames(x@siteCovs)], collapse=", "),
        " -- due to zero-variance or misspecification - dropping from PCA"
      ))
  }
  # subset our site-level covariates specified by the user
  covs <- covs[covs %in% colnames(x@siteCovs)]
  # fit a PCA to the site-level covs
  pca <- prcomp(x@siteCovs[,covs], scale.=scale, center=center)
  # figure out the final component to include that satisfies our
  # a priori variance threshold
  last_component <- find_min_variance_explained(pca)
  # sanity-check: can we meaningfully drop any input variables?
  if(last_component == length(covs)){
    warning(paste("we needed all of our components to satisfy the",
        " user-specified variance threshold (p=", var_threshold,"). Are",
        " you sure you want to do a dimensional reduction?",
        sep=""
      ))
    if(!force) return(NULL)
  }
  # fetch our non-focal (metadata) covs
  meta_vars <- x@siteCovs[ , !(colnames(x@siteCovs) %in% covs) ]
  # drop our original site-level covariates so that we only include those
  # that maximize the variance of our 1:n components
  covs <- names(
      pca$rotation[
        apply(
          abs(pca$rotation[,1:last_component]),
          MARGIN=2,
          FUN=which.max
        ),
        1 # columns are arbitrary here
      ]
    )
  # now re-fit a new PCA with an optimal subset of site covs taken from
  # the first n components from our initial PCA
  pca <- prcomp(x@siteCovs[,covs], scale.=scale, center=center)
  # predict() here will simply export the scores matrix ($x) for the training
  # dataset over the "keeper" components from our analysis
  x@siteCovs <- x@siteCovs[,covs]
  x@siteCovs <- as.data.frame(
      predict(
        pca,
        x@siteCovs
      ),
      stringsAsFactors = FALSE
    )
  # now bind our metadata covs back-in
  x@siteCovs <- cbind(meta_vars, x@siteCovs)
  return(list(x, pca))
}
#' testing : build a full PCA (for all covariates) and train a model to a subset
#' of components that explain some threshold of variance
quantile_pcr <- function(imbcr_df=NULL, siteCovs=NULL, threshold=0.7, K=NULL){
  if(is.null(K)){
    K <- calc_k(imbcr_df)
  }
  pca_m <- pca_dim_reduction(
      x=imbcr_df,
      covs=siteCovs,
      var_threshold=threshold,
      force=T # force a PCA, even if we can't satisfy the var_threshold
    )

  imbcr_df <- pca_m[[1]] # contains our PCA scores matrix and our model obj

  p_70_pcr_m <- unmarked::gdistsamp(
      as.formula(paste(
        "~",
        paste(colnames(pca_m[[2]]$rotation), collapse="+"),
        "+offset(log(effort))",
        sep=""
      )),
      ~1,
      as.formula(paste("~",paste(allDetCovs, collapse="+"))),
      data=imbcr_df,
      keyfun="halfnorm",
      mixture="P",
      se=T,
      K=K
    )
}
#' shorthand vector extraction function that performs a spatial join attributes
#' vector features in x with overlapping features in y. Will automatically
#' reproject to a consistent CRS.
#'
spatial_join <- function(x=NULL, y=NULL){
  over <- sp::over(
      x = sp::spTransform(
          x,
          sp::CRS(raster::projection(y))
        ),
      y = y
    )
  x@data <- cbind(x@data, over)
  return(x)
}
#' negative heurisitc filter against potential habitat variables -- dig through
#' everything in an unmarked data.frame's siteCovs slot and drop those variables
#' that we know are not habitat related. Return filtered vector to user for
#' consideration.
get_habitat_covs <- function(x){
  allHabitatCovs <- colnames(x@siteCovs)
  # heuristic -- drop anything wonky from the IMBCR dataset
  allHabitatCovs <- allHabitatCovs[!(
    allHabitatCovs %in%
    c("starttime","bcr","doy","effort","id","eightyeight","year",
      "date","stratum","observer","common.name","birdcode","sex","mgmtentity",
      "mgmtregion","mgmtunit","county","state","primaryhabitat","transectnum")
  )]
  return(allHabitatCovs)
}
#' positive heurisitc filter against potential detection variables -- dig through
#' everything in an unmarked data.frame's siteCovs slot and drop those variables
#' that we know are not detection related. Return filtered vector to user for
#' consideration.
get_detection_covs <- function(x){
  allDetCovs <- colnames(x@siteCovs)
  allDetCovs <- allDetCovs[(
    allDetCovs %in%
    c("starttime","bcr","doy")
  )]
  return(allDetCovs)
}
calc_table_summary_statistics <- function(x=NULL, vars=NULL){
  training_dataset_variable_ranges <- t(x[, vars])
  training_dataset_variable_ranges <- data.frame(
    mean=apply(training_dataset_variable_ranges, MARGIN=1, FUN=mean, na.rm=T),
    sd=apply(training_dataset_variable_ranges, MARGIN=1, FUN=sd, na.rm=T),
    min=apply(training_dataset_variable_ranges, MARGIN=1, FUN=min, na.rm=T),
    max=apply(training_dataset_variable_ranges, MARGIN=1, FUN=max, na.rm=T)
  )
  return(training_dataset_variable_ranges)
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
#' estimates a K parameter (upper-bound of integration) as a pre-cursor for
#' various count-based mixture models (e.g., poisson and negative binomial)
calc_k <- function(df=NULL, multiplier=1){
  return(floor(max(
      rowSums(unmarked:::getY(df)))*multiplier
    ))
}
#' there is some ambiguity in selecting a good upper-bounds for integration
#' in gdistamp (the K parameter). This builds an intercept model around a range
#' of K values and selects the smallest K value for which a decrease in AIC is no
#' longer observed. Might also consider checking for stability in the beta parameters
#' for this as well.
gdistsamp_find_optimal_k <- function(df=NULL, allHabitatCovs=NULL, allDetCovs=NULL, mixture=NULL, keyfun="halfnorm", K=NULL, multiple=1){
  if(is.null(K)){
    K <- unlist(lapply(
        seq(1, 2, by=0.15),
        FUN=function(x) calc_k(df, multiplier=x)
      ))
  }
  all_covs_m <- unlist(lapply(
     K,
     FUN=function(x){
       m <- try(suppressMessages(unmarked::gdistsamp(
          # abundance
          as.formula(paste(
            "~",
            paste(allHabitatCovs, collapse="+"),
            "+offset(log(effort))",
            sep=""
          )),
          ~1, # availability
          as.formula(paste("~",paste(allDetCovs, collapse="+"))), # detection
          data=imbcr_df,
          keyfun="halfnorm",
          mixture=mixture,
          unitsOut="kmsq",
          se=T,
          K=x
         )))
       if(class(m) == "try-error"){
           return(NA)
       } else {
           return(OpenIMBCR:::AIC(m))
       }
     }
  ))
  # bug-fix : did we not find a single valid K parameter?
  if(all(is.na(all_covs_m))){
    stop(
        "couldn't find convergence for any K parameters in our sequence -- this shouldn't happen. Maybe there's",
        "zero-inflation of our input data set?"
      )
  # bug-fix : did we only find one working value? Use it
  } else if( sum(!is.na(all_covs_m)) == 1 ){
    return(
      list(
        K=K[which(!is.na(all_covs_m))],
        AIC=all_covs_m[which(!is.na(all_covs_m))]
      )
    )
  }
  # bug-fix swap out our NA values with noise that
  # our second derivative won't settle on
  if(any(is.na(all_covs_m))){
    all_covs_m[is.na(all_covs_m)] <-
      which(is.na(all_covs_m)) * 2 * max(all_covs_m, na.rm=T)
  }
  # at what array position does the change in AIC start to bottom out?
  tail <- quantile(abs(round(diff(diff(all_covs_m)),3)), p=0.6) # median, but more tail-ey
  min_k <- min(which(abs(round(diff(diff(all_covs_m)),2)) < tail))
  return(
    list(
      K=K[min_k],
      AIC=all_covs_m[min_k]
    )
  )
}

check_correlation_matrix <- function(
  var=NULL,
  x_correlation_matrix=NULL,
  imbcr_df=NULL,
  correlation_threshold=0.5)
{
  correlation <- x_correlation_matrix[var,'PC1']
  # is this positive correlation? It shouldn't be
  if (correlation > 0) {
    warning(
      paste(
        "positive correlation observed between",var,"and our fragmentation metric. They should be the",
        "inverse of each other. Consider -1*PCA transformation")
    )
  }
  if (abs(correlation) > correlation_threshold){
      warning(paste("dropping",var,"-- it's strongly correlated with our fragmentation PCA"))
      imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern=var) ]
  }
  return(imbcr_df)
}
balance_zero_transects <- function(imbcr_df=NULL, multiple=2){
  nonzero_transects <- which(rowSums(imbcr_df@y)!=0)
  zero_transects <- which(rowSums(imbcr_df@y)==0)
  if(length(zero_transects)/length(nonzero_transects) > multiple){
    downsample_to <- floor(length(nonzero_transects)*multiple)
    zero_transects <- sample(zero_transects, size=downsample_to)
    # resample our input table by throwing out some of our zero transects
    imbcr_df@y <- imbcr_df@y[c(nonzero_transects,zero_transects) , ]
    imbcr_df@siteCovs <- imbcr_df@siteCovs[c(nonzero_transects,zero_transects) , ]
    imbcr_df@tlength <- rep(1, nrow(imbcr_df@y))
  }
  return(imbcr_df)
}
check_optimal_mixture_dist <- function(imbcr_df=NULL, allHabitatCovs=NULL, allDetCovs=NULL, dist=NULL, K=NULL){
  mixture <- gdistsamp_find_optimal_k(
    df=imbcr_df,
    allHabitatCovs=allHabitatCovs,
    allDetCovs=allDetCovs,
    mixture=dist,
    K=K
  )

  if(class(mixture) == "try-error"){
    imbcr_df <- balance_zero_transects(imbcr_df, multiple=0.5)
    mixture <- try(gdistsamp_find_optimal_k(
        df=imbcr_df,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        mixture=dist,
        K=K
      ))
    if(class(pois_aic)=="try-error"){
      warning("poisson:couldn't find an optimal K value, even after downsampling over-abundant zeros")
    } else {
      K_pois <- pois_aic$K
      pois_aic <- pois_aic$AIC
    }
  } else {
    K_pois <- pois_aic$K
    pois_aic <- pois_aic$AIC
  }
}
#
# MAIN
#

argv <- commandArgs(trailingOnly=T)
stopifnot(length(argv)>1)

# try and load a training dataset with spatial data we can
# join our IMBCR station points with
if(!file.exists(argv[1])){
  argv[1] <- "/global_workspace/thornburg/vector/units_attributed_training.shp"
  stopifnot(file.exists(argv[1]))
}

# what's our four-letter bird code?
if (nchar(argv[2])!=4){
  stop("expected first argument to be a four-letter bird code")
} else {
  argv[2] <- toupper(argv[2])
}

cat(" -- reading habitat training data and mean-centering\n")

units <- OpenIMBCR:::readOGRfromPath(argv[1])

# the raw IMBCR DataFrame will have a mess of variables in
# it that we don't need. Let's assume anything with the suffix
# _ar, _dst, and _ct define actual habitat variables.
habitat_covs <- colnames(units@data)
  habitat_covs <- habitat_covs[grepl(
      habitat_covs, pattern= c("_ar$|_dst$|_ct$")
    )]

# bug-fix : back-fill any units where no patches occur with NA values
units$inp_dst[units$inp_dst == 9999] <- NA

# calculate a total_area metric for a potential PCA (before mean-centering!)
# note that this isn't actually included in our habitat_covs
units$total_area <- calc_total_area(units, total_area_filter="_rd_")

# calculate some summary statistics for our habitat variables
# so that we can go back from scale() in the future when we
# project the model into novel conditions
habitat_vars_summary_statistics <- calc_table_summary_statistics(
    units@data,
    vars=habitat_covs
  )

# now mean-variance center our data
units@data[,habitat_covs] <- scale(units@data[, habitat_covs])

cat(" -- reading IMBCR data and parsing focal species observations\n")
# this returns an IMBCR SpatialPointsDataFrame
imbcr_observations <-
  scrub_imbcr_df(OpenIMBCR::imbcrTableToShapefile(
    list.files("/global_workspace/imbcr_number_crunching/",
         pattern="RawData_PLJV_IMBCR_20161201.csv$",
         recursive=T,
         full.names=T
       )[1]
    ),
    four_letter_code=argv[2],
    drop_na="none"  # keep aLL NA values
  )

# sanity-check: do we have enough observations of our bird to make a model?
if( sum(!is.na(imbcr_observations$radialdistance)) < 80){
  stop("we have fewer than 80 distance observations for ", argv[2])
}

cat(" -- calculating distance bins\n")

# define an arbitrary 10 breaks that will be used
# to construct distance bins
breaks <- append(0,as.numeric(quantile(as.numeric(
    imbcr_observations$radialdistance),
    na.rm=T,
    probs=seq(0.05,0.90,length.out=9))
  ))
imbcr_observations <- calc_dist_bins(
    imbcr_observations,
    breaks=breaks
  )[[2]]

cat(" -- calculating detection covariates\n")

imbcr_observations <- calc_day_of_year(imbcr_observations)
imbcr_observations <- calc_transect_effort(imbcr_observations)

# append detection covariate summary statistics to our
# habitat summary table so we can go-back from mean-variance
# scaling some point in the future when we go to predict()

habitat_vars_summary_statistics <- rbind(
  habitat_vars_summary_statistics,
  calc_table_summary_statistics(
    imbcr_observations@data,
    vars=c("doy","starttime","bcr","lat","lon")
  )
)

cat(" -- prepping input unmarked data.frame and performing PCA\n")

imbcr_observations@data[,c('starttime','doy','lat','lon')] <-
  scale(imbcr_observations@data[,c('starttime','doy','lat','lon')])

cat(" -- performing spatial join with our training units dataset\n")

imbcr_df <- spatial_join(
    imbcr_observations,
    units
  )

cat(
    " -- pooling IMBCR station observations -> transect and prepping for",
    "'unmarked'\n"
  )

# this will take us from IMBCR SpatialPointsDataFrame
# to an unmarkedFrameGDS so we can fit our model with
# the unmarked package.

imbcr_df <- scrub_unmarked_dataframe(
      build_unmarked_gds(
        df=imbcr_df,
        distance_breaks=breaks,
        drop_na_values=T # here we are dropping all NA values within transects and
      ),
      normalize=F,      # we already applied scale() to our input data
      prune_cutoff=0.1  # drop variables with low variance?
    )


# figure out an initial list of habitat covariates
allHabitatCovs <- get_habitat_covs(imbcr_df)

# drop the luke george version of habitat covariates
# for our initial round of testing
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="lg_")]

allDetCovs <- get_detection_covs(imbcr_df)

# explicitly define our metadata covariates
metaDataCovs <- c("effort", "id")

# drop anything lurking in the unmarked dataframe that isn't cogent
# to the PCA or modeling
imbcr_df@siteCovs <- imbcr_df@siteCovs[,c(metaDataCovs,allDetCovs,allHabitatCovs)]

# test (4) : First PCA axis after performing PCA reconstruction after dropping
# "total area" from a four-component PCA of our configuration statistics
if ( sum(grepl(colnames(imbcr_df@siteCovs), pattern=c("mn_p_ar|pat_ct|inp_dst"))) !=3 ) {
  warning("we dropped one or more of our fragmentation metrics due to ",
  "missingness while pruning the input dataset -- skipping PCA calculation")
} else {
  pca_m <- pca_reconstruction(imbcr_df, test=4, drop_total_area=F)
  imbcr_df_original <- imbcr_df
  imbcr_df <- pca_m[[1]]
}

#
# Build a correlation matrix for our retained covariates
#

x_correlation_matrix <- round(cor(imbcr_df@siteCovs),2)

#
# Check for positive correlation between SGP and MGP
#


if ( abs(x_correlation_matrix['ag_mgp_ar','ag_sgp_ar']) > 0.5){
    warning("found a relatively high level of correlation between spg and mgp area -- favoring sgp and dropping mgp")
    imbcr_df@siteCovs[,!grepl(colnames(imbcr_df@siteCovs), pattern="_sgp")]
    #imbcr_df@siteCovs[,'grass_ar'] <- rowSums(imbcr_df@siteCovs[, c('ag_sgp_ar','lg_sgp_ar')])
}

#
# Check for correlation between SGP and PC1
#
if (exists('pca_m')){
  imbcr_df <- check_correlation_matrix(var='ag_sgp_ar', x_correlation_matrix=x_correlation_matrix, imbcr_df=imbcr_df)
  imbcr_df <- check_correlation_matrix(var='ag_mgp_ar', x_correlation_matrix=x_correlation_matrix, imbcr_df=imbcr_df)
}
# drop roads from consideration
imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern="rd_ar")]

# update our detection covariates in-case they were dropped from the PCA
# due to missingness and our habitat covariates to account for the PCA
# fragmentation metric calculation

allHabitatCovs <- get_habitat_covs(imbcr_df)
allDetCovs <- get_detection_covs(imbcr_df)

#
# Determine a reasonable K from our input table and find
# minimum AIC values for both the Poisson and Negative Binomial
# that we can compare against to select an optimal mixture
# distribution
#


cat(
   " -- estimating a good 'K' parameter and building null (intercept-only) and ",
   " alternative (habitat PCA) models\n"
  )

K <- unlist(lapply(
    seq(1, 2.5, by=0.15),
    FUN=function(x) calc_k(imbcr_df, multiplier=x)
  ))

# if we have too many zeros in our dataset, our model will never
# converge -- if this happens, try randomly downsampling our number
# over zero transects to some 'multiple' of the number of non-zero
# transects

pois_aic <- try(gdistsamp_find_optimal_k(
    df=imbcr_df,
    allHabitatCovs=allHabitatCovs,
    allDetCovs=allDetCovs,
    mixture="P",
    K=K
  ))

if (class(pois_aic) == "try-error"){
    imbcr_df <- balance_zero_transects(imbcr_df, multiple=0.75)
    pois_aic <- try(gdistsamp_find_optimal_k(
        df=imbcr_df,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        mixture="P",
        K=K
      ))
    # if downsampling didn't help -- quit with an error. We don't
    # want to try to interpret a model with questionable data.
    if (class(pois_aic)=="try-error"){
      stop("poisson : couldn't find an optimal K value,",
           "even after downsampling over-abundant zeros")
    } else {
      K_pois <- pois_aic$K
      pois_aic <- pois_aic$AIC
    }
} else {
  K_pois <- pois_aic$K
  pois_aic <- pois_aic$AIC
}

negbin_aic <- try(gdistsamp_find_optimal_k(
    df=imbcr_df,
    allHabitatCovs=allHabitatCovs,
    allDetCovs=allDetCovs,
    mixture="NB",
    K=K
  ))

if (class(negbin_aic) == "try-error"){
    # this will do nothing if we've already downsampled to 'multiple'
    # for the Poisson mixture
    imbcr_df <- balance_zero_transects(imbcr_df, multiple=0.75)
    pois_aic <- try(gdistsamp_find_optimal_k(
        df=imbcr_df,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        mixture="NB",
        K=K
      ))
    # if downsampling didn't help -- quit with an error. We don't
    # want to try to interpret a model with questionable data.
    if (class(negbin_aic)=="try-error"){
      stop("poisson : couldn't find an optimal K value,",
           "even after downsampling over-abundant zeros")
    } else {
      K_negbin <- negbin_aic$K
      pois_aic <- negbin_aic$AIC
    }
} else {
  K_negbin <- negbin_aic$K
  negbin_aic <- negbin_aic$AIC
}

if(diff(c(negbin_aic, pois_aic)) > 4){
  K <- K_negbin
  mixture_dist <- "NB"
  # if we don't see a large improvement in AIC from using the
  # negative binomial, favor the simpler Poisson mixture
} else {
  K <- K_pois
  mixture_dist <- "P"
}

# test : build a model where fragmentation metrics were collapsed into a
# single covariate and all of our remaining habitat variables were not
# included in a PCA. Then use AIC to optimize variable inclusion.

all_covs_m <- unmarked::gdistsamp(
    as.formula(paste(
      "~",
      paste(allHabitatCovs, collapse="+"),
      "+offset(log(effort))",
      sep=""
    )),
    ~1,
    as.formula(paste("~",paste(allDetCovs, collapse="+"))),
    data=imbcr_df,
    keyfun="halfnorm",
    mixture=mixture_dist,
    unitsOut="kmsq",
    se=T,
    K=K
  )

#
# now for some model selection
#

# we are going to test for lat/lon with polynomials and log transformations
# explicitly, so let's drop them here

allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="lat|lon")]

model_selection_table <- OpenIMBCR:::allCombinations_dAIC(
  siteCovs=c(allHabitatCovs, "poly(lat,2)", "poly(lon,2)", "log(lat)", "log(lon)"),
  detCovs=c("doy","starttime"),
  step=100,
  umdf=imbcr_df,
  umFunction=unmarked::gdistsamp,
  mixture=mixture_dist,
  unitsOut="kmsq",
  K=K,
  se=T,
  keyfun="halfnorm",
  offset="offset(log(effort))"
)

#
# now calculate some akaike weights from our run table
#

model_selection_table$weight <- OpenIMBCR:::akaike_weights(
    model_selection_table$AIC
  )

cat("\n")
cat(" -- species:", argv[2], "\n")
cat(" -- dAIC (null - habitat):", intercept_m@AIC-all_covs_m@AIC, "\n")
cat("\n")

save(
    compress=T,
    list=c("argv",
           "habitat_vars_summary_statistics",
           "model_selection_table",
           "imbcr_df_original",
           "imbcr_df",
           "negbin_aic",
           "pois_aic",
           "mixture_dist",
           "K",
           "x_correlation_matrix",
           "intercept_m",
           "pca_m",
           "all_covs_m"),
    file=tolower(paste(
      tolower(argv[2]),
      "_imbcr_gdistsamp_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      ".rdata",
      sep=""))
    )
