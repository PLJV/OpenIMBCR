#
# accepts two arguments at runtime -- (1) a full path to the attributed
# USNG units training dataset and (2) the four-letter bird code for the
# species we are fitting our model to.
#

require(raster)
require(rgdal)
require(rgeos)
require(OpenIMBCR)

stopifnot(grepl(
    tolower(.Platform$OS.type), pattern = "unix"
  ))

system("clear")

#' hidden local function that will do a partial reconstruction of habitat area
#' and habitat fragmentation, providing a total area metric and a fragmentation
#' metric that are uncorrelated from each other.
pca_partial_reconstruction <- function(
  x=NULL,                   # source data.frame to use for our PCA
  frag_covs=NULL,           # explicitly define frag metrics, otherwise will pick defaults
  total_area_filter=NULL,   # i.e., drop any variables that match (grep)
  total_area_suffix="_ar$", # all habitat vars have this suffix
  drop_total_area=T,        # when we are done, should we keep a metric?
  scale=T,                  # by default, should we scale our input table?
  center=F,                 # by default, should we center our input table?
  test=4                    # debug: pick an implementation (1-4)
) {
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
    # Re-calculate a PCA from (above) partial reconstruction of our fragmentation metrics
    pca_2 <- prcomp(x_hat, scale.=T, center=T)
    # subset the scores matrix ($x) for the first principal component
    scores_matrix <- as.matrix(pca_2$x[,1])
    colnames(scores_matrix) <- "PC1"
    # bind our 'fragmentation' component and then drop
    # drop our lurking configuration metrics
    x@siteCovs <- cbind(
        x@siteCovs,
        scores_matrix
      )
    x@siteCovs <- x@siteCovs[ ,
        !grepl(
          colnames(
            x@siteCovs),
            pattern=ifelse(
                drop_total_area,
                "mn_p_ar$|pat_ct$|inp_dst$|total_area",
                "mn_p_ar$|pat_ct$|inp_dst$"
            )
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
quantile_pcr <- function(imbcr_df=NULL, siteCovs=NULL, detCovs=NULL, threshold=0.7, K=NULL){
  if(is.null(K)){
    K <- OpenIMBCR:::calc_k(imbcr_df)
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
      as.formula(paste("~",paste(detCovs, collapse="+"))),
      data=imbcr_df,
      keyfun="halfnorm",
      mixture="P",
      se=T,
      K=K
    )
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
    c(
      "starttime","bcr","doy",
      "endtime","sky_st","sky_end",
      "wind_st","wind_end","temp_st",
      "temp_end","effort","id",
      "eightyeight","year","date",
      "stratum","observer","common.name",
      "birdcode","sex","mgmtentity",
      "mgmtregion","mgmtunit","county",
      "state","primaryhabitat","transectnum"
    )
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
    c(
      "starttime","bcr","doy",
      "endtime","sky_st","sky_end",
      "wind_st","wind_end","temp_st",
      "temp_end"
    )
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
#' calculate a combined total area metric from a series of input total area
#' metrics
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

cat(" -- reading IMBCR data and parsing focal species observations\n")
# this returns an IMBCR SpatialPointsDataFrame
imbcr_observations <-
  OpenIMBCR:::scrub_imbcr_df(OpenIMBCR::imbcrTableToShapefile(
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

cat(" -- reading habitat training data and mean-centering\n")

units <- OpenIMBCR:::drop_overlapping_units(OpenIMBCR:::readOGRfromPath(
    argv[1]
  ))

# the raw IMBCR DataFrame will have a mess of variables in
# it that we don't need. Let's assume anything with the suffix
# _ar, _dst, and _ct define actual habitat variables.
habitat_covs <- colnames(units@data)
  habitat_covs <- habitat_covs[grepl(
      habitat_covs, pattern= c("_ar$|_ct$")
    )]

# back-fill any units where no patches occur with NA values
# units$inp_dst[units$inp_dst == 9999] <- NA

# calculate a total_area metric for a potential PCA (before mean-centering!)
# note that this isn't actually included in our habitat_covs
#units$total_area <- calc_total_area(units, total_area_filter="_rd_")

# calculate some summary statistics for our habitat variables
# so that we can go back from scale() in the future when we
# project the model into novel conditions
habitat_vars_summary_statistics <- calc_table_summary_statistics(
    units@data,
    vars=habitat_covs
  )

# now mean-variance center our data
#units@data[,habitat_covs] <- scale(units@data[, habitat_covs])

cat(" -- calculating distance bins\n")

# define an arbitrary 10 breaks that will be used
# to construct distance bins
#breaks <- append(
#    0,
#    as.numeric(quantile(as.numeric(
#      imbcr_observations$radialdistance),
#      na.rm=T,
#      probs=seq(0.05,0.90,length.out=5)
#    )
#  ))

breaks <- round(seq(
    from=0,
    to=quantile(imbcr_observations$radialdistance, 0.90, na.rm=T),
    length.out=7
  ))

# testing : drop observations that are greater than the 0.95 quantile
cutoff <- quantile(as.numeric(imbcr_observations$radialdistance), 0.95, na.rm=T)
 if(cutoff > 300){ # luke says most detections are within 300 meters of a transect
   cat(" -- dropping observations more than",cutoff,"from a station obs\n")
   distances <- as.numeric(imbcr_observations$radialdistance)
   distances[is.na(distances)] <- 0
   imbcr_observations <- imbcr_observations[
       distances < cutoff,
     ]
 }
imbcr_observations <- OpenIMBCR:::calc_dist_bins(
    imbcr_observations,
    breaks=breaks
  )[[2]]

cat(" -- calculating detection covariates\n")

# merge-in supplemental detection covariates
detection_metadata <- read.csv(
    "/global_workspace/imbcr_number_crunching/results/PLJV_IMBCR_SiteData_byPoint_2016-2017.csv"
  )
colnames(detection_metadata) <- tolower(colnames(detection_metadata))
imbcr_observations@data <- merge(
    imbcr_observations@data,
    detection_metadata,
    by=c("transectnum","point","year")
  )

# expand our ordinal 'sky' and 'wind' variables
imbcr_observations$sky_st <- (imbcr_observations$sky_st+1)^2
imbcr_observations$sky_end <- (imbcr_observations$sky_end+1)^2
imbcr_observations$wind_st <- (imbcr_observations$wind_st+1)^2
imbcr_observations$wind_end <- (imbcr_observations$wind_end+1)^2

imbcr_observations <- OpenIMBCR::calc_day_of_year(imbcr_observations)
imbcr_observations <- OpenIMBCR::calc_transect_effort(imbcr_observations)

# append detection covariate summary statistics to our
# habitat summary table so we can go-back from mean-variance
# scaling some point in the future when we go to predict()

habitat_vars_summary_statistics <- rbind(
  habitat_vars_summary_statistics,
  calc_table_summary_statistics(
    imbcr_observations@data,
    vars=
    c(
      "starttime","bcr","doy",
      "endtime","sky_st","sky_end",
      "wind_st","wind_end","temp_st",
      "temp_end"
      #"lat","lon",
      #"lat_2","lon_2","ln_lat","ln_lon"
    )
  )
)

cat(" -- prepping input unmarked data.frame and performing PCA\n")

imbcr_observations@data[,
  c(
    'starttime',
    'endtime',
    'sky_st',
    "sky_end",
    'wind_st',
    'wind_end',
    'temp_st',
    'temp_end',
    'doy')] <-
    #'lat',
    #'lon',
    #"lat_2",
    #"lon_2",
    #"ln_lat",
    #"ln_lon")] <-
  scale(imbcr_observations@data[,
    c(
    'starttime',
    'endtime',
    'sky_st',
    "sky_end",
    'wind_st',
    'wind_end',
    'temp_st',
    'temp_end',
    'doy'
    #'lat',
    #'lon',
    #"lat_2",
    #"lon_2",
    #"ln_lat",
    # "ln_lon")])
    )])

cat(" -- performing spatial join with our training units dataset\n")

imbcr_df <- OpenIMBCR:::spatial_join(
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

imbcr_df <- OpenIMBCR:::scrub_unmarked_dataframe(
      OpenIMBCR:::build_unmarked_gds(
        df=imbcr_df,
        distance_breaks=breaks,
        drop_na_values=T # drop all data with NA values in covs (for PCA)
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
#if ( sum(grepl(colnames(imbcr_df@siteCovs), pattern=c("mn_p_ar|pat_ct|inp_dst"))) !=3 ) {
#  warning("we dropped one or more of our fragmentation metrics due to ",
#  "missingness while pruning the input dataset -- skipping PCA calculation")
#} else {
#  pca_m <- pca_partial_reconstruction(imbcr_df, test=4, drop_total_area=F)
#  imbcr_df_original <- imbcr_df
#  imbcr_df <- pca_m[[1]]
#}

#
# Build a correlation matrix for our retained covariates
#

x_correlation_matrix <- round(cor(imbcr_df@siteCovs),2)

#
# Check for positive correlation between SGP and MGP
#
if(sum(grepl(colnames(x_correlation_matrix), pattern="sgp_ar|mgp_ar")) > 1) {
  if ( abs(x_correlation_matrix['ag_mgp_ar','ag_sgp_ar']) > 0.5){
      warning("found a relatively high level of correlation between spg and mgp area -- favoring sgp and dropping mgp")
      imbcr_df@siteCovs[,!grepl(colnames(imbcr_df@siteCovs), pattern="_sgp")]
      #imbcr_df@siteCovs[,'grass_ar'] <- rowSums(imbcr_df@siteCovs[, c('ag_sgp_ar','lg_sgp_ar')])
  }
}
#
# if we have a "grass" covariate, say from NASS, check for correlation with PC1
#
#if(sum(grepl(colnames(x_correlation_matrix), pattern="grass_ar")) > 0) {
#  imbcr_df <- OpenIMBCR:::check_correlation_matrix(
#    var='grass_ar',
#    x_correlation_matrix=x_correlation_matrix,
#    imbcr_df=imbcr_df
#  )
#}
#
# Check for correlation between SGP and PC1
#
#if(sum(grepl(colnames(x_correlation_matrix), pattern="sgp_ar|mgp_ar")) > 0) {
#  if (exists('pca_m')){
#    imbcr_df <- OpenIMBCR:::check_correlation_matrix(
#      var='ag_sgp_ar',
#      x_correlation_matrix=x_correlation_matrix,
#      imbcr_df=imbcr_df
#    )
#    imbcr_df <- OpenIMBCR:::check_correlation_matrix(
#      var='ag_mgp_ar',
#      x_correlation_matrix=x_correlation_matrix,
#      imbcr_df=imbcr_df
#    )
#  }
#}
# drop roads from consideration, if it's in the input table
#imbcr_df@siteCovs <-
#  imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern="rd_ar")]

# update our detection covariates in-case they were dropped from the PCA
# due to missingness and our habitat covariates to account for the PCA
# fragmentation metric calculation

allHabitatCovs <- get_habitat_covs(imbcr_df)
allDetCovs <- get_detection_covs(imbcr_df)
# post-hoc drop some of our detection covariates to keep from overfitting
allDetCovs <- allDetCovs[
    !grepl(allDetCovs, pattern="bcr|starttime|endtime|temp_end|sky_st|wind_st")
  ]

#
# Testing : only include a subset of our spatial covariates
# build a seperate model for spatial covs that we will overlay
# later
#
#allSpatialCovs <- allHabitatCovs[grepl(allHabitatCovs, pattern="lat|lon")]
# testing : add our first order and log terms to the habitat models
# for model selection -- but not our polynomial terms. Leave those for the
# spatial models only
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="lat|lon")]

# hack -- drop superfluous fragmentation metrics
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="inp_dst|mn_p_ar")]

# hack -- calculate first and second order terms for our habitat covariates
grass_ar <- as.data.frame(matrix(poly(imbcr_df@siteCovs[ , 'grass_ar'], degree=2), ncol=2))
  imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern='grass_ar') ]
colnames(grass_ar) <- c("first_grass_ar", "second_grass_ar")
  imbcr_df@siteCovs <- cbind(imbcr_df@siteCovs, grass_ar)
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="grass_ar")]
  allHabitatCovs <- c(allHabitatCovs, c("first_grass_ar", "second_grass_ar"))

shrub_ar <- as.data.frame(matrix(poly(imbcr_df@siteCovs[ , 'shrub_ar'], degree=2), ncol=2))
  imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern='shrub_ar') ]
colnames(shrub_ar) <- c("first_shrub_ar", "second_shrub_ar")
  imbcr_df@siteCovs <- cbind(imbcr_df@siteCovs, shrub_ar)
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="shrub_ar")]
  allHabitatCovs <- c(allHabitatCovs, c("first_shrub_ar", "second_shrub_ar"))

wetland_ar <- as.data.frame(matrix(poly(imbcr_df@siteCovs[ , 'wetland_ar'], degree=2), ncol=2))
  imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern='wetland_ar') ]
colnames(wetland_ar) <- c("first_wetland_ar", "second_wetland_ar")
  imbcr_df@siteCovs <- cbind(imbcr_df@siteCovs, wetland_ar)
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="wetland_ar")]
  allHabitatCovs <- c(allHabitatCovs, c("first_wetland_ar", "second_wetland_ar"))

pat_ct <- as.data.frame(matrix(poly(imbcr_df@siteCovs[ , 'pat_ct'], degree=2), ncol=2))
  imbcr_df@siteCovs <- imbcr_df@siteCovs[ , !grepl(colnames(imbcr_df@siteCovs), pattern='pat_ct') ]
colnames(pat_ct) <- c("first_pat_ct", "second_pat_ct")
  imbcr_df@siteCovs <- cbind(imbcr_df@siteCovs, pat_ct)
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="pat_ct")]
  allHabitatCovs <- c(allHabitatCovs, c("first_pat_ct", "second_pat_ct"))

# now mean-variance center our data
#units@data[,habitat_covs] <- scale(units@data[, habitat_covs])

#  allHabitatCovs <- append(allHabitatCovs, c("ln_lat","ln_lon"))
#
# Testing : select an optimal detection function
#
key_function <- OpenIMBCR:::gdistsamp_find_optimal_key_func(
    imbcr_df,
    allDetCovs
  )
#
# Determine a reasonable K from our input table and find
# minimum AIC values for both the Poisson and Negative Binomial
# that we can compare against to select an optimal mixture
# distribution
#

cat(
   " -- estimating a good 'K' parameter and building null ",
   "(intercept-only) and alternative (habitat) models\n"
  )

K <- unlist(lapply(
    seq(1, 2, by=0.25),
    FUN=function(x) OpenIMBCR:::calc_k(imbcr_df, multiplier=x)
  ))

# if we have too many zero transects in our dataset, our model will never
# converge -- if this happens, try randomly downsampling our number
# over zero transects to some 'multiple' of the number of non-zero
# transects

# This needs some debugging because it isn't working --
# probably something to do with scope
#
# test <- OpenIMBCR:::gdistsamp_find_optimal_k_with_transect_downsampling(
#     imbcr_df,
#     allHabitatCovs,
#     allDetCovs,
#     "P",
#     K
#   )

pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
    df=imbcr_df,
    allHabitatCovs=allHabitatCovs,
    allDetCovs=allDetCovs,
    mixture="P",
    keyfun=key_function,
    K=K
  ))

if (class(pois_aic) == "try-error"){
    imbcr_df <- OpenIMBCR:::balance_zero_transects(imbcr_df, multiple=0.5)
    pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
        df=imbcr_df,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        keyfun=key_function,
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

negbin_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
    df=imbcr_df,
    allHabitatCovs=allHabitatCovs,
    allDetCovs=allDetCovs,
    keyfun=key_function,
    mixture="NB",
    K=K
  ))

if (class(negbin_aic) == "try-error"){
    # this will do nothing if we've already downsampled to 'multiple'
    # for the Poisson mixture
    imbcr_df <- OpenIMBCR:::balance_zero_transects(imbcr_df, multiple=0.75)
    pois_aic <- try(OpenIMBCR:::gdistsamp_find_optimal_k(
        df=imbcr_df,
        allHabitatCovs=allHabitatCovs,
        allDetCovs=allDetCovs,
        keyfun=key_function,
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

# test : build an over-fit (all covariates) model and an intercept
# only model for debugging


intercept_m <- unmarked::gdistsamp(
    as.formula(paste(
      "~",
      "1",
      "+offset(log(effort))",
      sep=""
    )),
    ~1,
    as.formula(paste("~",paste(allDetCovs, collapse="+"))),
    data=imbcr_df,
    keyfun=key_function,
    mixture=mixture_dist,
    unitsOut="kmsq",
    se=T,
    K=K
  )

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
    keyfun=key_function,
    mixture=mixture_dist,
    unitsOut="kmsq",
    se=T,
    K=K
  )

#
# now for some model selection
#

model_selection_table <- OpenIMBCR:::allCombinations_dAIC(
  siteCovs=allHabitatCovs,
  detCovs=allDetCovs,
  step=500,
  umdf=imbcr_df,
  umFunction=unmarked::gdistsamp,
  mixture=mixture_dist,
  unitsOut="kmsq",
  K=K,
  se=T,
  keyfun=key_function,
  offset="offset(log(effort))"
)

# testing : do our akaike weights change when we subset using only those variables
# selected for inclusion in our top models?
all_variables_within_2aic <-
  model_selection_table$AIC < min(model_selection_table$AIC)+2
all_variables_within_2aic <- as.character(
    model_selection_table[all_variables_within_2aic,]$formula
  )
all_variables_within_2aic <- unique(unlist(lapply(
    all_variables_within_2aic,
    FUN=function(x) strsplit(strsplit(x, split="~")[[1]][2], split="[+]") )
  ))

all_variables_within_2aic <- all_variables_within_2aic[
    !grepl(all_variables_within_2aic, pattern="offset")
  ]

if(length(all_variables_within_2aic)<length(allHabitatCovs)){
  cat(
      " -- some of our variables were not found in models within 2 dAIC of the",
      " top model. Dropping unimportant variables and re-calculating a model selection",
      " table\n"
    )
  model_selection_table <- OpenIMBCR:::allCombinations_dAIC(
    siteCovs=all_variables_within_2aic,
    detCovs=allDetCovs,
    step=500,
    umdf=imbcr_df,
    umFunction=unmarked::gdistsamp,
    mixture=mixture_dist,
    unitsOut="kmsq",
    K=K,
    se=T,
    keyfun=key_function,
    offset="offset(log(effort))"
  )
}


#
# now calculate some akaike weights from our run table
#

model_selection_table$weight <- OpenIMBCR:::akaike_weights(
    model_selection_table$AIC
  )

save(
    compress=T,
    list=c("argv",
           "habitat_vars_summary_statistics",
           "model_selection_table",
           #"imbcr_df_original",
           "imbcr_df",
           "negbin_aic",
           "pois_aic",
           "mixture_dist",
           "key_function",
           "K",
           "x_correlation_matrix",
           "intercept_m",
           #"pca_m",
           "all_covs_m"),
    file=tolower(paste(
      tolower(argv[2]),
      "_imbcr_gdistsamp_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      ".rdata",
      sep=""))
    )
