require(raster)
require(rgdal)
require(rgeos)

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
                           drop_all_na=F){
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
  # zero-inflation fix (1) : drop transects without at least one detection
  if(!back_fill_all_na){
    valid_transects <- unique(detected@data[,transect_fieldname(df)])
    df <- df[df@data[,transect_fieldname(df)] %in% valid_transects,]
  }
  # zero-inflation fix (2) : drop all NA values
  if(drop_all_na){
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
#'
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
        to=quantile(df[,distance_fieldname(df)], p=0.90, na.rm=T),
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
#' accepts a formatted IMBCR data.frame and builds an unmarkedFrameGDS
#' data.frame from it. Will add latitude and longitude attributes (WGS84)
#' @export
build_unmarked_gds <- function(df=NULL,
                               numPrimary=1,
                               distance_breaks=NULL,
                               covs=NULL,
                               summary_fun=median
                               ){
  # if we have a SpatialPointsDataFrame, calculate spatial covariates
  # and tag them on to our data.frame for our summary calculations
  if(inherits(df,"Spatial")){
     s <- df
    df <- df@data
    # calculate lat/lon covariates in WGS84
    coords <- spTransform(s,"+init=epsg:4326")@coords
      colnames(coords) <- c("lon","lat")
    df <- cbind(df,coords)
      rm(coords)
    covs <- append(covs,c("lon","lat"))
  }
  # determine distance breaks / classes, if needed
  if(is.null(distance_breaks)){
    distance_breaks  = df$distance_breaks
    distance_classes = append(sort(as.numeric(unique(
                            df$processed_data$dist_class))),
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
  # build our unmarked frame and return to user
  return(unmarked::unmarkedFrameGDS(
      # distance bins
      y=transects[,grepl(names(transects),pattern="distance_")],
      # covariates that vary at the site (transect) level
      siteCovs=transects[,!grepl(colnames(transects),pattern="distance_")],
      # not used (covariates at the site-year level)
      yearlySiteCovs=NULL,
      survey="point",
      unitsIn="m",
      dist.breaks=distance_breaks,
      numPrimary=numPrimary # should be kept at 1 (no within-season visits)
    ))
}
#'
#'
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
  # extract a projection matrix representing our pca loadings (rotations)
  # for n components
  extractProjection <- function(n, princ) {
    # pull off the rotation.
    proj <- princ$rotation[,1:n]
    # sign was arbitrary, so flip in convenient form
    for(i in seq_len(n)) {
      si <- sign(mean(proj[,i]))
      if(si!=0) {
        proj[,i] <- proj[,i]*si
      }
    }
    return(proj)
  }
  if (is.null(covs)){
    stop("covs= argument must specify input covariates names for our PCA")
  }
  # by default, assume that the user has not mean-centered siteCovs
  # with scrub_unmarked_dataframe(). Duplicate re-scaling here shouldn't
  # change anything.
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
  } else {
    # fetch our non-focal (metadata) covs
    meta_vars <- x@siteCovs[ , !(colnames(x@siteCovs) %in% covs) ]
    # drop our original covariates to only those included in the PCA
    covs <- covs[1:last_component]
    pca <- prcomp(x@siteCovs[,covs], scale.=scale, center=center)
    # predict() here will export the PCA projection matrix for the training
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

#
# MAIN
#
# Accepts two arguments at runtime -- (1) a full path to the attributed
# USNG units training dataset and (2) the four-letter bird code for the
# species we are fitting our model to.
#

argv <- commandArgs(trailingOnly=T)
stopifnot(length(argv)>1)

if(!file.exists(argv[1])){
  argv[1] <- "/global_workspace/thornburg/vector/units_attributed_training.shp"
  stopifnot(file.exists(argv[1]))
}

cat(" -- reading habitat training data and mean-centering\n")

units <- OpenIMBCR:::readOGRfromPath(argv[1])

habitat_covs <- colnames(units@data)
  habitat_covs <- habitat_covs[grepl(
      habitat_covs, pattern= c("_ar$|_dst$|_ct$")
    )]

units@data[,habitat_covs] <- scale(units@data[, habitat_covs])

if(nchar(argv[2])!=4){
  stop("expected first argument to be a four-letter bird code")
} else {
  argv[2] <- toupper(argv[2])
}

cat(" -- reading IMBCR data and parsing focal species observations\n")

imbcr_observations <-
  scrub_imbcr_df(OpenIMBCR::imbcrTableToShapefile(
    list.files("/global_workspace/imbcr_number_crunching/",
         pattern="RawData_PLJV_IMBCR_20161201.csv$",
         recursive=T,
         full.names=T
       )[1]
    ),
    four_letter_code=argv[2],
    back_fill_all_na=F  # keep only NA values for transects with >= 1 spp det
    #back_fill_all_na=T # keep all NA values
  )

cat(" -- calculating distance bins\n")

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

cat(" -- performing spatial join with our training units dataset\n")

imbcr_df <- spatial_join(
    imbcr_observations,
    units
  )

cat(" -- pooling IMBCR station observations -> transect and",
    "prepping for 'unmarked'\n")

imbcr_df <- scrub_unmarked_dataframe(
      build_unmarked_gds(
        df=imbcr_df,
        distance_breaks=breaks
      ),
      normalize=F,      # we already applied scale() to our input data
      prune_cutoff=0.10 # drop variables with low variance
    )

cat(" -- mean-centering our spatial covariates and detection covariates\n")

imbcr_df@siteCovs[,c('lat','lon','starttime','doy')] <- scale(
    imbcr_df@siteCovs[,c('lat','lon','starttime','doy')]
  )

cat(" -- prepping input unmarked data.frame and performing PCA\n")

allHabitatCovs <- colnames(imbcr_df@siteCovs)
allHabitatCovs <- allHabitatCovs[!(
  allHabitatCovs %in%
  c("starttime","bcr","doy","effort","id","eightyeight","year",
    "date","stratum","observer","common.name","birdcode","sex","mgmtentity",
    "mgmtregion","mgmtunit","county","state","primaryhabitat","transectnum")
)]

# drop the luke george version of habitat covariates
# for our initial round of testing
allHabitatCovs <- allHabitatCovs[!grepl(allHabitatCovs, pattern="lg_")]

allDetCovs <- colnames(imbcr_df@siteCovs)
allDetCovs <- allDetCovs[(
  allDetCovs %in%
  c("starttime","bcr","doy")
)]

# explicitly define our metadata covariates
metaDataCovs <- c("effort", "id")

# drop anything lurking in the unmarked dataframe that isn't cogent
# to the PCA or modeling
imbcr_df@siteCovs <- imbcr_df@siteCovs[,c(metaDataCovs,allDetCovs,allHabitatCovs)]

# this is a traditional PCA reduction where all habitat covariates
# are considered as components (not just fragmentation)

pca_m <- pca_dim_reduction(
    x=imbcr_df,
    covs=allHabitatCovs,
    var_threshold=0.7
  )

imbcr_df <- pca_m[[1]] # contains our PCA projection matrix and our model obj

cat(" -- building null (intercept-only) and alternative (habitat PCA) models\n")

intercept_m <- unmarked::gdistsamp(
    ~1+offset(log(effort)), # abundance
    ~1,                     # availability
    as.formula(paste("~",paste(allDetCovs, collapse="+"))), # detection
    data=imbcr_df,
    keyfun="halfnorm",
    mixture="P",
    se=T,
    K=500,
  )

# poly_space_m <- unmarked::gdistsamp(
#     ~poly(lat,3)+poly(lon,3)+offset(log(effort)),
#     ~1,
#     as.formula(paste("~",paste(allDetCovs, collapse="+"))),
#     data=imbcr_df,
#     keyfun="halfnorm",
#     mixture="P",
#     se=T,
#     K=500,
#   )

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
    K=500,
  )

cat("\n")
cat(" -- species:", argv[2], "\n")
cat(" -- dAIC (null - habitat):", intercept_m@AIC-p_70_pcr_m@AIC, "\n")
cat("\n")

save(
    compress=T,
    list=c("argv","imbcr_df","intercept_m","pca_m","p_70_pcr_m"),
    file=paste(
      tolower(argv[2]),
      "_imbcr_gdistsamp_workflow_",
      paste(tolower(unlist(strsplit(date(), split=" "))[c(2,3)]), collapse="_"),
      ".rdata",
      sep=""
    )
  )
