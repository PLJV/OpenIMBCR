require(OpenIMBCR)

#' build a data.frame with columns that are individual strata
#' and rows that are 1.) # of transects 2.) # of species detected
calc_strata_summaries <- function(imbcr_df=NULL){
  strata <- unique(as.vector(imbcr_df$stratum))
  # number of transects in strata
  num_transects <- unlist(lapply(
      X=strata,
      FUN=function(x){
        length(unique(as.vector(
            imbcr_df[imbcr_df$stratum == x , ]$transectnum
          )))
    }))
  # number of species detected in strata
  num_spp <- unlist(lapply(
    X=strata,
    FUN=function(x){
      length(unique(as.vector(
          imbcr_df[imbcr_df$stratum == x , ]$birdcode
        )))
  }))
  # build a summary dataframe
  imbcr_df <- rbind(num_transects, num_spp)
  colnames(imbcr_df) <- strata
  rownames(imbcr_df) <- c("number of transects", "species detected")
  return(data.frame(imbcr_df))
}
#' build a data.frame with number of detections for birds, organized
#' in rows for all four-letter bird codes and columns for each stratum
calc_all_bird_detections_by_stratum <- function(imbcr_df=NULL){
  all_birds_observed <- unique(as.vector(
      imbcr_df[imbcr_df$radialdistance >= 0, ]$birdcode
    ))
  strata <- unique(as.vector(imbcr_df$stratum))
  # bird detections by stratum
  num_detections <- do.call(cbind, lapply(
      X=strata,
      FUN=function(x){
        focal <- imbcr_df[imbcr_df$stratum == x , ]
        focal <- data.frame(table(focal@data[,"birdcode"]))
        focal <- focal[as.character(focal[,1]) %in% all_birds_observed,]
    }))
  # subset our table of frequencies and drop duplicate name columns
  bird_codes <- num_detections[,1]
  num_detections <- num_detections[,1:ncol(num_detections) %% 2 == 0]
  rownames(num_detections) <- bird_codes
  colnames(num_detections) <- strata
  return(num_detections)
}

#
# MAIN
#

# accepts a single argument (two-letter code for state)
argv <- toupper(commandArgs(trailingOnly=T))


pljv_transects <- OpenIMBCR::imbcrTableToShapefile(
    filename=list.files(
        pattern="RawData", recursive=T, full.names=T
      )[1]
  )

# subset to our focal state
pljv_transects <- pljv_transects[pljv_transects$state == argv, ]

write.csv(
    x=calc_strata_summaries(pljv_transects),
    file=paste(argv,"_strata_summaries.csv", sep="")
  )

write.csv(
    x=calc_all_bird_detections_by_stratum(pljv_transects),
    file=paste(argv,"_bird_detections_by_strata.csv", sep="")
  )

