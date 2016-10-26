#
# A first-run at analyzing PLJV-IMBCR dataset -- mining for summary statistics with A. Bartuszevige for PLJV/PIF/SGCN
# species lists.
#
# Author: Kyle Taylor [2016]
#

require(raster)
require(rgdal)
require(landscapeAnalysis)

#
# local functions
#
recursiveFindFile <- function(name=NULL,root=Sys.getenv("HOME")){
  if(is.null(name)){
    return(NULL)
  } else {
    return(list.files(root,pattern=name,recursive=T,full.names=T))
  }
}

imbcrTableToShapefile <- function(filename=NULL,outfile=NULL,write=F){
  if(is.null(outfile) && write && !is.character(filename)){
    stop("cannot write shapefile to disk without an outfile= or filename.ext to parse")
  }
  # sanity-check to see if our output shapefile already exists
  s <- recursiveFindFile(outfile)[1]
  if(!is.null(s)){
    s <- landscapeAnalysis:::readOGRfromPath(s)
  } else {
    #
    # Parse ALL BCR raw data tables into a single table
    #
    if(class(filename) == "data.frame"){
      t <- filename
    } else {
      if(file.exists(filename)){
        t <- read.csv(filename)
      } else {
        t <- recursiveFindFile(name=filename)
          t <- lapply(t,read.csv)
            t <- do.call(what=rbind,t)
      }
    }
    names(t) <- tolower(names(t))
    #
    # iterate over each UTM zone in the table, creating SpatialPoints
    # projected to a focal UTM.  Then merge all of the zones together into
    # a single shapfile with an arbitrary CRS.
    #
    s <- list()
    for (zone in unique(na.omit(t$ptvisitzone))){
      s[[length(s) + 1]] <- na.omit(t[t$ptvisitzone == zone,])
      s[[length(s)]] <- SpatialPointsDataFrame(
                                  coords = data.frame(x = s[[length(s)]]$ptvisiteasting,
                                  y = s[[length(s)]]$ptvisitnorthing),
                                  data = s[[length(s)]],
                                  proj4string = CRS(projection(paste("+init=epsg:269", zone, sep = "")))
                        )
    }

    s <- lapply(s,FUN=spTransform,CRS(projection(s[[1]])))
      s <- suppressMessages(do.call(rbind,s))
        s$FID <- 1:nrow(s)
    # write to disk -- and allow some wiggle-room on filename conventions
    if(write){
      writeOGR(s,".",ifelse(is.null(outfile),gsub(filename,pattern=".csv",replacement=""),outfile),driver="ESRI Shapefile",overwrite=T)
    }
  }
  return(s)
}

#
# PLJV Region Priority Species (Breeding Season)
#

# PIF species trend is "unknown / slightly declining"
pif_waterbirds_shorebirds_trend_3 <-
c(
  "American Coot",
  "American White Pelican",
  "American Wigeon",
  "Black Tern",
  "Black-crowned Night-Heron",
  "Blue-winged Teal",
  "Canvasback",
  "Cattle Egret",
  "Cinnamon Teal",
  "Double-crested Cormorant",
  "Eared Grebe",
  "Forster's Tern",
  "Gadwall",
  "Great Blue Heron",
  "Green-winged Teal",
  "Lesser Scaup",
  "Mallard",
  "Northern Pintail",
  "Northern Shoveler",
  "Pied-billed Grebe",
  "Redhead",
  "Ruddy Duck",
  "Western Grebe",
  "White-faced Ibis",
  "American Avocet",
  "Spotted Sandpiper",
  "Willet"
)

# PIF species trend is "unknown / slightly declining"
pif_landbirds_trend_3 <-
c(
  "Burrowing Owl",
  "Great Horned Owl",
  "Greater Prairie-Chicken",
  "Greater Roadrunner",
  "Northern Bobwhite",
  "Orchard Oriole",
  "Painted Bunting",
  "Prairie Falcon",
  "Red-headed Woodpecker",
  "Wild Turkey"
)

sgcn_tier_2_swap_waterbirds <-
c(
  "Clark's Grebe",
  "King Rail"
)

sgcn_tier_2_swap_landbirds <-
c(
  "Bobolink",
  "Common Poorwill",
  "Lewis's Woodpecker",
  "Black-and-White Warbler",
  "Black-billed Magpie",
  "Brewer's Blackbird",
  "Brown Creeper",
  "Carolina Wren",
  "Cassin's Finch",
  "Cassin's Kingbird",
  "Clark's Nutcracker",
  "Cordilleran Flycatcher",
  "Dark-eyed Junco",
  "Gray Vireo",
  "Juniper Titmouse",
  "Kentucky Warbler",
  "Lazuli Bunting",
  "Louisiana Waterthrush",
  "Merlin",
  "Olive-sided Flycatcher",
  "Pileated Woodpecker",
  "Pine Siskin",
  "Pinyon Jay",
  "Plumbeous Vireo",
  "Prothonotary Warbler",
  "Purple Martin",
  "Pygmy Nuthatch",
  "Red-tailed Hawk",
  "Ruby-throated Hummingbird",
  "Savannah Sparrow",
  "Sedge Wren",
  "Sharp-shinned Hawk",
  "Spotted Towhee",
  "Townsend's Solitaire",
  "Tufted Titmouse",
  "Violet-green Swallow",
  "Virginia's Warbler",
  "White-eyed Vireo",
  "White-throated Swift",
  "Yellow-throated Vireo"
)

sgcn_tier_1_swap_waterbirds_shorebirds <-
c(
  "Upland Sandpiper"
)

esa_pif_declining_spp_landbirds <-
c(
  "Baltimore Oriole",
  "Bell's Vireo",
  "Black-capped Vireo",
  "Bewick's Wren",
  "Brown-headed Cowbird",
  "Bullock's Oriole",
  "Chichuahuan Raven",
  "Cassin's Sparrow",
  "Common Grackle",
  "Common Nighthawk",
  "Eastern Kingbird",
  "Eastern Meadowlark",
  "Ferruginous Hawk",
  "Grasshopper Sparrow",
  "Horned Lark",
  "Lark Bunting",
  "Lark Sparrow",
  "Lesser Prairie-Chicken",
  "Loggerhead Shrike",
  "Mississippi Kite",
  "Mourning Dove",
  "Northern Mockingbird",
  "Ring-necked Pheasant",
  "Scaled Quail",
  "Scissor-tailed Flycatcher",
  "Swainson's Hawk",
  "Western Kingbird",
  "Western Meadowlark"
)

esa_pif_declining_spp_waterbirds_shorebirds <-
c(
  "Killdeer",
  "Least Tern",
  "Long-billed Curlew",
  "Mountain Plover",
  "Piping Plover",
  "Snowy Plover"
)

sgcn_tier_1_swap_landbirds <-
c(
  "McCown's Longspur",
  "Bald Eagle",
  "Barn Owl",
  "Brewer's Sparrow",
  "Chuck-will's-widow",
  "Dickcissel",
  "Golden Eagle",
  "Harris's Hawk",
  "Henslow's Sparrow",
  "Northern Harrier",
  "Peregrine Falcon",
  "Short-eared Owl",
  "Sharp-tailed Grouse",
  "Carolina Chickadee",
  "Great-tailed Grackle",
  "Red-naped Sapsucker",
  "Rufous-crowned Sparrow",
  "Summer Tanager"
)



s <- imbcrTableToShapefile(filename=recursiveFindFile(name="RawData_PLJV_IMBCR_20161024.csv",root="/home/ktaylora/Incoming")[1])
