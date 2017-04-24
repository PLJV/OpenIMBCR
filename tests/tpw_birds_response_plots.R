find_max_ha <- function(x){
  x <- unlist(strsplit(x,split="x"))
    x <- as.numeric(x[length(x)]) # last char will always have window size
      x <- x*(30^2)*(10^-4) # maximum area possible : square meters -> ha
  return(x)
}

inset_plot <- function(m, var){
  max <- find_max_ha(var)
    max <- if (is.na(max)) NULL else c(0, max)
  OpenIMBCR::partialPredict(m, var=var, plot=T, new=F,
              ylim = c(0, 100),
              xlim = max,
              #xlab="Total Area Trees (Ha) [window size=5,012 ha]",
              ylab="Density (birds/ha)",
              xTransform="(x*(30^2))*(10^-4)", # square meters -> ha
              yTransform="y/100" # square kilometer -> ha
              )
}

par_plot <- function(m){
  # parse our focal explanatory vars from the model object
  vars <- strsplit(as.character(m@formula),"~")
    vars <- vars[[length(vars)]] # state co-vars are always last
      vars <- gsub(unlist(strsplit(vars, split="[+]")),
                          pattern=" ", replacement="")
  print(" -- plot variables considered:")
  print(vars)
  # determine how many panels are needed in our pairs plot
  nPanels <- length(vars) + length(vars)%%2
  par(mfrow = c(nPanels / 2, 2))
  # go crazy
  lapply(vars, FUN=inset_plot, m=m)
}

files <- list.files(pattern="rdata$")
for (file in files){
  load(file)
  par_plot(m_final)
}
