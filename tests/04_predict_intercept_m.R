require(unmarked)

pljv_area <- floor(rgeos::gArea(OpenIMBCR:::readOGRfromPath(
    "/gis_data/PLJV/PLJV_Boundary.shp"
  )) * 1e-6 )

population_size_estimates <- do.call(
    rbind,
    lapply(
    X=list.files(pattern="rdata"),
    FUN=function(f) {
      load(f);
      density <- mean(predict(intercept_m, type="lambda")[,1])
      return(data.frame(
          spp=substr(f,1,4),
          lambda=exp(coef(intercept_m)[1]),
          density=density,
          est_pop_sz=round(density*pljv_area),
          row.names=NULL
        ))
    })
  )

write.csv(
    population_size_estimates, 
    "population_size_estimates.csv"
  )
