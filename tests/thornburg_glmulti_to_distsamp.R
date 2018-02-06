load(commandArgs(trailingOnly = T)[1])
# consider using options(error=traceback)
options(warn = -1, error=traceback)

if (file.exists(gsub(
    r_data_file,
    pattern = "pois_glm",
    replacement = "hds_pois"
  ))){
  cat("previous workspace for this bird found in CWD (skipping)\n")
  q("no")
}

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65

#
# Local accessory functions, some of which may overload what's
# saved in the loaded Rdata file
#

plot_hn_det <- function(x=NULL, breaks=NULL){
  param <- exp(coef(x, type = "det"))
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability"
    )
  grid(); grid();
}

plot_model_pi <- function(tests = NULL, unmarked_m = NULL){
  hinge_m_pred <- density(as.vector(floor(
      predict(tests@objects[[1]], type = "response")))
    )
  hds_m_pred <- density(as.vector(floor(
      unmarked::predict(unmarked_m, type = "state")[, 1]))
    )

  xlim <- range(c(hinge_m_pred$x, hds_m_pred$x))
  ylim <- range(c(hinge_m_pred$y, hds_m_pred$y))

  dev.new();

  plot(
      hinge_m_pred,
      xlim = c(xlim[1] - (diff(xlim) * 0.1), xlim[2] + (diff(xlim) * 0.1)),
      ylim = c(ylim[1], ylim[2] + 0.01),
      col = "red",
      main = ""
    )

  lines(
      hds_m_pred,
      col = "blue"
    )

  grid(); grid();
}

quadratics_to_keep <-function(m){

     linear_terms <- grepl(names(unmarked::coef(m)), pattern = ")1")
  quadratic_terms <- grepl(names(unmarked::coef(m)), pattern = ")2")

  # test : are we negative and are we a quadratic term?
  keep <- (unmarked::coef(m) < 0) * quadratic_terms
    keep <- names(unmarked::coef(m))[keep==1]
      keep <- gsub(keep, pattern = ")2", replacement = ")")
  if(length(keep)>0){
    keep <- gsub(keep, pattern = "lam[(]", replacement = "")
      keep <- gsub(keep, pattern = "[)][)]", replacement = ")")
    return(keep)
  } else {
    return(NULL)
  }
}

quadratics_to_keep <- function(m){
  quadratic_terms <- grepl(names(coef(m)), pattern = "[)]2")
  # test : are we negative and are we a quadratic term?
  keep <- (coef(m) < 0) * quadratic_terms
  if(length(keep)==0){
    return(NULL)
  }
  quadratic_terms <- names(coef(m))[keep == 1]
  quadratic_terms <- gsub(quadratic_terms, pattern = "[)]2", replacement = ")")
  # test: are we a positive linear term and a negative quadratic
  direction_of_coeffs <- coef(m)/abs(coef(m))
  steps <- seq(2, length(direction_of_coeffs), by = 2)
  keep <- names(which(direction_of_coeffs[steps] + direction_of_coeffs[steps+1]  == 0))
  quadratic_terms <- gsub(keep, pattern = "[)]1", replacement = ")")
  # test are both our linear and quadratic terms negative? drop if so
  if (length(quadratic_terms) > 0){
    return(quadratic_terms)
  } else {
    return(NULL)
  }
}

calc_all_distsamp_combinations <- function(vars = NULL, poly_order=T){
  formulas <- gsub(OpenIMBCR:::mCombinations(
          siteCovs = paste("poly(",vars,",2,raw=T)",sep=""),
          availCovs = NULL,
          detCovs = NULL,
          offset = "offset(log(effort))")[,1],
        pattern = " ~1 ~1",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = " ",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = "[+]offset[(]log[(]effort[)][)]",
      replacement = ""
    )
  formulas <- gsub(formulas, pattern = "~", replacement = "")
  return(formulas)
}

calc_emp_dispersion_statistic <- function(x = NULL){
  observed <- sd(x)/round(mean(x))
  predicted <- rpois(n = length(x), round(mean(x)))
    predicted <- sd(predicted)/round(mean(predicted))
  return(observed/predicted)
}

glm_to_distsamp <- function(m=NULL, umdf=NULL){
  umdf@siteCovs <- m$data

  to_use <- names(coefficients(m))
    to_use <- to_use[2:length(to_use)]
  # drop any )1 or )2 poly() postfixes
  to_use <- unique(gsub(to_use, pattern="[)][0-9]", replacement=")"))

  to_use <- paste(
      "~1~",
      paste(to_use, collapse="+"),
      "+offset(log(effort))",
      sep=""
    )

  return(unmarked::distsamp(
      formula=as.formula(to_use),
      data=umdf,
      se=T,
      keyfun="halfnorm",
      unitsOut="kmsq",
      output="abund"
    ))
}

fit_distsamp <- function(lambdas=NULL, umdf=NULL){
  unmarked_models <- lapply(
      X=unmarked_models,
      FUN=function(m){
        unmarked::distsamp(
          formula=as.formula(paste("~1~",
            lambdas,
            sep=""
          )),
          data=umdf,
          se=T,
          keyfun="halfnorm",
          unitsOut="kmsq",
          output="abund"
        )
    })
}

fit_gdistsamp <- function(lambdas=NULL, umdf=NULL){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, varlist=c("umdf"))
  unmarked_models <- parallel::parLapply(
      cl=cl,
      X=lambdas,
      fun=function(m){
        unmarked::gdistsamp(
          pformula=as.formula("~1"),
          lambdaformula=as.formula(paste("~", m, sep="")),
          phiformula=as.formula("~1"),
          data=umdf,
          se=T,
          K=max(rowSums(umdf@y)),
          keyfun="halfnorm",
          unitsOut="kmsq",
          mixture="NB",
          output="abund"
        )
    })
  parallel::stopCluster(cl);
  return(unmarked_models)
}

aic_test_quadratic_terms_gdistsamp <- function(unmarked_models=NULL, original_formulas=NULL, umdf=NULL){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  # determine our run-time parameters
  quadratics <- lapply(unmarked_models, FUN=quadratics_to_keep)
        vars <- original_formulas
  # set-up our run and parallelize across our cores
  parallel::clusterExport(cl, varlist=c("umdf", "quadratics", "vars", "AIC_RESCALE_CONST", "AIC_SUBSTANTIAL_THRESHOLD"))
  vars <- parallel::parLapply(
    cl=cl,
    X=1:length(original_formulas),
    fun=function(i){
      # drop the lam() prefix
      quads <- gsub(gsub(quadratics[[i]], pattern="lambda[(]", replacement=""), pattern="[)][)]", replacement=")")
      v <- vars[i]
      if(length(quads)>0){
        # drop poly() notation from the list of all covariates using in this model
        v <- gsub(gsub(vars[i], pattern="poly[(]", replacement=""), pattern="*.[0-2].*..*", replacement="")
        v <- v[!as.vector(sapply(v, FUN=function(p){ sum(grepl(x=quads, pattern=p))>0  }))]
        # use AIC to justify our proposed quadratic terms
        for(q in quads){
          lin_var <- gsub(
            q,
            pattern=", 2,",
            replacement=", 1,"
          )
          lambda_formula <- ifelse(
            length(v)==0,
            # empty v?
            paste(
              "~",
              paste(c(lin_var, quads[!(quads %in% q)]), collapse="+"),
              "+offset(log(effort))",
              sep=""
            ),
            # valid v?
            paste(
              "~",
              paste(
                c(paste("poly(", paste(v, ", 1, raw=T)", sep=""), sep=""),
                lin_var, quads[!(quads %in% q)]),
                collapse="+"
              ),
              "+offset(log(effort))",
              sep=""
            )
          )
          m_lin_var <- OpenIMBCR:::AIC(unmarked::gdistsamp(
            pformula=as.formula("~1"),
            lambdaformula=as.formula(lambda_formula),
            phiformula=as.formula("~1"),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            mixture="NB",
            unitsOut="kmsq",
            output="abund"
          )) + AIC_RESCALE_CONST
          lambda_formula <- ifelse(
            length(v)==0,
            # empty v?
            paste(
              "~",
              paste(c(lin_var, quads),
              "offset(log(effort))",
              sep="+"),
              sep=""
            ),
            # valid v?
            paste(
              paste(
                "~",
                c(paste("poly(", paste(v, ", 1, raw=T)", sep=""), sep=""),
                lin_var, quads,
                collapse="+"
              ),
              "+offset(log(effort))",
              sep=""
            ))
          )
          m_quad_var <- OpenIMBCR:::AIC(unmarked::gdistsamp(
            pformula=as.formula("~1"),
            lambdaformula=as.formula(lambda_formula),
            phiformula=as.formula("~1"),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            mixture="NB",
            unitsOut="kmsq",
            output="abund"
          )) + AIC_RESCALE_CONST
          # if we don't improve our AIC with the quadratic by at-least 8 aic units
          # (pretty substatial support), keep the linear version
          if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD ){
            quads <- quads[!(quads %in% q)]
            v <- c(
              v,
              gsub(
                gsub(lin_var, pattern="poly[(]", replacement=""),
                pattern=", [0-9], raw.*=*.T[)]",
                replacement="")
            )
          }
        }
        v <- c(paste("poly(", paste(v, ", 1, raw=T)", sep=""), sep=""), quads)
      }
      return(v)
  })
  parallel::stopCluster(cl);
  return(vars);
}

aic_test_quadratic_terms_distsamp <- function(unmarked_model=NULL, original_formulas=NULL, umdf=NULL){
  quadratic_terms <- lapply(
      unmarked_models,
      FUN=quadratics_to_keep
    )
  # poisson version
  unmarked_models <- mapply(
    FUN=function(...){
      # drop the lam() prefix
      quadratics <- gsub(gsub(quadratics, pattern="lam[(]", replacement=""), pattern="[)][)]", replacement=")")
      # drop poly() notation from the list of all covariates using in this model
      vars <- gsub(gsub(vars, pattern="poly[(]", replacement=""), pattern="*.[0-2].*..*", replacement="")
      if(length(quadratics)>0){
        vars <- vars[!as.vector(sapply(vars, FUN=function(p){ sum(grepl(x=quadratics, pattern=p))>0  }))]
        # use AIC to justify our proposed quadratic terms
        for(q in quadratics){
          lin_var <- gsub(
            q,
            pattern=", 2,",
            replacement=", 1,"
          )

          m_lin_var <- OpenIMBCR:::AIC(unmarked::distsamp(
            formula=as.formula(paste("~1~",
                          paste(
                            c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""),lin_var,quadratics[!(quadratics %in% q)]), collapse="+"),
                          "+offset(log(effort))",
                          sep=""
            )),
            data=umdf,
            se=F,
            keyfun="halfnorm",
            unitsOut="kmsq",
            output="abund"
          ))+AIC_RESCALE_CONST
          m_quad_var <- OpenIMBCR:::AIC(unmarked::distsamp(
            formula=as.formula(paste("~1~",
                          paste(c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics), collapse="+"),
                          "+offset(log(effort))",
                          sep=""
            )),
            data=umdf,
            se=F,
            keyfun="halfnorm",
            unitsOut="kmsq",
            output="abund"
          ))+AIC_RESCALE_CONST
          # if we don't improve our AIC with the quadratic by at-least 8 aic units
          # (pretty substatial support), keep the linear version
          if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD){
            quadratics <- quadratics[!(quadratics %in% q)]
            vars <- c(vars,gsub(gsub(lin_var, pattern="poly[(]", replacement=""), pattern=", [0-9], raw.*=*.T[)]", replacement=""))
          }
        }

        vars <- c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics)
      }
      return(vars)
    },
    quadratics=quadratic_terms,
    vars=original_formulas
  )
  return(unmarked_models)
}
#
# MAIN
#

# define the vars we are going to use
vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct", "lat", "lon")

original_formulas <- unmarked_models <- calc_all_distsamp_combinations(vars)
  unmarked_models <- paste(unmarked_models, "+offset(log(effort))", sep="")

umdf <- unmarked::unmarkedFrameGDS(
    y=as.matrix(detections$y),
    siteCovs=s@data,
    dist.breaks=detections$breaks,
    numPrimary=1,
    survey="point",
    unitsIn="m"
  )

cat(" -- fitting a first round of negbin models and testing quadratic terms\n")

unmarked_models <- aic_test_quadratic_terms_gdistsamp(
  fit_gdistsamp(unmarked_models, umdf),
  original_formulas,
  umdf
)

# refit our models using the specification justified by our AIC
# quadratics test
unmarked_models <- fit_gdistsamp(
  X=lapply(unmarked_models, FUN=function(x){
    paste(paste(x, collapse="+"), "+offset(log(effort))", sep="")
  })
)

fit_gdistsamp
# how does our dispersion look?
dispersion <- sapply(
  X=unmarked_models,
  FUN=function(m){
    k <- sum(
      grepl(
        unlist(strsplit(as.character(m@formula)[3], split="[+]")),
        pattern="poly")
    )
    observed <- unmarked::getY(m@data)
    df <- (length(observed)-sum(is.na(observed)))-k
    return(round(AICcmodavg::Nmix.chisq(m)$chi.square / df, 2))
  }
)

# make a fitList
model_selection_table <- unmarked::modSel(unmarked::fitList(
  fits=unmarked_models
))

MOD_SEL_THRESHOLD <- max(which(model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD))
# select the top models (by array position) that satisfy our threshold
MOD_SEL_THRESHOLD <- as.numeric(model_selection_table@Full$model)[1:MOD_SEL_THRESHOLD]
# read from units attribute table and drop anything that isn't numeric
predict_df <- units@data

# standard model averaging from unmarked
predicted <- unmarked::predict(
  unmarked::fitList(unmarked_models[MOD_SEL_THRESHOLD]),
  type="state",
  dispersion=dispersion,
  newdata=predict_df
)



# fancy model averaging from AICcmodavg
predict_df <- AICcmodavg:::modavgPred.AICunmarkedFitDS(
    cand.set = unmarked_models[MOD_SEL_THRESHOLD],
    parm.type = "lambda",
    c.hat=2,
    newdata = predict_df
  )


# predict_df <- lapply(
#     X=1:length(aic_weights),
#     FUN=function(x){ predicted[[x]]$Predicted }
#   )

if (length(predict_df) > 1){
  # join our predictions across models into a single
  # matrix that we can apply a weight across
  predict_df <- do.call(cbind, predict_df)
  predict_df <- sapply(
      1:nrow(predict_df),
      FUN=function(i){
        weighted.mean(x=predict_df[i, ], w = aic_weights)
      }
    )
} else {
  predict_df <- unlist(predict_df)
}

predicted@data <- data.frame(pred = predicted);

cat(" -- writing to disk\n")

rgdal::writeOGR(
  predicted,
  ".",
  tolower(paste(argv[2],
      "_imbcr_hds_pois_prediction_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      sep="")
  ),
  driver="ESRI Shapefile",
  overwrite=T
)

r_data_file <- gsub(
    r_data_file,
    pattern="pois_glm",
    replacement="hds_pois"
  )

save(
    compress=T,
    list=ls(pattern="[a-z]"),
    file=r_data_file
  )
