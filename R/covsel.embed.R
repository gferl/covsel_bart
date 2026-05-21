#' covsel.embed
#'
#' Covariate selection with model-specific embedding (Step-2)
#'
#' @param covdata A data.frame containing continuous covariate values extracted at presence–absence ('pa') locations.
#' @param pa A numeric vector indicating species presences (1) and absences (0).
#' @param weights A numeric vector of weights corresponding to each value in 'pa' (same length as 'pa').
#' @param force An optional character vector specifying the name(s) of covariate(s) to be forced into the final set.
#' @param algorithms A character vector specifying the algorithm(s) to be used for the embedding procedure (options: "glm", "gam", "rf", "bart").
#' @param ncov An integer specifying the target number of covariates to include in the final set.
#' @param maxncov An integer specifying the maximum number of covariates allowed in the final set.
#' @param nthreads An integer specifying the number of cores to be used during parallel operations.
#' @param seed An integer specifying the random seed for reproducibility.
#' @param subbart Logical specifying whether an automatic BART variable selection process should start (FALSE by default).

#'
#' @return A list with three components:
#' \enumerate{
#'   \item \code{covdata}: A data.frame containing the covariates selected after the regularization, penalization, and ranking procedures.
#'   \item \code{ranks_1}: A data.frame containing the individual ranks of all covariates for each target algorithm.
#'   \item \code{ranks_2}: A data.frame containing the final average ranks of the selected covariates.
#' }
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @examples
#' library(covsel)
#' covdata<-data_covfilter
#' dim(covdata)
#' covdata_embed<-covsel.embed(covdata, pa=data_covsel$pa, algorithms=c('glm','gam','rf', "bart"), seed=12345)
#' dim(covdata_embed$covdata)
#' @export

giocovsel <- function (covdata, pa, weights = NULL, force = NULL, algorithms = c("glm", 
                                                                    "gam", "rf", "bart"), ncov = ceiling(log2(length(which(pa == 
                                                                                                                     1)))), maxncov = 12, nthreads = detectCores()/2, subbart=FALSE){
  ranks_1 <- data.frame()
  if (!is.numeric(weights)) 
    weights <- rep(1, length(pa))
  if ("glm" %in% algorithms) {
    form <- as.formula(paste0("as.factor(pa) ~ ", paste(paste0("poly(", 
                                                               names(covdata), ",2)"), collapse = " + "), 
                              "-1"))
    x <- model.matrix(form, covdata)
    mdl.glm <- suppressWarnings(cv.glmnet(x, as.factor(pa), 
                                          alpha = 0.5, weights = weights, family = "binomial", 
                                          type.measure = "deviance", parallel = TRUE))
    glm.beta <- as.data.frame(as.matrix(coef(mdl.glm, s = mdl.glm$lambda.1se)))
    glm.beta <- data.frame(covariate = row.names(glm.beta), 
                           coef = as.numeric(abs(glm.beta[, 1])))[which(glm.beta != 
                                                                          0), ][-1, ]
    if (nrow(glm.beta) < 1) {
      glm.beta <- as.data.frame(as.matrix(coef(mdl.glm, 
                                               s = mdl.glm$lambda.min)))
      glm.beta <- data.frame(covariate = row.names(glm.beta), 
                             coef = as.numeric(abs(glm.beta[, 1])))[which(glm.beta != 
                                                                            0), ][-1, ]
    }
    if (nrow(glm.beta) < 1) {
      print("No covariate selected after elastic-net regularization, skipping to next algorithm")
    }
    else {
      glm.beta <- data.frame(glm.beta[order(glm.beta$coef, 
                                            decreasing = TRUE), ], model = "glm")
      glm.beta$covariate <- stri_sub(glm.beta$covariate, 
                                     6, -6)
      glm.beta <- data.frame(setDT(glm.beta)[, .SD[which.max(coef)], 
                                             by = covariate])
      glm.beta$rank <- 1:nrow(glm.beta)
      ranks_1 <- rbind(ranks_1, glm.beta[, c("covariate", 
                                             "rank", "model")])
    }
  }
  if ("gam" %in% algorithms) {
    if (is.character(force)) {
      pointless10 <- integer(1)
      names(pointless10) <- "pointless10"
      df_force <- data.frame(covdata[, force])
      names(df_force) <- force
      pointless10 <- which(apply(df_force, 2, function(x) length(unique(x))) < 
                             10)
      if (length(pointless10) > 0) {
        form <- as.formula(paste0("pa ~ ", paste(paste0("s(", 
                                                        names(covdata)[names(covdata) != names(pointless10)], 
                                                        ",bs='cr')"), collapse = " + ")))
      }
      else {
        form <- as.formula(paste0("pa ~ ", paste(paste0("s(", 
                                                        names(covdata), ",bs='cr')"), collapse = " + ")))
      }
    }
    else {
      form <- as.formula(paste0("pa ~ ", paste(paste0("s(", 
                                                      names(covdata), ",bs='cr')"), collapse = " + ")))
    }
    mdl.gam <- suppressWarnings(mgcv::bam(form, data = cbind(covdata, 
                                                             as.factor(pa)), weights = weights, family = "binomial", 
                                          method = "fREML", select = TRUE, discrete = TRUE, 
                                          control = list(nthreads = nthreads)))
    t <- try(summary(mdl.gam), TRUE)
    if (class(t) == "try-error") {
      if (is.character(force)) {
        if (length(pointless10) > 0) {
          form <- as.formula(paste0("pa ~ ", paste(paste0("s(", 
                                                          names(covdata)[names(covdata) != names(pointless10)], 
                                                          ",bs='ts')"), collapse = " + ")))
        }
        else {
          form <- as.formula(paste0("pa ~ ", paste(paste0("s(", 
                                                          names(covdata), ",bs='ts')"), collapse = " + ")))
        }
      }
      mdl.gam <- suppressWarnings(mgcv::bam(form, data = cbind(covdata, 
                                                               as.factor(pa)), weights = weights, family = "binomial", 
                                            method = "fREML", select = TRUE, discrete = TRUE, 
                                            control = list(nthreads = nthreads)))
    }
    gam.beta <- data.frame(covariate = names(mdl.gam$model)[!names(mdl.gam$model) %in% 
                                                              c("(weights)", "pa")], summary(mdl.gam)$s.table, 
                           row.names = NULL)
    gam.beta <- gam.beta[gam.beta$p.value < 0.9, ]
    if (nrow(gam.beta) < 1) {
      print("No covariate selected after GAM (null-space penalization), skipping to next algorithm")
    }
    else {
      gam.beta <- data.frame(gam.beta[order(abs(gam.beta$Chi.sq), 
                                            decreasing = TRUE), ], rank = 1:nrow(gam.beta), 
                             model = "gam")
      ranks_1 <- rbind(ranks_1, gam.beta[, c("covariate", 
                                             "rank", "model")])
    }
  }
  if ("rf" %in% algorithms) {
    rf <- RRF(covdata, as.factor(pa), flagReg = 0)
    impRF <- rf$importance[, "MeanDecreaseGini"]
    imp <- impRF/(max(impRF))
    gamma <- 0.5
    coefReg <- (1 - gamma) + gamma * imp
    mdl.rf <- RRF(covdata, as.factor(pa), classwt = c(`0` = min(weights), 
                                                      `1` = max(weights)), coefReg = coefReg, flagReg = 1)
    rf.beta <- data.frame(covariate = row.names(mdl.rf$importance), 
                          mdl.rf$importance, row.names = NULL)
    
    
    if (nrow(rf.beta) < 1) {
      print("No covariate selected after RF (guided regularized random forest)")
    }
    else {
      rf.beta <- data.frame(rf.beta[order(rf.beta$MeanDecreaseGini, 
                                          decreasing = TRUE), ], rank = 1:nrow(rf.beta), 
                            model = "rf")
      
      rf.beta2 <- rf.beta
      rf.beta2$coefReg <- data.frame(covariate=rownames(as.data.frame(coefReg)),
                                     coefReg = as.data.frame(coefReg)$coefReg)$coefReg[match(rf.beta$covariate, data.frame(covariate=rownames(as.data.frame(coefReg)),                                                                     coefReg = as.data.frame(coefReg)$coefReg)$covariate)]
      
      rf.beta <- rbind(rf.beta2[rf.beta2$MeanDecreaseGini > 0, ], 
                       rf.beta2[rf.beta2$MeanDecreaseGini == 0, ][order(-rf.beta2[rf.beta2$MeanDecreaseGini == 0, ]$coefReg), ])[-5]
      
      ranks_1 <- rbind(ranks_1, rf.beta[, c("covariate", 
                                            "rank", "model")])
    }
  }
  
  
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  if ("bart" %in% algorithms) {
    library(embarcadero)
    mod_bart <- bart(x.train = covdata, y.train = pa, keeptrees = TRUE)
    
    ###################
    gio.varimp <- function(model) {
      
      varimps <- if (class(model) == "rbart") {
        rowMeans(model$varcount / colSums(model$varcount))
      } else {
        colMeans(model$varcount / rowSums(model$varcount))
      }
      
      df <- data.frame(
        names = names(varimps),
        varimps = as.numeric(varimps)
      )
      
      df$names <- reorder(df$names, -df$varimps)
      
      df
    }
    
    ###################
    bart.beta <- gio.varimp(mod_bart)
    colnames(bart.beta) <- c("covariate", "varimps")

    if (nrow(bart.beta) < 1) {
      print("No covariate selected after BART")
    }
    else {
      bart.beta <- data.frame(bart.beta[order(bart.beta$varimps, 
                                              decreasing = TRUE), ], rank = 1:nrow(bart.beta), 
                              model = "bart")
      
      ranks_1 <- rbind(ranks_1, bart.beta[, c("covariate", 
                                            "rank", "model")])
    }
    if(subbart==TRUE){
      varsel_bart <- variable.step(x.data = covdata, y.data = pa)
    } else {varsel_bart=NA}
    
  }

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  if (nrow(ranks_1) < 1) {
    print("No covariate selected after the embedding procedure ...")
    return(NULL)
  }
  else {
    intersect.tmp <- ranks_1[ranks_1$covariate %in% names(which(table(ranks_1$covariate) == 
                                                                  length(unique(ranks_1$model)))), ]
    intersect.tmp <- aggregate(intersect.tmp[, c("rank")], 
                               list(intersect.tmp$covariate), sum)
    colnames(intersect.tmp) <- c("covariate", "rank")
    intersect.sel <- data.frame(intersect.tmp[order(intersect.tmp$rank, 
                                                    decreasing = FALSE), ], rank.f = 1:nrow(intersect.tmp))
    union.tmp <- ranks_1[ranks_1$covariate %in% names(which(table(ranks_1$covariate) < 
                                                              length(unique(ranks_1$model)))), ]
    if (nrow(union.tmp) > 0) {
      union.tmp <- aggregate(union.tmp[, c("rank")], 
                             list(union.tmp$covariate), sum)
      colnames(union.tmp) <- c("covariate", "rank")
      union.sel.tmp <- data.frame(union.tmp[order(union.tmp$rank, 
                                                  decreasing = FALSE), ], rank.f = (max(intersect.sel$rank.f + 
                                                                                          1)):(max(intersect.sel$rank.f) + nrow(union.tmp)))
      ranks_2 <- rbind(intersect.sel, union.sel.tmp)
    }
    else {
      ranks_2 <- intersect.sel
    }
    if (ncov > maxncov) 
      ncov <- maxncov
    if (ncov > nrow(ranks_2)) 
      ncov <- nrow(ranks_2)
    ranks_2 <- ranks_2[1:ncov, ]
    if (is.character(force)) {
      tf <- force[which(!(force %in% ranks_2$covariate))]
      if (length(tf > 1)) {
        toforce <- data.frame(covariate = tf, rank = "forced", 
                              rank.f = "forced")
        ranks_2[c(nrow(ranks_2) - nrow(toforce) + 1):c(nrow(ranks_2)), 
        ] <- toforce
      }
    }
    ranks_2 <- ranks_2[, c("covariate", "rank.f")]
    covdata <- covdata[sub(".*\\.", "", unlist(ranks_2["covariate"]))]
    return(list(covdata = covdata, ranks_1 = ranks_1, ranks_2 = ranks_2, bart_subset=varsel_bart))
  }
}
