###########################################
## Separation detection function for     ##
## mixcure.penal.est under both ML/FT-PL ##
###########################################

separation.detection <- function(obj, pl, nsteps) {

  design.matrix <- model.frame(formula = obj$formula, data = obj$data, na.action = na.omit)
  out.cure <- matrix(0, nrow = nsteps, ncol = ncol(design.matrix))
  out.surv <- matrix(0, nrow = nsteps, ncol = ncol(design.matrix))

  for (no.i in 1:nsteps) {
 mix.out <- mixcure.penal.est(formula = obj$formula, data = obj$data, init = obj$init, pl=pl, iterlim = no.i)
 if (no.i>1) {
   out.cure[no.i, ] <- c(mix.out$coefficients$cure[,3])/out.cure[1, ]
   out.surv[no.i, ] <- c(mix.out$coefficients$surv[,3])/out.surv[1, ]

 } else{
  out.cure[no.i, ] <- c(mix.out$coefficients$cure[,3])
  out.surv[no.i, ] <- c(mix.out$coefficients$surv[,3])
 }
  }
  colnames(out.cure) <- c("Intercept",colnames(design.matrix)[-1])
  colnames(out.surv) <- c("Intercept",colnames(design.matrix)[-1])

  out.tab <- list(cure = out.cure,
                  surv = out.surv)
  return(out.tab)
}
