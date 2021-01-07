#' ANN breast cancer patients data with relatively low event rate
#'
#' Data contains the time-to-observation or disease free survival time of 5000 patients
#' whose HER2 biomarker data were measured as a binary variable,
#' as well as their censoring showing a rate of ~5% recurrence of distant metastasis
#'
#' @docType data
#'
#' @usage data(ANNbcBMdat2)
#'
#' @format A list containing a matrix `ANNbcBMdat2` with the vectors of `Time` and
#' `CENS` as the outcome and a vector of biomarker `Her2`..
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(ANNbcBMdat2)
#' \donttest{
#' mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat2,init=c(0.5,-0.1,-5,1,0.1), pl=T)
#' }
"ANNbcBMdat2"
