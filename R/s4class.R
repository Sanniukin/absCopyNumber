#' Class absCopyNumber
#' @description S4 class for storing input and results
#' @slot data data.table of Segmentation file containing all filtered segments
#' @slot snv.data table containing all filtered single nucletide variants
#' @slot params a list containing parameters control data filter and model construction etc.
#' @slot origin a list contain original segmentation and snv data
#' @slot summary table with basic input data stats
#' @slot result table containing estimated results ranked by their fitting errors
#' @slot absCN table containing the absolute copy number estimates of tumor cells for the top first solution (purity and ploidy pair)
#' @slot absSNV table containing the absolute multiplicity estimates of SNVs for the top first solution (purity and ploidy pair)
#'
#' @exportClass absCopyNumber
#' @import methods
absCopyNumber =  setClass(
    Class = "absCopyNumber",
    slots = c(
        data = "data.table",
        snv.data = "data.table",
        params = "list",
        origin = "list",
        summary = "data.table",
        result = "data.table",
        absCN = "data.table",
        absSNV = "data.table"
    )
)

setMethod(
    f = "show",
    signature = "absCopyNumber",
    definition = function(object) {
        cat(paste("An object of class ", class(object), "\nSummary information:\n"))
        print(object@summary)
    }
)
