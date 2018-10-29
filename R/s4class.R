#' Class absCopyNumber
#' @description S4 class for storing absCopyNumber input and results
#' @slot data a \code{data.table} of Segmentation file containing all filtered segments
#' @slot SNV a \code{data.table} containing input single nucletide variants.
#' @slot params a list containing parameters control data filter and model construction etc.
#' @slot estimation table containing estimated results ranked by their fitting errors
#' @slot TopResult a list containing result with minimal fitting errors, include 2 \code{data.table}:
#' \item{absCN}{table containing the absolute copy number estimates of tumor cells based on purity and ploidy pair}
#' \item{absSNV}{table containing the absolute multiplicity estimates of SNVs based on purity and ploidy pair}
#' @slot origin a data.table contain original segmentation data
#' @exportClass absCopyNumber
#' @import methods
absCopyNumber =  setClass(
    Class = "absCopyNumber",
    slots = c(
        data = "data.table",
        SNV = "data.table",
        params = "list",
        estimation = "data.table",
        TopResult = "list",  # purity, ploidy, absCN, absSNV
        origin = "data.table"
        # absCN = "data.table",
        # absSNV = "data.table"
    ), prototype = list(
        SNV = data.table::data.table(),
        params = list(),
        estimation = data.table::data.table(),
        TopResult = list(),
        origin = data.table::data.table()
    )
)


