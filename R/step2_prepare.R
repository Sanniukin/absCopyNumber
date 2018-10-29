


# prepare arguments
abs_prepare = function(absCopyNumber,
                       alpha.min = 0.2,
                       alpha.max = 1.0,
                       tau.min = 0.5,
                       tau.max = 8.0,
                       copyratio.min = -3,
                       copyratio.max = 3,
                       qmax = 10,
                       lamda = 0.5,
                       method = c("Grid Search", "Bayesian Optimization"),
                       min.sol.freq = 0.05,
                       snv.type = c("somatic", "germline"),
                       verbose = FALSE) {
    stopifnot(alpha.min>=0.1, alpha.max<=1, tau.min>=0.0, tau.max<=20, min.sol.freq>=0, min.sol.freq<=1,
              qmax<=100, lamda>0, lamda<=1)

    if(!inherits(absCopyNumber, "absCopyNumber")) {
        stop("Wrong input object, please check!")
    }

    method = match.arg(method)
    snv.type = match.arg(snv.type)

    params = list(
        alpha.min = alpha.min,
        alpha.max = alpha.max,
        tau.min = tau.min,
        tau.max = tau.max,
        copyratio.min = copyratio.min,
        copyratio.max = copyratio.max,
        qmax = qmax,
        lamda = lamda,
        method = method,
        min.sol.freq = min.sol.freq,
        snv.type = snv.type
    )

    absCopyNumber@params = params
    return(absCopyNumber)
}


subset.absCopyNumber = function(object, samples = NULL){
    if(is.null(samples)) {
        message("no sample name provided, return original object")
        return(object)
    }else{
        sample_exist = unique(object@data$sample)
        if(!all(samples %in% sample_exist)) {
            stop("detect invalid sample names:\n  ", paste(samples[!samples %in% sample_exist], collapse = " "))
        }else{
            object@data = object@data[sample %in% samples]
            if(!identical(object@SNV, data.table::data.table())){
                object@SNV = object@SNV[sample %in% samples]
            }
        }
        return(object)
    }
}

# setMethod(f = "abs_subset", signature = "absCopyNumber", definition = function(object, samples = NULL){
#     if(is.null(samples)) {
#         message("no sample name provided, return original object")
#         return(object)
#     }else{
#         sample_exist = unique(object@data$sample)
#         if(!all(samples %in% sample_exist)) {
#             stop("detect invalid sample names:\n  ", paste(samples[!samples %in% sample_exist], collapse = " "))
#         }else{
#             object@data = object@data[sample %in% samples]
#             if(!identical(object@SNV, data.table::data.table())){
#                 object@SNV = object@SNV[sample %in% samples]
#             }
#         }
#         return(object)
#     }
# })
