#' absCopyNumber calling with user specified parameters
#' @description This is the third step for absCopyNumber pipeline, absCopyNumber compute
#' solutions of purity and ploidy, then generate absolute copynumber and (optional) multiplicity
#' for top 1 solution.
#' @param absCopyNumber a \code{absCopyNumber} object generate from \code{abs_prepare} function.
#' @param samples character vector to subset samples for calling, default is \code{NULL}.
#' @param verbose if \code{TRUE}, print extra information.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a \code{absCopyNumber} object with calling result.
#' @import data.table
#' @export
#'
#' @examples
#' file_cn = system.file("extdata/example.cn.txt.gz", package = "absCopyNumber")
#' file_snv = system.file("extdata/example.snv.txt.gz", package = "absCopyNumber")
#'
#' \donttest{
#' res1 =  abs_initialize(seg = file_cn, snv = file_snv, verbose = T)
#' res1 =  abs_prepare(res1)
#' res1 =  abs_calling(res1, verbose = TRUE)
#' }
#'
#' @seealso \code{\link[absCopyNumber]{absCopyNumber}}, \code{\link[absCopyNumber]{abs_initialize}}, , \code{\link[absCopyNumber]{abs_prepare}}
abs_calling = function(absCopyNumber, samples = NULL, verbose = FALSE){

    stopifnot(is.character(samples) | is.null(samples))
    if(!inherits(absCopyNumber, "absCopyNumber")) {
        stop("Wrong input object, please check!")
    }

    object = absCopyNumber

    if(!is.null(samples)){
        object = abs_subset(object, samples = samples)
    }

    # method calling with expression
    inner_calling = expression({
        if (verbose) cat("Loading data...\n")
        seg = .SD
        sampleN = .BY[[1]]
        snv = object@SNV[sample == sampleN]

        r <- seg$normalized.ratio
        params = object@params
        max.r.cutoff <- params$copyratio.max
        min.r.cutoff <- params$copyratio.min
        #----- assign SNP to each segment
        if (!identical(snv, data.table::data.table()) & nrow(snv) != 0) {
            # make chrom as char vector
            if(verbose){
                cat(paste0("Detect snv data for sample '", sampleN, "'...\n"))
            }

            if (class(snv$chrom) != "character") {
                snv[, chrom := as.character(chrom)]
            } else{
                if (grepl(pattern = "chr",
                          snv$chrom[1],
                          ignore.case = T)) {
                    snv$chrom = gsub(
                        pattern = "^chr",
                        replacement = "",
                        snv$chrom,
                        ignore.case = T
                    )
                }
            }

            snp2seg <- NULL
            fb <- NULL
            het.ind <- NULL
            for (i in 1:nrow(snv)) {
                ind <-
                    which(
                        seg$chrom == as.character(snv[i, "chrom"]) &
                            seg$loc.start <= as.integer(snv[i, "position"]) &
                            seg$loc.end   >= as.integer(snv[i, "position"])
                    )
                if (length(ind) == 1) {
                    if (r[ind] <= max.r.cutoff & r[ind] >= min.r.cutoff) {
                        snp2seg <- c(snp2seg, ind)
                        fb <- c(fb, as.numeric(snv[i, "tumor_var_freq"]))
                        het.ind <- c(het.ind, i)
                    }
                }
            }

            if (verbose) {
                cat(
                    "Assign SNP to each segment: ========================================\n"
                )
                cat("# of SNVs used:", length(het.ind), "\n")
                cat(
                    "====================================================================\n"
                )
            }
        } else {
            if (verbose){
                cat(paste0("No snv data find for sample '", sampleN, "'... \n"))
                cat("Skip variant assignment...\n")
            }
        }

        #---- Process parameters
        alpha.min = params$alpha.min
        alpha.max = params$alpha.max
        tau.min   = params$tau.min
        tau.max   = params$tau.max
        copyratio.min = params$copyratio.min
        copyratio.max = params$copyratio.max
        qmax      = params$qmax
        lamda     = params$lamda
        method    = params$method
        min.sol.freq = params$min.sol.freq
        snv.type  = params$snv.type

        if (verbose) {
            cat("Loading arguments: ================================================\n")
            cat("alpha.min =", alpha.min, "\n")
            cat("alpha.max =", alpha.max, "\n")
            cat("tau.min   =", tau.min, "\n")
            cat("tau.max   =", tau.max, "\n")
            cat("copyratio.min =", copyratio.min, "\n")
            cat("copyratio.max =", copyratio.max, "\n")
            cat("qmax      =", qmax, "\n")
            cat("lamda     =", lamda, "\n")			# only CN when lamda equals 1.0
            cat("method    =", method, "\n")
            cat("min.sol.freq   =", min.sol.freq, "\n")
            cat("snv.type  =", snv.type, "\n")
            cat(
                "====================================================================\n"
            )
        }

        outlier.frac.1 <-
            length(which(r > max.r.cutoff)) / length(r)
        if (verbose)
            cat(
                100 * outlier.frac.1,
                paste0(
                    "% segments with copy ratio r > ",
                    max.r.cutoff,
                    "\n"
                )
            )

        outlier2.frac.1 <-
            length(which(r < min.r.cutoff)) / length(r)
        if (verbose) {
            cat(
                100 * outlier2.frac.1,
                paste0(
                    "% segments with copy ratio r < ",
                    min.r.cutoff,
                    "\n"
                )
            )
            cat(
                "====================================================================\n"
            )
        }

        #------- calculate weight for CN data
        gf <- seg$loc.end - seg$loc.start + 1  # vector of genome fraction
        gf <- gf / sum(gf)

        if(verbose) {
            cat("Assign weights:\n")
            cat("Assign to copy-number-based ratio estimator...\n")
            cat("Copy ratio outside range (copyratio.min, copyration.max) will get weight value 0...\n")
        }

        wts <- rep(1, length(r))
        wts[r > max.r.cutoff] <- 0.0
        wts[r < min.r.cutoff] <- 0.0
        wts.cn <- wts / sum(wts)

        if(exists("snp2seg")){
            if (verbose) {
                cat("Assign to SNV-frequency-based estimator...\n")
            }
            wts.het <- 1.0 / length(snp2seg) * (1.0 - lamda)
        }

        #------ Assign Method
        # if (!any(method %in% c("Grid Search", "Bayesian Optimization"))) {
        #     stop("Invalid method, can only be one of \"Grid Search\", \"Bayesian Optimization\".")
        # }else
        if (method == "Grid Search") {
            #-- construct solution grid
            alpha.grid <-
                seq(from = alpha.min, to = alpha.max, by = 0.05)
            tau.grid <- seq(from = tau.min, to = tau.max, by = 0.05)
            search.grid <-
                expand.grid(alpha.grid = alpha.grid, tau.grid = tau.grid)
            n.grid <- nrow(search.grid)

            #-- construct grid search result data.frame
            search.res <-
                data.frame(
                    alpha = numeric(),
                    tau = numeric(),
                    tau.def = numeric(),
                    mse = numeric(),
                    alpha0 = numeric(),
                    tau0 = numeric()
                )

            #-- set text progress bar for Grid Search
            if (verbose) {
                cat("Grid searching:\n")
                total <- n.grid
                # create progress bar
                pb <-
                    txtProgressBar(
                        min = 0,
                        max = total,
                        style = 3,
                        width = 60,
                        char = ">"
                    )
            }

            #-- choose run function and then run
            if (exists("snp2seg")) {
                for (k in 1:n.grid) {
                    if (verbose)
                        setTxtProgressBar(pb, k)

                    alpha0 <- search.grid[k, 1]
                    tau0 <- search.grid[k, 2]

                    a.res <-
                        optimize.alpha(
                            alpha0,
                            alpha.min,
                            alpha.max,
                            tau0,
                            wts.cn,
                            r,
                            gf,
                            eta = 1.0,
                            wts.het,
                            snp2seg,
                            fb,
                            qmax,
                            snv.type
                        )
                    if (!is.null(a.res)) {
                        search.res <- rbind(search.res, a.res)
                    }
                }

                if (verbose) close(pb)
            } else {
                for (k in 1:n.grid) {
                    # update progress bar
                    if (verbose) setTxtProgressBar(pb, k)

                    alpha0 <- search.grid[k, 1]
                    tau0 <- search.grid[k, 2]

                    a.res <-
                        optimize.alpha.simple(alpha0,
                                              alpha.min,
                                              alpha.max,
                                              tau0,
                                              wts.cn,
                                              r,
                                              gf,
                                              eta = 1.0,
                                              qmax)
                    if (!is.null(a.res)) {
                        search.res <- rbind(search.res, a.res)
                    }
                }

                if (verbose)
                    close(pb)
            }

            #--- check if have no any result
            res_check <-
                data.frame(
                    alpha = numeric(),
                    tau = numeric(),
                    tau.def = numeric(),
                    mse = numeric(),
                    alpha0 = numeric(),
                    tau0 = numeric()
                )

            if (identical(search.res, res_check)) {
                if (verbose){
                    cat("Find no result, suspect a normal sample...\n")
                    cat("====================\n")
                }
               res = list(alpha = 1,
                           tau = 2,
                           mse = 0,
                           count = integer(),
                           rank = 1L)
                res
            } else {
                #--- tidy result data.frame
                colnames(search.res) <-
                    c("alpha", "tau", "tau.def", "mse", "alpha0", "tau0")

                search.res[, c("alpha", "tau", "tau.def")] <-
                    round(search.res[, c("alpha", "tau", "tau.def")], 2)
                tmp <-
                    aggregate(search.res[, c("tau.def", "mse")], search.res[, c("alpha", "tau")], mean)
                tmp2 <-
                    aggregate(search.res[, 1], search.res[, c("alpha", "tau")], length)

                if (verbose) {
                    cat("Total", nrow(search.res), "results\n")
                }

                search.res <- data.frame(tmp, count = tmp2[, 3])
                if (verbose)
                    cat(nrow(search.res), "unique results\n")

                oo <- order(search.res$mse)
                search.res <- search.res[oo,]
                rownames(search.res) <- NULL

                # only cluster when result number more than 1
                if (nrow(search.res) > 1) {
                    if (verbose)
                        cat("Clustering results...\n")
                    # do clustering
                    search.res <-
                        cluster.solution(search.res, alpha.cut = 0.10, tau.cut = 0.15)
                }

                if (verbose)
                    cat("Filtering impossible solutions...\n")
                # filtering out impossible solutions

                min.sol.freq <- min.sol.freq * sum(search.res[, "count"])
                proper.ind <-
                    which(
                        search.res[, "alpha"] >= alpha.min &
                            search.res[, "alpha"] <= alpha.max &
                            search.res[, "tau.def"] >= tau.min &
                            search.res[, "tau.def"] <= tau.max &
                            search.res[, "count"] >= min.sol.freq
                    )
                if (length(proper.ind) > 0) {
                    search.res <- search.res[proper.ind, , drop = FALSE]
                }
                if (verbose) cat("Final solution number:", nrow(search.res), "\n=========\n")

                res = list(alpha = search.res$alpha,
                           tau = search.res$tau,
                           mse = search.res$mse,
                           count = search.res$count,
                           rank = seq_along(search.res$alpha))
                res
            }
        } else if (method == "Bayesian Optimization") {
            # TODO
            NULL
        } else {
                stop("invalid method.")
            }



    })

    object@estimation = object@data[, eval(inner_calling), by = sample]


    #--- evaluate top result
    if (verbose)
        cat("Evaluating the best result...\n")
    object@TopResult = abs_obtain(object = object, rank = 1, onlyCN = FALSE, verbose = FALSE)

    cat("Done.\n")
    object
}


#' Obtain absolute copy number and optinal multiplicity (if SNV data provided)
#'
#' @description According to the user specified \code{rank} value, this function retrieve corresponding
#' purity and ploidy solution pair, then compute absolute copy number, optional muliplicity.
#' @param object a \code{absCopyNumber} object after absCopyNumber calling
#' @param rank integer, rank value of solution, i.e. \code{estimation} slot.
#' @param samples character vector, use if you wanna just get solution for specified samples.
#' @param onlyCN if \code{TRUE}, just get absolute copy number result (a \code{data.table}), otherwise
#' return a list include both absolute copy number and snv data with multiplicity.
#' @param verbose if \code{TRUE}, print extra information.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a \code{list} or a \code{data.table}, details see \code{onlyCN} option.
#' @import data.table
#' @export
abs_obtain = function(object,
                      rank = 1,
                      samples = NULL,
                      onlyCN = FALSE,
                      verbose = FALSE) {
    stopifnot(is.numeric(rank), length(rank) == 1)
    if (!inherits(object, "absCopyNumber")) {
        stop(
            "Input should be a absCopyNumber object after calling with abs_calling function."
        )
    } else {
        if (verbose)
            cat("Detect absCopyNumber object as input.\n")
        if (identical(object@params, list())) {
            stop("Should be a absCopyNumber object after calling with abs_calling function.")
        }


        if (!is.null(samples)) {
            object = abs_subset(object, samples = samples)
        }

        # assign snp to seg
        inner_calling = expression({
            if (verbose)
                cat("Loading data...\n")
            seg = .SD
            sampleN = .BY[[1]]
            snv = object@SNV[sample == sampleN]

            r <- seg$normalized.ratio
            params = object@params
            max.r.cutoff <- params$copyratio.max
            min.r.cutoff <- params$copyratio.min
            #----- assign SNP to each segment
            if (!identical(snv, data.table::data.table()) & nrow(snv) != 0) {
                # make chrom as char vector
                if (verbose) {
                    cat(paste0(
                        "Detect snv data for sample '",
                        sampleN,
                        "'...\n"
                    ))
                }

                if (class(snv$chrom) != "character") {
                    snv[, chrom := as.character(chrom)]
                } else{
                    if (grepl(pattern = "chr",
                              snv$chrom[1],
                              ignore.case = T)) {
                        snv$chrom = gsub(
                            pattern = "^chr",
                            replacement = "",
                            snv$chrom,
                            ignore.case = T
                        )
                    }
                }

                snp2seg <- NULL
                fb <- NULL
                het.ind <- NULL
                for (i in 1:nrow(snv)) {
                    ind <-
                        which(
                            seg$chrom == as.character(snv[i, "chrom"]) &
                                seg$loc.start <= as.integer(snv[i, "position"]) &
                                seg$loc.end   >= as.integer(snv[i, "position"])
                        )
                    if (length(ind) == 1) {
                        if (r[ind] <= max.r.cutoff & r[ind] >= min.r.cutoff) {
                            snp2seg <- c(snp2seg, ind)
                            fb <-
                                c(fb, as.numeric(snv[i, "tumor_var_freq"]))
                            het.ind <- c(het.ind, i)
                        }
                    }
                }

                if (verbose) {
                    cat(
                        "Assign SNP to each segment: ========================================\n"
                    )
                    cat("# of SNVs used:", length(het.ind), "\n")
                    cat(
                        "====================================================================\n"
                    )
                }
                list(snp2seg = snp2seg, het.ind = het.ind, fb = fb)
            } else {
                if (verbose) {
                    cat(paste0(
                        "No snv data find for sample '",
                        sampleN,
                        "'... \n"
                    ))
                    cat("Skip variant assignment...\n")

                }
                list(snp2seg = integer(), het.ind = integer(), fb = double())
            }
        })

        if (!onlyCN){
            snv_assign = object@data[, eval(inner_calling), by = sample]
        }

        res_cn = object@data
        res_cn$r.estimate = double()
        res_cn$copynumber = integer()

        if (!onlyCN){
            res_snv = object@SNV
            if (!identical(res_snv, data.table::data.table())){
                res_snv$multiplicity = integer()
            }
        }

        samples = unique(object@estimation$sample)
        qmax = object@params$qmax
        snv.type = object@params$snv.type

        for (i in samples) {

            #--- absolute copy number
            data = object@estimation[sample == i]
            r = object@data[sample == i]$normalized.ratio

            alpha.opt <- as.numeric(data[rank, "alpha"])
            tau.opt   <- as.numeric(data[rank, "tau"])
            tmp <-
                (r * (alpha.opt * tau.opt + (1.0 - alpha.opt) * 2) - (1.0 - alpha.opt) *
                     2) / alpha.opt
            qs <- round(tmp)
            qs[qs < 0] <- 0
            qs[qs > qmax] <- qmax
            rhat <-
                (alpha.opt * qs + (1.0 - alpha.opt) * 2) / (alpha.opt * tau.opt + (1.0 -
                                                                                       alpha.opt) * 2)

            res_cn[sample == i]$r.estimate = rhat
            res_cn[sample == i]$copynumber = as.integer(qs)

            #--- multiplicity

            if (exists("snv_assign") & !onlyCN & !identical(res_snv, data.table::data.table())) {
                snp2seg = snv_assign[sample == i]$snp2seg
                if (is.null(snp2seg)) {
                    res_snv[sample == i]$multiplicity = integer()
                    break
                }
                het.ind = snv_assign[sample == i]$het.ind
                fb = snv_assign[sample == i]$fb
                q.het <- qs[snp2seg]
                if (snv.type == "somatic") {
                    tmp <-
                        ((alpha.opt * q.het + (1.0 - alpha.opt) * 2) * fb) / alpha.opt
                } else if (snv.type == "germline") {
                    tmp <-
                        ((alpha.opt * q.het + (1.0 - alpha.opt) * 2) * fb - (1.0 - alpha.opt)) /
                        alpha.opt
                }
                mg <- round(tmp)
                mg[mg < 0] <- 0
                tmp.ind <- which(mg > q.het)
                if (length(tmp.ind) > 0) {
                    mg[tmp.ind] <- q.het[tmp.ind]
                }

                res_snv[sample == i][het.ind]$multiplicity = as.integer(mg)
            }
        }
    }

    if (onlyCN) {
        return(res_cn)
    } else {
        return(list(absCN = res_cn, absSNV = res_snv))
    }
}


# call_gridsearch = function(){
#
# }
#
# call_bayesianOP = function(){
#
# }
