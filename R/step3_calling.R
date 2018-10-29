# call algorithm
abs_calling = function(absCopyNumber, samples = NULL){

    stopifnot(is.character(samples) | is.null(samples))
    if(!inherits(absCopyNumber, "absCopyNumber")) {
        stop("Wrong input object, please check!")
    }

    object = absCopyNumber

    if(!is.null(samples)){
        object = subset.absCopyNumber(object, samples = samples)
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
        if (!identical(snv, data.table::data.table())) {
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
                NULL
            }

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
            if (verbose) cat("Final solution number:", nrow(search.res), "\n")

            res = list(alpha = search.res$alpha,
                 tau = search.res$tau,
                 mse = search.res$mse,
                 count = search.res$count)
            res

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
        cat("Outputing the best result...\n")


    cat("Done.\n")
    object
}


abs_obtain = function(data, alpha = NULL, tau = NULL, qmax = 10, snv.type = "somatic"){

    if (inherits(data, "absCopyNumber")) {
        data = data@estimation
    }

    min.ind <- which(dasta$mse == data$mse[1])
    #cat("# of optimal parameters:",length(min.ind),"\n")

    alpha.opt <- data[1, "alpha"]
    tau.opt   <- data[1, "tau"]

    tmp <-
        (r * (alpha.opt * tau.opt + (1.0 - alpha.opt) * 2) - (1.0 - alpha.opt) *
             2) / alpha.opt
    qs <- round(tmp)
    qs[qs < 0] <- 0
    qs[qs > qmax] <- qmax
    rhat <-
        (alpha.opt * qs + (1.0 - alpha.opt) * 2) / (alpha.opt * tau.opt + (1.0 -
                                                                               alpha.opt) * 2)
    abs.cn <- data.frame(data[, 1:4],
                         r = r,
                         rhat = rhat,
                         CN = qs)

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
    abs.mg <-
        data.frame(snv.data[het.ind,], multiplicity = mg)

}

# call_gridsearch = function(){
#
# }
#
# call_bayesianOP = function(){
#
# }
