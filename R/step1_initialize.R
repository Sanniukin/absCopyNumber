#' Initialize absCopyNumber object from file or data.frame
#'
#' This is the first step for user to use a flexible way
#' @import data.table methods
abs_initialize = function(
    seg,
    min.seg.len = NULL,
    snv = NULL,
    isMaf = FALSE,
    sample_seg = NULL,
    sample_snv = NULL,
    ratio.min = -3,
    ratio.max = 3,
    platform = c("WES", "WGS", "MicroArray"),
    #alpha.min = 0.2,
    #alpha.max = 1.0,
    #tau.min = 0.5,
    #tau.max = 8.0,
    #min.sol.freq = 0.05,
    #qmax = 10,
    #lamda = 0.5,
    verbose = FALSE
) {
    stopifnot(is.logical(isMaf),
              is.logical(verbose))

    # check platform
    platform <- match.arg(platform)
    if (platform == "WES") {
        if (is.null(min.seg.len))
            min.seg.len <- 200
    } else if (platform == "WGS") {
        if (is.null(min.seg.len))
            min.seg.len <- 3000
    } else if (platform == "MicroArray") {
        if (is.null(min.seg.len))
            min.seg.len <- 0
    } else {
        stop("ERROR: you must specify the platform: WES or WGS or MicroArray.")
    }

    # 1: Read segmentation file if it’s a file or convert to data.table i
    if (is.data.frame(x = seg)) {
        seg = data.table::setDT(seg)
    } else{
        if (verbose)
            message("reading segmentation file...")

        if (as.logical(length(grep(
            pattern = 'gz$',
            x = seg,
            fixed = FALSE
        )))) {
            # If system is Linux use fread, else use gz connection to read gz file.
            if (Sys.info()[['sysname']] == 'Windows') {
                seg.gz = gzfile(description = seg, open = 'r')
                suppressWarnings(seg <-
                                     data.table::as.data.table(
                                         read.csv(
                                             file = seg.gz,
                                             header = TRUE,
                                             sep = '\t',
                                             stringsAsFactors = FALSE,
                                             comment.char = "#"
                                         )
                                     ))
                close(seg.gz)
            } else if (Sys.info()[['sysname']] == 'Darwin') {
                seg = suppressWarnings(
                    data.table::fread(
                        cmd = paste('gunzip -c ', seg),
                        sep = '\t',
                        stringsAsFactors = FALSE,
                        verbose = FALSE,
                        data.table = TRUE,
                        showProgress = TRUE,
                        header = TRUE,
                        fill = TRUE
                    )
                )
            } else{
                seg = suppressWarnings(
                    data.table::fread(
                        cmd = paste('zcat <', seg),
                        sep = '\t',
                        stringsAsFactors = FALSE,
                        verbose = FALSE,
                        data.table = TRUE,
                        showProgress = TRUE,
                        header = TRUE,
                        fill = TRUE
                    )
                )
            }
        } else {
            suppressWarnings(
                seg <-
                    data.table::fread(
                        input = seg,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        verbose = FALSE,
                        data.table = TRUE,
                        showProgress = TRUE,
                        header = TRUE,
                        fill = TRUE
                    )
            )
        }
    }

    # 2: same for snv if snv is not NULL
    if (!is.null(snv)) {
        if (is.data.frame(x = snv)) {
            snv = data.table::setDT(snv)
        } else{
            if (verbose)
                message("reading somatic mutation file...")

            if (as.logical(length(grep(
                pattern = 'gz$',
                x = snv,
                fixed = FALSE
            )))) {
                # If system is Linux use fread, else use gz connection to read gz file.
                if (Sys.info()[['sysname']] == 'Windows') {
                    snv.gz = gzfile(description = snv, open = 'r')
                    suppressWarnings(snv <-
                                         data.table::as.data.table(
                                             read.csv(
                                                 file = snv.gz,
                                                 header = TRUE,
                                                 sep = '\t',
                                                 stringsAsFactors = FALSE,
                                                 comment.char = "#"
                                             )
                                         ))
                    close(snv.gz)
                } else if (Sys.info()[['sysname']] == 'Darwin') {
                    snv = suppressWarnings(
                        data.table::fread(
                            input = paste('gunzip -c ', snv),
                            sep = '\t',
                            stringsAsFactors = FALSE,
                            verbose = FALSE,
                            data.table = TRUE,
                            showProgress = TRUE,
                            header = TRUE,
                            fill = TRUE
                        )
                    )
                } else {
                    snv = suppressWarnings(
                        data.table::fread(
                            input = paste('zcat <', snv),
                            sep = '\t',
                            stringsAsFactors = FALSE,
                            verbose = FALSE,
                            data.table = TRUE,
                            showProgress = TRUE,
                            header = TRUE,
                            fill = TRUE
                        )
                    )
                }
            } else {
                suppressWarnings(
                    snv <-
                        data.table::fread(
                            input = snv,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            verbose = FALSE,
                            data.table = TRUE,
                            showProgress = TRUE,
                            header = TRUE,
                            fill = TRUE
                        )
                )
            }
        }

    }

    # 3: validata seg data and transform data as standard form
    status_seg = validate_seg(seg)

    # function to standardize segmentation input
    sd_seg_input = function(seg, sample_seg, verbose = FALSE) {
        # check if 'sample' column exists
        if ("sample" %in% colnames(seg)) {
            seg = seg[, c(
                "sample",
                "chrom",
                "loc.start",
                "loc.end",
                "eff.seg.len",
                "normalized.ratio"
            )]
        } else if (!is.null(sample_seg)) {
            sample_seg = as.character(sample_seg)
            if (sample_seg %in% colnames(seg)) {
                if (verbose)
                    message("treat sample_seg argument as column name")

                cols = c(
                    sample_seg,
                    "chrom",
                    "loc.start",
                    "loc.end",
                    "eff.seg.len",
                    "normalized.ratio"
                )
                seg = seg[, ..cols]
                colnames(seg)[1] = "sample"
            } else {
                if (verbose)
                    message("treat sample_seg argument as sample name (vector)")
                seg = data.table::data.table(sample = sample_seg,
                                             seg[, c("chrom",
                                                     "loc.start",
                                                     "loc.end",
                                                     "eff.seg.len",
                                                     "normalized.ratio")])
            }

        } else {
            if (verbose)
                message(
                    "no sample name find and sample_seg argument is not provided, fill with 'sample'."
                )
            seg = data.table::data.table(sample = "sample",
                                         seg[, c("chrom",
                                                 "loc.start",
                                                 "loc.end",
                                                 "eff.seg.len",
                                                 "normalized.ratio")])
        }

        seg
    }

    if (status_seg == "standard ngs") {
        seg = sd_seg_input(seg, sample_seg, verbose)
    } else {
        if (verbose)
            message(
                "transforming standard segmentation format to standard input of absCopyNumber..."
            )
        seg[, ':='(
            chrom = Chromosome,
            loc.start = Start,
            loc.end = End,
            eff.seg.len = End - Start + 1,
            normalized.ratio = Segment_Mean
        )]
        seg = sd_seg_input(seg, sample_seg, verbose)

    }

    # 4: validata snv data and transform data as standard form
    if (isMaf) {
        if (verbose)
            message("checking required columns in Maf...")

        # https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
        maf_cols = c(
            "Chromosome",
            "Start_position",
            "Variant_Type",
            "t_ref_count",
            "t_alt_count",
            "Tumor_Sample_Barcode"
        )
        if (!all(maf_cols %in% colnames(snv))) {
            stop(
                "Minimal required columns \"Tumor_Sample_Barcode\", \"Chromosome\", \"Start_position\", \"Variant_Type\", \"t_ref_count\", \"t_alt_count\" in Maf is not satisfied. \nPlease check your input or not provide snv using Maf."
            )
        }
        else {
            snv = snv[Variant_Type == "SNP", ..maf_cols]
            snv[, ':='(
                sample = Tumor_Sample_Barcode,
                chrom = Chromosome,
                position = Start_position,
                tumor_var_freq = t_alt_count / (t_ref_count + t_alt_count)
            )]
            snv = snv[, c("sample", "chrom",  "position", "tumor_var_freq")]
        }
    } else if (is.null(snv)) {
        snv = data.table::data.table()
    } else {
        status_snv = validata_snv(snv, verbose)
        if (status_snv == "standard snv") {
            if ("sample" %in% colnames(snv)) {
                if (verbose) message("treat sample_snv argument as sample column")
                snv = snv[, c("sample",
                              "chrom",
                              "position",
                              "tumor_var_freq")]
            } else if (!is.null(sample_snv)) {
                if (sample_snv %in% colnames(snv)) {
                    if (verbose)
                        message("treat sample_snv argument as column name")
                    cols = c(sample_snv,
                             "chrom",
                             "position",
                             "tumor_var_freq")
                    snv = snv[, ..cols]
                    colnames(snv)[1] = "sample"
                } else {
                    if (verbose)
                        message("treat sample_snv argument as sample name (vector)")
                    snv = data.table::data.table(sample = sample_snv,
                                                 snv[, c("chrom",  "position", "tumor_var_freq")])
                }
            } else {
                if (verbose)
                    message(
                        "no sample name find and sample_snv argument is not provided, fill with 'sample'."
                    )
                snv = data.table::data.table(sample = 'sample',
                                             snv[, c("chrom",  "position", "tumor_var_freq")])
            }
        }

    }

    # 5: read parameters
    if (verbose) {
        cat("Reading arguments: ================================================\n")
        cat("min.seg.len =", min.seg.len, "\n")
    }

    # 6: filtering based on eff.seg.len
    seg_old = data.table::copy(seg)

    n.seg.1 <- nrow(seg)
    seg <- na.omit(seg)
    seg <- seg[eff.seg.len >= min.seg.len]
    n.seg.2 <- nrow(seg)
    if (verbose){
        cat("Filtering segments with NAs and length less than eff.seg.len ...\n")
        cat("All segments: ", n.seg.1, "\n")
        cat("Retained segments: ", n.seg.2, "\n")
        cat("====================================================================\n")
    }

    gf <-
        seg$loc.end - seg$loc.start + 1  # vector of genome fraction: spaced length
    gf <- gf / sum(gf)
    eta <- 1.0

    r <- seg$normalized.ratio

    max.r.cutoff <- ratio.max
    min.r.cutoff <- ratio.min

    outlier.frac.1 <-
        length(which(r > max.r.cutoff)) / length(r)
    outlier.gf.1 <- sum(gf[r > max.r.cutoff])
    if (verbose)
        cat(
            100 * length(which(r > max.r.cutoff)) / length(r),
            paste0("% segments with copy ratio r>", max.r.cutoff, "  before rescaling\n")
        )

    outlier2.frac.1 <-
        length(which(r < min.r.cutoff)) / length(r)
    outlier2.gf.1 <- sum(gf[r < min.r.cutoff])
    if (verbose) {
        cat(
            100 * length(which(r < min.r.cutoff)) / length(r),
            paste0("% segments with copy ratio r<", min.r.cutoff, "  before rescaling\n")
        )
        cat("Filtering them...\n")
        cat(
            "====================================================================\n"
        )
    }

    # make chrom as char vector
    if(class(seg$chrom) != "character") {
        seg[, chrom:=as.character(chrom)]
    }else{
        if (grepl(pattern = "chr", seg$chrom[1], ignore.case = T)) {
            seg$chrom = gsub(pattern = "^chr", replacement = "", seg$chrom, ignore.case = T)
        }
    }

    # 7: assign SNP to each segment after filtering
    if(!identical(snv, data.table::data.table())){

        # make chrom as char vector
        if(class(snv$chrom) != "character") {
            snv[, chrom:=as.character(chrom)]
        }else{
            if (grepl(pattern = "chr", snv$chrom[1], ignore.case = T)) {
                snv$chrom = gsub(pattern = "^chr", replacement = "", snv$chrom, ignore.case = T)
            }
        }

        snp2seg <- NULL
        fb <- NULL
        het.ind <- NULL
        for (i in 1:nrow(snv)) {
            ind <-
                which(seg$chrom == as.character(snv[i, "chrom"]) &
                          seg$loc.start <= as.integer(snv[i, "position"]) &
                          seg$loc.end   >= as.integer(snv[i, "position"]))
            if (length(ind) == 1) {
                if (r[ind] <= max.r.cutoff & r[ind] >= min.r.cutoff) {
                    snp2seg <- c(snp2seg, ind)
                    fb <- c(fb, snv[i, "tumor_var_freq"])
                    het.ind <- c(het.ind, i)
                }
            }
        }

        if (verbose) {
            cat("Assign SNP to each segment: ========================================\n")
            cat("# of SNVs used:", length(het.ind), "\n")
            cat(
                "====================================================================\n"
            )
        }
    }


    absSummary = data.table::data.table(
        nsample = length(unique(seg$sample)),
        nchrom = length(unique(seg$chrom)),
        max.seg = max(seg$eff.seg.len),
        min.seg = min(seg$eff.seg.len),
        max.ratio = max(seg$normalized.ratio),
        min.ratio = min(seg$normalized.ratio)
    )
    res = absCopyNumber(data = seg,
                        snv.data = snv[if(exists("het.ind")) return(het.ind) else TRUE,],
                        params=list(),
                        origin=list(
                            seg = seg_old,
                            snv = snv
                        ),
                        summary=absSummary,
                        result=data.table::data.table(),
                        absCN=data.table::data.table(),
                        absSNV=data.table::data.table())

    if (verbose) {
        message("Done !")
    }

    return(res)
}

validate_seg = function(seg, verbose = FALSE) {
    seg_col1 = c("chrom",
                 "loc.start",
                 "loc.end",
                 "eff.seg.len",
                 "normalized.ratio")
    seg_col2 = c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
    if (all(seg_col2 %in% colnames(seg))) {
        if (verbose)
            message("standard segmentation columns are detected.")
        return("standard seg")
    } else if (all(seg_col1 %in% colnames(seg))) {
        if (verbose)
            message("standard NGS segmentation required columns are detected")
        return("standard ngs")
    } else {
        stop("invalid seg input, please check!")
    }
}

validata_snv = function(snv, verbose = FALSE) {
    snv_col = c("chrom",  "position", "tumor_var_freq")
    if (all(snv_col %in% colnames(snv))) {
        if (verbose)
            message("standard snv required columns are detected.")
        return("standard snv")
    } else {
        stop("invalid snv input, please check!")
    }
}
