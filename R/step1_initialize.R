
#' Initialize absCopyNumber object from file or data.frame
#'
#' @description  This is the first step for user to use a flexible way for copynumber data input,
#' absolute copy number (with purity and ploidy pair) calling. This function help
#' user input their copy segmentation data as file or \code{data.frame}, input
#' their SNV data as file or \code{data.frame}.
#'
#' @details when process data with multiple samples, absCopyNumber require data.frame/file with one
#' column for the sample variable, see \code{sample_seg} option.
#' Support input format by absCopyNumber package please run \code{abs_supportfiles} function.
#'
#' @param seg The name of file or a \code{data.frame} containing the segmentation data.
#' @param min.seg.len The minimum
#' length of a segment to be included in computation. The default value is 200
#' bp for WES, 3000 bp for WGS and 0 bp for MicroArray.
#' @param snv optinal. The name of file or a \code{data.frame} containing a set of somatic single
#' nucleotide variants (SNVs).
#' @param isMaf logical. Whether the specified \code{snv} argument is provide as MAF file/data.frame.
#' @param sample_seg character, default is \code{NULL}. At default function will look for 'sample' column in
#' data (the best way to input multiple samples) of CNV or SNV. If 'sample' column is not find
#' and \code{sample_seg} is specified,
#' \code{sample_seg} will be treated as column name storing sample name or sample names (character vector),
#' otherwise data will be treated as single sample data and sample name is 'sample'.
#' @param sample_snv same as \code{sample_seg} but for snv data.
#' @param platform character, default is 'WES'. Must be one of 'WES', 'WGS' or 'MicroArray'.
#' @param verbose if \code{TRUE}, print extra information.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a \code{absCopyNumber} object
#' @import data.table methods
#' @export
#' @examples
#' file_cn = system.file("extdata/example.cn.txt.gz", package = "absCopyNumber")
#' file_snv = system.file("extdata/example.snv.txt.gz", package = "absCopyNumber")
#'
#' \donttest{
#' res1 =  abs_initialize(seg = file_cn, snv = file_snv, verbose = T)
#' res2 =  abs_initialize(seg = file_cn, snv = file_snv, sample_seg = "test1",
#'                        sample_snv = "test1", verbose = T)
#' }
#'
#' @seealso \code{\link[absCopyNumber]{absCopyNumber}}, \code{\link[absCopyNumber]{abs_prepare}}, , \code{\link[absCopyNumber]{abs_calling}}
abs_initialize = function(
    seg,
    min.seg.len = NULL,
    snv = NULL,
    isMaf = FALSE,
    sample_seg = NULL,
    sample_snv = NULL,
    platform = c("WES", "WGS", "MicroArray"),
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

    origin = data.table::copy(seg)

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
                            cmd = paste('gunzip -c ', snv),
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
                            cmd = paste('zcat <', snv),
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
            normalized.ratio = 2 ^ Segment_Mean
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

    # make chrom as char vector
    if(class(seg$chrom) != "character") {
        seg[, chrom:=as.character(chrom)]
    }else{
        if (grepl(pattern = "chr", seg$chrom[1], ignore.case = T)) {
            seg$chrom = gsub(pattern = "^chr", replacement = "", seg$chrom, ignore.case = T)
        }
    }

    res = absCopyNumber(data = seg, SNV = snv, origin = origin, params = list(platform=platform))
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


# batch initializ : TODO?
# abs_binitailize = function(batch_seg,
#                            min.seg.len = NULL,
#                            batch_snv = NULL,
#                            isMaf = FALSE,
#                            samplenames = NULL,
#                            platform = c("WES", "WGS", "MicroArray")
# ){
#
# }

# Hack global variable  ---------------------------------------------------
base::suppressMessages(utils::globalVariables(c("..cols", "Chromosome",
                                                "Start", "End", "Segment_Mean",
                                                "Variant_Type", "..maf_cols",
                                                "Tumor_Sample_Barcode",
                                                "Start_position", "t_alt_count",
                                                "t_ref_count", "eff.seg.len")))
