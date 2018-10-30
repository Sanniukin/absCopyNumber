
#' Print information about input of absCopyNumber and extra message
#'
#' @return Nothing
#' @export
#'
#' @examples
#' abs_supportfiles()
abs_supportfiles = function(){
    cat("===================================================================\n")
    cat("================= Welcome to use absCopyNumber ====================\n")
    cat("===================================================================\n")
    cat("Files metioned below are TAB separated, comment with # is allowed.\n")
    cat("\n")
    cat("===> 2 type of support CNV file format:\n")
    cat("#1 --- Standard input (muliple samples best has a 'sample' column)\n")
    cat("chrom loc.start   loc.end eff.seg.len normalized.ratio\n")
    cat("   10     92986  39138142      814538        0.9074532\n")
    cat("   10  42356124 135491096     2462313        0.9013927")
    cat("#2 --- Standard segmentation file\n")
    cat("sample Chromosome    Start      End Num_Probes Segment_Mean\n")
    cat(" V-01           1  3218610 14449771       5840      -0.1329\n")
    cat(" V-01           1 14450545 14453003          2      -1.9378\n")
    cat("\n")
    cat("===> 2 type of support SNV file format:\n")
    cat("#1 --- Standard input (muliple samples best has a 'sample' column)\n")
    cat("chrom  position tumor_var_freq\n")
    cat("   14 106452833         0.4349\n")
    cat("   12 113440906         0.3132\n")
    cat("#2 --- Standard MAF file\n")
    cat("detail please see https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/\n")
    cat("Minimal required columns:\n")
    cat("\"Tumor_Sample_Barcode\", \"Chromosome\", \"Start_position\", \"Variant_Type\", \"t_ref_count\", \"t_alt_count\"\n")

    cat("==============\n")
    cat("Notes:\n")
    cat("1) chrom       -- The chromosome number of a segment, integer number or character is allowed.\n")
    cat("2) loc.start   -- The start position of a segment.\n")
    cat("3) loc.end     -- The end position of a segment.\n")
    cat("4) eff.seg.len -- For exome sequencing, due to the nature of highly uneven coverage\n")
    cat("                  (zero coverage for introns), this column gives the number of base pairs\n")
    cat("                  with actual observed coverage. can be derived by concatenate all the VARSCAN bins\n")
    cat("                  within a segment.\n")
    cat("5) normalized.ratio -- The mean copy ratio (tumor DNA vs. germline DNA) of a segment.\n")
    cat("                       Note that the copy ratio should be normalized to eliminate any sequencing\n")
    cat("                       throughput difference between tumor and germline DNA.\n")
    cat("                       For example, samtools can be used to count total reads that were properly\n")
    cat("                                    paired/aligned for tumor and germline DNA. The difference then\n")
    cat("                                   need to be adjusted accordingly.\n")
    cat("6) tumor_var_freq  --  The proportion of reads supporting the somatic SNV allele. Must be a fraction number.\n")

    cat("==============\n")
    cat("END.")
}
