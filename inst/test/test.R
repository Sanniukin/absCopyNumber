test_seg = data.table::fread(input = "gunzip -c inst/extdata/SNP6_solid_tumor.seg.txt.gz")
test_maf = data.table::fread(input = "gunzip -c inst/extdata/SNP6_solid_tumor.maf.txt.gz")
test_cn = data.table::fread(input = "gunzip -c inst/extdata/example.cn.txt.gz")
test_snv = data.table::fread(input = "gunzip -c inst/extdata/example.snv.txt.gz")

# test basic input as data.frame
abs_initialize(seg = test_seg, verbose = T)
abs_initialize(seg = test_seg, verbose = F)
abs_initialize(seg = test_cn, verbose = T)

# test basic input (add snv file)
abs_initialize(seg = test_cn, snv = test_snv, verbose = T)
abs_initialize(seg = test_seg, snv = test_maf, isMaf = T, verbose = T)

# test basic input (set sample name)
res1 = abs_initialize(seg = test_cn, snv = test_snv, sample_seg = "test1", sample_snv = "test2", verbose = T)

# test basic input (set sample column)
test_cn2 = data.table::copy(test_cn)
test_cn2[, yyy:="yyy"]
test_snv2 = data.table::copy(test_snv)
test_snv2[, zzz:="zzz"]

res2 = abs_initialize(seg = test_cn, snv = test_snv, sample_seg = "yyy", sample_snv = "zzz", verbose = T)

#------ test file path
file_seg = "inst/extdata/SNP6_solid_tumor.seg.txt.gz"
file_maf = "inst/extdata/SNP6_solid_tumor.maf.txt.gz"
file_cn =  "inst/extdata/example.cn.txt.gz"
file_snv = "inst/extdata/example.snv.txt.gz"

# test basic input as file
abs_initialize(seg = file_seg, verbose = T)
abs_initialize(seg = file_seg, verbose = F)
abs_initialize(seg = file_cn, verbose = T)

# test basic input (add snv file)
abs_initialize(seg = file_cn, snv = file_snv, verbose = T)
abs_initialize(seg = file_seg, snv = file_maf, isMaf = T, verbose = T)

# test basic input (set sample name)
res3 = abs_initialize(seg = file_cn, snv = file_snv, sample_seg = "test1", sample_snv = "test2", verbose = T)

identical(res1, res3)

res4 = abs_initialize(seg = file_seg, snv = file_maf, isMaf = TRUE, sample_seg = "Sample",
                      verbose = T)

# test abs_calling
test_cn = data.table::fread(input = "gunzip -c inst/extdata/example.cn.txt.gz")
test_snv = data.table::fread(input = "gunzip -c inst/extdata/example.snv.txt.gz")
test_cn2 = data.table::copy(test_cn)
test_snv2 = data.table::copy(test_snv)
test_cn$sample = test_snv$sample = "sample1"
test_cn2$sample = test_snv2$sample = "sample2"

d2sample = abs_initialize(seg = rbind(test_cn, test_cn2), snv = rbind(test_snv, test_snv2))
object = data.table::copy(d2sample)
object = abs_prepare(object)

unique(subset(object, samples = "sample1")@data$sample)
unique(subset(object, samples = c("sample1", "sample2"))@data$sample)
