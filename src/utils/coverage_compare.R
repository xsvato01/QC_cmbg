library(data.table)
library(formattable)

args = commandArgs(trailingOnly=TRUE)
path = args[1]
path_stats = args[2]
# path = "/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/mutliqc/AteroHemo_gene_coverage.txt"
# path_stats = "/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/mutliqc/Seq_stats.txt"
SeqStats = fread(paste0(path_stats), sep="\t")

table = fread(paste0(path), sep="\t")


final = table %>% arrange(Sample, Exon) %>%
  mutate(zscore_mad = round(((Coverage_median-SeqStats$median_median)/SeqStats$mad), digits = 2)) %>% 
  filter(zscore_mad < -3)

ft = formattable(final, list(
  Coverage_mean = FALSE,
  zscore_mad = color_tile("transparent", "lightpink")
))

htmlwidgets::saveWidget(as.htmlwidget(ft), file="PerExon_mqc.html")