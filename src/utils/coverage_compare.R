library(data.table)
library(formattable)
library(dplyr)
library(DT)

print(list.files(pattern = "*_AteroHemo_gene_coverage.txt"))
args = commandArgs(trailingOnly=TRUE)
#path = args[1]
path_stats = args[2]
# path = "/storage/home/450402/000000-My_Documents/QC_test_run/work/3a/40fa8b18f655ffac59b23eb79d647b/Atero_test_AteroHemo_gene_coverage.txt"
# path_stats = "/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/mutliqc/Seq_stats.txt"
SeqStats = fread(paste0(path_stats), sep="\t")

table = fread(paste0(args[1],"_AteroHemo_gene_coverage.txt"), sep="\t")
exon_start = stringr::str_split(table$Exon, "_")


final = table %>% arrange(Sample, Exon) %>%
 # mutate(zscore_mad = round(((Coverage_median-SeqStats$median_median)/SeqStats$mad), digits = 2)) %>%
#  filter(zscore_mad < -1.5) %>%
  filter(Coverage_min < 30) %>%
  arrange(Coverage_min) %>%
  mutate(Start = sapply(stringr::str_split(Exon, "_"), "[", 2)) %>%
  mutate(Exon = sapply(stringr::str_split(Exon, "_"), "[", 1) )%>%
  relocate(Chrom, .before = Exon) %>%
  relocate(Start, .after = Exon) 


ft = formattable(final, list(
  Coverage_mean = FALSE,
  Coverage_min = color_tile("red", "lightpink")
))

htmlwidgets::saveWidget(as.datatable(ft), file="PerExon.html", title = "Statistika exonů")

#htmlwidgets::saveWidget(as.htmlwidget(ft), file="PerExon.html", title = "Statistika exonů")
