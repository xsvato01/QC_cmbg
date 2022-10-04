library(formattable)
install.packages("https://cran.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_2.1.tar.gz", repos=NULL, type="source")
install.packages("htmlwidgets")

table = read.table('/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/mutliqc/AteroHemo_gene_coverage.txt',
                   header = TRUE)

tbl_summary = table %>% group_by(Exon) %>%
  summarise(Cov_min_total = min(Coverage_min),
            Cov_max_total = max(Coverage_max),
            std_median = sd(Coverage_median),
            mean_median = mean(Coverage_median),
            median_median = median(Coverage_median),
            mad = mad(Coverage_median)
            )

write.table(tbl_summary,"/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/mutliqc/Seq_stats.txt",
            sep="\t", quote = F, row.names = F)

table_one = table %>% arrange(Exon) %>%
  mutate(zcore_mad = (Coverage_median-tbl_summary$median_median)/tbl_summary$mad)

final = table %>% arrange(Sample, Exon) %>%
  mutate(zscore_mad = round(((Coverage_median-tbl_summary$median_median)/tbl_summary$mad), digits = 2)) %>% 
  relocate(Sample, .before = Exon) %>% filter(zscore_mad < -3)

ft = formattable(final, list(
  Exon.1 = FALSE,
  Coverage_mean = FALSE,
  zscore_mad = color_tile("transparent", "lightpink")
  ))

htmlwidgets::saveWidget(as.htmlwidget(ft), file="PerExon_mqc.html")

