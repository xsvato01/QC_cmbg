library(data.table)
args = commandArgs(trailingOnly=TRUE)
path<-args[1]
#path<-"/storage/home/450402/000000-My_Documents/QC_k8s_testing/launch/000000-My_Documents/coverage"

samples_exons = list()
print(list.files(path, pattern = "*PBcoverage.txt"))
for(i in list.files(path, pattern = "*PBcoverage.txt")){
	print(i)  
  file <- fread(paste0(path,"/",i), sep="\t")
  colnames(file) <- c("Chrom", "Start", "End", "Exon", "Base_number", "Coverage")
  file_name_stat <- gsub(".PBcoverage.txt", "", i)
  file$sample <- file_name_stat
  exon_sample_stat <- file[sample==file_name_stat,.(Sample = file_name_stat, Chrom = Chrom[1],
                                                             Coverage_min = min(Coverage), Coverage_max = max(Coverage),
                                                            Coverage_mean = round(mean(Coverage), digit = 2),
                                                            Coverage_median = median(as.double(Coverage))),
                                                            by = Exon] %>% relocate(Sample, .before = Exon) 
  samples_exons = rbind(samples_exons,exon_sample_stat)
}

# write.table(samples_exons, file=paste0(args[2], ".perexon_samples.txt"),
#             sep="\t", quote = F, row.names = F)


write.table(samples_exons, file=paste0(path,"/",args[2], "_AteroHemo_gene_coverage.txt"),
            sep="\t", quote = F, row.names = F)
