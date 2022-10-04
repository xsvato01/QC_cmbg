run = "${params.datain}".split("/")
run = run[run.size()-2]
launchDir = "${launchDir}/${run}"


process FASTQC {
	tag "FASTQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/fastQC", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	path '*'
	//path "fastqc_out/*"

	script:
	"""
	fastqc $reads -o ./
	rm -r ?
	"""
}

process ALIGN_CPU {
	tag "CPU align on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4*${task.cpus} ${params.refindex} $reads > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}


process ALIGN_GPU {
  tag "GPU align on $name using $task.cpus CPUs and $task.memory memory"
  publishDir "${launchDir}/mapped/", mode:'copy'
	accelerator 1
	container "cerit.io/ceitec/clara-parabricks:3.8.0-1.ampere"
	memory 20.GB

  input:
   tuple val(name), path(reads)

  output:
   tuple val(name), path("${name}.bam")

  script:
   """
   echo $PATH
	 /parabricks/run_pipeline.py fq2bam --ref ${params.refgpu} \
	 --in-fq $reads \
	 --out-bam ${name}.bam \
	 --tmp-dir .
   """
}

process SORT_INDEX {
  tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
  publishDir "${launchDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(bam)

	output:
  tuple val(name), path("${name}.sorted.bam")
  tuple val(name), path("${name}.sorted.bai")

	script:
	"""
  samtools sort ${name}.bam -o ${name}.sorted.bam
  samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process COVERAGE_STATS {
	tag "Creating coverage stats on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy'

	input:
	tuple val(name), path(bam)

	output:
	path "*"
	//path "*.{txt,interval_list,aln_metrics,flagstat,samstats}"

	script:
	"""
	bedtools coverage -abam ${params.covbed} -b $bam -d > ${name}.PBcoverage.txt
	qualimap bamqc -bam $bam -gff ${params.covbed} -outdir ${name} -outfile ${name}.qualimap -outformat HTML
	#mv ${name}/genome_results.txt ${name}.genome_results.txt
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
  picard BedToIntervalList -I ${params.covbed} -O ${name}.interval_list -SD ${params.ref}/GRCh38-p10.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}/GRCh38-p10.fa -O ${name}.aln_metrics
	rm -r ?
	"""
}


process MULTIQC {
	tag "first QC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/mutliqc/", mode:'copy'

	input:
	path ("QC_results/*")

	output:
	path "*"

	script:
	"""
	Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $run
	Rscript --vanilla ${params.covcompare} "AteroHemo_gene_coverage.txt" ${params.seqstats}
	multiqc -c ${params.QCconfig} . -n report.html
	"""

}


workflow {
 rawfastq =	channel.fromFilePairs("${params.datain}/*R{1,2}*", checkIfExists: true)
 fastqced =	FASTQC(rawfastq)
 bam =		ALIGN_CPU(rawfastq)
 sortedbam =	SORT_INDEX(bam)
 stats =	COVERAGE_STATS(sortedbam[0])
//stats[0].collect().view()

  MULTIQC(stats[0].mix(fastqced).collect())
// MULTIQC(stats[0].collect())
}
