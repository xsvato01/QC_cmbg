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

process FIRST_ALIGN_BAM {
	tag "first align on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/mapped/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
//	path "*" 
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4*${task.cpus} ${params.refindex} $reads > ${name}.sam 
	samtools view -Sb ${name}.sam -o ${name}.bam
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
	multiqc -c ${params.QCconfig} . -n report.html
	"""

}

 
workflow {
 rawfastq =	channel.fromFilePairs("${params.datain}/aB*R{1,2}*", checkIfExists: true)
 fastqced =	FASTQC(rawfastq)
 sortedbam =	FIRST_ALIGN_BAM(rawfastq)
 stats =	COVERAGE_STATS(sortedbam[0])			
//stats[0].collect().view()
		
  MULTIQC(stats[0].mix(fastqced).collect())	
// MULTIQC(stats[0].collect())	
}
