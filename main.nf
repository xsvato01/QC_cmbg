run = "${params.datain}".split("/")
run = run[run.size()-2]
launchDir = "${launchDir}/${run}"


process FIRST_ALIGN_BAM {
	tag "first align on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/mapped/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
    tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4 ${params.refindex} $reads \
	| samtools view -Sb -o - -| samtools sort -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai	
	"""
}



process COVERAGE {
	tag "Creating coverage on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	tuple val(name), path("*.txt")
	
	script:
	"""
	bedtools coverage -abam ${params.covbed} -b $bam -d > ${name}.PBcoverage.txt   
	"""
}

process COVERAGE_STATS {
	tag "Creating coverage stats on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(annotated_normed)
	
	script:
	"""
	Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $run  
	"""
}


process MULTIQC {
	tag "first QC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	path "*"

	script:
	"""
	qualimap bamqc -bam $bam -outdir ./
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
    picard BedToIntervalList -I ${params.covbedpicard} -O ${name}.interval_list -SD ${params.ref}.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}.fa -O ${name}.aln_metrics
	multiqc . -n report.html
	"""

}

 
workflow {
 rawfastq = channel.fromFilePairs("${params.datain}/BTK*R{1,2}*", checkIfExists: true)
 
	sortedbam	= FIRST_ALIGN_BAM(rawfastq)
	covered		= COVERAGE(sortedbam[0])
	COVERAGE_STATS(covered)					
	MULTIQC(sortedbam[0])	
}
