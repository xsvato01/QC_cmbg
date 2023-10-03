run = "${params.datain}".split("/")
run = run[run.size()-1]
launchDir = "${launchDir}/${run}"


process FASTQC {
	tag "FASTQC on $name using $task.cpus CPUs and $task.memory memory"
	//publishDir  "${launchDir}/fastQC", mode:'copy'

	input:
	tuple val(type), val(name), path(reads)

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
	// publishDir "${launchDir}/mapped/", mode:'copy'

	input:
	tuple val(type), val(name), path(reads)

	output:
	tuple val(type), val(name), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	#minimap2 -a ${params.refgpu} $reads > ${name}.sam
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
	# echo $PATH
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
	tuple val(type), val(name), path(bam)

	output:
	tuple val(type), val(name), path("${name}.sorted.bam")
	tuple val(type), val(name), path("${name}.sorted.bai")

	script:
	"""
	samtools sort ${name}.bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process COVERAGE_STATS {
	tag "Creating coverage stats on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy', pattern: "*PBcoverage.txt"

	input:
	tuple val(type), val(name), path(bam), val(specific_bed)
	//val ('specific_bed')

	output:
	path "*"
	
	script:
	"""
	bedtools coverage -abam ${specific_bed} -b $bam -d > ${name}.PBcoverage.txt
	qualimap bamqc -bam $bam -gff ${params.covbed} -outdir ${name} -outfile ${name}.qualimap -outformat HTML
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard BedToIntervalList -I ${params.covbed} -O ${name}.interval_list -SD ${params.ref}/GRCh38-p10.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}/GRCh38-p10.fa -O ${name}.aln_metrics
	rm -r ?
	"""
}


process MULTIQC {
	tag "first QC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/multiqc/", mode:'copy'

	input:
	path ("QC_results/*")

	output:
	path "*"

	script:
	"""
	Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $run
	Rscript --vanilla ${params.covcompare} $run ${params.seqstats}
	cat PerExon.html | perl  -pe 's/[^[:ascii:]]//g;' > PerExon_mqc.html
	echo "log_filesize_limit: 30000000" > multiqc_config.yaml
	multiqc . -n report_coverage.html
	cp report_coverage.html ${params.datain}
	"""

}


workflow {
	Atero = Channel
    .fromFilePairs("${params.datain}/raw_fastq/a*R{1,2}*" )
    .map {it -> [ "Atero", it[0], it[1]]}
 	Hemato =	Channel
    .fromFilePairs("${params.datain}/raw_fastq/h*R{1,2}*" )
    .map {it -> [ "Hemato", it[0], it[1]]}
		Fialka =	Channel
    .fromFilePairs("${params.datain}/raw_fastq/f*R{1,2}*" )
    .map {it -> [ "Fialka", it[0], it[1]]}
		Pank =	Channel
    .fromFilePairs("${params.datain}/raw_fastq/p*R{1,2}*" )
    .map {it -> [ "Pank", it[0], it[1]]}

 	Concat = Atero.mix(Hemato,Pank,Fialka)
		Concat.view()

 fastqced =	FASTQC(Concat)
 bam =		ALIGN_CPU(Concat)
 sortedbam =	SORT_INDEX(bam)
	sortedbam[0].view()

	sortedbam[0].branch {
        Atero: it[0] == "Atero"
									return [it[0], it[1], it[2], params.Atero23_re]
        Hemato: it[0] == "Hemato"
									return [it[0], it[1], it[2],  params.Hemato23_re]
							Pank: it[0] == "Pank"
									return [it[0], it[1], it[2],  params.Pank23]
							Fialka: it[0] == "Fialka"
									return [it[0], it[1], it[2],  params.Fialka23]
    }
    .set{sorted}

		// sortedbam[0].map({
  //       if it[0] == "Atero"
		// 							return [it, params.Atero23_re]
  //       if it[0] == "Hemato"
		// 							return [it, params.Hemato23_re]
  //   }).view()

	 sorted.Atero.view { "$it is Atero" }
  sorted.Hemato.view {"$it is Hemato"}
		sorted.Pank.view {"$it is Pank"}
		sorted.Fialka.view {"$it is Fialka"}
	 sorted = sorted.Atero.mix(sorted.Hemato,sorted.Pank,sorted.Fialka)
										//.view{"$it is merged"}
  stats =	COVERAGE_STATS(sorted)
  MULTIQC(stats[0].mix(fastqced).collect())
}
