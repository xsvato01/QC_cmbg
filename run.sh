nextflow kuberun xsvato01/QC_cmbg -r main -pod-image 'cerit.io/nextflow/nextflow:22.11.1' \
	-resume -process.echo -with-report nf-report.html -c zaloha_nextflow.config  --datain /mnt/shared/MedGen/sequencing_results/primary_data/230922_Next_oI_NHA
