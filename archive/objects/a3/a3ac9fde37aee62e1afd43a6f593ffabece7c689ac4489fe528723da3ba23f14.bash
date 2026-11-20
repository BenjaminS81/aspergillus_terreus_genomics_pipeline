for dir in */; do
	cd "$dir" || continue
	echo "$dir"

	if [[ -f "contigs.fasta" && ! -d ater ]]; then
		conda run -n seqkit seqkit replace -p ".+" -r "contig_{nr}" contigs.fasta > contigs_renamed.fasta
		conda run -n repeats BuildDatabase -name repeats_db contigs_renamed.fasta
		conda run -n repeats RepeatModeler -database repeats_db -threads 16
		conda run -n repeats RepeatMasker -pa 16 -lib repeats_db-families.fa -xsmall contigs_renamed.fasta
		conda run -n funannotate2 funannotate2 clean -f contigs_renamed.fasta.masked -o cleaned_masked_contigs.fasta
		conda run -n funannotate2 funannotate2 train -f cleaned_masked_contigs.fasta -s "Aspergillus alabamensis" -o aalb
		conda run -n funannotate2 funannotate2 predict -i aalb
		conda run -n funannotate2 funannotate2 annotate -i aalb
	fi
	cd ..
done
