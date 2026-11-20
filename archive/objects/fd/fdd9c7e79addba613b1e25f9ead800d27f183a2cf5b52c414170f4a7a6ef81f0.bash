for file in *.fasta; do
	echo "$file"

	if [[ ! -d ater ]]; then
		#conda run -n seqkit seqkit replace -p ".+" -r "contig_{nr}" contigs.fasta > contigs_renamed.fasta
		conda run -n repeats BuildDatabase -name repeats_db $file
		conda run -n repeats RepeatModeler -database repeats_db_$file -threads 16
		conda run -n repeats RepeatMasker -pa 16 -lib repeats_db-families.fa -xsmall $file_renamed.fasta
		conda run -n funannotate2 funannotate2 clean -f $file_renamed.fasta -o cleaned_$file_renamed.fasta
		conda run -n funannotate2 funannotate2 train -f $file_renamed.fasta -s "Aspergillus terreus" -o $file_ater
		conda run -n funannotate2 funannotate2 predict -i $file_ater
		conda run -n funannotate2 funannotate2 annotate -i $file_ater
	fi
done
