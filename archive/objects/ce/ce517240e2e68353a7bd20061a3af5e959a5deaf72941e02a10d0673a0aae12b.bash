for file in $(cat core_orthogroups.txt); do
	conda run -n mafft mafft --inputorder --adjustdirectionaccurately Orthogroup_Sequences/"$file".fa > mafft/"$file"_mafft.fa
	echo "$file".fa
done
