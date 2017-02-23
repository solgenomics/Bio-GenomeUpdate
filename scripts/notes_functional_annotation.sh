#Notes for ITAG3.0 Functional annotation 
#=======================================

#Author: Mirella Flores
#Description: For ITAG3.0 Functional annotation. Uses: compare_seq.pl, get_functional_annot.pl


#1. Run and filter InterproScan

#-Trim asterisks
	sed --i 's/*//g' ITAG3.0_proteins.fasta >ITAG3.0_proteins_woast.fasta 

#-Split fasta
	split_multifasta.pl --in ITAG3.0_proteins_woast.fasta  --output_dir=./ --f=ITAG3.0_proteins_ --seqs_per_file=1000

#-Run InterproScan
	i=1
	while [ $i -le 35 ]
	do
	sh interproscan-5.22-61.0/interproscan.sh -goterms -iprlookup -pa -i ITAG3.0_proteins_$i.fsa -f tsv,gff3
	echo "Round $i finished"
	i=$(( i+5 ))
	done

#-Filter Interpro

	grep "GO:" ITAG3.0_proteins.fsa.tsv | awk -F'\t' '{if (($8-$7)/$3 > 0.5){ print $1"\t"$12"\t"$14"\t"($8-$7)/$3"\t"$9 } }' | uniq > pre.txt
	
	awk '!seen[$1]++' pre.txt > ITAG3.0_proteins.interpro_filtered.txt
	
	awk -F'\t' '{if (($8-$7)/$3 > 0.5){ print $1"\t"$12"\t"$4"\t"$5 }}' ITAG3.0_proteins.fsa.tsv | sort -k1 |  grep "Pfam" | awk '!a[$1]++' >ITAG3.0_proteins.pfam_interpro.txt
	awk -F'\t' '{if (($8-$7)/$3 > 0.5){ print $1"\t"$12"\t"$14 }}' ITAG3.0_proteins.fsa.tsv | sort -k1 |  grep "GO:" | awk '!a[$1]++' >ITAG3.0_proteins.go_interpro.txt


#2. To get Solyc without version and match NCBI 

	paste  <( cut -f1,2 /export/projects/SL3.0/ITAG3.0/v1_ITAG3.0_without_functional_annotation/ncbiID_Solyc.txt | sed 's/mRNA://g') <( sed 's/\.[0-9]//g' /export/projects/SL3.0/ITAG3.0/v1_ITAG3.0_without_functional_annotation/ncbiID_Solyc.txt | cut -f2 | sed 's/mRNA://g') <( cut -d " " -f2 /export/projects/SL3.0/ITAG3.0/v1_ITAG3.0_without_functional_annotation/ncbiID_Solyc.txt) | column -s $'\t' -t >ncbi_match.txt


#3. To get NCBI description with match (2)

	join -1 1 -2 3 <(sort -k1 ncbi_match.txt) <(sort -k3 /export/projects/SL3.0/ITAG3.0/v1_ITAG3.0_without_functional_annotation/curated_proteins_NCBI.txt) | tac | awk '{out=$2"\t"$1"\t"$3"\t"$4"\t"; for(i=7;i<=NF;i++){out=out" "$i}; print out}' > ncbi_annot.tmp

	awk '!a[$1]++' ncbi_annot.tmp > ncbi_annot.txt


#4. Get new version 

#	- Run script to compare sequences and get a list 

	./compare_seq.pl

	cat <(grep "0\.2\.1" listSolyc.txt) <(grep "0\.1\.1" listSolyc.txt) | cut -f1 > newversion.txt

	grep -Fxv -f list.txt newversion.txt >newversion2.txt 
	
	paste <(sed 's/\.[0-9]//g' newversion2.txt) <(awk -F'.' '{ print $0" "$1"."$2+1}' newversion2.txt) | awk -F' ' '{print $1"\t"$2"\t"$3}' > list4versioning.txt

#5. Get loci file

#6. merge loci - ncbi

#	- Looks for in second file which ones don't match first file
	
	grep -Fxv -f <(cut -f1 solyc_curated_loci.txt | sed 's/ //g' ) <(cut -f3 ncbi_annot.txt | sed 's/ //g' ) >list_wo_loci.txt

#	-Filter no loci annotation
	
	awk 'FNR==NR{a[$1];next} !($3 in a)' list_wo_loci.txt ncbi_annot.txt > ncbi_woloci_annot.txt


#7. Filter AHRD

	grep -v "Unknown protein" ahrd_ITAG_tomato_output.csv | awk -F"\t" '{ print $1"\t"$4" (AHRD V3.3 "$3"\t"$2 }' >ahrd_ITAG_filtered.txt
	awk -F"\t" '{ print $1"\t"$4" (AHRD V3.3 "$3"\t"$2 }' ahrd_ITAG_tomato_output.csv >ahrd_ITAG_filtered.txt

#	And delete header

#8. To create a list 

	cut -f9 /export/projects/SL3.0/ITAG3.0/v1_ITAG3.0_without_functional_annotation/ITAG3.0_gene_models.gff | grep "ID=mRNA" | cut -d";" -f1,4 | sed 's/ID=mRNA://g' | sed 's/_AED=//g' | awk -F";" '{ if ($2>=1){ print $1"\t"$2"\tNote=LOW QUALITY:";} else {print $1"\t"$2"\tNote=";}  }' >listSolyc.txt    

#8. Main script to merge all previous results
 
#	Need all these files: (change names as needed in the script)
#	ahrd_ITAG_filtered.txt
#	solyc_curated_loci.txt
#	ncbi_woloci_annot.txt
#	ITAG3.0_proteins.go_interpro.txt
#	ITAG3.0_proteins.pfam_interpro.txt

	./get_functional_annot.pl

#9. Create dropped files from ITAG2.4

	grep -Fv -f oldgenes.txt ITAG2.4_gene_models.gff3

