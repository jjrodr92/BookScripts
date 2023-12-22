#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jlrodr92@gmail.com
#SBATCH --mem=140GB
#SBATCH -n 40   
#SBATCH -J VIRALDB
#SBATCH -p medium

DIA=`date +"%d/%m/%Y"`
HORA=`date +"%H:%M"`
echo "Start $DIA $HORA"
echo "Sample $1"

db="/home/jrodriguez/Julio/DB/"
dir="/home/jrodriguez/Julio/PMS_files"
dir1="/home/jrodriguez/Julio/PMS_files/fastq"
dir3="/home/jrodriguez/Julio/sambam/NEW"
results="/home/jrodriguez/Julio/blast"
results1="/home/jrodriguez/Julio/fasta"
results2="/home/jrodriguez/Julio/tables"
results3="/home/jrodriguez/Julio/blast/$1"
results4="/home/jrodriguez/Julio/blast/$1/aftercap"
results5="/home/jrodriguez/Julio/blast/$1/aftercap/NCBI"
mkdir $results3
mkdir $results4
mkdir $results5
mkdir $dir3

##Step1
#module load gcc/6.4.0 diamond/0.9.22
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1

# Blast

diamond blastx -d $db/marcoDB1 -q $dir/Trinity_$1.fasta -o $results/$1-ViralaliNEW.blast -f 0 -k 2 -e 0.00001 --unal 0 --more-sensitive

diamond blastx -d $db/marcoDB1 -q $dir/Trinity_$1.fasta -o $results/$1-ViralNEW.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 2 -e 0.00001 --unal 0 --more-sensitive

#Table of results
mkdir $results3/temp
grep -e ">" $results/$1-ViralaliNEW.blast > $results3/list.txt
cut -d "[" -f 2 $results3/list.txt > $results3/list1.txt
sed -i 's/]/ /g' $results3/list1.txt
cut -d "[" -f 1 $results3/list.txt > $results3/list2.txt
sed -i 's/>//g' $results3/list2.txt
paste $results/$1-ViralNEW.blast $results3/list2.txt > $results3/NCBItable1.blast
paste $results3/NCBItable1.blast $results3/list1.txt > $results3/NCBItable2.blast
cat <(head -1 header) $results3/NCBItable2.blast > $results2/$1_ViraltablefinalNEW.blast

#Extract list1
grep -e "Query=" $results/$1-ViralaliNEW.blast > $results3/temp/Scaffoldslistnew.txt
cat $results3/temp/Scaffoldslistnew.txt | sort | uniq > $results3/temp/Scaffoldslistnewuni.txt
cut  -d " " -f 2 $results3/temp/Scaffoldslistnewuni.txt > $results3/$1-list1.txt
echo "Number of list1"
grep -o "TRINITY" $results3/$1-list1.txt | wc -l

#Extract fasta from Trinity
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $dir/Trinity_$1.fasta | perl -pe "s/>/\n>/g" > $dir/Trinity_$1_oneline.fasta
grep -A1 --no-group-separator -f $results3/$1-list1.txt $dir/Trinity_$1_oneline.fasta > $dir/$1-significant1NEW.fasta
echo "Number in fasta file"
grep -o ">" $dir/$1-significant1NEW.fasta | wc -l

# CAP3
#module load cap3/10.2011-1-intel-x86_64
source activate jrodriguez 
an
echo "CAP3 processing"
mkdir $dir/cap3_data_$1
mv $dir/$1-significant1NEW.fasta $dir/cap3_data_$1/
cap3 $dir/cap3_data_$1/$1-significant1NEW.fasta > $dir/cap3_data_$1/$1-significant1NEW.cap3
#srun mv $dir/*cap.ace $dir/*cap.contigs $dir/*cap.contigs.links $dir/*cap.contigs.qual $dir/*cap.info $dir/*cap.singlets $dir/cap3_data_$1/
cat $dir/cap3_data_$1/*cap.singlets $dir/cap3_data_$1/*cap.contigs > $dir/$1_mapNEW.fasta

module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
# Blast after Cap3
diamond blastx -d $db/marcoDB1 -q $dir/$1_mapNEW.fasta -o $results/$1-ViralaliNEWafterCAP.blast -f 0 -k 2 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -d $db/marcoDB1 -q $dir/$1_mapNEW.fasta -o $results/$1-ViralNEWafterCAP.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 2 -e 0.00001 --unal 0 --more-sensitive

#Table of results

grep -e ">" $results/$1-ViralaliNEWafterCAP.blast > $results4/list.txt
cut -d "[" -f 2 $results4/list.txt > $results4/list1.txt
sed -i 's/]/ /g' $results4/list1.txt
cut -d "[" -f 1 $results4/list.txt > $results4/list2.txt
sed -i 's/>//g' $results4/list2.txt
paste $results/$1-ViralNEWafterCAP.blast $results4/list2.txt > $results4/NCBItable1.blast
paste $results4/NCBItable1.blast $results4/list1.txt > $results4/NCBItable2.blast
cat <(head -1 header) $results4/NCBItable2.blast > $results2/$1_ViraltablefinalNEWafterCAP.blast

#Extract list2
grep -e "Query=" $results/$1-ViralaliNEWafterCAP.blast > $results4/$1-list2.txt
cat $results4/$1-list2.txt | sort | uniq > $results4/$1-list2uni.txt
cut -d " " -f 2 $results4/$1-list2uni.txt > $results4/$1-list2.txt

echo "Number of list2"
grep -o -e "TRINITY" -e "Contig" $results4/$1-list2.txt | wc -l

#Extract fasta from CAP3 results
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $dir/$1_mapNEW.fasta | perl -pe "s/>/\n>/g" > $dir/$1_mapNEW_oneline.fasta
grep -A1 --no-group-separator -f $results4/$1-list2.txt $dir/$1_mapNEW_oneline.fasta > $dir/$1-significant2NEW.fasta
echo "Number fasta file"
grep -o ">" $dir/$1-significant2NEW.fasta | wc -l

##Step2 BLAST TO NCBI
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
diamond blastx -d $db/NEWNR/nr -q $dir/$1-significant2NEW.fasta -o $results/$1-NCBINEW.blast -f 0 -k 2 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -d $db/NEWNR/nr -q $dir/$1-significant2NEW.fasta -o $results/$1-NCBItab1NEW.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 2 -e 0.00001 --unal 0 --more-sensitive

#Table of results
grep -e ">" $results/$1-NCBINEW.blast > $results5/list.txt
cut -d "[" -f 2 $results5/list.txt > $results5/list1.txt
sed -i 's/]/ /g' $results5/list1.txt
cut -d "[" -f 1 $results5/list.txt > $results5/list2.txt
sed -i 's/>//g' $results5/list2.txt
paste $results/$1-NCBItab1NEW.blast $results5/list2.txt > $results5/$1-NCBItable1.blast
paste $results5/$1-NCBItable1.blast $results5/list1.txt > $results5/$1-NCBItable2.blast
cat <(head -1 header) $results5/$1-NCBItable2.blast > $results2/$1-NCBItablefinalNEW.blast


# Mapping
module load BWA/0.7.17-GCC-8.2.0-2.31.1

bwa index $dir/$1-significant2NEW.fasta
bwa mem -t 40 $dir/$1-significant2NEW.fasta $dir1/$1_map.fastq.gz > $dir3/$1-significant2NEW.sam

module load SAMtools/1.9-GCC-8.2.0-2.31.1

samtools view -bS $dir3/$1-significant2NEW.sam  -o $dir3/$1-significant2NEW.bam
samtools sort -@ 40 $dir3/$1-significant2NEW.bam -o $dir3/$1-significant2NEW.sort.bam
samtools index $dir3/$1-significant2NEW.sort.bam
samtools idxstats $dir3/$1-significant2NEW.sort.bam > $results2/$1_idxstats.txt
#rm $results1/$1-significant2.fasta.*
rm $dir3/$1-significant2NEW.sam

DIA=`date +"%d/%m/%Y"`
HORA=`date +"%H:%M"`
echo "End $DIA $HORA"
echo "Sample $1"