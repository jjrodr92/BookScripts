#!/bin/sh 
resources="./bbmap-39.01-0/resources/"

#Sket
sendsketch.sh in=./sampleID/sampleID_1.fastq.gz out=./sampleID/sketch_1.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus nt ow
sendsketch.sh in=./sampleID/sampleID_2.fastq.gz out=./sampleID/sketch_2.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus nt ow
sendsketch.sh in=./sampleID/sampleID_1.fastq.gz out=./sampleID/sketch_1.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus refseq append
sendsketch.sh in=./sampleID/sampleID_2.fastq.gz out=./sampleID/sketch_2.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus refseq append
sendsketch.sh in=./sampleID/sampleID_1.fastq.gz out=./sampleID/sketch_1.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus silva minkeycount=2 append
sendsketch.sh in=./sampleID/sampleID_2.fastq.gz out=./sampleID/sketch_2.txt reads=200000 fname=sampleID.fastq.gz minprob=0.2 samplerate=1.0 merge printname0=f records=20 overwrite=true color=false depth depth2 unique2 volume sortbyvolume contam2=genus silva minkeycount=2 append

#Clumpify
clumpify.sh pigz=t unpigz=t zl=4 reorder in1=./sampleID/sampleID_1.fastq.gz out1=./sampleID/TEMP1_1.fastq.gz in2=./sampleID/sampleID_2.fastq.gz out2=./sampleID/TEMP1_2.fastq.gz passes=1
#CLEAN Step1 - adapter removal
bbduk.sh -threads=8 ktrim=r ordered minlen=49 minlenfraction=0.33 mink=11 tbo tpe rcomp=t overwrite=true k=23 hdist=1 hdist2=1 ftm=5 pigz=t unpigz=t zl=4 ow=true in1=./sampleID/TEMP1_1.fastq.gz out1=./sampleID/TEMP2_1.fastq.gz in2=./sampleID/TEMP1_2.fastq.gz out2=./sampleID/TEMP2_2.fastq.gz rqc=hashmap outduk=./sampleID/ktrim_kmerStats1.txt stats=./sampleID/ktrim_scaffoldStats1.txt loglog ref=resources/adapters.fa

#CLEAN Step2 - artefacts removal
bbduk.sh -threads=8 maq=10,0 trimq=25 qtrim=r ordered overwrite=true maxns=1 minlen=49 minlenfraction=0.33 k=25 hdist=1 pigz=t unpigz=t zl=6 cf=t barcodefilter=crash ow=true in1=./sampleID/TEMP2_1.fastq.gz out1=./sampleID/TEMP3_1.fastq.gz in2=./sampleID/TEMP2_2.fastq.gz out2=./sampleID/TEMP3_2.fastq.gz outm=./sampleID/synth1.fq.gz rqc=hashmap outduk=./sampleID/kmerStats1.txt stats=./sampleID/scaffoldStats1.txt loglog ref=$resources/phix_adapters.fa.gz,$resources/lambda.fa.gz,$resources/sequencing_artifacts.fa.gz

#CLEAN Step3 - short removal
bbduk.sh -threads=8 ordered overwrite=true k=20 hdist=1 pigz=t unpigz=t zl=6 ow=true in1=./sampleID/TEMP3_1.fastq.gz out1=./sampleID/TEMP4_1.fastq.gz in2=./sampleID/TEMP3_2.fastq.gz out2=./sampleID/TEMP4_2.fastq.gz outm=./sampleID/synth2.fq.gz outduk=./sampleID/kmerStats2.txt stats=./sampleID/scaffoldStats2.txt loglog ref=$resources/short.fa

#CLEAN Step4 – ribosomal removal
bbduk.sh -threads=8 ordered k=31 ref=$resources/ribokmers.fa.gz ow=true in1=./sampleID/TEMP4_1.fastq.gz out1=./sampleID/TEMP5_1.fastq.gz in2=./sampleID/TEMP4_2.fastq.gz out2=./sampleID/TEMP5_2.fastq.gz outm=./sampleID/ribo.fq.gz outduk=./sampleID/ribo_Stats1.txt stats=./sampleID/ribo_Stats2.txt
 
#MERGE PAIRED-END
bbmerge.sh -threads=8 loose overwrite=true in1=./sampleID/TEMP5_1.fastq.gz in2=./sampleID/TEMP5_2.fastq.gz ihist=./sampleID/ihist_merge.txt outc=./sampleID/cardinality.txt pigz=t unpigz=t zl=9 adapters=$resources/adapters.fa
reformat.sh in=./sampleID/TEMP5_#.fastq.gz out=./sampleID/sampleID_map.fastq.gz
