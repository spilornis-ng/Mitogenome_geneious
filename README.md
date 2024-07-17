##A crude pipeline for mitogenome recovery from low-coverage WGS using reference mitogenome
##This pipeline uses Geneious prime
##You need trimmed reads for running this pipeline
##Here the example name Tepo refers to Tephrodornis pondicerianus. Use whatever code or sample name you want to.

##Checklist of softwares required
1. bwa
2. bedtools
3. samtools
4. qualimap
5. Geneious or any open source software

##Steps 1-5 will require a linux OS or WSL, and Step 6 requires geneious in this particular case, hence can be done on any computer with the software


##########################


##Step1 - Download the reference mitogenome from NCBI or other online sources in fasta format. Closer the species to your taxa, the better.

##Step2 - Index the reference mitogenome (BWA)

```
bwa index ~/path/to/your/reference.fasta/
```
##Step3 - Map the trimmed reads of low-coverage WGS to the indexed mitogenome reference (BWA and samtools)
```
bwa mem -t 20 -M -R "@RG\tID:Tepo_5x\tSM:Tepo_\tLB:IlluminaWGS\tPL:ILLUMINA" \ #here I have used -t flag to specify cores, -R to add read-group to bam
~/path/to/your/indexed/reference.fasta/ \
~/path/to/your/R1.fq.gz \
~/path/to/your/R2.fq.gz | \
samtools view -bh - | \
samtools sort -@20 -T tmp -n -o Tepo_mapped_sorted.bam
```
##Step4 - Run BAMQC to get an idea of mapping coverage and quality (qualimap)
```
qualimap bamqc -bam Tepo_mapped_sorted.bam -outfile results.pdf
```
##Inspect the output to understand the coverage and quality of the reads mapped to the reference mitogenome

##Step5 - This maybe a bit complex, but it just seemed easier for me to think at the time.
##We extract only the reads that map to the mitogenome reference. This will lead to fairly small and manageable R1 and R2 fastq files.
##Using this it makes our downstream work less intensive as we will discard all the reads that do not map to the mitogenome (>99 percent of total reads)
```
bedtools bamtofastq -i Tepo_mapped_sorted.bam -fq mitomapped_R1.fq.gz -fq2 mitomapped_R2.fq.gz
```
##Step6 - Using Geneious prime. This is definitely redundant but since we are dealing with a small data size now, Geneious makes visualizations convinient.

##Load your reference mitogenome and mitomapped_R1.fq.gz and mitomapped_R2.fq.gz files into geneious

##Map the R1 and R2 files to the reference in geneious. For this I use default settings and use geneious aligner.

##You can now visualize where all the reads map on the reference mitogenome.

##Once you get a good idea, you can generate a consensus sequence.

##After generating the consensus sequence, transfer the annotations from the reference to your consensus sequence.

##Make sure to use a similarity of 40-50% to ensure that annotations are transfered from your reference to the consensus sequence.

##Your annotated consensus sequence will be the mitogenome you have extracted from your WGS dataset. 
