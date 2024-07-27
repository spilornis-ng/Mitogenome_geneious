## Introduction
1. A crude pipeline for mitogenome recovery from low-coverage WGS using reference nuclear+mitogenome or just mitogenome
2. This pipeline uses Geneious prime
3. You need trimmed reads for running this pipeline

Checklist of softwares required
1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [qualimap](http://qualimap.conesalab.org/)
4. [Geneious Prime](https://www.geneious.com/)

While recovering mitogenome from whole genome resequencing data we need to address potential [NUMTs](https://en.wikipedia.org/wiki/Nuclear_mitochondrial_DNA_segment#:~:text=Nuclear%20mitochondrial%20DNA%20(NUMT)%20segments,nuclear%20genome%20of%20eukaryotic%20organisms) in our dataset. To do this we ideally require both nuclear and mitochondrial reference genomes for species of our interest or the next closest sister taxa. However, in most cases reference genomes for non-model species are not available (although this is changing rapidly). In such a case one can recover mitogenome using only mitochondrial reference of the target species or its closest sister, but this would lead to inclusion of NUMTs in the mitogenome assembly and any results of downstream analyses should be interpreted carefully. More or less the pipeline is similar for both the cases, and has been expanded upon below.

Although I have used geneious to generate a consensus and annotate, one can use any open source software to do the same.

Please note that Step 5 is common for both the case scenarios.

Scripts to run the pipeline on linux and HPC will be uploaded soon.

## Mitogenome recovery using reference nuclear and mitgenome

### Step 1 - Download the reference genomes and indexing (bwa) 
From NCBI or other online sources, download the nuclear and mitochondrial reference genome of your species of interest. In case no reference is available we can use a close sister species to our taxa.
Merge the nuclear reference and mitochondrial reference. This will allow us to reduce the proportion of NUMTs in our mitochondrial assembly downstream. 

```
cat nuclear_ref.fna mito_ref.fasta > final_ref.fna
```

Now we index the reference using bwa

```
bwa index ~/path/to/final_ref.fna
```


### Step 2 - Map the trimmed reads to the indexed reference (BWA and samtools) and extract only reads mapping to the mitogenome region

Now we will map the trimmed reads of our data to the indexed reference we generated in the previous step. To do that we will use bwa.

```
bwa mem -t 20 -M -R "@RG\tID:\tSM:\tLB:IlluminaWGS\tPL:ILLUMINA" \ #here I have used -t flag to specify cores, -R to add read-group to bam, do not forget to add details to the bam-header
~/path/to/indexed/final_ref.fna \
~/path/to/R1.fq.gz \
~/path/to/R2.fq.gz | \
samtools view -bh - | \
samtools sort -@20 -T tmp -o output_mapped_sorted.bam
```
We will now extract the reads mapping to only the mitochondrial region. To do that we will first check the header of our mitochondrial reference fasta to know what to look for in the generated bam file.

```
head /path/to/mito_ref.fasta
```
```
>NC_042191.1 Culicicapa ceylonensis mitochondrion, complete genome
GTCCCTGTAGCTTACAAAAAGCATGACACTGAAGATGTCAAGACGGCTGCCATATGCACCCAAGGACAAA
AGACTTAGTCCTAACCTTACTGTTGGTTTTTGCTAGACATATACATGCAAGTATCCGCGCGCCAGTGTAG
```
The head command will print out the first few lines of the fasta file. Note down the accession number afte **>** (in this case **NC_042191.1**). We are now going to use samtools to extract the reads mapping to this region.

```
samtools view output_mapped_sorted.bam "NC_042191.1" > output_mitomapped_sorted.bam
```
We now have the bam file consisting of only reads mapped to the mitogenome region. Do note that the file size reduces considerably as < 98% of the reads will map to the mitogenome.


### Step 3 - Run BAMQC to get an idea of mapping coverage and quality (qualimap)
Using qualimap bamqc we will assess the mapping coverage and quality
```
qualimap bamqc -bam output_mitomapped_sorted.bam -outfile results.pdf
```
Inspect the output to understand the coverage and quality of the reads mapped to the mitogenome region. We can now move to Step 5 to generate a consensus mitogenome and annotate.


## Mitogenome recovery using mitogenome reference only
### Step 1 - Download the mitogenome reference and index (bwa) 
From NCBI or other online sources, download the mitochondrial reference genome of your species of interest. In case no reference is available we can use a close sister species to our taxa.

Now we index the reference using bwa

```
bwa index ~/path/to/mito_ref.fasta
```

### Step 2 - Map the trimmed reads to the indexed reference (BWA and samtools)

Now we will map the trimmed reads of our data to the indexed reference we generated in the previous step. To do that we will use bwa.

```
bwa mem -t 20 -M -R "@RG\tID:\tSM:\tLB:IlluminaWGS\tPL:ILLUMINA" \ #here I have used -t flag to specify cores, -R to add read-group to bam, do not forget to add details to the bam-header
~/path/to/indexed/mito_ref.fasta/ \
~/path/to/R1.fq.gz \
~/path/to/R2.fq.gz | \
samtools view -bh - | \
samtools sort -@20 -T tmp -o output_mapped_sorted.bam
```

### Step 3 - Run BAMQC to get an idea of mapping coverage and quality (qualimap)
Using qualimap bamqc we will assess the mapping coverage and quality
```
qualimap bamqc -bam output_mapped_sorted.bam -outfile results.pdf
```
Inspect the output to understand the coverage and quality of the reads mapped to the reference mitogenome. Note that < 98% of our reads have mapped.

### Step 4 - Extract only mapped reads (samtools)
As we just used mitogenome reference, we can proceed to extract only mapped reads. To do this we will use samtools
```
samtools view -b -F 4 -o output_mitomapped_sorted.bam output_mapped_sorted.bam
```

We are now ready to proceed to the Step 5 to generate consensus mitogenome and annotate.

## Step 5 - Using Geneious prime
Import the final bam file into geneious using **File > Import > Files**

Download the same reference mitogenome (only mitogenome reference) that was used in generating the bam file using NCBI database in the geneious software.
![download (6)](https://github.com/user-attachments/assets/1b84ac18-7dfa-45a5-b4ed-f0ad85e9e4c2)
![download (7)](https://github.com/user-attachments/assets/1dc1ea70-b614-4c79-a4ae-d797b100ef9e)


You can now visualize where all the reads map on the reference mitogenome.

![download](https://github.com/user-attachments/assets/c54d43c2-d9e0-498d-9c9c-cf6191eeb01a)


Once you get a good idea, you can generate a consensus sequence. To do this select **Tools > Generate Consensus Seqeunce**

![download (1)](https://github.com/user-attachments/assets/d0af7eb9-1ba3-4a7c-a493-dbb41b5f12eb)
![download (2)](https://github.com/user-attachments/assets/57cc17b3-cc44-4f93-9b7d-a579ad9fed83)



After generating the consensus sequence, transfer the annotations from the reference mitogenome sequence by checking the checkbox under Annotations and Tracks tab.

![download (3)](https://github.com/user-attachments/assets/9b1f43f6-ec2e-4892-9de6-416655b32ded)

Make sure to use a similarity of 40-50% to ensure that annotations are transfered from your reference to the consensus sequence.

In case we are interested in circularizing the mitogenome we can do it using **Sequence > Circularize Sequence**
![download (4)](https://github.com/user-attachments/assets/288b0874-6912-4ca1-a480-dad0002623b8)
![download (5)](https://github.com/user-attachments/assets/19b67558-2d9c-4fad-ae98-7f5a8dc64b8d)

This will give us the final circularized annotated mitogenome which has beene recovered from the whole-genome dataset

