# Genomic Islands in *Wolbachia* Prophages Drive Amplification and Diversification of Cytoplasmic Incompatibility Genes in *Culex pipiens*

This repository accompanies the manuscript investigating the genomic environment of the cytoplasmic incompatibility factor (*cif*) genes in several *w*Pip strains of *Wolbachia*, the endosymbiotic bacterium infecting mosquitoes of the *Culex* complex.

The study aimed to characterize the genomic context and identify key genes flanking the *cif* loci in *w*Pip strains by analyzing Nanopore long reads polished with Illumina short reads. As *de novo* assembly attempts were unsuccessful for circularization of chromosome, we performed a manual analysis of polished Nanopore reads mapped against the *cif* genes.

This GitHub repository provides the main commands, scripts, and key intermediate files used throughout the analyses.

## 1. *De novo* assembly

## 2. Selection of Nanopore polished reads by mapping *cif* genes

## 3. Identification of *cif* variants in the Nanopore polished reads: *cinA-cinB* monomorphy and *cidA-cidB* variants

## 4. Annotation of Nanopore polished reads and assembly of partial prophage contigs

## 5. Analyse of specific genes : identity, phylogeny, recombination

### 5.2. Phylogenies and network analyses

#### Genes forming the *w*Pip *cid* genomic island: phylogenies and network analyses

#### Prophage key genes: phylogenies and network analyses

##### srWO phylogeny

Recombinase sequences used by Bordenstein & Bordenstein (2022) (<https://doi.org/10.1371/journal.pgen.1010227>) were recovered to determine which srWO type were our new wPip prophage sequences. We aligned knwown-srWO type recombinases with wPip Pel-and-JHB recombinases, as well as our Tunis, Slab and Harash prophage recombinases, and included it in a `.faa` file named `srWO_alignment.faa`.
The best ML substitution model and the phylogenetic tree were then created:

```
modeltest-ng -i srWO_alignment.faa -p 12 -T raxml -d aa

raxml-ng --all --msa srWO_alignment.faa --model FLU+G4 --prefix srWO_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

#### wPip phylogeny

We used sequences of five MLST genes from Atyame et al. (2011) (<https://doi.org/10.1093/molbev/msr083>) which we compared to sequences from the wPip Pel-JHB-and-Mol reference genomes, as well as our Harash sequences (the Slab and Tunis sequences were those from Atyame et al.). First, sequences of the pk1, pk2, MutL, GP12, and GP15 genes were independently aligned using `Clustal Omega` (Sievers and Higgins, 2017, <https://doi.org/10.1002/pro.3290>) implemented in `Unipro UGENE` (Okonechnikov et al, 2012, <https://doi.org/10.1093/bioinformatics/bts091>), and concatenated into a .fasta file named `wPip_strain_alignment.fasta`. Next, substitution models were evaluated using `modeltest (v.0.1.7)` (Darriba et al, 2019, <https://doi.org/10.1093/molbev/msz189>) to determine the most appropriate ML substitution model (based on the AICc criterion) for phylogenetic tree construction with `raxml-ng (v.1.1.0)` (Kozlov et al, 2019, <https://doi.org/10.1093/bioinformatics/btz305>):

```
####bash####

modeltest-ng -i wPip_strain_alignment.fasta -p 12 -T raxml -d nt

raxml-ng --all --msa wPip_strain_alignment.fasta --model HKY+I+G4 --prefix wPip_strain_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

Finally, the phylogenetic tree was visualized and modified using figtree (<https://github.com/rambaut/figtree/>) and MEGA7 (<https://megasoftware.net/>)


## 6. Check-quality of newly assembled prophage contigs by mapping raw Nanopore longreads

#### Step 1. Remove the reads of the host *Culex pipiens* 

Prepare the files usefull for mapping:

```
####bash####
# Pool all the sequencing .fastq.gz files together
cat *.fastq.gz > FBF16085_pass_all.fastq.gz
# importe a reference genome of Culex pipiens. Here: GCA_041146695.1, renamed "ref_culpip.fasta"
# importe a reference genome of Wolbachia wPip : Here: NC_010981.1, renamed "wPipPel.fasta"
```

To best prepare for mapping the raw sequencing reads onto a clean host Culex pipiens genome, we create a subset of sequences from `ref_culpip.fasta` that are not also found in Wolbachia wPip. Here is the script used to generate the subset dataset `culpip_masked.fa` using `bbmap (v.38.87)` (<https://bbmap.org/>):

```
####bash####
# Cut wPipPel.fasta into small artificial fragments of 80 bp
shred.sh in=wPipPel.fasta out=wPipPel_shredded.fa length=80 minlength=70 overlap=40

# Map wPipPel_shredded.fa onto the genome ref_culpip.fasta
bbmap.sh ref=ref_culpip.fasta in=wPipPel_shredded.fa outm=culpip_to_wPipPel_mapped.sam minid=0.85 maxindel=2 threads=6

# Mask in the genome ref_culpip.fasta the regions that were aligned with wPipPel
bbmask.sh in=ref_culpip.fasta out=culpip_masked.fa entropy=0.7 sam=culpip_to_wPipPel_mapped.sam threads=10
```

Then, we used `Minimap2 (v.2.24)` (Li, 2018, <https://doi.org/10.1093/bioinformatics/bty191>, <https://github.com/lh3/minimap2>), `samtools (v.1.9)` (<https://github.com/samtools/samtools>) and `anvi'o (v.8)` (Eren et al, 2020, <https://doi.org/10.1038/s41564-020-00834-3>, <https://anvio.org/>) to map the raw sequencing reads onto `culpip_masked.fa` and remove the mapped reads for downstream analyses:

```
####bash####
# Map the reads
minimap2 -ax map-ont -a culpip_masked.fa --sam-hit-only fastq_pass/FBF16085_pass_all.fastq.gz > FBF16085_to_culpip_masked.sam

# Extract the IDs of reads mapped to to_culpip_masked so they can be removed later
samtools view FBF16085_to_culpip_masked.sam | cut -f 1 > FBF16085_IDs_to_remove_culpip.txt

# Remove from the FASTQ all reads whose ID is in the list (here, those mapped to culpip_masked)
iu-remove-ids-from-fastq -i fastq_pass/FBF16085_pass_all.fastq.gz -l FBF16085_IDs_to_remove_culpip.txt -d " "

# Rename the FASTQ containing the reads remaining after removing culpip_masked
mv fastq_pass/FBF16085_pass_all.fastq.gz.survived fastq_pass/FBF16085_no_culpip.fastq.gz
```

Here is a synthesis of the raw sequencing reads' nature of the Harash line as example:

```
# Total number of reads processed: 2,475,160 (Total number of raw sequencing reads)
# Number of reads removed        : 2,132,771 (Total number of reads mapping on "culpip_masked.fa")
# Remaining reads                : 342,389
```


#### Step 2. Select the reads mapped to *Wolbachia* *w*Pip

We mapped the remaining reads longer than 10 kb onto our prophage contigs Tunisp-a–d, Slabp-a–c, and Harashp-a–d. Here is the example for Harash reads:

```
####bash####
# Pool the four harashp_*.fasta contigs into a single fasta file, then reformate it for further analyse
awk 'FNR==1{print ""}1' Harashp_*.fasta > Harashp_ALL.fasta
awk '/^>/ {print; next} {gsub(/[RYSWKMBDHV-]/,"N"); print}' Harashp_ALL.fasta > Harashp_ALL_1.fasta
awk '/^>/{skip=seen[$0]++} !skip' Harashp_ALL_1.fasta > Harashp_ALL.fasta
rm Harashp_ALL_1.fasta

# Subselect reads longer than 10kb
seqkit seq -m 10000 fastq_pass/FBF16085_no_culpip.fastq.gz > FBF16085_no_culpip_10kb.fastq

# Map the subselect reads longer than 10kb to the Harashp_* contigs
minimap2 -ax map-ont -a Harashp_ALL.fasta FBF16085_no_culpip_10kb.fastq --sam-hit-only > FBF16085_no_culpip_10kb_on_Harashp.sam
# .sam to .bam
samtools view -bS FBF16085_no_culpip_10kb_on_Harashp.sam -o FBF16085_no_culpip_10kb_on_Harash-raw.bam

# Filter the BAM file to keep only reads with mapping quality ≥ 2
samtools view -bq 2 FBF16085_no_culpip_10kb_on_Harash-raw.bam > FBF16085_no_culpip_10kb_on_Harash_quality.bam

# Extract read IDs mapped to reference "Harashp_*" from the filtered BAM file. Keep only unique read IDs and save them into a text file
samtools view FBF16085_no_culpip_10kb_on_Harash_quality.bam | awk '$3=="Harashp_a" {print $1}' | sort -u > Harashp_a_IDs.txt
seqkit grep -f Harashp_a_IDs.txt FBF16085_no_culpip_10kb.fastq | seqkit fq2fa > Harashp_a_reads.fasta
grep -o ">" Harashp_a_reads.fasta | wc -l

samtools view FBF16085_no_culpip_10kb_on_Harash_quality.bam | awk '$3=="Harashp_b" {print $1}' | sort -u > Harashp_b_IDs.txt
seqkit grep -f Harashp_b_IDs.txt FBF16085_no_culpip_10kb.fastq | seqkit fq2fa > Harashp_b_reads.fasta
grep -o ">" Harashp_b_reads.fasta | wc -l

samtools view FBF16085_no_culpip_10kb_on_Harash_quality.bam | awk '$3=="Harashp_c" {print $1}' | sort -u > Harashp_c_IDs.txt
seqkit grep -f Harashp_c_IDs.txt FBF16085_no_culpip_10kb.fastq | seqkit fq2fa > Harashp_c_reads.fasta
grep -o ">" Harashp_c_reads.fasta | wc -l

samtools view FBF16085_no_culpip_10kb_on_Harash_quality.bam | awk '$3=="Harashp_d" {print $1}' | sort -u > Harashp_d_IDs.txt
seqkit grep -f Harashp_d_IDs.txt FBF16085_no_culpip_10kb.fastq | seqkit fq2fa > Harashp_d_reads.fasta
grep -o ">" Harashp_d_reads.fasta | wc -l
```

Here are the commands and synthesis of the remaining reads longer than 10 kb matching with the Harashp contigs:

```
####bash####
samtools view FBF16085_no_culpip_10kb_on_Harash_quality.bam | awk '{print $1"\t"$3}' | sort -u | cut -f2 | sort | uniq -c

# Mapped on Harashp_a                                    : 504
# Mapped on Harashp_b                                    : 92
# Mapped on Harashp_c                                    : 475
# Mapped on Harashp_d                                    : 407
# Total reads longer than 10kb mapping on Harashp contigs: 1,478
```

Global read coverage were vizualised on the respective Harashp contigs using `Integrative Genomics Viewer (v.2.19.7)` (Robinson et al. (2011, <https://doi.org/10.1038/nbt.1754>, <https://igv.org/>) by opening `-sorted.bam` files generated as follow:

```
####bash####
anvi-init-bam FBF16085_no_culpip_10kb_on_Harash_quality.bam -o FBF16085_no_culpip_10kb_on_Harash_quality-sorted.bam
```
