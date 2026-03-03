# Genomic Islands in *Wolbachia* Prophages Drive Amplification and Diversification of Cytoplasmic Incompatibility Genes in *Culex pipiens*

This repository accompanies the manuscript investigating the genomic environment of the cytoplasmic incompatibility factor (*cif*) genes in several *w*Pip strains of *Wolbachia*, the endosymbiotic bacterium infecting mosquitoes of the *Culex* complex.

The study aimed to characterize the genomic context and identify key genes flanking the *cif* loci in *w*Pip strains by analyzing Nanopore long reads polished with Illumina short reads. As *de novo* assembly attempts were unsuccessful for circularization of chromosome, we performed a manual analysis of polished Nanopore reads mapped against the *cif* genes.

This GitHub repository provides the main commands, scripts, and key intermediate files used throughout the analyses.

## 1. Polishing and *De novo* assembly

Here is an example of the workflow used for Harash.

#### Step 1. Remove the reads of the host *Culex pipiens* 

Prepare the files usefull for mapping:

```
####bash####
# Pool all the sequencing .fastq.gz files together
cat *.fastq.gz > FBF16085_pass_all.fastq.gz
# importe a reference genome of Culex pipiens. Here: GCA_041146695.1, renamed "ref_culpip.fasta"
# importe a reference genome of Wolbachia wPip : Here: AM999887.1, renamed "wPipPel.fasta"
```

To best prepare for mapping the raw sequencing reads onto a clean host *Culex pipiens* genome, we create a subset of sequences from `ref_culpip.fasta` that are not also found in Wolbachia wPip. Here is the script used to generate the subset dataset `culpip_masked.fa` using `bbmap (v.38.87)` (<https://bbmap.org/>):

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

#### Step 2. Polishing with Illumina short reads

Finir avec un fichier "FBF16085_no_culpip_polished_Illumina.fastq.gz"

#### Step 3. Assembly attempts



## 2. Selection of Nanopore polished reads by mapping *cif* genes

We used `Minimap2 (v.2.24)` and `samtools (v.1.9)` to map the Nanopore polished reads longer than 10 kb onto *cidA* (*cifA* type I), *cidB* (*cifB* type I), *cinA* (*cifA* type IV) and *cinB* (*cifB* type IV) extracted from the *w*Pip Pel refence genome (accession number: `AM999887.1`), so called: `WP0282`, `WP0283`, `WP0294`, `WP0295` and the ramining cidB fragment `WP1291`. `Minimap2 (v.2.24)` and `samtools (v.1.9)` (<https://github.com/samtools/samtools>) were used. Here is the script for Harash as example:

```
####bash####
# Pool the five cif genes *.fasta into a single fasta file, then reformate it for further analyse
awk 'FNR==1{print ""}1' cidA.fa cidB.fa cinA.fa cinB.fa WP1291.fa > cif_ALL.fasta

# Subselect reads longer than 10kb
seqkit seq -m 10000 fastq_pass/FBF16085_no_culpip_polished_Illumina.fastq.gz > FBF16085_no_culpip_polished_Illumina_10kb.fastq

# Map the subselect reads longer than 10kb with cif genes
minimap2 -ax map-ont -a cif_ALL.fasta FBF16085_no_culpip_polished_Illumina_10kb.fastq --sam-hit-only > FBF16085_no_culpip_polished_Illumina_10kb_with_cif.sam
# .sam to .bam
samtools view -bS FBF16085_no_culpip_polished_Illumina_10kb_with_cif.sam -o FBF16085_no_culpip_polished_Illumina_10kb_with_cif.bam

# Filter the BAM file to keep only reads with mapping quality ≥ 2
samtools view -bq 2 FBF16085_no_culpip_polished_Illumina_10kb_with_cif.bam > FBF16085_no_culpip_polished_Illumina_10kb_with_cif_quality.bam

# Extract read IDs mapped to one at least one cif gene from the filtered BAM file. Keep only unique read IDs and save them into a text file
samtools view -F 4 FBF16085_no_culpip_polished_Illumina_10kb_with_cif_quality.bam | cut -f1 | sort -u > FBF16085_cif_IDs.txt
seqkit grep -f FBF16085_cif_IDs.txt FBF16085_no_culpip_polished_Illumina_10kb.fastq | seqkit fq2fa > FBF16085_cif_reads.fasta
```

## 3. Annotation of Nanopore polished reads, identification of *cid* variants and assembly of partial prophage contigs

The longest Nanopore polished reads containing *cif* genes were annotated using manual BLASTn (megablast program) against the *w*Pip Pel reference genome (accession number: `AM999887.1`). This genome was already annotated with the `NCBI Prokaryotic Genome Annotation Pipeline v6.9`, so the correspondign matches could be attrubuted to predictive genes. A match was considered positive if the nucleotide identity with *w*Pip Pel was above 75%. For genes with undetermined functions in *w*Pip Pel, we renseigne the new predicted functions by Bordenstein and Bordenstein of the corresponding genes (<https://doi.org/10.1371/journal.pgen.1010227>). For transposases, we conducted BLAST searches against `ISfinder` to classify them into homolog groups (Siguier et al. 2006, <https://doi.org/10.1093/nar/gkj014>, <https://isfinder.biotoul.fr/blast.php>) and use `HHpred` (Zimmermann et al. 2018, <https://doi.org/10.1016/j.jmb.2017.12.007>, <https://toolkit.tuebingen.mpg.de/tools/hhpred>) with the `SCOPe70 v.2.08`, `Pfam-A v37.0`, `SMART v6.0`, and `COG/KOG v1.0 databases` to characterized the PDDEXK2 tranposase flanking *cidA* genes. For the group II intron retrotransposons (maturases), we used Zbase for their identification (Candales et al. 2012, <https://doi.org/10.1093/nar/gkr1043>, <http://webapps2.ucalgary.ca/~groupii/cgi-bin/main/blastusr.php>).

Nanopore polished reads containing *cid* genes were individually screened against databases of *cidA* and *cidB* variants previously identified and extensively characterized in more than 40 *w*Pip strains collected worldwide (Namias et al., 2023, <https://doi.org/10.1016/j.csbj.2023.07.012>, Namias et al., 2025 <https://doi.org/10.1093/molbev/msaf200>). These databases are available in the corresponding GitHub repository under the names `CIDA_FULL_database.fasta` and `CIDB_FULL_database.fasta`. For this analysis, each Nanopore polished read was aligned to the *cid* databases using `MAFFT` implemented on `Geneious` (<https://www.geneious.com/>) and manually curated to precisely identify *cid* variants. 
For Nanopore polished reads containing *cin* genes, alignments were performed against (*cinA*) and `WP0295.fasta` (*cinB*) to confirm the monomorphy of this gene pair in the *w*Pip strains Tunis, JHB, Slab, Harash, and Mol.

Nanopore polished reads sharing identical *cidA-cidB* variants and exhibiting the same annotated gene composition were aligned together using `MAFFT` (treshold: 99% nucleotide identity) implemented on `UGENE v52.0` (Okonechnikov et al. (2012), <https://doi.org/10.1093/bioinformatics/bts091>, <https://ugene.net/>) to create contigs. The consensus sequences of the final contigs are available on NBCI under accession numbers: `PX571960-PX571970`.

## 4. Phylogenies, networks and whole identities of target genes

### 4.1. Identities

To estimate similarity percentages for selected target genes (for instance, comparing structural prophage genes with genes from the *cid* genomic island), each gene was aligned independently using MAFFT implemented on UGENE. By comparing the sequences from our wPip prophage contigs and those of the reference genomes Pel and JHB, we then calculated the percentage of similarity as the number of variable positions over the total number of aligned positions fo each gene:

```
####R####
library(ape)
alignment <- read.dna("Terminase_AL.fa", format = "fasta") # where Terminase_Al.fa is an alignement file of the target gene
aln_matrix <- as.matrix(alignment)
is_polymorphic <- function(column) {
  length(unique(column)) > 1
}
polymorphic_sites <- apply(aln_matrix, 2, is_polymorphic)
n_poly <- sum(polymorphic_sites)
alignment_length <- ncol(aln_matrix)
percent_poly <- (n_poly / alignment_length) * 100
cat("Total sites:", alignment_length, "\n")
cat("Polymorphic sites:", n_poly, "\n")
```

### 4.2. Phylogenies

All the phylogenetic trees generated in this study were calculated using the `Maximum Likelihood` method. For each gene, operon, or concatenated sequence of interest targeted for phylogenetic analysis, an alignment file was required. All the alignment files were generated by aligning the desired sequences using `Clustal Omega` (Sievers and Higgins (2017),<https://doi.org/10.1002/pro.3290>, <https://github.com/GSLBiotech/clustal-omega>) or `MAFFT` (Katoh and Standley (2013), <https://doi.org/10.1093/molbev/mst010>, <https://github.com/GSLBiotech/mafft>) implemented on `UGENE v52.0`, and then the positions with gaps '-' were removed.

Then, for each ALIGNMENT.fa files, substitution models were evaluated using `modeltest v0.1.7` (Darriba et al. (2019), <https://doi.org/10.1093/molbev/msz189>, <https://github.com/ddarriba/modeltest>) to determine the most appropriate ML substitution model (based on the AICc criterion), followed by phylogenetic tree construction with `raxml-ng v1.1.0` (Kozlov et al. (2019) <https://doi.org/10.1093/bioinformatics/btz305>, <https://github.com/amkozlov/raxml-ng>):

```
####bash####
# predict the best ML model
modeltest-ng -i ALIGNMENT.fa -p 12 -T raxml -d aa # for example: FLU+G4

# calculate the ML tree
raxml-ng --all --msa ALIGNMENT.fa --model FLU+G4 --prefix Your_Tree-raxmlng --seed 5 --threads 4 --bs-trees 1000
raxml-ng --support --tree Your_Tree-raxmlng.raxml.bestTree --bs-trees 1000 --prefix Your_Tree-boot --threads 2
```

Finally, the phylogenetic tree was visualized and adapted using `figtree v.1.4.4` (<https://github.com/rambaut/figtree/>) and `MEGA11` (<https://megasoftware.net/>)

For the wPip phylogeny, we used sequences of five MLST genes (pk1, pk2, MutL, GP12, and GP15) from Atyame et al. (2011) (<https://doi.org/10.1093/molbev/msr083>) which we compared to sequences from the wPip Pel,JHB, and Mol reference genomes, as well as our Harash sequences (the Slab and Tunis sequences were those from Atyame et al.).

For phage gene phylogenies (srWO, Minor capsid, GpA, Tail tape measure, Tail fiber, GpV, RepA and methylase), we compared sequences found on wPip prophages (Pel, JHB, Mol, Tunis, Slab and Tunis) with few key examples of Wolbachia prophages characterized in Bordenstein & Bordenstein (2022) (<https://doi.org/10.1371/journal.pgen.1010227>).

### 4.3. Network analyses

To vizualise recombination events for some specific genes, we used the `SplitsTree App v.4.19.2` (Huson and Bryant (2024), <https://doi.org/10.1038/s41592-024-02406-3>) using the Uncorrected_P distance method and the NeighborNet network method.
For each alignement, a formal test of recombination (PHI test) were performed using SplitsTree.

## 5. Check-quality of newly assembled prophage contigs by mapping Nanopore raw reads

### 5.1. Global coverage of the assembled prophage contigs with Nanopore raw reads

To check the quality of our newly assembled prophage contigs, we first mapped the raw Nanopore reads longer than 10 kb onto our prophage contigs Tunisp-a–d, Slabp-a–c, and Harashp-a–d to visualize the global read coverages. Here is the example for Harash reads with `Minimap2 (v.2.24)` and `samtools (v.1.9)`:

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

### 5.2. Nanopore raw read analysis of genomic island–forming genes to assess potential Illumina polishing artifacts

`Minimap2 (v.2.24)` and `Bedtools (v.31.1)` (Quinland (2015), <https://doi.org/10.1002/0471250953.bi1112s47>, <https://github.com/arq5x/bedtools2>) were used to retrieve the sequences of three genes (Transposase PDDEXK2, Intron group II and *rnhA*) onto both the Nanopore polished reads and the Nanopore raw reads to assess potential Illumina polishing artifacts:

```
####bash####
# preparation: put in a folder the Nanopore raw reads which map the prophage contigs Tunisp*, Slabp* and Harashp* (these raw reads were called Harash_reads.fasta, Slab_reads.fasta and Tunis_reads.fasta), the Nanopore polished reads used to create the contigs (Harash_polished.fasta, Slab_polished.fasta and Tunis_polished.fasta) and the three target gene queries (querie_TrpPDDEX.fa, querie_intII.fa and querie_rnhA.fa)

for query in querie_*.fa; do
    query_name=$(basename "$query" .fa)
    mkdir -p results_2_${query_name}
    echo "=== Processing $query ==="
    for fasta in *.fasta; do
        sample="${fasta%.*}"
        echo "Mapping $fasta vs $query"
        minimap2 -x map-ont -k13 -w5 -A1 -B2 -O2,24 -E2,1 -t 8 "$query" "$fasta" > tmp.paf
        if [ -s tmp.paf ]; then
            awk '$11 > 300 {print $1"\t"$3"\t"$4}' tmp.paf > coords.bed
            bedtools getfasta -fi "$fasta" -bed coords.bed -fo results_2_${query_name}/${sample}_mapped_regions.fasta
            echo "Regions extracted"
        else
            echo "No mapping"
        fi  
        rm -f tmp.paf coords.bed
    done
done
```

Raw read sequences of three target genes (Transposase PDDEXK2, Intron group II and *rnhA*) extracted from Nanopore raw reads and Nanopore polished reads where aligned using `MAFFT`, reverse-complemented if necessary, and visulized on `UGENE`. A synthesis file of theses alignment is available on this GitHub repository, called `Polishing_error_verification.pdf`. Point mutations are detected throughout the reads and in the majority of them. However, none of these mutations are shared between reads; all reads exhibit a common structural pattern, which is more or less clear depending on sequencing quality. These results indicate that only mutations attributable to sequencing errors were identified, demonstrating that Illumina polishing did not introduce artificial artifacts but rather corrected isolated single-nucleotide discrepancies.
