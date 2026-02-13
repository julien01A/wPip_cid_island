# wPip_cid_island

### Mapping of Nanopore longreads on contigs

#### Step 1. Remove the reads of the host *Culex pipiens*

Prepare the files usefull for mapping:

```
####bash####
# Pool all the sequencing .fastq.gz files together
cat *.fastq.gz > FBF16085_pass_all.fastq.gz
# importe a reference genome of Culex pipiens. Here: GCA_041146695.1, renamed "ref_culpip.fasta"
#  importe a reference genome of Wolbachia wPip : Here: NC_010981.1, renamed "wPipPel.fasta"
```

To best prepare for mapping the raw sequencing reads onto a clean host Culex pipiens genome, we create a subset of sequences from `ref_culpip.fasta` that are not also found in Wolbachia wPip. Here is the script used to generate the subset dataset `culpip_masked.fa`:

```
####bash####
# Cut wPipPel.fasta into small artificial fragments of 80 bp
shred.sh in=wPipPel.fasta out=wPipPel_shredded.fa length=80 minlength=70 overlap=40

# Map wPipPel_shredded.fa onto the genome ref_culpip.fasta
bbmap.sh ref=ref_culpip.fasta in=wPipPel_shredded.fa outm=culpip_to_wPipPel_mapped.sam minid=0.85 maxindel=2 threads=6

# Mask in the genome ref_culpip.fasta the regions that were aligned with wPipPel
bbmask.sh in=ref_culpip.fasta out=culpip_masked.fa entropy=0.7 sam=culpip_to_wPipPel_mapped.sam threads=10
```

Then, we used `Minimap2 (v.2.24)` (Li, 2018, <https://doi.org/10.1093/bioinformatics/bty191>, <https://github.com/lh3/minimap2>) and samtools (v.1.9) (<https://github.com/samtools/samtools>) to map the raw sequencing reads onto `culpip_masked.fa` and remove the mapped reads for downstream analyses:

```
####bash####
# Map the reads
minimap2 -ax map-ont -a culpip_masked.fa --sam-hit-only fastq_pass/FBF16085_pass_all.fastq.gz > FBF16085_to_culpip_masked.sam

# Extract the IDs of reads mapped to to_culpip_masked so they can be removed later
samtools view FBF16085_to_culpip_masked.sam | cut -f 1 > FBF16085_IDs_to_remove_culpip.txt

# Remove from the FASTQ all reads whose ID is in the list (here, those mapped to culpip)
iu-remove-ids-from-fastq -i fastq_pass/FBF16085_pass_all.fastq.gz -l FBF16085_IDs_to_remove_culpip.txt -d " "

# Rename the FASTQ containing the reads remaining after removing culpip
mv fastq_pass/FBF16085_pass_all.fastq.gz.survived fastq_pass/FBF16085_no_culpip.fastq.gz
```

Here is a synthesis of the raw sequencing reads' nature of the Harash line as example:

```
# Total number of reads processed: XXXX (Total number of raw sequencing reads)
# Number of reads removed        : XXXX (Total number of reads mapping on "culpip_masked.fa")
# Remaining reads: XXXX
```


#### Step 2. Select the reads mapped to Wolbachia wPip

We mapped the remaining reads longer than 10 kb onto our proephage contigs Tunisp-a–d, Slabp-a–c, and Harashp-a–d.

```
####bash####
```

### wPip phylogeny

We used sequences of five MLST genes from Atyame et al. (2011) (<https://doi.org/10.1093/molbev/msr083>) which we compared to sequences from the wPip Pel-JHB-and-Molestus reference genomes, as well as our Harash sequences (the Slab and Tunis sequences were those from Atyame et al.). First, sequences of the pk1, pk2, MutL, GP12, and GP15 genes were independently aligned using `Clustal Omega` (Sievers and Higgins, 2017, <https://doi.org/10.1002/pro.3290>) implemented in `Unipro UGENE` (Okonechnikov et al, 2012, <https://doi.org/10.1093/bioinformatics/bts091>), and concatenated into a .fasta file named `wPip_strain_alignment.fasta`. Next, substitution models were evaluated using `modeltest (v.0.1.7)` (Darriba et al, 2019, <https://doi.org/10.1093/molbev/msz189>) to determine the most appropriate ML substitution model (based on the AICc criterion) for phylogenetic tree construction with `raxml-ng (v.1.1.0)` (Kozlov et al, 2019, <https://doi.org/10.1093/bioinformatics/btz305>):

```
modeltest-ng -i wPip_strain_alignment.fasta -p 12 -T raxml -d nt

raxml-ng --all --msa wPip_strain_alignment.fasta --model HKY+I+G4 --prefix wPip_strain_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

Finally, the phylogenetic tree was visualized and modified using figtree (<https://github.com/rambaut/figtree/>) and MEGA7 (<https://megasoftware.net/>)

### srWO phylogeny

Recombinase sequences used by Bordenstein & Bordenstein (2022) (<https://doi.org/10.1371/journal.pgen.1010227>) were recovered to determine which srWO type were our new wPip prophage sequences. We aligned knwown-srWO type recombinases with wPip Pel-and-JHB recombinases, as well as our Tunis, Slab and Harash prophage recombinases, and included it in a `.faa` file named ```srWO_alignment.faa```.
The best ML substitution model and the phylogenetic tree were then created:

```
modeltest-ng -i srWO_alignment.faa -p 12 -T raxml -d aa

raxml-ng --all --msa srWO_alignment.faa --model FLU+G4 --prefix srWO_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

