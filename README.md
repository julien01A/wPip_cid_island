# WO_cid_island

### wPip phylogeny
We used sequences of five MLST genes from Atyame et al. (2011) (<doi:10.1093/molbev/msr083>) which we compared to sequences from the wPip Pel-JHB-and-Molestus reference genomes, as well as our Harash sequences (the Slab and Tunis sequences were those from Atyame et al.). First, sequences of the pk1, pk2, MutL, GP12, and GP15 genes were independently aligned using Clustal Omega (<doi:10.1002/pro.3290>) implemented in Unipro UGENE (<doi:10.1093/bioinformatics/bts091>), and concatenated into a .fasta file named ```wPip_strain_alignment.fasta```. 
Next, substitution models were evaluated using modeltest v0.1.7 (<doi:10.1093/molbev/msz189>) to determine the most appropriate ML substitution model (based on the AICc criterion) for phylogenetic tree construction with raxml-ng v1.1.0 (<doi:10.1093/bioinformatics/btz305>):
```
modeltest-ng -i wPip_strain_alignment.fasta -p 12 -T raxml -d nt

raxml-ng --all --msa wPip_strain_alignment.fasta --model HKY+I+G4 --prefix wPip_strain_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

Finally, the phylogenetic tree was visualized and modified using figtree (<https://github.com/rambaut/figtree/>) and MEGA11 (<https://megasoftware.net/>)

### srWO phylogeny
Recombinase sequences used by Bordenstein & Bordenstein (2022) (<doi:10.1371/journal.pgen.1010227>) were recovered to determine which srWO type were our new wPip prophage sequences. We aligned knwown-srWO type recombinases with wPip Pel-and-JHB recombinases, as well as our Harash and Tunis prophage recombinases, and included it in a .faa file named ```srWO_alignment.faa```.
The best ML substitution model and the phylogenetic tree were then created:
```
modeltest-ng -i srWO_alignment.faa -p 12 -T raxml -d aa

raxml-ng --all --msa srWO_alignment.faa --model FLU+G4 --prefix srWO_alignment-raxmlng --seed 5 --threads 4 --bs-trees 1000

```

