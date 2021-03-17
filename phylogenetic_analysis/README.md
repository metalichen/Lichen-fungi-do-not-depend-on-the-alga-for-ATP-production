## Details on the phylogenetic analysis

### ATP9 tree. The full list of commands:
* Align the protein sequences and trim the alignement to remove positions with >90% gaps
```
mafft --genafpair --maxiterate 10000 atp9_renamed.faa  > atp9_renamed_aligned.faa
trimal -in atp9_renamed_aligned.faa  -out atp9_renamed_gpmph_intact.trimal.phylip -gt 0.1 -phylip
```
* Estimate the best substitution model
```
iqtree -s atp9_renamed_gpmph_intact.trimal.phylip -m MF
```
Estimated best substitution model LG+F+G4

* Run ML phylogenetic analysis
```
iqtree -s tp9_renamed_gpmph_intact.trimal.phylip -bb 50000 -pre atp9 -seed 12345 -m LG+F+G4 -nt 12
```
Output of the analysis and the alignement is in /tree

### dN/dS Analysis
Done following [this protocol](https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e?step=4)
* Produce alignments
```
mafft --genafpair --maxiterate 10000 lecanoro_atp9_nuclear_transcript.fna > lecanoro_atp9_nuclear_transcript_aligned.fna
mafft --genafpair --maxiterate 10000 lecanoro_atp9_nuclear.faa > lecanoro_atp9_nuclear_aligned.faa
```
* Get dN/dS ratios
```
~/bin/pal2nal.v14/pal2nal.pl ../protein_seqs/lecanoro_atp9_nuclear_aligned.faa ../transcript_seqs/lecanoro_atp9_nuclear_transcript.fna -output paml -nogap > lecanoro_atp9_nuclear.pal2nal
```
