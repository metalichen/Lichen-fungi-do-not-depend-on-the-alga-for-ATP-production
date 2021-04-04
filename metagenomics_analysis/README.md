# Detailed description of the metagenomics analysis:

1. scripts/Snakefile contatins the full pipeline from raw reads to binned metagenome. We used metaWRAP to filter from human contamination, assemble them and bin the metagenomic assembly

2. scripts/run_busco.sh contains a wrapper script we used to run BUSCO4 on each bin in each assembly. In each metagenome we identified a bin containing the ascomycete MAG:
	* TS1974/short_summary.specific.ascomycota_odb10.bin.22.txt C:96.8%[S:96.4%,D:0.4%],F:0.2%,M:3.0%,n:1706 
	* GTX0158/short_summary.specific.ascomycota_odb10.bin.70.txt C:95.4%[S:95.3%,D:0.1%],F:0.3%,M:4.3%,n:1706  
	* GTX0161/short_summary.specific.ascomycota_odb10.bin.25.txt C:75.5%[S:75.0%,D:0.5%],F:0.4%,M:24.1%,n:1706  
	* GTX0163/short_summary.specific.ascomycota_odb10.bin.55.txt C:96.9%[S:96.6%,D:0.3%],F:0.2%,M:2.9%,n:1706 

3. In case of three metagenomes (TS1974, GTX0161, GTX0163), the ascomycete MAG was split between multiple bins. 
To get the full MAG we used gc-cov plots (../gc_cov_plots/gc_cov_plots.R):
	* We located bin identified as the ascomycete MAG by BUSCO
	* In case of TS1974, GTX0161, and GTX0163, this bin was part of a linear-shaped cloud
	* We merged the bin with other bins constituting the cloud
	* We ran BUSCO4 on the merged bins to confirm that the merging increased completeness while maintaining low contamination rate:
		- TS1974/bin.8_39_22.fa C:96.9%[S:96.5%,D:0.4%],F:0.2%,M:2.9%,n:1706 
		- GTX0161/bins.25_18_31.fa C:97.2%[S:96.4%,D:0.8%],F:0.2%,M:2.6%,n:1706  
		- GTX0163/bin.55_67_82_19.fa C:96.9%[S:96.7%,D:0.2%],F:0.2%,M:2.9%,n:1706
		
4. scripts/funannotate_annotation.sh contains a wrapper file we used to annotate the ascomycete MAGs using the funannotate pipeline

5. Folder output_MAG_annotations contains the final genome annotations for the four ascomycete MAGs
