########################################################################################################################################

                                                ####### MEME SUITE ########

########################################################################################################################################

Online tool
https://meme-suite.org/meme/tools/meme-chip

BiocManager::install("memes")
library(memes)
https://snystrom.github.io/memes-manual/
https://github.com/snystrom/memes  

#### Downloading databases
  https://meme-suite.org/meme/db/motifs
JASPAR/JASPAR2024_CORE_vertebrates_non-redundant.meme
  
# human genome fasta (required for meme-chip)
wget https://meme-suite.org/meme/meme-software/Databases/sequences/fasta_human.tgz



#### list of FUNCTIONS in R 
Function      Name 	Use 	                    Sequence-Input 	Motif Input 	Output
runStreme() 	Motif Discovery (short motifs) 	Yes 	          No 	          universalmotif_df
runDreme()   	Motif Discovery (short motifs) 	Yes 	          No 	          universalmotif_df
runAme() 	    Motif Enrichment 	              Yes 	          Yes 	        data.frame (optional: sequences column)
runFimo()   	Motif Scanning 	                Yes 	          Yes 	        GRanges of motif positions
runTomTom() 	Motif Comparison              	No 	            Yes 	        universalmotif_df w/ best_match_motif and tomtom columns*
runMeme()   	Motif Discovery (long motifs) 	Yes     	      No 	          universalmotif_df

meme-chip #(in terminal only)



# Difference
  

  | Tool          | Type                  | Goal                                     | Input                              | Output                           |
  | --------------| --------------------- | ---------------------------------------- | ---------------------------------- | -------------------------------- |
  | meme-chip     | Combo tool            | Discover & annotate motifs               | FASTA sequences                    | Known + novel motifs             |
  
  | runMeme()     | De novo discovery     | Discover longer, de novo motifs          | FASTA sequences                    | PWM motifs (MEME format)         |
  | runDreme()    | De novo discovery     | Find short enriched motifs               | Positives vs. background sequences | Short novel motifs               |
  | runAme()      | Enrichment test       | Test known motifs for enrichment         | Known motif DB + sequences         | Ranked list of enriched motifs   |
  | runTomTom()   | Motif comparison      | Match discovered motifs to known         | Query motifs + target DB           | Similar known motifs with scores |
  | runCentriMo() | Positional enrichment | Detect central motif enrichment in peaks | Motifs + peak-centered sequences   | Motifs with central enrichment   |
  
##### 
  

##### Running meme-chip in terminal

# Extracting the FASTA files from the bed defined regions
bedtools getfasta -fi hg38.fa -bed upreg_diffbind.bed -fo upreg_diffbind.fa # conversion vers fasta pour le meme suite 

-fi 	# Input FASTA file (your reference genome) — here, hg38 (human)
-bed 	# Input BED file with genomic coordinates (regions of interest, e.g., peaks)
-fo.  # Output FASTA file to write extracted sequences



#### Running meme-chip
#loading as variable
JASPAR2024=/Users/administrateur/Desktop/Bio_info/reference_genomes/MEME/JASPAR2024_CORE_vertebrates_non-redundant.meme
hg38=/Users/administrateur/Desktop/Bio_info/reference_genomes/ENSEMBL/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
mm10=/Users/administrateur/Desktop/Bio_info/reference_genomes/UCSC/mm10.fa.gz


### Getting the FASTA seq in the peak intervals
cd /Users/administrateur/Desktop/Bio_info/ChiPseq/ChIPSeq_Yap_ActiveMotiv/data/working_data/meme/
  
  ls *.bed | while read file

bedtools getfasta -fi $mm10 -bed Up-regulated_diffbind.bed -fo Up-regulated_diffbind.fa # conversion vers fasta pour le meme suite



### # Actually Running meme
meme-chip -oc meme_chip_output -db $JASPAR2024 -meme-maxw 15 upreg_diffbind.fa



| Flag                        | Description                                                                     |
  | --------------------------- | ------------------------------------------------------------------------------- |
  | -oc <dir>                | Output directory (required)                                                     |
  | -db <motif.meme>         | Motif database to compare discovered motifs (e.g., JASPAR)                      |
  | -meme-mod <mode>         | Motif distribution model for MEME: oops, zoops, or anr(default: zoops) |
  | -meme-minw <int>         | Minimum motif width (default: 6)                                                |
  | -meme-maxw <int>         | Maximum motif width (default: 30)                                               |
  | -meme-nmotifs <int>      | Number of motifs to find with MEME (default: 3)                                 |
  | -dreme/ -nodreme         | Enable or disable DREME (short motif discovery)                                 |
  | -centrimo/ -nocentrimo   | Enable or disable CentriMo (positional enrichment)                              |
  | -neg <background.fa>     | Background sequences for DREME/CentriMo (optional)                              |
  | -noecho                  | Do not copy input FASTA to output directory                                     |
  | -verbosity <level>       | Control verbosity (1–5, default 2)                                              |
  