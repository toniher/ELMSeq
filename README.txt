[What is ELMSeq?]
ELM-seq (expression level monitoring by DNA methylation) uses DamID (Escherichia coli DNA adenine methylase as a reporter coupled with methylation-sensitive restriction enzyme digestion and high-throughput sequencing) to enable in vivo quantitative analyses of upstream regulatory sequences.

[ELMSeq.py]
This program is taking ELMSeq sequencing files and generating DAMRatios to correponding sequences.

[How to use]
We assumed that the sequencing format is the same as we described in Nature Communication 2017 Paper.

# 1. Need to set proper parameter in ELMSeq.py file 

# User needs to setup output file folder path
OUTPUT_FOLDER = "./output"

# User needs to setup corresponding fastq file path
# Example sequencing data for transcription study (promoter)
UNCUT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_CGTGAT.fastq.gz"            # No Enzyme treated sample
DPNI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_ACATCG.fastq.gz"             # DpnI treated sample
MBOI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_GCCTAA.fastq.gz"             # MboI treated sample

# User needs to specify which experiment it is among followings: ( promoter / utr_with_strong_promoter / utr_with_weak_promoter )
EXPERIMENT_TYPE = "promoter" # "promoter", "utr_with_strong_promoter", "utr_with_weak_promoter" 

# 2. After properly setting the parameters, User can run the program as following command

# python ELMSeq.py

# 3. This is the typical command line output

[ Starting ELMSeq Analysis ] It might take a long time (more than hours) depending on sequencing size
[ Read No Enzyme Treated Sequence File ] ./data/dam_screen_lib_2_12466_CGTGAT.fastq.gz
[ Read DpnI Treated Sequence File ] ./data/dam_screen_lib_2_12466_ACATCG.fastq.gz
[ Read MboI Treated Sequence File ] ./data/dam_screen_lib_2_12466_GCCTAA.fastq.gz
[ Calculating DAMRatios ]
[ Finishing ELMSeq Analysis ]

# 4. In the output folder, you can see DAMRatios for all the sequences


[Output File Format]
# 1. The first three line is commented. It gives number of processed reads.
# 2. The following results show the DAMRatio and corresponding sequence. Results were sorted by DAMRatio. 

# uncut_all_cnt	uncut_match_cnt	DpnI_all_cnt	DpnI_match_cnt	MboI_all_cnt	MboI_match_cnt
# 2000000	363302	2000000	353958	2000000	350541
# DAMRatio	motif	uncut_cnt	DpnI_cnt	MboI_cnt	uncut_CPM	DpnI_CPM	MboI_CPM
1.575207e+02	TTGTTACGGCGTCGTTCTAAAACATTATAATTATCGAT	17	0	155	49.545557	2.825194	445.026402
1.322770e+02	TCGCTTACGTCAAACAGATCTTTAGTATAATTGACCCT	14	0	130	41.287964	2.825194	373.708068
...
