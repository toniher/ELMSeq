[What is ELMSeq?]


This program is taking ELMSeq sequencing files and generating DAMRatios to correponding sequences.

We assumed that the sequencing format is the same as we described in Nature Communication 2017 Paper.

# After setting proper parameter in ELMSeq.py file 
# User can run the program as following command

# python ELMSeq.py

# Typical command line output

[ Starting ELMSeq Analysis ] It might take a long time (more than hours) depending on sequencing size
[ Read No Enzyme Treated Sequence File ] ./data/dam_screen_lib_2_12466_CGTGAT.fastq.gz
[ Read DpnI Treated Sequence File ] ./data/dam_screen_lib_2_12466_ACATCG.fastq.gz
[ Read MboI Treated Sequence File ] ./data/dam_screen_lib_2_12466_GCCTAA.fastq.gz
[ Calculating DAMRatios ]
[ Finishing ELMSeq Analysis ]

# It will gives DAMRatios for all the sequences


[Output Files]

