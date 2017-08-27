'''
 *
 *  ELMSeq CommandLine v1.0
 *
 *  Created by Jae-Seong Yang on 04/05/17.
 *  Copyright 2017 CRG. All rights reserved.
 *
 *   
'''

#======================================================================================
# This program calculate DAMRatio from the ELMSeq data
# If it is necessary, user should change corresponding file paths or prefix sequence
#======================================================================================

# User needs to setup output file folder path
OUTPUT_FOLDER = "./output"

# User needs to setup corresponding fastq file path
# Example sequencing data for transcription study (promoter)
# For EXPERIMENT_TYPE = "promoter"  
UNCUT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_CGTGAT.fastq.gz"
DPNI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_ACATCG.fastq.gz"
MBOI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_GCCTAA.fastq.gz"

# Example sequencing data for translation study (with strong promoter)
# For EXPERIMENT_TYPE = "utr_with_strong_promoter"  
#UNCUT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_TGGTCA.fastq.gz"
#DPNI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_CACTGT.fastq.gz"
#MBOI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_ATTGGC.fastq.gz"

# Example sequencing data for translation study (with weak promoter)
# For EXPERIMENT_TYPE = "utr_with_weak_promoter"  
#UNCUT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_GATCTG.fastq.gz"
#DPNI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_TCAAGT.fastq.gz"
#MBOI_FASTQ_PATH = "./data/dam_screen_lib_2_12466_CTGATC.fastq.gz"

# User needs to specify which experiment it is among followings: ( promoter / utr_with_strong_promoter / utr_with_weak_promoter )
EXPERIMENT_TYPE = "promoter" # "promoter", "utr_with_strong_promoter", "utr_with_weak_promoter" 
READ_COUNT_CUT_OFF = 100
ALPHA = 1.0 # Used when calculating CPM

# If it is necessary, user needs to change prefix sequence filter
PROMOTER_RANDOM_SEQ_LEN = 38
PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)GACCGGAACTTCTATGATCGAGATCGAGATCGAGATCGCGGCCGCAAC\w*"
USE_PROMOTER_TATTAT_FILTER = True  # True or False

UTR_RANDOM_SEQ_LEN = 25
UTR_STRONG_PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)AGTTTATATTATAACACTTTAACCTATGGC\w+" 
UTR_WEAK_PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)TGCAATTATTCTAACAAACCCCAAACTTATTTCAA\w+"

WRITE_ALL_DAMRATIOS = False  # True or False

#======================================================================================

import os
import sh
import subprocess
import time
import re
import psutil


def read_summary_file( input_filepath="./output/exact_match/summary/Translation_Weak_Promoter_llmp200.result.txt" ):
    '''
    # uncut_all_cnt	uncut_match_cnt	DpnI_all_cnt	DpnI_match_cnt	MboI_all_cnt	MboI_match_cnt
    # 86070864	14926650	95232908	16485908	81139744	13929837
    # DAMRatio	motif	uncut_cnt	DpnI_cnt	MboI_cnt	uncut_CPM	DpnI_CPM	MboI_CPM
    5.562432e+02	AGTAGATTTTTCGATTACTCAATTTTATAATCATTTAA	51	0	469	3.483702	0.060658	33.740524
    4.852334e+02	TCATGCAACAGATAAAGGCTGTTGTTATAATTAACAAT	42	0	409	2.880754	0.060658	29.433223
    4.497285e+02	GATGAAAGAAAAATTGGTTTTATAATATAATATTTTTA	42	0	379	2.880754	0.060658	27.279573
    '''
    data = []
    f = open( input_filepath )
    title = f.next()
    total_cnt = f.next()
    [ uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt ] = total_cnt[2:].split()

    cnt_infos = [ int(uncut_all_cnt), int(uncut_match_cnt), int(DpnI_all_cnt), int(DpnI_match_cnt), int(MboI_all_cnt), int(MboI_match_cnt) ]

    for line in f.xreadlines():
        if line[0] == "#": continue
        [ DAMRatio, motif, uncut_cnt, DpnI_cnt, MboI_cnt, uncut_CPM, DpnI_CPM, MboI_CPM ] = line[:-1].split( "\t" )
        if "N" in motif: continue
        if motif[25:25+6] != "TATAAT": continue
        data.append( [ float(DAMRatio), motif, int(uncut_cnt), int(DpnI_cnt), int(MboI_cnt), float(uncut_CPM), float(DpnI_CPM), float(MboI_CPM) ] )
    f.close()

    print len( data )
    return title, cnt_infos, data


def MonitorMemory():
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss/1024.0/1024.0)


def getCommandOutput(command, exedir = None):
    if exedir != None:
        os.chdir( exedir )
    task = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    return task

def ReadSequenceFile( filepath ):
    if filepath[-3:] == ".gz":
        cmd = "zcat %s" % ( filepath )
    else:
        cmd = "cat %s" % ( filepath )
    task = getCommandOutput( cmd )
    return task


def revcomp(dna, reverse):
    bases = 'ATGCNTACGN'
    complement_dict = {bases[i]:bases[i+5] for i in range(5)}
    if reverse:
        dna = reversed(dna)
    result = [complement_dict[base] for base in dna]
    return ''.join(result)

def sort_by_cnt( aDic, _reverse ):
    output = []
    for key in aDic:
        output.append( [ aDic[ key ], key ] )
    output.sort( reverse=_reverse )
    return output


def match_cnt( ref, query ):
    match = 0
    for i in range(len(ref)):
        if ref[i] == query[i]: match += 1
    return match

def motif_cnt( f, pattern = "\w+TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)AGTTTATATTATAACACTTTAACCTATGGC\w+" ):
    motif_cnt_dic = {}
    all_cnt = 0
    match_cnt = 0

    for seq in f.xreadlines():
        all_cnt += 1
        m = re.match( pattern, seq )
        if m == None: continue
        match_cnt += 1
        motif_rc = m.group( "motif" )
        motif = revcomp( motif_rc, True )
        motif_cnt_dic[ motif ] = motif_cnt_dic.get( motif, 0 ) + 1

    return motif_cnt_dic, match_cnt, all_cnt


def ReadSequenceFiles( output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path, pattern = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)GACCGGAACTTCTATGATCGAGATCGAGATCGAGATCGCGGCCGCAAC\w*" ):
    uncut_motif_cnt_dic = {}
    DpnI_motif_cnt_dic = {}
    MboI_motif_cnt_dic = {}

    # 1. Read No enzyme treated file
    print "[ Read No Enzyme Treated Sequence File ]", uncut_fastq_path
    fastq = ReadSequenceFile( uncut_fastq_path )
    [ uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    # 2. Read DpnI treated file
    print "[ Read DpnI Treated Sequence File ]", dpnI_fastq_path
    fastq = ReadSequenceFile( dpnI_fastq_path )
    [ DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    # 3. Read MboI treated file
    print "[ Read MboI Treated Sequence File ]", mboI_fastq_path
    fastq = ReadSequenceFile( mboI_fastq_path )
    [ MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    return uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt, DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt, MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt


def CalculateDAMRatio( experiment_type, output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path ):
    # User needs to specify which experiment it is among followings: ( promoter / utr_with_strong_promoter / utr_with_weak_promoter )
    print "[ Starting ELMSeq Analysis ]", "It might take a long time (more than hours) depending on sequencing size"

    assert( experiment_type in [ "promoter", "utr_with_strong_promoter", "utr_with_weak_promoter" ] )
    if experiment_type == "promoter": 
        pattern = PROMOTER_SEQ_PATTERN
    if experiment_type == "utr_with_strong_promoter": 
        pattern = UTR_STRONG_PROMOTER_SEQ_PATTERN
    if experiment_type == "utr_with_weak_promoter": 
        pattern = UTR_WEAK_PROMOTER_SEQ_PATTERN

    [ uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt, DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt, MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt ] = ReadSequenceFiles( output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path, pattern )
    _motifs = set( uncut_motif_cnt_dic.keys() )
    _motifs.update( DpnI_motif_cnt_dic.keys() )
    _motifs.update( MboI_motif_cnt_dic.keys() )

    print "[ Calculating DAMRatios ]"

    ################################################
    # Random sequence length and TATAAT motif filter
    ################################################
    motifs = set()
    if experiment_type == "promoter": 
        for motif in _motifs:
            if len( motif ) != PROMOTER_RANDOM_SEQ_LEN: continue
            if USE_PROMOTER_TATTAT_FILTER and motif[25:25+6] != "TATAAT": continue
            motifs.add( motif )
    else:
        for motif in _motifs:
            if len( motif ) != UTR_RANDOM_SEQ_LEN: continue
            motifs.add( motif )

    ################################################
    # Calculate CPM and DAMRatio
    ################################################
    data_all = []
    for motif in motifs:
        uncut_cnt = uncut_motif_cnt_dic.get( motif, 0 )
        DpnI_cnt = DpnI_motif_cnt_dic.get( motif, 0 )
        MboI_cnt = MboI_motif_cnt_dic.get( motif, 0 )
        CPM_uncut = ( float( uncut_cnt ) + ALPHA ) / uncut_match_cnt * (10**6)
        CPM_DpnI = ( float( DpnI_cnt ) + ALPHA) / DpnI_match_cnt * (10**6)
        CPM_MboI = ( float( MboI_cnt ) + ALPHA ) / MboI_match_cnt * (10**6)
        data_all.append( [ CPM_MboI / CPM_DpnI, motif, uncut_cnt, DpnI_cnt, MboI_cnt, CPM_uncut, CPM_DpnI, CPM_MboI ] )
    data_all.sort( reverse = True )

    ################################################
    # Write output files 
    ################################################
    output_all_filepath = "%s/%s.all.txt" % ( output_folder, experiment_type )
    output_filepath = "%s/%s.txt" % ( output_folder, experiment_type )
    fout = open( output_filepath, "w" )
    if WRITE_ALL_DAMRATIOS: fout_all = open( output_all_filepath, "w" )
    print >> fout, "# uncut_all_cnt\tuncut_match_cnt\tDpnI_all_cnt\tDpnI_match_cnt\tMboI_all_cnt\tMboI_match_cnt"
    print >> fout, "# %d\t%d\t%d\t%d\t%d\t%d" % ( uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt )
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# uncut_all_cnt\tuncut_match_cnt\tDpnI_all_cnt\tDpnI_match_cnt\tMboI_all_cnt\tMboI_match_cnt"
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# %d\t%d\t%d\t%d\t%d\t%d" % ( uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt )
    print >> fout, "# DAMRatio\tmotif\tuncut_cnt\tDpnI_cnt\tMboI_cnt\tuncut_CPM\tDpnI_CPM\tMboI_CPM"
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# DAMRatio\tmotif\tuncut_cnt\tDpnI_cnt\tMboI_cnt\tuncut_CPM\tDpnI_CPM\tMboI_CPM"
    for DAMRatio, motif, uncut_cnt, DpnI_cnt, MboI_cnt, CPM_uncut, CPM_DpnI, CPM_MboI in data_all:
        if WRITE_ALL_DAMRATIOS: print >> fout_all, "\t".join( [ "%e"%DAMRatio, motif, "%d"%uncut_cnt, "%d"%DpnI_cnt, "%d"%MboI_cnt, "%f"%CPM_uncut, "%f"%CPM_DpnI, "%f"%CPM_MboI ] )
        if uncut_cnt >= READ_COUNT_CUT_OFF or DpnI_cnt >= READ_COUNT_CUT_OFF or MboI_cnt >= READ_COUNT_CUT_OFF:
            print >> fout, "\t".join( [ "%e"%DAMRatio, motif, "%d"%uncut_cnt, "%d"%DpnI_cnt, "%d"%MboI_cnt, "%f"%CPM_uncut, "%f"%CPM_DpnI, "%f"%CPM_MboI ] )
    fout.close()
    if WRITE_ALL_DAMRATIOS: fout_all.close()

    print "[ Finishing ELMSeq Analysis ]"


if __name__ == "__main__":
    CalculateDAMRatio( EXPERIMENT_TYPE, OUTPUT_FOLDER, UNCUT_FASTQ_PATH, DPNI_FASTQ_PATH, MBOI_FASTQ_PATH )

