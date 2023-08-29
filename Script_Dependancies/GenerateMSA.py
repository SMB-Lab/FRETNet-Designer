"""
Author: Frank Duffy
Affil: Sanabria SMB Lab
Last updated: 8/24/2023

def generate_MSA(settings)
    
    Takes a settings dictionary containing the file location of the seed 
    sequence and hmm, easel, and pfam filter commands in addition to other user 
    variables and builds an hmm profile, performs an hmm search against 
    the provided database with that profile, reformats the produced alignment 
    file as afa, and filters the alignment file by a max gap number.

"""

import subprocess as sp

def generate_msa(settings):
    
    #build hmm using local hmmbuild executable
    command = settings['hmmer_bin'] + 'hmmbuild'
    hmm_name = settings['outputHMM']
    seed_seq = settings['input_seed']
    cmd = [command, hmm_name, seed_seq] 
    sp.run(cmd)

    #get alignments using local HMMER executable
    command = settings['hmmer_bin'] + 'hmmsearch'
    flag = '-A'
    #hmm_name = hmm_name
    database = settings['database_path']
    cmd = [command, flag, "get_alignments_tempfile1", hmm_name, database] 
    sp.run(cmd)

    #reformat alignment file as afa
    command = settings['easel_miniapps'] + 'esl-reformat'
    flag = '-o'
    output_file_name = settings['outputMSA_filename']
    frmt = 'afa'
    cmd = [command, flag, output_file_name, frmt, "get_alignments_tempfile1"] 
    sp.run(cmd)

    #filter afa formatted alignment file with python script
    command = 'python'
    script = settings['filter_pfam_path']
    output_file_name = settings['outputMSA_filename']
    max_gaps = settings['max_align_gaps']
    cmd = [command, script, output_file_name, max_gaps] 
    sp.run(cmd)