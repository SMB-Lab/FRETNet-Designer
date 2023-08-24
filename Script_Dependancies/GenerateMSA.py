"""
Author: Frank Duffy
Affil: Sanabria SMB Lab
Last updated: 8/24/2023

def generate_MSA(settings)
    
    Takes a settings dictionary and makes a filtered MSA alignment file.

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
    hmm_name = hmm_name
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
    max_gaps = '30'
    cmd = [command, script, output_file_name, max_gaps] 
    sp.run(cmd)