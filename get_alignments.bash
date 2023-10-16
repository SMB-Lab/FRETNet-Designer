#This program is meant to read in an HMM Profile file and search it against the UniProtKB database.  
#It returns an AFA file containing a multiple sequence alignment of each hit.  It should be noted 
#that the searching algorithm logs a different entry for each domain match.  This means that if a 
#single protein seqeunce contains two seperate regions matching the HMM Profile, then there will be
#two seperate entries for that protein, each one containing one part of the total seqeunce that 
#contains the matching sequence.

#!/bin/bash

#Asks for seed sequence or HMM?
printf "\n\n\nDo you have an HMM Profile?  Otherwise it will ask for a file containing seed sequences. (yes/no)\n"
read hmmOrSeed

if [ "$hmmOrSeed" = "yes" ]; then
	#Gets name of hmm profile file
	printf "Please enter the name of the hmm profile file.\n"
	read hmmProfile

	if [ -f "$hmmProfile" ]; then
	   	printf ""
	else
	   	printf "\nFile $hmmProfile does not exist, closing program.\n"
		exit
	fi
elif [ "$hmmOrSeed" = "no" ]; then
	#Gets name of seed seqeunces file
	printf "Please enter the name of the seed seqeunces file.\n"
	read seedSeq

	if [ -f "$seedSeq" ]; then
	   	#Gets name of hmm profile file
		printf "What would you like your HMM Profile file to be called?.\n"
		read hmmProfile
	else
	   	printf "\nFile $seedSeq does not exist, closing program.\n"
		exit
	fi
else
	printf "\nYou did not enter either "yes" or "no", closing program.\n"
	exit
fi

#Gets name of hmm profile file
printf "Do you want to limit the search to just the annotated sequences database?  Otherwise it will search the total database. (yes/no)\n"
read whichDatabase

if [ "$whichDatabase" = "yes" ]; then
	database="/zfs/smblab/group_software/HMMER/db/uniprot_sprot.fasta.gz"
elif [ "$whichDatabase" = "no" ]; then
	database="/zfs/smblab/group_software/HMMER/db/uniprot_trembl.fasta.gz"
elif [ "$whichDatabase" = "pfam" ]; then
	database="/zfs/smblab/group_software/HMMER/db/Pfam-A.fasta.gz"
elif [ "$whichDatabase" = "pdz" ]; then
        database="/zfs/smblab/group_software/HMMER/db/PDZ12_complete.fasta"
else
	printf "\nYou did not enter either "yes" or "no", closing program.\n"
	exit
fi

#Gets maximum number of gaps to be cut out
printf "Please enter the number of maximum number of continuous gaps per sequence.\n"
read maxGaps

reg='^-?[0-9]+([.][0-9]+)?$'
if ! [[ $maxGaps =~ $reg ]] ; then
   	printf "\nInput was not a positive number. Exiting program.\n"
	exit
else
   	printf ""
fi

#Gets name of output
printf "What would you like your output file to be called?\n"
read outputFileName

#Builds HMM profile
if [ "$hmmOrSeed" = "no" ]; then
	printf "\n\nBuilding HMM Profile...\n"

	/zfs/smblab/group_software/HMMER/install/bin/hmmbuild $hmmProfile $seedSeq > /dev/null 2>&1

	printf "\nFinished building HMM Profile...\n"
else
   	printf "\n"
fi

printf "\nSearching the databases using hmmsearch...\n"

#runs hmmsearch
/zfs/smblab/group_software/HMMER/install/bin/hmmsearch -A get_alignments_tempfile1 $hmmProfile $database > get_alignments_tempfile2

printf "\nhmmearch completed!\n"

printf "\nConverting output to AFA.\n"

cp get_alignments_tempfile1 get_alignments_tempfile1_preConverted

#runs esl-reformat
/zfs/smblab/group_software/HMMER/hmmer-3.3.2/easel/miniapps/esl-reformat -o $outputFileName afa get_alignments_tempfile1

printf "\nConversion complete!\n"

printf "\nFiltering output...\n"

python /zfs/smblab/fduffy/source/repos/FRETNet-Designer/filter_pfam_args.py $outputFileName $maxGaps

printf "\nDone filtering output!\n"

#rm get_alignments_tempfile1
#rm get_alignments_tempfile2

