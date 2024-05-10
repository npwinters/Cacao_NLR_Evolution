#!/bin/bash
# This script aligns NB-ARC domains from classified NLRs for phylogenetic tree construction
############################################################
# Functions                                                #
############################################################
# Display help function
Help()
{
   # Display Help
   echo "This script generates the files necessary create NLR alignments and phylogenies from NB-ARC domains."
   echo "Dependencies: MAFFT v7.407, seqtk v1.3-r106, R v3.6+ with tidyverse and argparse"
   echo
   echo "Usage: generate_nbarc_alignments.sh [options] ORTHO_FILE NBARC_PROT_FILE"
   echo
   echo "options:"
   echo "p    File prefix"
   echo "o    Logical. If true, creates a directory for orthogroup fasta files and alignments"
   echo "a    Logical. If true, alignments are made for the extracted NB-ARC domains."
   echo 
   echo
   echo "arguments:"
   echo "ORTHO_FILE        Orthogroup file from PlantTribes' GeneFamilyClassifier."
   echo "NBARC_PROT_FILE   FASTA file containing protein sequences for only the NB-ARC domains, from classify_nlr.sh."
   echo "NBARC_POS_FILE    Tab-delimited file containing NB-ARC positions for each NLR, from classify_nlr.sh"
   echo 
   echo 
}

# Inner join function
MERGE () 
{
  		awk 'BEGIN {
        		FS = OFS = "\t"
        		}
      		NR == FNR {
        		# while reading the 1st file
        		# store its records in the array f
        		f[$1] = $0
        		next
        		}
      		$1 in f {
        		# when match is found
        		# print all values
        		print f[$1], $0
        		}' $1 $2
 }

############################################################
############################################################
# Main program                                             #
############################################################
############################################################
if [ $# -eq 0 ]; then
    Help
else
  # Set variable defaults
  prefix="fileOut"
  orthogroups=false
  alignments=false
  
  ############################################################
  # Process the input options. Add options as needed.        #
  ############################################################
  # Get the options
  while getopts ":h:p:o:a:" option; do
     case $option in
        h) # display Help
           Help
           exit;;
        p) # Genotype
           prefix=$OPTARG;;
        o) # Genotype
           orthogroups=$OPTARG;;
        a) # Alignments
           alignments=$OPTARG;;
       \?) # Invalid option
           echo "Error: Invalid option"
           exit;;
     esac
  done
  
  ###########################################
  # Process the positional arguments        #
  ###########################################
  # Ortho group file
  ORTHO=${@:$OPTIND:1}
  
  # FASTA file
  FASTA=${@:$OPTIND+1:1}
  
  # File of NB-ARC domain positions
  NBARC=${@:$OPTIND+2:1}
  
  echo
  # Identify genes with multiple NB-ARC domains 
  echo "Identifying the best NB-ARC domains for alignment and extracting protein sequences."
  Rscript choose_longest_nbarc.R --domain_position_file $NBARC --orthogroup_file $ORTHO > ${NBARC%.txt}_longest.txt
  cut -f 1 ${NBARC%.txt}_longest.txt > tmp
  seqtk subseq $FASTA tmp > ${NBARC%.txt}_longest.faa
  rm tmp
  echo
  
  # Conditional--should orthogroup fasta files and alignments be created
  if [[ "$orthogroups" = true ]]
  then
  
    echo "Creating directory for orthogroup fasta files."
    mkdir ${prefix}_orthoFastaFiles
    echo
    
    echo "Subsetting NB-ARC protein file according to orthogroup classification."
    cat $ORTHO | cut -f 2 | tail -n +2 | sort | uniq > unique_orthos.txt
    echo
    
    while read ortho; do 
      awk -v var=$ortho '$2 == var' ${NBARC%.txt}_longest.txt | cut -f 1 > ${prefix}_orthoFastaFiles/${ortho}_geneID.txt
    done < unique_orthos.txt
    
    for file in $(ls ${prefix}_orthoFastaFiles/*_geneID.txt); do 
      seqtk subseq ${NBARC%.txt}_longest.faa $file > ${file%.txt}.faa
    done
    
  fi
  
  if [[ "$alignments" = true ]]
  then 
    echo "Creating NB-ARC alignments (this may take awhile)"
    echo
    
    # Align all NB-ARC domains using the linsi algorithm
    mafft --maxiterate 1000 --localpair ${NBARC%.txt}_longest.faa > ${NBARC%.txt}_longest.aln
    
    # Align NB-ARC domains in each orthogroup using the linsi algorithm
    for file in $(ls ${prefix}_orthoFastaFiles/*.faa); do 
      mafft --maxiterate 1000 --localpair $file > ${file%.faa}.aln
    done
  fi

fi 
