#!/bin/bash
############################################################
# Functions                                                #
############################################################
Help()
{
   # Display Help
   echo "This script parses the MCScanX HTML file to identify tandem arrays.."
   echo
   echo "Usage: parse_tandem_dups.sh [options] HTML_File"
   echo
   echo "options:"
   echo "o    Output suffix for file name."
   echo 
   echo 
   echo
   echo "arguments:"
   echo "HTML_File    HTML file from MCScanX that highlights tandem arrays in red."
   echo
   echo 
   echo
}

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
  output="tandemDups"
  
  ############################################################
  # Process the input options. Add options as needed.        #
  ############################################################
  # Get the options
  while getopts ":h:o:" option; do
     case $option in
        h) # display Help
           Help
           exit;;
        o) # Output
           output=$OPTARG;;
       \?) # Invalid option
           echo "Error: Invalid option"
           exit;;
     esac
  done
  
  ###########################################
  # Process the positional arguments        #
  ###########################################
  # FASTA file
  HTML=${@:$OPTIND:1}
  
  cat $HTML |\
  grep "^<tr" |\
  cut -f 3 -d ' ' |\
  sed -e 's/\<td\>.*//g' -e 's/\<\/td\>//g'|\
  sed 's/>/ /g' |\
  tr -d "'" |\
  sed 's/bgcolor=//g' |\
  sed 's/\#ee0000/Tandem/g' |\
  sed 's/\#dddddd/Non-Tandem/g' |\
  awk '{print $1"\t"$2}' |\
  nl |\
  sed 's/^[[:space:]]*//' |\
  ( echo -e 'Position\tDup_Type\tGene'; cat -- )> ${HTML%.html}_${output}.txt
  
fi
