#!/bin/bash
# Noah Winters, August 2021
# This script classifies NLRs from a fasta file using cacao-specific HMM classifiers and NLR-parser
# Dependencies: HMMER v3.3, MEME Suite v4.9.1, JRE v1.6+, seqtk v1.3-r106

if [ $# -eq 0 ]
  then
    echo "----------------------------------------------------"
    echo "This script classifies NLRs from a fasta file using cacao-specific HMM classifiers and NLR-parser"
    echo "Usage: bash classify_nlr.sh <arguments>"
    echo "Arguments: PROT.FASTA
           CDS.FASTA
           HMM_CUTOFF--default 0.0001
           PATH_TO_MEME_XML_FILE
           PATH_TO_NLRPARSER"
    echo "----------------------------------------------------"
else

        PROT=$1
        CDS=$2
        HMM_CUTOFF=$3
        MEME_XML=$4
        NLRPARSER=$5
        
        ##########################
        ### HMM Classification ###
        ##########################
        
        # NB-ARC domains
        hmmsearch\
          --domtblout ${PROT%.fasta}_nbarc_domTable.out\
          --incE $HMM_CUTOFF\
          cacao_nbarc.hmm\
          $PROT > ${PROT%.fasta}_nbarc.out
        
        # RPW8 domains
        hmmsearch\
          --domtblout ${PROT%.fasta}_rpw8_domTable.out\
          --incE $HMM_CUTOFF\
          cacao_rpw8.hmm\
          $PROT > ${PROT%.fasta}_rpw8.out
        
        # TIR domains
        hmmsearch\
          --domtblout ${PROT%.fasta}_tir_domTable.out\
          --incE $HMM_CUTOFF\
          cacao_tir.hmm\
          $PROT > ${PROT%.fasta}_tir.out

        # Parse domTable output from HMM for each NBARC, RPW8, and TIR domains
          for file in $(ls *_domTable.out)
          do 
            cat $file |\
            sed '/^#,*/d' |\
            tr -s " " |\
            cut -f 1,7,20,21 -d " " |\
            sed 's/ / /g' |\
            awk -v var=$HMM_CUTOFF ' $2 < var  ' > tmp
            echo -e 'Gene\tE-Value\tDom_Start\tDom_End' | cat - tmp  > ${file%_domTable.out}_domPos.txt
            rm tmp
          done
        
        #################################
        ### NLR-Parser Classification ###
        #################################
        
        # Run proteins through mast -- need NLR-Parser XML file
        mast $MEME_XML $PROT\
          -o ${PROT%.fasta}_mastOut
          
        mv ${PROT%.fasta}_mastOut/mast.xml ${PROT%.fasta}_mastOut/${PROT%.fasta}_mastOut.xml
        
        # Run NLR-Parser -- path specific as argument
        java -jar $NLRPARSER\
          -i ${PROT%.fasta}_mastOut/${PROT%.fasta}_mastOut.xml\
          -o ${PROT%.fasta}_nlrparser.txt
          
          
        ###########################
        ### Combine Annotations ###
        ###########################
        
        # Combine NB-ARC, RPW8, and TIR domain annotations via HMM
        cat ${PROT%.fasta}_nbarc_domPos.txt\
          ${PROT%.fasta}_rpw8_domPos.txt\
          ${PROT%.fasta}_tir_domPos.txt |\
          cut -f 1 -d " "|\
          sort |\
          uniq > ${PROT%.fasta}_nbarc_rpw8_tir_hmm.txt
          
          
        # Combine NB-ARC, RPW8, and TIR domain annotations via HMM
        cat ${PROT%.fasta}_nlrparser.txt | cut -f 1,2,8 | awk '$3 ~ /16/ || $3 ~ /17/' > tmp
        echo -e 'Gene\tnlr_parser_class\tmotiflist' | cat - tmp > ${PROT%.fasta}_nlrparser_ccOnly.txt
        cat ${PROT%.fasta}_nlrparser_ccOnly.txt\
          ${PROT%.fasta}_nbarc_rpw8_tir_hmm.txt |\
          cut -f 1 |\
          sort |\
          uniq |\
          grep -v -e "SequenceName" -e "Gene" > ${PROT%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.txt
        rm tmp
        
      # Extract the protein sequences for all NLRs
      seqtk subseq $PROT ${PROT%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.txt > ${PROT%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.faa
        
      # Extract the CDS sequences for all NLRs
      seqtk subseq $CDS ${PROT%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.txt > ${CDS%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.fna
        
      # Extract NBARC domains from protein sequences
      cat ${PROT%.fasta}_nbarc_domPos.txt | tr -s " " | cut -f 1,3,4 -d " " | sed 's/ / /g' > tmp
      seqtk subseq ${PROT%.fasta}_nbarc_rpw8_tir_hmm_nlrparser.faa tmp > ${PROT%.fasta}_nbarc_domPos.faa
      rm tmp
          
fi
