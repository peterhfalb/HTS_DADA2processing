#!/bin/bash
Help(){
# Display Help
	echo "This script checks what 16S, ITS, or 18S primers you have"
	echo "primercheck.sh sample_R1_L001.fastq.gz sample_R2_L001.fastq.gz"
}	
while getopts ":h" OPTION;
do
    case $OPTION in
        h | * )
        Help >&2
        exit 1
        ;;
    esac
done

        while IFS=' ' read -r PRIME LOC || [ -n "${PRIME}" ]; do
        check=$(zgrep --max-count=100 -c ${PRIME} $1);
        if [ $check == 100 ]; then
          echo "found: $PRIME"
          R1=$(echo $PRIME | sed 's/\.\*/N/g')
          break 2
        exit
fi
done < /users/4/goul0109/dada2processing/list_to_search_R1.txt

        while IFS=' ' read -r PRIME LOC || [ -n "${PRIME}" ]; do
        check=$(zgrep --max-count=100 -c ${PRIME} $2);
        if [ $check == 100 ]; then
            R2=$(echo $PRIME | sed 's/\.\*/N/g')
            echo "found: $PRIME"
            break 2
         exit
    fi
    done < /users/4/goul0109/dada2processing/list_to_search_R2.txt
echo "recommended primer trimming command:"
echo "../dada2processing/primer_trim.sh -p ${R1} -q ${R2}"

# usage primercheck3.sh forwardread reverseread
