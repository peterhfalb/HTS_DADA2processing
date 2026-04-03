#!/bin/bash
set -euo pipefail 
usage() { echo "Usage: $0 [-m < 16S|ITS|18S>] [-q <good|bad>]" 1>&2; exit 1; }
Help(){
# Display Help
	echo "This script starts a dada2 script"
	echo "options:"
	echo "-m Method valid options: [16S|ITS|18S]"
	echo "-q quality valid options: [good|bad] default(good)"
}	

quality=good
while getopts "m:q:" OPTION; 
do
    case $OPTION in
        m)
            if [[ $OPTARG =~ ^(16S|ITS|18S)$ ]]; then
            METHOD=${OPTARG}
            else 
            echo "-m must be one of 16S|ITS|18S"
                exit 2
            fi
            ;;
        q)
            if [[ $OPTARG =~ ^(good|bad)$ ]]; then
            QUALITY=${OPTARG}
            else 
            echo "-q must be one of good|bad"
            	exit 2
            fi
            ;;
        \?)
        printf "Unknown option: -%s\n" $OPTARG
        Help >&2
        exit 1
        ;;
        :) printf "missing argument for -%s\n" $OPTARG >&2
        Help >&2
        exit 1
        ;;
        h | * )
        Help >&2
        exit 1
        ;;
    esac
done

echo "running: $METHOD"
BASEDIR=$(dirname "$0")
if [[ $METHOD == "16S" ]] ; then Rscript $BASEDIR/run_16S_dada2.R $QUALITY
fi
if [[ $METHOD == "ITS" ]] ; then Rscript $BASEDIR/run_ITS_dada2.R $QUALITY
fi
if [[ $METHOD == "18S" ]] ; then Rscript $BASEDIR/run_18S_dada2.R $QUALITY
fi
