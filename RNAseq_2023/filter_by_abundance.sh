#!/usr/bin/bash

# Check if input file and threshold are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input-file> <threshold>"
    exit 1
fi

INFILE=$1
THRESHOLD=$2
OUTFILE="${INFILE}.${THRESHOLD}up.tsv"

# AWK command to process the file
awk -v thresh="${THRESHOLD}" '{
    sum=0;
    for(i=2; i<=NF; i++) {
        sum+=$i
    };
    if(sum/(NF-1)>=thresh){
        print
    }
}' "$INFILE" > $OUTFILE

