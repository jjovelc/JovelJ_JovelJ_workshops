#!/bin/bash

infile=$1
threshold=$2

awk -v threshold="$threshold" '{
    for (i = 2; i <= NF; i++) {
        sum += $i
    } 
    if ((sum / (NF - 1) >= threshold) || NF == 1) {
        print
    } 
    sum = 0
}' < "$infile"
