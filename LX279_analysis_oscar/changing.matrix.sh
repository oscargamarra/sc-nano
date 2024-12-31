#!/bin/bash

#This script changes the MMF created with R
#In 'transform_to_mmf.R', we changed 0s to -1
#We set 0 values to -1, now we change them back to 0.
#We do this because Souporcell MMF includes 0s, but traditional MMFs do not.

# enter input and output files
input_file="market.ref.mtx" #Output of 'transform_to_mmr.R'
output_file="market.ref.sop.mtx"

# extract first line
head -n 1 "$input_file" > "$output_file"

# add an extra line as souporcell does
echo "% written by sprs" >> "$output_file"

# change rest of mtx, substituting -1 for 0
tail -n +3 "$input_file" | sed 's/\(^[0-9]\+\s\+[0-9]\+\s\+\)-1$/\10/' >> "$output_file"

echo "File $output_file is ready for Souporcell"
