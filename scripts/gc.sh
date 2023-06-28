#!/bin/bash

# Calculates global GC value.
# Requires fasta file as argument

full=`grep -v '>' $1 | tr -d '\n' | wc -c `
gc=`grep -v '>' $1 | tr -d '\n' | tr -d 'aAtT' | wc -c`
echo "scale=2; $gc / $full" | bc
