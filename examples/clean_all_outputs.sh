#!/usr/bin/env sh

# Small helper script to clean output directories after running examples.
for DIRECTORY in $(ls -d *_*/)
do
    echo "Cleaning:" ${DIRECTORY}output/\*
    rm -f ${DIRECTORY}output/*
done
                 
