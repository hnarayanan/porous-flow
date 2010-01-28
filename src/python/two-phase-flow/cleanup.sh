#!/usr/bin/env bash

for file in *.vtu
do
    gawk -f truncate.gawk $file  > tmp_$file
    mv tmp_$file $file
done