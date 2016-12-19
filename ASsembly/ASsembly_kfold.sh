#!/bin/bash

for file in test1 training1 test2 training2 test3 training3; do
cp -r "/media/carlos/ADATA NH13/SimulatedBam_gsnap_ALL/${file}/bam_files" /media/carlos/Data/LocalData/bam_files
sudo bash ASsembly.sh "/media/carlos/Data/LocalData"
mkdir -p "/media/carlos/ADATA NH13/SimulatedBam_gsnap_ALL/${file}/int"
cp -rp /media/carlos/Data/LocalData/int/* "/media/carlos/ADATA NH13/SimulatedBam_gsnap_ALL/${file}/int/" && rm -r /media/carlos/Data/LocalData/int/*
rm -rf /media/carlos/Data/LocalData/bam_files
done
