#!/bin/bash -i
eval "$(conda shell.bash hook)"
conda activate sina

sina -i output/tmp_fasta.fasta -r ~/PhD/Data_Storage/AlignmentDBs/SILVA138.1/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb -o output/tmp_alignment.fasta

conda deactivate
