#!/usr/bin/env python3

import os
import gzip
from pathlib import Path
import string

data_dir = os.path.join("..", "data")
process_data_dir = os.path.join("..", "processed_data")

file_1 = "fragments_use.B1.bed.gz"
file_2 = "ATAC_all.hg19.fragments.B2.tsv.gz"
file_3 = "ATAC.hg19.fragments.B3.tsv.gz"

file_output_1 = "ATAC.hg19.fragments.all.updated.B1.tsv.gz"
file_output_2 = "ATAC.hg19.fragments.all.updated.B2.tsv.gz"
file_output_3 = "ATAC.hg19.fragments.all.updated.B3.tsv.gz"

Path(process_data_dir).mkdir(parents=True, exist_ok=True)


# for file_3, add ",B3" at the end of each cell ID (4th column)
with gzip.open(os.path.join(data_dir, file_3),'rt') as input_f, \
    gzip.open(os.path.join(process_data_dir, file_output_3),'wt') as output_f:
    for line in input_f:
        tmp = str.split(line, "\t")
        tmp[3] = tmp[3] + ",B3"
        tmp = "\t".join(tmp)
        
        output_f.write(tmp)

# for file_2, add ",B2" at the end of each cell ID (4th column)
with gzip.open(os.path.join(data_dir, file_2),'rt') as input_f, \
    gzip.open(os.path.join(process_data_dir, file_output_2),'wt') as output_f:
    for line in input_f:
        tmp = str.split(line, "\t")
        tmp[3] = tmp[3] + ",B2"
        tmp = "\t".join(tmp)
        
        output_f.write(tmp)

# for file_1, 
# first add an extra digit for R1/R2/R3 to keep it consistent with the previous two files
# then add ",B1" at the end of each cell ID (4th column)
with gzip.open(os.path.join(data_dir, file_1),'rt') as input_f, \
    gzip.open(os.path.join(process_data_dir, file_output_1),'wt') as output_f:
    for line in input_f:
        tmp = str.split(line, "\t")
        tmp[3] = tmp[3][:3] + '0' + tmp[3][3:9] + '0' + tmp[3][9:15] + '0' + tmp[3][15:]
        tmp[3] = tmp[3] + ",B1"
        tmp = "\t".join(tmp)
        
        output_f.write(tmp)

