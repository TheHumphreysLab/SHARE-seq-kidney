#at the shell interface, combine the three files together:
cat ATAC.hg19.fragments.all.updated.B1.tsv.gz \
ATAC.hg19.fragments.all.updated.B2.tsv.gz \
ATAC.hg19.fragments.all.updated.B3.tsv.gz > ATAC.hg19.fragments.B123.tsv.gz

gunzip ATAC.hg19.fragments.B123.tsv.gz

#sort the file
TMPDIR=../tmp_sort/ sort -k1,1 -k2,2n -k3,3n -k4,4 --parallel=8 -S 40G ATAC.hg19.fragments.B123.tsv > ATAC.hg19.fragments.B123.sorted.tsv

bgzip ATAC.hg19.fragments.B123.sorted.tsv

tabix -p bed ATAC.hg19.fragments.B123.sorted.tsv.gz