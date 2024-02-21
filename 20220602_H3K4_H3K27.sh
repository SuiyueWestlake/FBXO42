#!/bin/sh
#Sample 1(Jurkat WT H3K4m1) VS 7(Jurkat Igg)
chilin simple --threads 8 -i Sample1vs7 -o Sample1vs7 -u suiyue -s hg38 -t "1_FKDL210222437-1a_1.fq.gz,1_FKDL210222437-1a_2.fq.gz" -c "7_FKDL210222443-1a_1.fq.gz,7_FKDL210222443-1a_2.fq.gz" -p both --pe

#Sample 2(Jurkat FBXO42KO H3K4m1) VS 8(Jurkat FBXO42KO Igg)
chilin simple --threads 8 -i Sample2vs8 -o Sample2vs8 -u suiyue -s hg38 -t "2_FKDL210222438-1a_1.fq.gz,2_FKDL210222438-1a_2.fq.gz" -c "8_FKDL210222444-1a_1.fq.gz,8_FKDL210222444-1a_2.fq.gz" -p both --pe

#Sample 3(Jurkat WT H3K27AC) VS 7(Jurkat Igg)
chilin simple --threads 8 -i Sample3vs7 -o Sample3vs7 -u suiyue -s hg38 -t "3_FKDL210222439-1a_1.fq.gz,3_FKDL210222439-1a_2.fq.gz" -c "7_FKDL210222443-1a_1.fq.gz,7_FKDL210222443-1a_2.fq.gz" -p both --pe

#Sample 4(Jurkat FBXO42KO H3K27AC) VS 8(Jurkat FBXO42KO Igg)
chilin simple --threads 8 -i Sample4vs8 -o Sample4vs8 -u suiyue -s hg38 -t "4_FKDL210222440-1a_1.fq.gz,4_FKDL210222440-1a_2.fq.gz" -c "8_FKDL210222444-1a_1.fq.gz,8_FKDL210222444-1a_2.fq.gz" -p both --pe
