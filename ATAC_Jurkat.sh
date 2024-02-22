#!/bin/sh
#ATAC seq analysis
#Sample 1 Jurkat_WT
chilin simple --threads 8 -i Jurkat_WT -o Jurkat_WT -u suiyue -s hg38 -t "Jurkat_WT_1.fq.gz,Jurkat_WT_2.fq.gz" -r dnase -p both --pe

#Sample 2 Jurkat_FBXO42KO
chilin simple --threads 8 -i Jurkat_FBXO42KO -o Jurkat_FBXO42KO -u suiyue -s hg38 -t "Jurkat_FBXO42KO_1.fq.gz,Jurkat_FBXO42KO_2.fq.gz" -r dnase -p both --pe
