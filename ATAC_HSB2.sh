#!/bin/sh
#ATAC seq analysis
#Sample 1 HSB2_WT
chilin simple --threads 8 -i HSB2_WT -o HSB2_WT -u suiyue -s hg38 -t "HSB2_WT_1.fq.gz,HSB2_WT_2.fq.gz" -r dnase -p both --pe

#Sample 2 HSB2_FBXO42KO
chilin simple --threads 8 -i HSB2_FBXO42KO -o HSB2_FBXO42KO -u suiyue -s hg38 -t "HSB2_FBXO42KO_1.fq.gz,HSB2_FBXO42KO_2.fq.gz" -r dnase -p both --pe
