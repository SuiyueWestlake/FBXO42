#!/bin/sh
#ATAC seq analysis
#Sample 1 293T_WT
chilin simple --threads 8 -i 293T_WT -o 293T_WT -u suiyue -s hg38 -t "293T_WT_1.fq.gz,293T_WT_2.fq.gz" -r dnase -p both --pe

#Sample 2 293T_FBXO42KO
chilin simple --threads 8 -i 293T_FBXO42KO -o 293T_FBXO42KO -u suiyue -s hg38 -t "293T_FBXO42KO_1.fq.gz,293T_FBXO42KO_2.fq.gz" -r dnase -p both --pe
