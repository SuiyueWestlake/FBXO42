dir_fq=~/projects/dux4/data/mgy/chipseq
dir_project=~/projects/dux4/analysis/mgy/chipseq

mkdir -p ${dir_project}/{logs,bowtie2,macs3,homer2}/{DUX4_IGH,Vector,DUX4_IGH_RAG1_RAG2,RAG1_RAG2,RAG1_RAG2_H3}/
mkdir -p ${dir_project}/fastqc/

conda activate sim
nproc=0
for s in DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3
do # DUX4_IGH Vector 
fastqc -t 20 ${dir_fq}/${s}_R1.fq.gz ${dir_fq}/${s}_R2.fq.gz -o ${dir_project}/fastqc/ 1>${dir_project}/logs/${s}/fastqc 2>&1 &

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done

conda activate assw
nproc=0
for s in DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3
do # DUX4_IGH Vector 
(bowtie2 -p 20 --rg-id "${s}" --rg "ID:${s}" --rg "SM:${s}" --rg "LB:DNA" --rg "PL:ILLUMINA" -x ~/doc/reference/bowtie2/index -1 ${dir_fq}/${s}_R1.fq.gz -2 ${dir_fq}/${s}_R2.fq.gz -S ${dir_project}/bowtie2/${s}/${s}.sam && echo "bowtie2 finished successfully" && \
samtools view -@ 20 -b -o ${dir_project}/bowtie2/${s}/${s}.bam ${dir_project}/bowtie2/${s}/${s}.sam && echo "samtools sam2bam finished successfully" && \
samtools sort -@ 20 -o ${dir_project}/bowtie2/${s}/${s}.sorted.bam -T ${dir_project}/bowtie2/${s}/tmp.${s} ${dir_project}/bowtie2/${s}/${s}.bam && echo "samtools sort finished successfully" && \
samtools index ${dir_project}/bowtie2/${s}/${s}.sorted.bam && echo "samtools index sorted bam finished successfully" && \
picard MarkDuplicates REMOVE_DUPLICATES=false I=${dir_project}/bowtie2/${s}/${s}.sorted.bam O=${dir_project}/bowtie2/${s}/${s}.sorted.mkdup.bam M=${dir_project}/bowtie2/${s}/marked_dup_metrics.txt && echo "picard markduplicates finished successfully" && \
samtools index ${dir_project}/bowtie2/${s}/${s}.sorted.mkdup.bam && echo "samtools index mkdup bam finished successfully" && \
picard MarkDuplicates REMOVE_DUPLICATES=true I=${dir_project}/bowtie2/${s}/${s}.sorted.bam O=${dir_project}/bowtie2/${s}/${s}.sorted.rmdup.bam M=${dir_project}/bowtie2/${s}/removed_dup_metrics.txt && echo "picard removeduplicates finished successfully" && \
samtools index ${dir_project}/bowtie2/${s}/${s}.sorted.rmdup.bam && echo "samtools index rmdup bam finished successfully" &) 1> ${dir_project}/logs/${s}/mapping 2>&1

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done


treats=(DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3) # DUX4_IGH 
control=Vector

source ~/MyPythonEnv/bin/activate
nproc=0
for s in ${treats[@]}
do
bam_treat=${dir_project}/bowtie2/${s}/${s}.sorted.rmdup.bam
bam_control=${dir_project}/bowtie2/${control}/${control}.sorted.rmdup.bam
macs3 callpeak --SPMR -B -f BAMPE -g hs -q 0.01 --keep-dup all -t ${bam_treat} -c ${bam_control} --outdir ${dir_project}/macs3/${s}/ -n ${s}_${control} --seed 1 --tempdir ~/tmp/ 1>${dir_project}/logs/${s}/macs3 2>&1 &

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done
deactivate


treats=(DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3) # DUX4_IGH 
control=Vector

conda activate chip
nproc=0
for s in ${treats[@]}
do
sort -r -g -k 5 ${dir_project}/macs3/${s}/${s}_${control}_summits.bed > ${dir_project}/macs3/${s}/${s}_${control}_sorted_summits.bed && \
head -n 5000 ${dir_project}/macs3/${s}/${s}_${control}_sorted_summits.bed | cut -f 1,2,3,4,9 > ${dir_project}/macs3/${s}/${s}_${control}_peaks_top_motif.bed && \
findMotifsGenome.pl ${dir_project}/macs3/${s}/${s}_${control}_peaks_top_motif.bed hg38 ${dir_project}/homer2/${s}/ -size 200 -p 20 -mis 2 -mask -S 25 1>${dir_project}/logs/${s}/homer2 2>&1 &

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done


conda activate assw
gtf=~/doc/reference/gtf/gencode.v32.annotation.gtf

for s in DUX4_IGH DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3
do
bam=${dir_project}/bowtie2/${s}/${s}.sorted.rmdup.bam

bamCoverage \
-p 10 \
-b ${bam} \
-o ${bam}.normalized.bw \
--binSize 20 \
--normalizeUsing RPKM \
--smoothLength 60 \
--extendReads \
--centerReads \
--ignoreForNormalization chrX chrM &
done
# --ignoreDuplicates
wait

computeMatrix \
reference-point \
-p 10 \
--referencePoint TSS \
-R ${gtf} \
-S ${dir_project}/bowtie2/DUX4_IGH/DUX4_IGH.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2/RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2_H3/RAG1_RAG2_H3.sorted.rmdup.bam.normalized.bw \
-b 2000 \
-a 2000 \
--skipZeros \
-o ${dir_project}/deeptools/merged_rmdup.matrix_TSS.gz \
--outFileSortedRegions ${dir_project}/deeptools/merged_rmdup.regions_TSS.bed

plotHeatmap \
-m ${dir_project}/deeptools/merged_rmdup.matrix_TSS.gz \
-out ${dir_project}/deeptools/merged_rmdup.TSS_heatmap2.pdf \
--colorMap RdBu \
--perGroup

computeMatrix \
scale-regions \
-p 30 \
-R ~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/target/DUX4_IGH_peaks_target_up.bed \
~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/target/DUX4_IGH_peaks_target_down.bed \
-S ${dir_project}/bowtie2/DUX4_IGH/DUX4_IGH.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2/RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2_H3/RAG1_RAG2_H3.sorted.rmdup.bam.normalized.bw \
-b 2000 \
-a 2000 \
--skipZeros \
-o ${dir_project}/deeptools/merged_rmdup_target.matrix_TSS.gz \
--outFileSortedRegions ${dir_project}/deeptools/merged_rmdup_target.regions_TSS.bed

plotHeatmap \
-m ${dir_project}/deeptools/merged_rmdup_target.matrix_TSS.gz \
-out ${dir_project}/deeptools/merged_rmdup_target.heatmap.pdf \
--colorMap RdBu \
--perGroup

computeMatrix \
reference-point \
--referencePoint center \
-p 30 \
-R ~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/target/DUX4_IGH_summits_target_up.bed \
~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/target/DUX4_IGH_summits_target_down.bed \
-S ${dir_project}/bowtie2/DUX4_IGH/DUX4_IGH.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2/RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
${dir_project}/bowtie2/RAG1_RAG2_H3/RAG1_RAG2_H3.sorted.rmdup.bam.normalized.bw \
-b 2000 \
-a 2000 \
--skipZeros \
-o ${dir_project}/deeptools/merged_rmdup_target_summit.matrix_TSS.gz \
--outFileSortedRegions ${dir_project}/deeptools/merged_rmdup_target_summit.regions_TSS.bed

plotHeatmap \
-m ${dir_project}/deeptools/merged_rmdup_target_summit.matrix_TSS.gz \
-out ${dir_project}/deeptools/merged_rmdup_target_summit.heatmap2.pdf \
--colorMap RdBu


conda activate chip

bam1=${dir_project}/bowtie2/DUX4_IGH/DUX4_IGH.sorted.rmdup.bam
bam2=${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.sorted.rmdup.bam
peaks1=${dir_project}/macs3/DUX4_IGH/DUX4_IGH_Vector_peaks.narrowPeak
peaks2=${dir_project}/macs3/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2_Vector_peaks.narrowPeak

mkdir -p ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/

run_spp.R -c=${bam1} -p=15 -odir=${dir_project}/bowtie2/DUX4_IGH -savp -rf -out=${dir_project}/bowtie2/DUX4_IGH/phantomPeakStatsReps.tab
run_spp.R -c=${bam2} -p=15 -odir=${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2 -savp -rf -out=${dir_project}/bowtie2/DUX4_IGH_RAG1_RAG2/phantomPeakStatsReps.tab

bamToBed -i ${bam1} > ${bam1%%.bam}.bed
bamToBed -i ${bam2} > ${bam2%%.bam}.bed

awk -v OFS="\t" '{if($6 ~ /-/){$2=$3-1;}else{$3=$2+1}; print $0}' ${bam1%%.bam}.bed | \
slopBed -s -l 0 -r 145 -g ~/doc/reference/star_2.7.5a/chrNameLength.txt -i - > ${bam1%%.bam}_ext.bed
awk -v OFS="\t" '{if($6 ~ /-/){$2=$3-1;}else{$3=$2+1}; print $0}' ${bam2%%.bam}.bed | \
slopBed -s -l 0 -r 145 -g ~/doc/reference/star_2.7.5a/chrNameLength.txt -i - > ${bam2%%.bam}_ext.bed

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${peaks1} > ${peaks1}.bed
awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${peaks2} > ${peaks2}.bed
cat ${peaks1}.bed ${peaks2}.bed > ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged.bed
awk -v FS="\t" '{if($5>=20){print $0}}' ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged.bed > ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged_0.01.bed
sort -k1,1 -k2,2n ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged_0.01.bed > ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged_0.01_sorted.bed

intersectBed -a ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged_0.01_sorted.bed -b ${bam1%%.bam}_ext.bed -c > ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/DUX4_IGH.count.dat
intersectBed -a ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/narrowpeaks_merged_0.01_sorted.bed -b ${bam2%%.bam}_ext.bed -c > ${dir_project}/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.count.dat


computeMatrix \
reference-point \
-p 30 \
--referencePoint TSS \
-R ~/projects/dux4/analysis/mgy/chipseq/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/up.bed \
~/projects/dux4/analysis/mgy/chipseq/macs3/DUX4_IGH_vs_DUX4_IGH_RAG1_RAG2/down.bed \
-S ~/projects/dux4/analysis/mgy/chipseq/bowtie2/DUX4_IGH/DUX4_IGH.sorted.rmdup.bam.normalized.bw \
~/projects/dux4/analysis/mgy/chipseq/bowtie2/DUX4_IGH_RAG1_RAG2/DUX4_IGH_RAG1_RAG2.sorted.rmdup.bam.normalized.bw \
-b 2000 \
-a 2000 \
--skipZeros \
-o ~/projects/dux4/analysis/mgy/chipseq/deeptools/up_down.matrix.gz \
--outFileSortedRegions ~/projects/dux4/analysis/mgy/chipseq/deeptools/up_down.regions.bed

plotHeatmap \
--refPointLabel peak \
-m ~/projects/dux4/analysis/mgy/chipseq/deeptools/up_down.matrix.gz \
-out ~/projects/dux4/analysis/mgy/chipseq/deeptools/up_down.heatmap.pdf \
--colorMap RdBu


for s in DUX4_IGH DUX4_IGH_RAG1_RAG2 RAG1_RAG2 RAG1_RAG2_H3
do
chilin simple --threads 20 --pe --skip 1 --dont_remove -u mhj -t "${dir_fq}/${s}_R1.fq.gz,${dir_fq}/${s}_R2.fq.gz" -c "${dir_fq}/Vector_R1.fq.gz,${dir_fq}/Vector_R2.fq.gz" -s hg38 -r tf -p narrow -o ${dir_project}/chilin/bwa_mem/${s} -i ${s} --mapper bwa_mem 1>${dir_project}/logs/${s}/chilin.bwa_mem 2>&1 &
done


BETA plus --gs ~/doc/reference/fa/hg38.fa -p ~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/DUX4_IGH_peaks.bed -e ~/projects/dux4/analysis/mgy/rnaseq/deg/DUX4_VEC.txt -k O -g hg38 --gname2 --info 1,3,7 -o ~/projects/dux4/analysis/mgy/chipseq/beta/distance/ --pn `wc -l ~/projects/dux4/analysis/mgy/chipseq/chilin/DUX4_IGH/DUX4_IGH_peaks.bed` --method distance -n dux4 -d 50000 --df 0.05 --da 1 -c 0.05





















