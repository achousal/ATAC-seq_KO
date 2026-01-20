
# bsub -P acc_BiNGS -q private -m hassod01-1 -n 4 -W 12:00 -R rusage[mem=10000] -Is /bin/bash



cd /sc/arion/projects/BiNGS/bings_omics/data/bings/2024/bings_teaching/icahn/bsr2402/atacsseq/raw
raw_dir="/sc/arion/projects/BiNGS/bings_omics/data/bings/2024/bings_teaching/icahn/bsr2402/atacsseq/raw"

bl_regions=$raw_dir/"hg38-blacklist.v2.bed"
bowtie_index=$raw_dir/"bowtie2_index"
homer_preparsed_dir=$raw_dir/"preparsed_dir"


# create directories for each packge in the data directory before running script
fastq_dir="../analysis/fastq"
fastqc_dir="../analysis/fastqc"
trimmgalore_dir="../analysis/trimmgalore"
bowtie2_dir="../analysis/bowtie2"
samtools_dir="../analysis/samtools"
picard_dir="../analysis/picard"
macs_dir="../analysis/macs"
deeptools_dir="../analysis/deeptools"
subread_dir="../analysis/subread"
homer_dir="../analysis/homer"
container_dir="../../../containers"

for file in ${raw_dir}/*_chr5_R1.fastq;  

  do 
  
  file="$(basename $file)"

  #get samples name for sam name
  samplename=$(echo "$file" | sed -e 's/_ATAC_chr5_R1.fastq//') ;

  #status
  echo "processing ${samplename}"




 #################
 #### run fastqc
 #################
echo "running fastqc: ${samplename}"

mkdir -p $fastqc_dir
ml fastqc/0.11.9

start_time=`date +%s`

fastqc --threads 8 \
       ${raw_dir}/${samplename}_ATAC_chr5_R1.fastq ${raw_dir}/${samplename}_ATAC_chr5_R2.fastq \
       --outdir ${fastqc_dir}

end_time=`date +%s`

echo execution time was `expr $end_time - $start_time` s.

ml unload fastqc/0.11.9

#######################################
#### run trimmgalore to remove adapters
#######################################
echo "trimming adaptors: ${samplename}" 


ml trim_galore/0.6.6
mkdi -p $trimmgalore_dir
start_time=`date +%s` 

trim_galore --fastqc \
             --paired \
             --cores 8 \
             --gzip \
             --output_dir ${trimmgalore_dir} \
              ${raw_dir}/${samplename}_ATAC_chr5_R1.fastq ${raw_dir}/${samplename}_ATAC_chr5_R2.fastq


end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload trim_galore/0.6.6

##########################
#### run bowtie2 to align
##########################
echo "aligning: ${samplename}"

mkdir -p $bowtie2_dir
ml bowtie2/2.4.4

start_time=`date +%s` 

bowtie2 -p 8 \
  -x $bowtie_index/"index" \
  -1 ${trimmgalore_dir}/${samplename}_ATAC_chr5_R1_val_1.fq.gz \
  -2 ${trimmgalore_dir}/${samplename}_ATAC_chr5_R2_val_2.fq.gz \
  -S ${bowtie2_dir}/${samplename}_noA.sam \
  -X 2000 \
  2> ${bowtie2_dir}/${samplename}.log 

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload bowtie2/2.4.4

###############################
#### run samtools to index/sort
###############################
echo "filtering/indexing/sorting:  ${samplename}"

ml use /hpc/packages/minerva-rocky9/modulefiles
ml samtools/1.17
mkdir -p $samtools_dir
start_time=`date +%s` 

#remove chrM and keep only "chr" chromosoes
samtools view -h ${bowtie2_dir}/${samplename}_noA.sam | awk  '($3 != "chrM" && $3 != "chrUn")' >  ${bowtie2_dir}/${samplename}_noA_noM_chr_q20.sam

#convert to bam  
samtools view -h -q 20 ${bowtie2_dir}/${samplename}_noA_noM_chr_q20.sam -o ${samtools_dir}/${samplename}_noA_noM_chr_q20.bam

#sort
samtools sort  -@ 8 ${samtools_dir}/${samplename}_noA_noM_chr_q20.bam  -o ${samtools_dir}/${samplename}_noA_noM_chr_q20_sorted.bam 

#index
samtools index -@ 8  -b ${samtools_dir}/${samplename}_noA_noM_chr_q20_sorted.bam 

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload samtools/1.17

############################### 
#### run picard to de-duplicate
###############################
echo "deduplicating: ${samplename}"

ml picard/3.1.1
mkdir -p $picard_dir
start_time=`date +%s`

java -jar $PICARD MarkDuplicates \
          I=${samtools_dir}/${samplename}_noA_noM_chr_q20_sorted.bam  \
          O=${picard_dir}/${samplename}_noA_noM_chr_q20_sorted_nd.bam  \
          REMOVE_DUPLICATES=true \
          VALIDATION_STRINGENCY=LENIENT \
          M=${picard_dir}/${samplename}_noA_noM_chr_q20_sorted_nd.txt 

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload picard/3.1.1

#######################################################
#### run samtools re-index/sort deduplicated bam file
#######################################################
echo "re-indexing/sorting: ${samplename}"

ml samtools/1.17
start_time=`date +%s`

samtools sort  -@ 8 ${picard_dir}/${samplename}_noA_noM_chr_q20_sorted_nd.bam -o ${samtools_dir}/${samplename}_final.bam 

samtools index -@ 8  -b ${samtools_dir}/${samplename}_final.bam  

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload samtools/1.17

#######################################################
#### create bigwigs of final bam files
#######################################################
echo "generating bigwig: ${samplename}"

ml macs/2.1.0
mkdir -p $deeptools_dir
start_time=`date +%s`

bamCoverage  --bam ${samtools_dir}/${samplename}_final.bam \
             --outFileName $deeptools_dir/${samplename}_final.bw \
             --outFileFormat bigwig  \
             --binSize=10 \
             --normalizeUsing RPKM \
             --extendReads=200  \
             --numberOfProcessors 8

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload macs/2.1.0

 done

################################
#### Merging bam files##########
################################

echo "Merging bam files"

ml samtools/1.17


samtools merge ${samtools_dir}/master_atac.bam ${samtools_dir}/*_final.bam
samtools index -@ 8 -b ${samtools_dir}/master_atac.bam 



############################  
## run macs2 to call peaks
############################
echo "Calling peaks"

ml macs/2.1.0
mkdir -p $macs_dir
start_time=`date +%s`

macs2 callpeak --nomodel \
      -t ${samtools_dir}/master_atac.bam  \
      --outdir ${macs_dir} \
      -n master_atac \
      -f BAMPE \
      -g hs \
      --keep-dup all \
      --slocal 1000 \
      2> ${macs_dir}/master_atac_macs2.log


end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload macs/2.1.0


######################
# remove bl regions and generating bed file for gb
#######################

echo "Peak filtering/formating"

ml bedtools/2.31.0


bedtools intersect -a ${macs_dir}/master_atac_peaks.narrowPeak -b $bl_regions -v > ${macs_dir}/master_atac_peaks_bl.narrowPeak

awk '{print $1,$2,$3,$4}' ${macs_dir}/master_atac_peaks_bl.narrowPeak > ${macs_dir}/temp.bed

echo "track name=\"SKmel147_ATAC_master_regions\" description=\"SKmel147_ATAC_master_regions\"" | cat - ${macs_dir}/temp.bed  > ${macs_dir}/master_atac_peaks_bl.bed
rm ${macs_dir}/temp.bed

ml unload bedtools/2.31.0

##################################################  
#### run deeptools bamcoverage to generate bigwig
##################################################
echo "generating bigwig: master bam"

ml macs/2.1.0
start_time=`date +%s`

bamCoverage  --bam ${samtools_dir}/master_atac.bam \
             --outFileName $deeptools_dir/master_atac.bw \
             --outFileFormat bigwig  \
             --binSize=10 \
             --normalizeUsing RPKM \
             --extendReads=200  \
             --numberOfProcessors 4

# for file in ${samtools_dir}/*_final.bam;
  

#   do 

#   file="$(basename $file)"

#   #get samples name for sam name
#   samplename=$(echo "$file" | sed -e 's/.bam//') ;

#   bamCoverage --bam ${samtools_dir}/$file \
#               --outFileName $deeptools_dir/${samplename}.bw \
#               --outFileFormat bigwig  \
#               --binSize=10 \
#               --normalizeUsing RPKM \
#               --extendReads=200  \
#               --numberOfProcessors 4

#  done

# end_time=`date +%s`
# echo execution time was `expr $end_time - $start_time` s.




# ##################################################  
# #### run featurecounts to generate coun matrix
# ##################################################

echo "generating count matrix"

ml subread/2.0.1
mkdir -p $subread_dir

start_time=`date +%s`
awk 'OFS="\t" {print $4, $1, $2, $3, "."}' ${macs_dir}/master_atac_peaks_bl.narrowPeak > ${macs_dir}/master_atac_peaks_bl.saf

featureCounts -p -a ${macs_dir}/master_atac_peaks_bl.saf -F SAF -o ${subread_dir}/master_atac_peaks_bl_subread.txt ${samtools_dir}/*final.bam -T 12

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.




##################################################  
#### run deeptools to generate heatmap
##################################################
echo "generating deeptool heatmap - individual sample"

ml macs/2.1.0

start_time=`date +%s`

#### run deeptools to generate heatmap (split by sample)

awk 'FNR==NR {a[$4]; next} FNR> 1 && $4 in a' ${macs_dir}/master_atac_peaks_bl.bed ${macs_dir}/master_atac_summits.bed > ${macs_dir}/master_atac_summits_bl.bed

computeMatrix reference-point \
    -R ${macs_dir}/master_atac_summits_bl.bed \
    --skipZeros \
    -S ${deeptools_dir}/SKMel147_ARID2WT_final.bw ${deeptools_dir}/SKMel147_ARID2KO_final.bw \
    -o ${deeptools_dir}/master_atac_summits_bl_2_samples.gz \
    -b 5000 -a 5000 \
    --referencePoint center \
    -p 4 \
    --samplesLabel ARID2WT ARID2KO

plotHeatmap -m ${deeptools_dir}/master_atac_summits_bl_2_samples.gz \
    --plotFileFormat pdf \
    -out ${deeptools_dir}/master_atac_summits_bl_2_samples.pdf \
    --outFileSortedRegions ${deeptools_dir}/master_atac_summits_bl_2_samples_sorted.bed \
    --dpi 720 \
    --missingDataColor White \
    --colorMap Reds \
    --regionsLabel peaks \
    --heatmapHeight 13 

#### run deeptools to generate heatmap (merged)

echo "generating deeptool heatmap - merged sample"

computeMatrix reference-point \
    -R ${macs_dir}/master_atac_summits_bl.bed \
    --skipZeros \
    -S ${deeptools_dir}/master_atac.bw  \
    -o ${deeptools_dir}/master_atac_summits_bl_merged.gz \
    -b 5000 -a 5000 \
    --referencePoint center \
    -p 4 \
    --samplesLabel master_atac 

plotHeatmap -m ${deeptools_dir}/master_atac_summits_bl_merged.gz \
    --plotFileFormat pdf \
    -out ${deeptools_dir}/master_atac_summits_bl_merged.pdf \
    --outFileSortedRegions ${deeptools_dir}/master_atac_summits_bl_merged_sorted.bed \
    --dpi 720 \
    --missingDataColor White \
    --colorMap Reds \
    --regionsLabel peaks \
    --heatmapHeight 13 


module unload macs/2.1.0

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


##################################################  
#### homer motif analysis
##################################################

echo "motif enrichment analysis using homer"
ml homer/4.10
mkdir -p $homer_dir



start_time=`date +%s`

bed_path=${macs_dir}/master_atac_summits_bl.bed

findMotifsGenome.pl  $bed_path  hg38 ${homer_dir} -preparsedDir $homer_preparsed_dir  -size 200  -p 6 


end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.




###############################################################  
#### Normalize and compute log2(KO/WT) add differential status
##############################################################

echo "Identify differential accessible regions "

tail -n +3 ${subread_dir}/master_atac_peaks_bl_subread.txt > ${subread_dir}/master_atac_peaks_bl_subread_no_header.txt

#calculate total reads per sample
total_reads_wt=$(awk '{ sum += $8 } END { print sum }' ${subread_dir}/master_atac_peaks_bl_subread_no_header.txt)
total_reads_ko=$(awk '{ sum += $7 } END { print sum }' ${subread_dir}/master_atac_peaks_bl_subread_no_header.txt)

#replace 0s
awk '{ if ($7 == 0) $7 = 1; print }' ${subread_dir}/master_atac_peaks_bl_subread_no_header.txt > ${subread_dir}/master_atac_peaks_bl_subread_no_header_nozero.txt

#calculate rpkm for wt 
awk -v total_reads_wt="$total_reads_wt" '{ printf "%s\t%s\t%.2f\t%.2f\n", $1, $6,  ($8 / ($6 / 1000)) / (total_reads_wt / 1000000), $7 }' ${subread_dir}/master_atac_peaks_bl_subread_no_header_nozero.txt > ${subread_dir}/master_atac_peaks_bl_subread_norm.txt
awk -v total_reads_ko="$total_reads_ko" '{ printf "%s\t%s\t%.2f\t%.2f\n", $1, $2, $3, ($4 / ($2 / 1000)) / (total_reads_ko / 1000000) }' ${subread_dir}/master_atac_peaks_bl_subread_norm.txt > ${subread_dir}/master_atac_peaks_bl_subread_norm_final.txt

#calculate log2fc
awk '{ {$5 = log($4/$3)/log(2)} print }' ${subread_dir}/master_atac_peaks_bl_subread_norm_final.txt > ${subread_dir}/master_atac_peaks_bl_subread_no_header_log2diff.txt


# save up in KO peaks
awk '$5 >= 1' ${subread_dir}/master_atac_peaks_bl_subread_no_header_log2diff.txt > ${subread_dir}/master_atac_peaks_bl_up_KO.txt
awk 'FNR==NR {a[$1]; next} FNR> 1 && $4 in a' ${subread_dir}/master_atac_peaks_bl_up_KO.txt ${macs_dir}/master_atac_summits.bed > ${subread_dir}/master_atac_peaks_bl_up_KO_summit.bed

# save down in KO peaks
awk '$5 <= -1' ${subread_dir}/master_atac_peaks_bl_subread_no_header_log2diff.txt > ${subread_dir}/master_atac_peaks_bl_down_KO.txt
awk 'FNR==NR {a[$1]; next} FNR> 1 && $4 in a' ${subread_dir}/master_atac_peaks_bl_down_KO.txt ${macs_dir}/master_atac_summits.bed > ${subread_dir}/master_atac_peaks_bl_down_KO_summit.bed

# save static peaks
awk '$5 < 1 && $5 > -1 ' ${subread_dir}/master_atac_peaks_bl_subread_no_header_log2diff.txt > ${subread_dir}/master_atac_peaks_bl_static.txt
awk 'FNR==NR {a[$1]; next} FNR> 1 && $4 in a' ${subread_dir}/master_atac_peaks_bl_static.txt ${macs_dir}/master_atac_summits.bed > ${subread_dir}/master_atac_peaks_bl_static_summit.bed


######################################################
#### plot heatmap on differential and static regions
######################################################


echo "generating deeptool heatmap - differential and static peaks"

ml macs/2.1.0
start_time=`date +%s`

for file in ${subread_dir}/*_summit.bed;  

  do 
    # file="$(basename $file)"
    # samplename=$(echo "$file" | sed 's/.*\(_bl_\)/\1/g')
    samplename=$(echo "$file" | sed 's/.*bl.//g')

computeMatrix reference-point \
        -R ${file} \
        --skipZeros \
        -S ${deeptools_dir}/SKMel147_ARID2WT_final.bw ${deeptools_dir}/SKMel147_ARID2KO_final.bw \
        -o ${deeptools_dir}/${samplename}.gz \
        -b 5000 -a 5000 \
        --referencePoint center \
        -p 4 \
        --samplesLabel ARID2WT ARID2KO

plotHeatmap -m ${deeptools_dir}/${samplename}.gz \
        --plotFileFormat pdf \
        -out ${deeptools_dir}/${samplename}.pdf \
        --dpi 720 \
        --missingDataColor White \
        --colorMap Reds \
        --regionsLabel ${samplename} \
        --heatmapHeight 13 
done



end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

ml unload macs/2.1.0

##################################################################
#### homer motif analysis - split by differential and static peaks 
##################################################################

echo "motif enrichment analysis using homer, split by differential status"

ml homer/4.10

mkdir -p $homer_dir

preparsedDir=${homer_dir}/preparsed_dir


start_time=`date +%s`


for file in ${subread_dir}/*_summit.bed;  

  do 

    samplename=$(echo "$file" | sed 's/.*bl.//g' | sed 's/_summit.bed//')
    output_dir=${homer_dir}/$samplename
    mkdir -p $output_dir

    bed_path=$file
    
    findMotifsGenome.pl  $bed_path  hg38 $output_dir -preparsedDir $homer_preparsed_dir  -size 200  -p 6 

done



end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


##################################################################
#### homer motif analysis - up/down vs static background
##################################################################

echo "motif enrichment analysis using homer"

ml homer/4.10

mkdir -p $homer_dir

preparsedDir=${homer_dir}/preparsed_dir


start_time=`date +%s`


for file in ${subread_dir}/*_KO_summit.bed;  

  do 

    samplename=$(echo "$file" | sed 's/.*bl.//g' | sed 's/_KO_summit.bed//')
    samplename=${samplename}_vs_static

    output_dir=${homer_dir}/$samplename
    mkdir -p $output_dir

    bed_path=$file
    
    findMotifsGenome.pl  $bed_path  hg38 $output_dir -preparsedDir $homer_preparsed_dir  -size 200  -p 4 -bg ${subread_dir}/master_atac_peaks_bl_static_summit.bed

done



end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


ml unload homer/4.10


##################################################################
#### associate peaks with nearest gene promoter
##################################################################

##Peak annotation
echo "annotate peaks with nearest gene and genomic localization"

ml homer/4.10


for file in ${subread_dir}/*_summit.bed;  

  do 

    samplename=$(echo "$file" | sed 's/.*bl.//g' | sed 's/_summit.bed//')
    
    annotatePeaks.pl  ${file} hg38 -size 200 > ${subread_dir}/${samplename}_annotated.txt
    
done


ml unload homer/4.10



#######################################
#Deeptools heatmap with ChIPseq data and kmeans clustering
########################################

echo "plotting deeptools heatmap with chipseq data"

ml macs/2.1.0
mkdir -p $deeptools_dir
start_time=`date +%s`


computeMatrix reference-point \
    -R ${macs_dir}/master_atac_summits_bl.bed \
    --skipZeros \
    -S ${deeptools_dir}/SKMel147_ARID2WT_final.bw  \
       ${deeptools_dir}/SKMel147_ARID2KO_final.bw  \
       ${raw_dir}/SKmel147_ARID2_final.bw \
       ${raw_dir}/SKmel147_FOSL2_final.bw \
    -o ${deeptools_dir}/atac_summits_bl_w_chipseq.gz \
    -b 5000 -a 5000 \
    --referencePoint center \
    -p 4 

plotHeatmap -m ${deeptools_dir}/atac_summits_bl_w_chipseq.gz \
        --plotFileFormat pdf \
        -out ${deeptools_dir}/atac_summits_bl_w_chipseq.pdf \
        --dpi 720 \
        --missingDataColor White \
        --colorMap Reds \
        --heatmapHeight 13  \
        --samplesLabel ATAC_WT ATAC_KO ARID2 FOSL2 \
        --zMax 1500 1500 500 1000 

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


ml unload macs/2.1.0