#!/bin/bash
#SBATCH --job-name=flair2
#SBATCH --output=flair2.out
#SBATCH --error=flair2.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4000


cd /home/yangyalan/scratch-midway2/flair/NvsP_rename
ref=/project2/xczhang/Yalan/reference

#conda activate flair

#for i in Neuron Progenitor
for i in Neuron
do
minimap2 -t 30 -ax splice -uf --secondary=no -C5 /project2/xczhang/Yalan/reference/GRCh38.primary_assembly.genome.fa ../$i\.rename.fa > $i\.mapped.sam 
samtools view -bS $i\.mapped.sam|samtools sort > $i\.mapped.bam
samtools index $i\.mapped.bam
bam2Bed12 -i $i\.mapped.bam > $i\.bed12
done

#cat Neuron.bed12 Progenitor.bed12 |sort|uniq > flair.aligned.bed

#flair correct -q Neuron.bed12 -g $ref/GRCh38.primary_assembly.genome.fa -f $ref/gencode.v40.annotation.gtf --output Neuron 

#flair correct -q Progenitor.bed12 -g $ref/GRCh38.primary_assembly.genome.fa -f $ref/gencode.v40.annotation.gtf --output Progenitor

#cat Neuron_all_corrected.bed Progenitor_all_corrected.bed > flair_all_corrected.bed

#flair collapse -g $ref/GRCh38.primary_assembly.genome.fa -r ../NsvP/Neuron.rename.fasta,../NsvP/Progenitor.rename.fasta -q flair_all_corrected.bed --gtf $ref/gencode.v40.annotation.gtf

#flair quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa 

#flair diffSplice -i flair.collapse.isoforms.bed -q flair.quantify.counts.tsv -o diff
#diff_iso_usage flair.quantify.counts.tsv Neuron_Neuron_batch1 Progentior_Progentior_batch1 neuron_progentior.diff_isos.txt

#diffsplice_fishers_exact diff/diffsplice.es.events.quant.tsv Neuron_Neuron_batch1 Progentior_Progentior_batch1 neuron_progentior.es.fishers.tsv


