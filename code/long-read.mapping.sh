#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --output=mapping.out
#SBATCH --error=mapping.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4000


ref=/project2/xczhang/Yalan/reference

#conda activate flair
cd /home/yangyalan/scratch-midway2/flair
for i in Immature.Neuron OPC Dividing.Progenitor Cajal–Retzius.Cells Neural.Progenitor Inhibitory.Neuron Excitatory.Neuron
for i in CajalRetzius.Cells
do
#seqkit replace -p '.+' -r $i'_{nr}' ../celltype/$i.fasta > $i\.rename.fa
minimap2 -t 30 -ax splice -uf --secondary=no -C5 /project2/xczhang/Yalan/reference/GRCh38.primary_assembly.genome.fa $i\.rename.fa > $i\.mapped.sam 
samtools view -bS $i\.mapped.sam|samtools sort > $i\.mapped.bam
samtools index $i\.mapped.bam
samtools view -h $i\.mapped.bam | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -bS - > $i\.filter.bam
done
