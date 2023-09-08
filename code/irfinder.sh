#!/bin/bash
#SBATCH --job-name=irfinder
#SBATCH --output=irfinder.out
#SBATCH --error=irfinder.err
#SBATCH --time=30:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4000


ref=/project2/xczhang/Yalan/reference
IRfinder=/project2/xczhang/Yalan/miniconda3/IRFinder-2.0-beta
fasta=/home/yangyalan/scratch-midway2/flair


#conda activate flair
cd pwd
ln -s $ref/GRCh38.primary_assembly.genome.fa Hg38_2/genome.fa 
ln -s $ref/gencode.v40.annotation.gtf Hg38_2/transcripts.gtf

$IRfinder/bin/IRFinder BuildRefProcess -r Hg38_2 -e $IRfinder/REF/extra-input-files/RNA.SpikeIn.ERCC.fasta.gz -R $IRfinder/REF/extra-input-files/Human_hg38_nonPolyA_ROI.bed -m $IRfinder/REF/Mapabilities/hg38/MapabilityExclusion.70bp.bed.gz

for i in C11_C21_Dividing C1_Excitatory C6_Excitatory CajalRetzius.Cells Progenitor C12_Excitatory C3_Neural C8_Neural Immature.Neuron Neuron C15_Neural C4_Excitatory C9_Excitatory Excitatory.Neuron Inhibitory.Neuron OPC
do
$IRfinder/bin/IRFinder Long -r Hg38 -d $i $fasta/$i\.rename.fa
done

for i in Neuron Progenitor OPC
do
for j in Progenitor OPC
do
$IRfinder/bin/analysisWithNoReplicates.pl -A $i/IRFinder-IR-dir.txt -B $j/IRFinder-IR-dir.txt > $i\-vs-$j\.txt
done
done


for i in CajalRetzius.Cells Excitatory.Neuron Immature.Neuron Inhibitory.Neuron
do
    for j in Excitatory.Neuron Immature.Neuron Inhibitory.Neuron
	do 
	$IRfinder/bin/analysisWithNoReplicates.pl -A $i/IRFinder-IR-dir.txt -B $j/IRFinder-IR-dir.txt > $i\-vs-$j\.txt
done
done

for i in C12_Excitatory C1_Excitatory C6_Excitatory C4_Excitatory C9_Excitatory 
  do
   for j in C1_Excitatory C6_Excitatory C4_Excitatory C9_Excitatory
	do
	$IRfinder/bin/analysisWithNoReplicates.pl -A $i/IRFinder-IR-dir.txt -B $j/IRFinder-IR-dir.txt > $i\-vs-$j\.txt
done
done

for i in C11_C21_Dividing C15_Neural C3_Neural C8_Neural
do
	for j in C15_Neural C3_Neural C8_Neural
	do 
        $IRfinder/bin/analysisWithNoReplicates.pl -A $i/IRFinder-IR-dir.txt -B $j/IRFinder-IR-dir.txt > $i\-vs-$j\.txt
	done
done
