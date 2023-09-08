#!/bin/bash
#SBATCH --job-name=flair2
#SBATCH --output=flair2.out
#SBATCH --error=flair2.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4000


cd /home/yangyalan/scratch-midway2/flair/flair_diff
ref=/project2/xczhang/Yalan/reference

#for i in Neuron Progentior OPC
#  do
#    for j in Progentior OPC
#        do
#                diff_iso_usage flair.quantify.counts.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\_diff_isos.txt
#                diffsplice_fishers_exact diffsplice/diffsplice.es.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.es.fishers.tsv
#                diffsplice_fishers_exact diffsplice/diffsplice.alt3.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.alt3.fishers.tsv
#                diffsplice_fishers_exact diffsplice/diffsplice.alt5.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.alt5.fishers.tsv
#                diffsplice_fishers_exact diffsplice/diffsplice.ir.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.ir.fishers.tsv
#        done
done

#for i in CajalRetzius Excitatory  Immature Inhibitor 
#  do
   # for j in Excitatory Immature Inhibitor
	do
		diff_iso_usage flair.quantify.counts.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\_diff_isos.txt
		diffsplice_fishers_exact diffsplice/diffsplice.es.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.es.fishers.tsv
		diffsplice_fishers_exact diffsplice/diffsplice.alt3.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.alt3.fishers.tsv
		diffsplice_fishers_exact diffsplice/diffsplice.alt5.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.alt5.fishers.tsv
		diffsplice_fishers_exact diffsplice/diffsplice.ir.events.quant.tsv $i\_$i\_batch1 $j\_$j\_batch1 $i\_$j\.ir.fishers.tsv
	done
done


