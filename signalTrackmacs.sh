

input2=/data/ltt/projects/LiTing/202105Chip/results/bam/Input20427/Input20427.bam.sorted.rmd.bam
ig2=/data/ltt/projects/LiTing/202105Chip/results/bam/IgG20427/IgG20427.bam.sorted.rmd.bam
ig4=/data/ltt/projects/LiTing/202105Chip/results/bam/IgG4/IgG4.bam.sorted.rmd.bam
input4=/data/ltt/projects/LiTing/202105Chip/results/bam/Input4/Input4.bam.sorted.rmd.bam
input3=/data/ltt/projects/LiTing/202105Chip/results/bam/Input30427/Input30427.bam.sorted.rmd.bam
ig3=/data/ltt/projects/LiTing/202105Chip/results/bam/IgG30427/IgG30427.bam.sorted.rmd.bam
k2=/data/ltt/projects/LiTing/202105Chip/results/bam/KAN20427/KAN20427.bam.sorted.rmd.bam
k3=/data/ltt/projects/LiTing/202105Chip/results/bam/KAN30427/KAN30427.bam.sorted.rmd.bam
k4=/data/ltt/projects/LiTing/202105Chip/results/bam/KAN4/KAN4.bam.sorted.rmd.bam


#cmd="macs2 callpeak -t $k2 -c $input2 -f BAM -g hs -n KANSLvsInput_k2 -B --SPMR --nomodel &"
#echo $cmd
#eval $cmd
#
#cmd="macs2 callpeak -t $k3 -c $input3 -f BAM -g hs -n KANSLvsInput_k3 -B --SPMR --nomodel &"
#echo $cmd
#eval $cmd
#
#cmd="macs2 callpeak -t $k4 -c $input4 -f BAM -g hs -n KANSLvsInput_k4 -B --SPMR --nomodel &"
#echo $cmd
#eval $cmd

#  macs2 bdgcmp -t kansl1vsInput_pileup.bdg -c kansl1vsInput


pileups=(`ls *pileup.bdg`)
lambdas=(`ls *lambda.bdg`)

#for i in 0 1 2
#do
#	echo ${pileups[i]}; echo ${lambdas[i]}
#	cmd="macs2  bdgcmp -t ${pileups[i]} -c ${lambdas[i]} -o kansl1VsInput_${i}_EE.bdg -m FE &"
#	echo $cmd
#	eval "$cmd"
#	echo
#done	

#for i in `ls *EE.bdg`
#do
#	echo $i
#	cmd="sort -k1,1 -k2,2n $i > ${i}.sorted.bdg &"
#	echo $cmd
#	eval "$cmd"
#done

#chromInfo=/data/ltt/ref/HiC/hicpro/doc/homo/chrom_hg19.sizes
chromInfo=/ref/fasta/hg19/hg19.chrom.sizes
#for i in `ls *.sorted.bdg`
#do
#	bedGraphToBigWig $i $chromInfo ${i}.bw &
#done	

#for i in `ls *sorted.bdg.bw`
#do
#	cmd="bigWigToBedGraph $i ${i}.STX17region.bdg -chrom=chr9 -start=100400000 -end=102900000 &"
#	echo $cmd; eval "$cmd"
#done

for i in `ls *STX17region.bdg`
do
	bedGraphToBigWig $i $chromInfo ${i}.bw 
done

