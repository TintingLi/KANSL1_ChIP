computeMatrix reference-point  --referencePoint TSS -a 5000 -b 5000 \
	-R kansg_regular_peaks.narrowPeak.tss.bed -S bigwig/kansg.sorted.rmd.bam.bw -o matrix_kansl1_tss.gz \
	--outFileSortedRegions regions_kansl1.bed &
	
plotHeatmap -m matrix_kansl1_tss.gz -out kansl1.png --heatmapHeight 15
plotHeatmap -m matrix_kansl1_tss.gz -out kansl1.png --heatmapHeight 15
plotHeatmap -m matrix_kansl1_tss.gz -out kansl1.png --heatmapHeight 15 --colorList 'blue,white,red' 
