#!/bin/bash
#
# Roberto Toro, 12 May 2015
#
# Script to compute GRMs
#
#

# Path to command line tols
gcta=/home1/Ghfc/rto/bin/gcta
genes2snps=/home1/Ghfc/rto/bin/genes2snps/genes2snps
plink=/home1/Ghfc/rto/bin/plink

# Directory to store results
imagen=/home1/Ghfc/rto/imagen

# Data directories
bed=/home1/Ghfc/rto/imagen/geno/all-pruned
hg18genes=/home1/Ghfc/rto/bin/genes2snps/hg18genes.txt
hg18snps=/home1/Ghfc/rto/bin/genes2snps/hg18snp

if [ true ]; then
# compute a GRM using all (filtered, R2 pruned) SNPs
#-------------------------------------------------------
	# Remove results from a previous run
	rm -r $imagen/grm-all
	
	# make directories for new results
	mkdir -p $imagen/grm-all/log
	
	# Make grm in parallel for 22 chromosomes
	for((i=1;i<=22;i++)); do
		qsub	-N			"grm-chr$i"\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-all/log \
				-o			$imagen/grm-all/log \
				<< EOF
		$gcta	--bfile		$bed\
				--chr		$i\
				--make-grm\
				--out		"$imagen/grm-all/grm-chr$i"
EOF
	done
	
	# Combine the 22 GRMs into one single GRM
	for((i=1;i<=22;i++)); do echo $imagen/grm-all/grm-chr$i >> $imagen/grm-all/all.multi.txt; done
	dep=$(for((i=1;i<=22;i++)); do echo -n "grm-chr$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
	qsub	-N			"grm-all"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-all/log \
			-o			$imagen/grm-all/log \
			-hold_jid	$dep \
			<< EOF
	$gcta	--mgrm	$imagen/grm-all/all.multi.txt \
			--make-grm \
			--out	"$imagen/grm-all/grm-all"
EOF
	
	# Filter for genetic relationship <0.025
	qsub	-N			"grm-all-0.025"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-all/log \
			-o			$imagen/grm-all/log \
			-hold_jid	"grm-all" \
			<< EOF
	$gcta	--grm			$imagen/grm-all/grm-all\
			--grm-cutoff	0.025 \
			--make-grm\
			--out			"$imagen/grm-all/grm-all-0.025"
EOF
	
	# Compute 10 principal components
	qsub	-N			"grm-all-0.025.pca"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-all/log \
			-o			$imagen/grm-all/log \
			-hold_jid	"grm-all-0.025" \
			<< EOF
	$gcta	--grm			$imagen/grm-all/grm-all-0.025\
			--pca			10 \
			--out			"$imagen/grm-all/grm-all-0.025"
EOF
fi

if [ true ]; then
# Compute 1 GRM per chromosome
#-------------------------------------------------------
	rm -r $imagen/grm-perchr
	mkdir $imagen/grm-perchr
	mkdir $imagen/grm-perchr/log
	grm025id=$imagen/grm-all/grm-all-0.025.grm.id
	for((i=1;i<=22;i++)); do
		qsub	-N			"grm-chr$i-0.025"\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-perchr/log \
				-o			$imagen/grm-perchr/log \
				-hold_jid	"grm-all-0.025"\
				<< EOF
		$gcta	--grm	$imagen/grm-all/grm-chr$i \
				--keep	$grm025id \
				--make-grm\
				--out	"$imagen/grm-perchr/grm-chr$i-0.025"
EOF
	done
fi

if [ true ]; then
# Compute GRMs from genic and nongenic SNPs
#-------------------------------------------------------	
	# remove results from a previous run
	rm -r $imagen/grm-genic
	
	# make directories for new results
	mkdir -p $imagen/grm-genic/log
	
	# Make a gene list
	cut -d' ' -f 4 $hg18genes|sort -u > $imagen/grm-genic/genic.txt
	
	# 3 definitions of genic: strict, +/- 20kbp, +/- 50kbp
	for margin in 0 20 50; do
		# make a list of SNPs, .snplist
		$genes2snps	$imagen/grm-genic/genic.txt\
					$hg18snps"$margin"Kbp.set\
					> $imagen/grm-genic/genic"$margin".snplist
		# make .bed files for genic and nongenic SNPs
		qsub	-N	genic"$margin".bed\
				-S	/bin/bash\
				-q	ghfc\
				-e	$imagen/grm-genic/log \
				-o	$imagen/grm-genic/log \
				<< EOF
		$plink	--bfile		$bed\
				--extract	$imagen/grm-genic/genic"$margin".snplist\
				--make-bed\
				--noweb\
				--out		$imagen/grm-genic/genic"$margin"
EOF
		qsub	-N	nongenic"$margin".bed\
				-S	/bin/bash\
				-q	ghfc\
				-e	$imagen/grm-genic/log \
				-o	$imagen/grm-genic/log \
				<< EOF
		$plink	--bfile		$bed\
				--exclude	$imagen/grm-genic/genic"$margin".snplist\
				--make-bed\
				--noweb\
				--out		$imagen/grm-genic/nongenic"$margin"
EOF
		
		# make genic GRM
		for((i=1;i<=22;i++)); do echo $imagen/grm-genic/"genic$margin-$i" >> $imagen/grm-genic/genic"$margin".multi.txt; done
		for((i=1;i<=22;i++)); do
			qsub	-N			"genic$margin-$i.grm"\
					-S			/bin/bash\
					-q			ghfc\
					-hold_jid	genic"$margin".bed\
					-e			$imagen/grm-genic/log \
					-o			$imagen/grm-genic/log \
					<< EOF
			$gcta	--bfile	$imagen/grm-genic/genic"$margin"\
					--chr	$i\
					--make-grm\
					--out	$imagen/grm-genic/"genic$margin-$i"
EOF
		done
		dep=$(for((i=1;i<=22;i++)); do echo -n "genic$margin-$i.grm"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
		qsub	-N			genic"$margin".grm\
				-S			/bin/bash\
				-q			ghfc\
				-hold_jid	$dep\
				-e			$imagen/grm-genic/log \
				-o			$imagen/grm-genic/log \
				<< EOF
		$gcta	--mgrm	$imagen/grm-genic/genic"$margin".multi.txt\
				--make-grm\
				--out	$imagen/grm-genic/genic"$margin"
EOF
		qsub	-N			"genic$margin-0.025"\
				-S			/bin/bash\
				-q			ghfc\
				-hold_jid	"genic$margin.grm","grm-all-0.025"\
				-e			$imagen/grm-genic/log \
				-o			$imagen/grm-genic/log \
				<< EOF
		$gcta	--grm	$imagen/grm-genic/genic"$margin"\
				--keep	$grm025id\
				--make-grm\
				--out	$imagen/grm-genic/"genic$margin-0.025"
EOF
		
		# make nongenic GRM
		for((i=1;i<=22;i++)); do echo $imagen/grm-genic/"nongenic$margin-$i" >> $imagen/grm-genic/nongenic"$margin".multi.txt; done
		for((i=1;i<=22;i++)); do
			qsub	-N			"nongenic$margin-$i.grm"\
					-S			/bin/bash\
					-q			ghfc\
					-hold_jid	nongenic"$margin".bed\
					-e			$imagen/grm-genic/log \
					-o			$imagen/grm-genic/log \
					<< EOF
			$gcta	--bfile	$imagen/grm-genic/nongenic"$margin"\
					--chr	$i\
					--make-grm\
					--out	$imagen/grm-genic/"nongenic$margin-$i"
EOF
		done
		dep=$(for((i=1;i<=22;i++)); do echo -n "nongenic$margin-$i.grm"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
		qsub	-N			nongenic"$margin".grm\
				-S			/bin/bash\
				-q			ghfc\
				-hold_jid	$dep\
				-e			$imagen/grm-genic/log \
				-o			$imagen/grm-genic/log \
				<< EOF
		$gcta	--mgrm	$imagen/grm-genic/nongenic"$margin".multi.txt\
				--make-grm\
				--out	$imagen/grm-genic/nongenic"$margin"
EOF
		qsub	-N			"nongenic$margin-0.025"\
				-S			/bin/bash\
				-q			ghfc\
				-hold_jid	nongenic"$margin".grm,"grm-all-0.025"\
				-e			$imagen/grm-genic/log \
				-o			$imagen/grm-genic/log \
				<< EOF
		$gcta	--grm	$imagen/grm-genic/nongenic"$margin"\
				--keep	$grm025id\
				--make-grm\
				--out	$imagen/grm-genic/"nongenic$margin-0.025"
EOF
	done
fi

if [ true ]; then
# Make GRMs for low, medium and high MAF
#-------------------------------------------------------
	rm -r $imagen/grm-maf
	mkdir -p $imagen/grm-maf/log
	# make 1 grm per chr, for the 3 maf levels
	for((i=1;i<=22;i++)); do
		echo $imagen/grm-maf/"maf.05-20.$i" >> $imagen/grm-maf/"maf.05-20.multi.txt"
		echo $imagen/grm-maf/"maf.20-35.$i" >> $imagen/grm-maf/"maf.20-35.multi.txt"
		echo $imagen/grm-maf/"maf.35-50.$i" >> $imagen/grm-maf/"maf.35-50.multi.txt"
	done
	for((i=1;i<=22;i++)); do
	qsub	-N			"maf.05-20.$i"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			<< EOF
	$gcta	--bfile		$bed\
			--maf		0.05\
			--max-maf	0.20\
			--chr		$i\
			--make-grm\
			--out		$imagen/grm-maf/"maf.05-20.$i"
EOF
	qsub	-N			"maf.20-35.$i"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			<< EOF
	$gcta	--bfile		$bed\
			--maf		0.20\
			--max-maf	0.35\
			--chr		$i\
			--make-grm\
			--out		$imagen/grm-maf/"maf.20-35.$i"
EOF
	qsub	-N			"maf.35-50.$i"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			<< EOF
	$gcta	--bfile		$bed\
			--maf		0.35\
			--max-maf	0.50\
			--chr		$i\
			--make-grm\
			--out		$imagen/grm-maf/"maf.35-50.$i"
EOF
	done
	
	# mix the chromosomic grm files into a single one, for the 3 maf levels
	dep=$(for((i=1;i<=22;i++)); do echo -n "maf.05-20.$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
	qsub	-N			"maf.05-20"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	$dep\
			<< EOF
	$gcta	--mgrm		$imagen/grm-maf/"maf.05-20.multi.txt"\
			--make-grm\
			--out		$imagen/grm-maf/"maf.05-20"
EOF
	dep=$(for((i=1;i<=22;i++)); do echo -n "maf.20-35.$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
	qsub	-N			"maf.20-35"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	$dep\
			<< EOF
	$gcta	--mgrm		$imagen/grm-maf/"maf.20-35.multi.txt"\
			--make-grm\
			--out		$imagen/grm-maf/"maf.20-35"
EOF
	dep=$(for((i=1;i<=22;i++)); do echo -n "maf.35-50.$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
	qsub	-N			"maf.35-50"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	$dep\
			<< EOF
	$gcta	--mgrm		$imagen/grm-maf/"maf.35-50.multi.txt"\
			--make-grm\
			--out		$imagen/grm-maf/"maf.35-50"
EOF
	
	# filter subjects closer than 0.025, for the 3 maf levels
	qsub	-N			"maf.05-20.0.025"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	"maf.05-20","grm-all-0.025"\
			<< EOF
	$gcta	--grm		$imagen/grm-maf/"maf.05-20"\
			--keep		$grm025id\
			--make-grm\
			--out		$imagen/grm-maf/"maf.05-20.0.025"
EOF
	qsub	-N			"maf.20-35.0.025"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	"maf.20-35","grm-all-0.025"\
			<< EOF
	$gcta	--grm		$imagen/grm-maf/"maf.20-35"\
			--keep		$grm025id\
			--make-grm\
			--out		$imagen/grm-maf/"maf.20-35.0.025"
EOF
	qsub	-N			"maf.35-50.0.025"\
			-S			/bin/bash\
			-q			ghfc\
			-e			$imagen/grm-maf/log \
			-o			$imagen/grm-maf/log \
			-hold_jid	"maf.35-50","grm-all-0.025"\
			<< EOF
	$gcta	--grm		$imagen/grm-maf/"maf.35-50"\
			--keep		$grm025id\
			--make-grm\
			--out		$imagen/grm-maf/maf.35-50.0.025
EOF
fi

if [ true ]; then
# Make GRMs for different functional gene partitions
# Currently: genes preferentially expressed in CNS and
# genes involved in neurodevelopment.
#-------------------------------------------------------	
	for A in cnsexpression neurodev; do
		rm -r $imagen/grm-$A
		mkdir -p $imagen/grm-$A/log
		
		B=no$A
		
		# make A.snplist
		$genes2snps	$imagen/grm/$A.txt\
					"$hg18snps"50Kbp.set\
					> $imagen/grm-$A/$A.snplist
		
		# make A, B .bed
		qsub	-N			$A.bed\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				<< EOF
		$plink	--bfile		$bed\
				--extract	$imagen/grm-$A/$A.snplist\
				--make-bed\
				--noweb\
				--out		$imagen/grm-$A/$A
EOF
		
		qsub	-N			$B.bed\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				-hold_jid	"genic50.bed" \
				<< EOF
		$plink	--bfile		$imagen/grm-genic/genic50 \
				--exclude	$imagen/grm-$A/$A.snplist\
				--make-bed\
				--noweb\
				--out		$imagen/grm-$A/$B
EOF
		
		# make A grm
		for((i=1;i<=22;i++)); do
			echo $imagen/grm-$A/"$A-$i" >> $imagen/grm-$A/"$A".multi.txt
		done
		for((i=1;i<=22;i++)); do
			qsub	-N			"$A-$i"\
					-S			/bin/bash\
					-q			ghfc\
					-e			$imagen/grm-$A/log \
					-o			$imagen/grm-$A/log \
					-hold_jid	"$A".bed\
					<< EOF
			$gcta	--bfile		$imagen/grm-$A/$A\
					--chr		$i\
					--make-grm\
					--out		$imagen/grm-$A/"$A-$i"
EOF
		done
		dep=$(for((i=1;i<=22;i++)); do echo -n "$A-$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
		qsub	-N			$A\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				-hold_jid	$dep\
				<< EOF
		$gcta	--mgrm		$imagen/grm-$A/"$A".multi.txt\
				--make-grm\
				--out		$imagen/grm-$A/$A
EOF
		qsub	-N			"$A-0.025"\
				-S			/bin/bash\
				-q			ghfc\
				-hold_jid	$A,"grm-all-0.025"\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				<< EOF
		$gcta	--grm		$imagen/grm-$A/$A\
				--keep		$grm025id\
				--make-grm\
				--out		$imagen/grm-$A/"$A-0.025"
EOF
		
		# make B grm
		for((i=1;i<=22;i++)); do
			echo $imagen/grm-$A/"$B-$i" >> $imagen/grm-$A/"$B".multi.txt
		done
		for((i=1;i<=22;i++)); do
			qsub	-N			"$B-$i"\
					-S			/bin/bash\
					-q			ghfc\
					-e			$imagen/grm-$A/log \
					-o			$imagen/grm-$A/log \
					-hold_jid	"$B".bed\
					<< EOF
			$gcta	--bfile		$imagen/grm-$A/$B\
					--chr		$i\
					--make-grm\
					--out		$imagen/grm-$A/"$B-$i"
EOF
		done
		dep=$(for((i=1;i<=22;i++)); do echo -n "$B-$i"; if [ $i -lt 22 ]; then echo -n ",";fi;done)
		qsub	-N			$B\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				-hold_jid	$dep\
				<< EOF
		$gcta	--mgrm		$imagen/grm-$A/"$B".multi.txt\
				--make-grm\
				--out		$imagen/grm-$A/$B
EOF
		qsub	-N			"$B-0.025"\
				-S			/bin/bash\
				-q			ghfc\
				-e			$imagen/grm-$A/log \
				-o			$imagen/grm-$A/log \
				-hold_jid	$B,"grm-all-0.025"\
				<< EOF
		$gcta	--grm		$imagen/grm-$A/$B\
				--keep		$grm025id\
				--make-grm\
				--out		$imagen/grm-$A/"$B-0.025"
EOF
	done
fi
