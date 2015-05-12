#!/bin/bash
#
# Roberto Toro, 12 May 2015
#
# Estimation of phenotypic variance captured by SNPs
# in different genomic partitions
#

# links to command line tools
mygcta=/home1/Ghfc/rto/bin/mygcta
gcta=/home1/Ghfc/rto/bin/gcta
pheno=/home1/Ghfc/rto/imagen/data
imagen=/home1/Ghfc/rto/imagen

# Loop through phenotypes
for x in l10icv l10bv l10hip l10th l10ca l10pa l10pu l10amy l10acc height viq piq; do

if [ true ]; then
# table 1: Using GRMs computed from all SNPs, including 10 PCs as covariates
rm -r $imagen/hsq-all
mkdir -p $imagen/hsq-all/log
qsub	-N			"all.$x" \
		-S			/bin/bash \
		-q			ghfc \
		-e			$imagen/hsq-all/log \
		-o			$imagen/hsq-all/log \
		-cwd \
		<< EOF
$mygcta $gcta	--pheno		$pheno/$x.txt \
				--grm		$imagen/grm-all/grm-all \
				--qcovar	$pheno/pred-age.txt \
				--qcovar	$imagen/grm-all/grm-all-0.025.eigenvec \
				--covar		$pheno/centre.txt \
				--covar		$pheno/sex.txt \
				--keep		$imagen/hsq/"keep.txt" \
				--out		$imagen/hsq-all/all.$x \
				--reml
EOF
fi

if [ true ]; then
# table 2: Using GRMs computed from all SNPs, not including 10 PCs as covariates
rm -r $imagen/hsq-nopca
mkdir -p $imagen/hsq-nopca/log
qsub	-N			"nopca.$x" \
		-S			/bin/bash \
		-q			ghfc \
		-e			$imagen/hsq-nopca/log \
		-o			$imagen/hsq-nopca/log \
		-cwd \
		<< EOF
$mygcta $gcta	--pheno		$pheno/$x.txt \
				--grm		$imagen/grm-all/grm-all \
				--qcovar	$pheno/pred-age.txt \
				--covar		$pheno/centre.txt \
				--covar		$pheno/sex.txt \
				--keep		$imagen/hsq/"keep.txt" \
				--out		$imagen/hsq-nopca/nopca.$x \
				--reml
EOF
fi

if [ true ]; then
# table 3: Using GRMs computed from genic and nongenic SNP partitions
	rm -r $imagen/hsq-genic
	mkdir $imagen/hsq-genic
	mkdir $imagen/hsq-genic/log
	for y in 0 20 50; do
		echo $imagen/grm-genic/genic$y-0.025 >> $imagen/hsq-genic/genic$y.test.txt
		echo $imagen/grm-genic/nongenic$y-0.025 >> $imagen/hsq-genic/genic$y.test.txt
		for z in 1 2; do
			qsub	-N			"genic$y.$z.$x" \
					-S			/bin/bash \
					-q			ghfc \
					-e			$imagen/hsq-genic/log \
					-o			$imagen/hsq-genic/log \
					-cwd \
					<< EOF
			$mygcta $gcta	--pheno		$pheno/$x.txt \
							--mgrm		$imagen/hsq-genic/genic$y.test.txt \
							--qcovar	$pheno/pred-age.txt \
							--qcovar	$imagen/grm-all/grm-all-0.025.eigenvec \
							--covar		$pheno/centre.txt \
							--covar		$pheno/sex.txt \
							--keep		$imagen/hsq/"keep.txt" \
							--out		$imagen/hsq-genic/genic$y.$z.$x \
							--reml \
							--reml-lrt	$z
EOF
		done
	done
fi

if [ true ]; then
# table 4: Using GRMs computed from genic SNPs preferentially expressed
# in the CNS, versus remaining genic SNPs, versus nongenic SNPs
	rm -r $imagen/hsq-cnsexpression
	mkdir -p $imagen/hsq-cnsexpression/log
	echo $imagen/grm-cnsexpression/cnsexpression-0.025 >> $imagen/hsq-cnsexpression/cnsexpression.test.txt
	echo $imagen/grm-cnsexpression/nocnsexpression-0.025 >> $imagen/hsq-cnsexpression/cnsexpression.test.txt
	echo $imagen/grm-genic/nongenic50-0.025 >> $imagen/hsq-cnsexpression/cnsexpression.test.txt
	for z in 1 2 3; do
		qsub	-N			"cnsexpression.$z.$x" \
				-S			/bin/bash \
				-q			ghfc \
				-e			$imagen/hsq-cnsexpression/log \
				-o			$imagen/hsq-cnsexpression/log \
				-cwd \
				<< EOF
		$mygcta $gcta	--pheno		$pheno/$x.txt \
						--mgrm		$imagen/hsq-cnsexpression/cnsexpression.test.txt \
						--qcovar	$pheno/pred-age.txt \
						--qcovar	$imagen/grm-all/grm-all-0.025.eigenvec \
						--covar		$pheno/centre.txt \
						--covar		$pheno/sex.txt \
						--keep		$imagen/hsq/"keep.txt" \
						--out		$imagen/hsq-cnsexpression/cnsexpression.$z.$x \
						--reml-lrt	$z \
						--reml
EOF
	done
fi

if [ true ]; then
# table 5: Using GRMs computed from genic SNPs involved in neurodev, versus remaining
# genic SNPs, versus nongenic SNPs
	rm -r $imagen/hsq-neurodev
	mkdir -p $imagen/hsq-neurodev/log
	echo $imagen/grm-neurodev/neurodev-0.025 >> $imagen/hsq-neurodev/neurodev.test.txt
	echo $imagen/grm-neurodev/noneurodev-0.025 >> $imagen/hsq-neurodev/neurodev.test.txt
	echo $imagen/grm-genic/nongenic50-0.025 >> $imagen/hsq-neurodev/neurodev.test.txt
	for z in 1 2 3; do
	qsub	-N			"neurodev.$z.$x" \
			-S			/bin/bash \
			-q			ghfc \
			-e			$imagen/hsq-neurodev/log \
			-o			$imagen/hsq-neurodev/log \
			-cwd \
			<< EOF
	$mygcta $gcta	--pheno		$pheno/$x.txt \
					--mgrm		$imagen/hsq-neurodev/neurodev.test.txt \
					--qcovar	$pheno/pred-age.txt \
					--qcovar	$imagen/grm-all/grm-all-0.025.eigenvec \
					--covar		$pheno/centre.txt \
					--covar		$pheno/sex.txt \
					--keep		$imagen/hsq/"keep.txt" \
					--out		$imagen/hsq-neurodev/neurodev.$z.$x \
					--reml-lrt $z \
					--reml
EOF
	done
fi

if [ true ]; then
# table 6: Using partitions based on MAF
rm -r $imagen/hsq-maf
mkdir -p $imagen/hsq-maf/log
echo $imagen/grm-maf/maf.05-20.0.025 >> $imagen/hsq-maf/maf.test.txt
echo $imagen/grm-maf/maf.20-35.0.025 >> $imagen/hsq-maf/maf.test.txt
echo $imagen/grm-maf/maf.35-50.0.025 >> $imagen/hsq-maf/maf.test.txt
for z in 1 2 3; do
qsub	-N			"maf.$z.$x" \
		-S			/bin/bash \
		-q			ghfc \
		-e			$imagen/hsq-maf/log \
		-o			$imagen/hsq-maf/log \
		-cwd \
		<< EOF
$mygcta $gcta	--pheno		$pheno/$x.txt \
				--mgrm		$imagen/hsq-maf/maf.test.txt \
				--qcovar	$pheno/pred-age.txt \
				--qcovar	$imagen/grm-all/grm-all-0.025.eigenvec \
				--covar		$pheno/centre.txt \
				--covar		$pheno/sex.txt \
				--keep		$imagen/hsq/"keep.txt" \
				--out		$imagen/hsq-maf/maf.$z.$x \
				--reml-lrt	$z \
				--reml
EOF
done
fi

if [ true ]; then
# table 7: Using one GRM per chromosome
rm -r $imagen/hsq-perchr
mkdir $imagen/hsq-perchr
mkdir $imagen/hsq-perchr/log
for ((i=1;i<=22;i++)); do
qsub	-N			"chr$i.$x" \
		-S			/bin/bash \
		-q			ghfc \
		-e			$imagen/hsq-perchr/log \
		-o			$imagen/hsq-perchr/log \
		-cwd \
		<< EOF
$mygcta $gcta	--pheno		$pheno/$x.txt \
				--grm		$imagen/perchr/grm-chr$i-0.025 \
				--qcovar	$pheno/pred-age.txt \
				--covar		$pheno/centre.txt \
				--covar		$pheno/sex.txt \
				--keep		$imagen/hsq/"keep.txt" \
				--out		$imagen/hsq-perchr/chr$i.$x \
				--reml
EOF
done
for((i=1;i<=22;i++)); do echo $imagen/grm-all/grm-chr$i-0.025 >> $imagen/hsq-perchr/perchr.test.txt; done
qsub	-N			"allchr.$x" \
		-S			/bin/bash \
		-q			ghfc \
		-e			$imagen/hsq-perchr/log \
		-o			$imagen/hsq-perchr/log \
		-cwd \
		<< EOF
$mygcta $gcta	--pheno		$pheno/$x.txt \
				--mgrm		$imagen/hsq-perchr/perchr.test.txt \
				--qcovar	$pheno/pred-age.txt \
				--covar		$pheno/centre.txt \
				--covar		$pheno/sex.txt \
				--keep		$imagen/hsq/"keep.txt" \
				--out		$imagen/hsq-perchr/allchr.$x \
				--reml
EOF
fi

done
