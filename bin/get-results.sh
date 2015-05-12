#!/bin/bash
#
# Roberto Toro, 12 May 2015
#
# Script to collect results from hsq-sge.sh
#
#

r=(l10icv l10bv l10hip l10th l10ca l10pa l10pu l10amy l10acc height viq piq)
n=(ICV BV Hip Th Ca Pa Pu Amy Acc Height VIQ PIQ)

sz=${#r[@]}

# all snps
total=0;
for ((i=1;i<=22;i++)); do
	snps=$(cat grm-all/log/grm-chr$i.o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
	total=$(($total+$snps))
done
echo "Total SNPs included in analysis: "$total

echo

# perchr
for ((i=1;i<=22;i++)); do
	echo "Total SNPs in chr$i: "$(cat grm-all/log/grm-chr$i.o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
done

echo

# genic/nongenic SNPs
for margin in 0 20 50 ;do
	echo "Genic = Ref. Seq. ± "$margin"Kbp"
	total=0;
	for ((i=1;i<=22;i++)); do
		snps=$(cat grm-genic/log/genic$margin-$i.grm.o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
		total=$(($total+$snps))
	done
	echo "Total genic SNPs: "$total
	total=0;
	for ((i=1;i<=22;i++)); do
		snps=$(cat grm-genic/log/nongenic$margin-$i.grm.o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
		total=$(($total+$snps))
	done
	echo "Total nongenic SNPs: "$total
done

echo

# maf
for maf in "05-20" "20-35" "35-50" ;do
	total=0;
	for ((i=1;i<=22;i++)); do
		snps=$(cat grm-maf/log/maf.$maf.$i.o*|awk 'BEGIN{i=0}{i++;if(i==30)print $10}');
		total=$(($total+$snps))
	done
	echo "Total MAF "$maf" SNPs: "$total
done

echo

# gene lists
for genelist in cnsexpression neurodev; do
	total=0;
	for ((i=1;i<=22;i++)); do
		snps=$(cat grm-$genelist/log/$genelist-$i.o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
		total=$(($total+$snps))
	done
	echo "Total "${genelist}"+ SNPs: "$total
	total=0
	for ((i=1;i<=22;i++)); do
		snps=$(cat "grm-$genelist"/log/no${genelist}"-$i".o*|awk 'BEGIN{i=0}{i++;if(i==21)print $1}');
		total=$(($total+$snps))
	done
	echo "Total "${genelist}"- SNPs: "$total
	echo
done

echo

# table 1: all
echo "Table 1. Estimates of variance explained by genetic factors VG/VP."
echo -e "Phenotype\tVG/VP ± s.e.\tP-value\tN"
for ((i=0;i<$sz;i++)); do
	echo -en ${n[$i]}"\t"
	cat hsq-all/all.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==5)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==9)printf("%.4f\t",$2);if(i==10)printf("%i\n",$2)}'
done

echo

# table 2: no pca
echo "Table 2. Estimates of variance explained by genotyped markers, without control for population stratification."
echo -e "Phenotype\tVG/VP ± s.e.\tP-value\tN"
for ((i=0;i<$sz;i++)); do
	echo -en ${n[$i]}"\t"
	cat hsq-nopca/nopca.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==5)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==9)printf("%.4f\t",$2);if(i==10)printf("%i\n",$2)}'
done

echo

# table 3: genic
echo "Table 3. Estimates of variance explained by genic subsets of genotyped markers, including neighbouring regulatory regions within 0, ±20 and ±50 Kbp."
for j in 0 20 50; do
	echo "Genic regions = Ref. Seq. ± "$j"Kbp"
	echo -e "Phenotype\tVgenic/VP ± s.e.\tVnongenic/VP ± s.e.\tPgenic\tPnongenic\tN"
	for ((i=0;i<$sz;i++)); do
		echo -en ${n[$i]}"\t"
		cat hsq-genic/genic$j.1.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==6)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==7)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==11)printf("%.4f\t",$2)}'
		cat hsq-genic/genic$j.2.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==11)printf("%.4f\t",$2);if(i==12)printf("%i\n",$2)}'
	done
done

echo

# table 4: cnsexpression
echo "Table 4. CNS Expression"
echo -e "Phenotype\tVCNS+/VP ± s.e.\tVCNS-/VP ± s.e.\tVnongenic/VP ± s.e.\tPCNS+\tPCNS-\tPnongenic\tN"
for ((i=0;i<$sz;i++)); do
	echo -en ${n[$i]}"\t"
	
	cat hsq-cnsexpression/cnsexpression.1.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==7)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==8)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==9)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==13)printf("%.4f\t",$2)}'
	cat hsq-cnsexpression/cnsexpression.2.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2)}'
	cat hsq-cnsexpression/cnsexpression.3.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2);if(i==14)printf("%i\n",$2)}'
done

echo

# table 5: neurodev
echo "Table 5. Neurodevelopment"
echo -e "Phenotype\tVneuro+/VP ± s.e.\tVneuro-/VP ± s.e.\tVnongenic/VP ± s.e.\tPneuro+\tPneuro-\tPnongenic\tN"
for ((i=0;i<$sz;i++)); do
	echo -en ${n[$i]}"\t"
	
	cat hsq-neurodev/neurodev.1.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==7)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==8)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==9)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==13)printf("%.4f\t",$2)}'
	cat hsq-neurodev/neurodev.2.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2)}'
	cat hsq-neurodev/neurodev.3.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2);if(i==14)printf("%i\n",$2)}'
done

echo

# table 6: maf
echo "Table 6. Allele frequency"
echo -e "Phenotype\tV5-20/VP ± s.e.\tV20-35/VP ± s.e.\tV35-50/VP ± s.e.\tP5-20\tP20-35\tP35-50\tN"
for ((i=0;i<$sz;i++)); do
	echo -en ${n[$i]}"\t"
	
	cat hsq-maf/maf.1.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==7)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==8)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==9)printf("%.0f ± %.0f\t",$2*100,$3*100);if(i==13)printf("%.4f\t",$2)}'
	cat hsq-maf/maf.2.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2)}'
	cat hsq-maf/maf.3.${r[$i]}.hsq|awk 'BEGIN{i=0}{i++;if(i==13)printf("%.4f\t",$2);if(i==14)printf("%i\n",$2)}'
done
