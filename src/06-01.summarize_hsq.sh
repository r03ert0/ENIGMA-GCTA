#!/usr/bin/env bash

set -e

shopt -s extglob

#Cut prefix - up to first prefix wildcard match: ${variable#prefix}
#Cut prefix - up to last prefix wildcard match: ${variable##prefix}
#Cut suffix - up to first suffix wildcard match: ${variable%suffix}
#Cut suffix - up to last suffix wildcard match: ${variable%%suffix}

## first arg is the directory containing the GCTA heritability results
dir=$1
## second arg is the directory containing the GRM
dirgrm=$2
## second arg contains prefix of output file 
out=$3


### === extract number of SNPs used to compute each GRM ===
#nbinfiles = `find $dirgrm -name *N.bin`


### === extract univariate heritability results

# get all univariate heritability estimates
if true; then

    groups=("hsq-all" "hsq-nopca" "hsq-perchr") #"hsq-nofilter"
    
    for group in "${groups[@]}"; do
        echo "$group"

        outfile=${out}_${group}.txt
        printf "directory\tsnpgroup\tPhenotype\tVG/VP ± s.e.\tnb_samples\n" > "$outfile" #\tLRT\tdf\tP-value

        for file in "${dir}"/"${group}"/*.hsq; do
            fn=${file##*/}
            fn=${fn%.hsq}
            groupan=${fn%.*}
            pheno=${fn##*.}
        printf "%s\t%s\t%s\t" "${group}" "${groupan}" "${pheno}" >> "$outfile"
        if [[ "$groupan" != "allchr" ]]; then
                awk 'BEGIN{i=0}{i++;if(i==5) {printf("%.2f ± %.2f\t",$2*100,$3*100)} else if (i == 11) printf(" %.0f\n",$2)}' "$file" >> "$outfile"
        else

            awk 'BEGIN{i=0}{i++;if(i==49) {printf("%.2f ± %.2f\t",$4*100,$5*100)} else if (i == 51) printf(" %.0f\n",$2)}' "$file" >> "$outfile"
        fi

        done
    done

fi

if true; then

    group="hsq-maf"
    outfile=${out}_$group.txt
    echo $group
    groupname=${group#hsq-}
    printf "directory\tsnpgroup\tPhenotype\tV%s.0.001-0.05/VP ± s.e.\tV%s.0.05-0.2/VP ± s.e.\tV%s.0.2-0.35/VP ± s.e.\tV%s.0.35-0.5/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" $groupname $groupname $groupname $groupname> "$outfile" #\tLRT5-20\tLRT20-35\tLRT35-50\tN"

    for file in "${dir}"/"$group"/*.hsq; do
    
        fn=${file##*/}
        fn=${fn%.hsq}
        groupan=${fn%.*}
        pheno=${fn##*.}
            printf "hsq-genic/%s\t%s\t%s\t" ${group} "${groupan}" "${pheno}" >> "$outfile"
            awk 'BEGIN{i=0}{i++;if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==10)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==13)printf("%.2f ± %.2f\t",$4*100,$5*100);} END{printf("%.0f\n",$2)};' "$file" >> "$outfile"

    done


    groups=("hsq-cnsexpression" "hsq-neurodev")
        
    for group in "${groups[@]}"; do
        outfile=${out}_${group}.txt
        echo "$group"
        groupname=${group#hsq-}
        if [[ $group = "hsq-maf" ]]; then
            printf "directory\tsnpgroup\tPhenotype\tV%s.0.001-0.05/VP ± s.e.\tV%s.0.05-0.2/VP ± s.e.\tV%s.0.2-0.35/VP ± s.e.\tV%s.0.35-0.5/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" "$groupname" "$groupname" "$groupname" "$groupname"> "$outfile" #\tLRT5-20\tLRT20-35\tLRT35-50\tN"
        else
            printf "directory\tsnpgroup\tPhenotype\tV%s/VP ± s.e.\tVnon%s/VP ± s.e.\tVnongenic/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" "$groupname" "$groupname" > "$outfile" #\tLRT${group}+\tLRT${group}-\${nongenic}\tN
        fi

        for file in "$dir"/"$group"*/*.hsq; do

            fn=${file##*/}
            fn=${fn%.hsq}
            groupan=${fn%.*}
            pheno=${fn##*.}
            printf "%s\t%s\t%s\t" "${group}" "${groupan}" "${pheno}" >> "$outfile"
            #awk 'BEGIN{i=0}{i++;if(i==7)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\n",$4*100,$5*100)};' $file >> $outfile
            awk 'BEGIN{i=0}{i++;if(i==7)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\t",$4*100,$5*100);} END{printf("%.0f\n",$2)};' "$file" >> "$outfile"
            #awk 'BEGIN{i=0}{i++;if(i==10)printf("%.4f\t",$2)}'  $file >> $outfile
            #awk 'BEGIN{i=0}{i++;if(i==10)printf("%.4f\t",$2);if(i==11)printf("%i\n",$2)}'  $file >> $outfile
            #echo '' >> $outfile
        done
    done


    group="updown-margin"
    echo $group
    tutu=("${dirgrm}"/grm-genic/updown-margin[0-9][0-9])
    groups=("${tutu[@]##*/}")
    #echo "${groups[0]}"

    outfile=${out}_hsq-genic_${group}.txt
    printf "directory\tsnpgroup\tPhenotype\tVgenic-margin0/VP ± s.e.\tVupdown-margin/VP ± s.e.\tVnongenic/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" > "$outfile"

    for group in "${groups[@]}"; do
        #echo $group
        groupname=$group

        for file in "${dir}"/hsq-genic/"${group}"*.hsq; do
        
            fn=${file##*/}
            fn=${fn%.hsq}
            groupan=${fn%.*}
            pheno=${fn##*.}

            printf "hsq-genic/%s\t%s\t%s\t" "$group" "$groupan" "$pheno" >> "$outfile"
            awk 'BEGIN{i=0}{i++;if(i==7)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\t",$4*100,$5*100);} END{printf("%.0f\n",$2)};' "$file" >> "$outfile"

        done
    done


    group="updown-margin-2bins"
    tutu=("${dirgrm}"/grm-genic/updown-margin[0-9][0-9]-[0-9][0-9])
    groups=("${tutu[@]##*/}")
    outfile=${out}_hsq-genic_${group}.txt
    printf "directory\tsnpgroup\tPhenotype\tVgenic-margin0/VP ± s.e.\tVupdown-margin-bin1/VP ± s.e.\tVupdown-margin-bin2/VP ± s.e.\tVnongenic/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" > "$outfile"

    for group in "${groups[@]}"; do
        #echo $group
        groupname=$group

        for file in "$dir"/hsq-genic/"$group"*.hsq; do
        
            fn=${file##*/}
            fn=${fn%.hsq}
            groupan=${fn%.*}
            pheno=${fn##*.}

            printf "hsq-genic/%s\t%s\t%s\t" "$group" "$groupan" "$pheno" >> "$outfile"
            awk 'BEGIN{i=0}{i++;if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==10)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\t",$2*100,$3*100); if(i==13)printf("%.2f ± %.2f\t",$4*100,$5*100);;} END{printf("%.0f\n",$2)};' "$file" >> "$outfile"

        done
    done


    # group=("genic-exon1intron1")
    # tutu=`for path in   ${dirgrm}/grm-genic/${group}; do path=${path##*/}; echo "${path%.*}"; done | sort -u`
    # groups=($(echo "$tutu" | tr ' ' '\n')) # string to array
    # outfile=${out}_hsq-genic_${group}.txt
    # rm $outfile
    # printf "directory\tsnpgroup\tPhenotype\tVgenic-exon1intron1/VP ± s.e.\tVgenic-notexon1intron1/VP ± s.e.\tVupdown-margin20/VP ± s.e.\tVnongenic-margin20/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" >> $outfile

    # for group in "${groups[@]}"; do
    #     #echo $group
    #     groupname=`basename $group | sed 's/hsq-//' | sed 's/\\*//'`

    #     for file in ${dir}/hsq-genic/${group}*.hsq; do

    #         fn=${file##*/}
    #         fn=${fn%.hsq}
    #         groupan=${fn%.*}
    #         pheno=${fn##*.}

    #         printf "hsq-genic/${group}\t${groupan}\t${pheno}\t" >> $outfile
    #         awk 'BEGIN{i=0}{i++;if(i==8)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==10)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==11)printf("%.2f ± %.2f\t",$2*100,$3*100); if(i==13)printf("%.2f ± %.2f\t",$4*100,$5*100);;} END{printf("%.0f\n",$2)};' $file >> $outfile

    #     done
    # done


    groups=("hsq-genic")
    
    for group in "${groups[@]}"; do
        outfile=${out}_${group}.txt
        echo "$group"
        groupname=${group#hsq-}
        printf "directory\tsnpgroup\tPhenotype\tV%s/VP ± s.e.\tVnon%s/VP ± s.e.\tSum_of_V(G)/Vp\tnb_samples\n" "${groupname}" "${groupname}" > "$outfile" #\tLRT${group}+\tLRT${group}-
        
        for file in "$dir"/"$group"/genic*.hsq; do #$(find $dir/ -maxdepth 1 -mindepth 1 -name "${snpgroup}*.hsq"); do egrep '^[a-z0-9_]+\.[a-z0-9_]+\.hsq'
            fn=${file##*/}
            fn=${fn%.hsq}
            groupan=${fn%.*}
            pheno=${fn##*.}

            printf "%s\t%s\t%s\t" "$group" "$groupan" "$pheno" >> "$outfile"
            #awk 'BEGIN{i=0}{i++;if(i==6)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==7)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==15)printf("%.0f\n", $2)}' $file >> $outfile #if(i==9)printf("%.2f ± %.2f",$2*100,$3*100);}
        awk 'BEGIN{i=0}{i++;if(i==6)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==7)printf("%.2f ± %.2f\t",$2*100,$3*100);if(i==9)printf("%.2f ± %.2f\t",$4*100,$5*100);} END{printf("%.0f\n", $2)}' "$file" >> "$outfile"
            #awk 'BEGIN{i=0}{i++;if(i==10)printf("%.4f\t",$2)}'  $file >> $outfile
            #awk 'BEGIN{i=0}{i++;if(i==10)printf("%.4f\t",$2);if(i==11)printf("%i\n",$2)}'  $file >> $outfile
            #echo '' >> $outfile

        done
    done



fi
 
       



### === extract bivariate heritability results
if true; then


    dirbiv=$dir/hsq-biv

    outfile=${out}_hsq-biv.txt
    #outfile=hsq_bivariate_out_3011.txt

    snpgroups=("all.rg=0" "all.rg=1" "")


    printf "snpgroup\tphenotype1\tphenotype2\tV(G)/Vp_phenotype1 ± s.e.\tV(G)/Vp_phenotype2 ± s.e.\trG ± s.e.\tpvalue\tnb_samples\n" > "$outfile" #\tLRT\tdf\tP-value\tN
    # printf "${snpgroup}\tV(G)/Vp_${pheno1} ± s.e.\tV(G)/Vp_${pheno2} ± s.e.\t" >> $outfile

    for snpgroup in "${snpgroups[@]}"   #"$@"
    do

        if [[ "${snpgroup}" = '' ]]; then
            for file in "$dirbiv"/!(*.*).!(*.*).hsq;
            do
                if [[ ! -f $file ]]; then
                    continue
                fi
                fn=${file##*/}
                fn=${fn%.hsq}
                phenos=$fn
                pheno1=${phenos%.*}
                pheno2=${phenos#*.}
                {
                    printf "all\t%s\t%s\t" "$pheno1" "$pheno2"
                    awk 'BEGIN{i=0}{i++;if(i==10)printf("%.2f ± %.2f\t",$2,$3);if(i==11)printf("%.2f ± %.2f\t",$2,$3);if(i==12)printf("%.2f ± %.2f\t",$2,$3);if(i==14)printf("NA\t%.0f",$2)};' "$file" #if(i==15)printf("%.4f\t",$2);if(i==16)printf("%i\t",$2);if(i==17)printf("%.4f\t",$2);if(i==18)printf("%i\n",$2)}' $file >> $outfile
                    echo ''
                } >> "$outfile"
                    #est=$(awk '$1=="rG"{print $2,$3}' $file)
            done


        else 
            for file in "$dirbiv"/"$snpgroup"*.hsq
            do
                if [[ ! -f $file ]]; then
                    continue
                fi
                fn=${file##*/}
                fn=${fn%.hsq}
                phenos=${fn#${snpgroup}[.-]}
                pheno1=${phenos%.*}
                pheno2=${phenos#*.}
                {
                    printf "%s\t%s\t%s\t" "$snpgroup" "$pheno1" "$pheno2"
                    awk 'BEGIN{i=0}{i++;if(i==10)printf("%.2f ± %.2f\t",$2,$3);if(i==11)printf("%.2f ± %.2f\t",$2,$3);if(i==12)printf("%.2f ± %.2f\t",$2,$3);if(i==17)printf("%s\t",$2);if(i==18)printf("%.0f",$2)}' "$file" #if(i==15)printf("%.4f\t",$2);if(i==16)printf("%i\t",$2);if(i==17)printf("%.4f\t",$2);if(i==18)printf("%i\n",$2)}' $file >> $outfile
                    echo ''
                } >> "$outfile"
                #est=$(awk '$1=="rG"{print $2,$3}' $file)
            done

        fi 


    done

fi


# get histogram of -log10(P-values) from GWAS using all SNPs
# echo "Histogram of -log10(P-values) from GWAS using all SNPs"
# if [ true ]; then
# for ((i=0;i<$sz;i++)); do
#         printf ${n[$i]}" "
#         awk 'NR>1{v=int(-log($9)/log(10));if(v>7)v=7;arr[v]+=1}END{for(v=0;v<7;v++)printf "%i ",arr[v];printf "%i\n",arr[7]}' $dir/gwas-all/all.${r[$i]}.assoc.linear
# done
# echo
# fi

