#!/usr/bin/env python3

"""GRM using SNPs located in genic regions only"""

import os
import preprocessing
# from configall import PLINK, GCTA, NBPROC, USE_SBATCH
import config_dataset
import sys
import pandas as pd


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)
    # in_prefix = 'all'
    in_file_pruned = os.path.join(config.pru_dir, 'all')
    # in_dir_grm_filtered = os.path.join(GRM_DIR, 'all-0.025')
    in_genicregions = os.path.join(config.external_dir, 'hg19genes.txt')
    in_geneids = os.path.join(config.external_dir, 'genic.txt')
    out_dir_grm_genic = os.path.join(config.grm_dir, 'grm-genic')
    out_dir_grm_nongenic = os.path.join(config.grm_dir, 'grm-genic')

    margins = [0, 10, 20, 30, 40, 50]

    # =========  Genic / non-genic  =============

    for margin in margins:
        # Extract SNPs located in the genic regions
        listsnp_margin = os.path.join(out_dir_grm_genic, 'genic-margin' + str(margin) + '.snplist')
        preprocessing.plink_make_set(bfile=in_file_pruned,
                                     border=margin,
                                     subset=in_geneids,
                                     gene_set=in_genicregions,
                                     out_file=os.path.splitext(listsnp_margin)[0],
                                     plink=config.plink
                                     )

        out_file_grm_genic_margin = os.path.join(out_dir_grm_genic, 'genic-margin' + str(margin),
                                                 'genic-margin' + str(margin))
        out_file_grm_nongenic_margin = os.path.join(out_dir_grm_nongenic, 'nongenic-margin'
                                                    + str(margin), 'nongenic-margin' + str(margin))

        # Compute genic GRM using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_genic_margin,
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--extract", str(listsnp_margin),
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch
                               )

        # Compute non genic GRM using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_nongenic_margin,
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--exclude", str(listsnp_margin),
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch
                               )

        # ========= Genic / xxk upstream and downstream / non-genic  =============


        if margin > 0:
            # in_50kupdownregions = os.path.join(config.external_dir, 'hg19genes_50kupstream_50kdownstream.txt')
            # out_dir_grm_50kupdown = os.path.join(config.grm_dir, 'grm-genic')

            # SNPs located in the genic regions, with margin = 0
            listsnp_margin0 = os.path.join(out_dir_grm_genic, 'genic-margin' + str(0) + '.snplist')
            readsnp_genic_margin0 = pd.read_table(listsnp_margin0)
            readsnp_genic_margin0 = readsnp_genic_margin0.set_index(readsnp_genic_margin0.iloc[:, 0])

            # SNPs located in the genic regions, with margin = margin of the iteration
            readsnp_genic_margin = pd.read_table(listsnp_margin)
            readsnp_genic_margin = readsnp_genic_margin.set_index(readsnp_genic_margin.iloc[:, 0])

            # difference between the two previous sets gives the upstream/downstream SNPs
            filesnp_updown_margin = os.path.join(out_dir_grm_genic, 'updown-margin' + str(margin) + '.snplist')
            readsnp_updown_margin = readsnp_genic_margin.index.difference(readsnp_genic_margin0.index)
            readsnp_updown_margin = pd.DataFrame(data=readsnp_updown_margin)
            readsnp_updown_margin.to_csv(filesnp_updown_margin, sep='\t', index=False)

            out_file_grm_updown = os.path.join(out_dir_grm_genic, 'updown-margin' + str(margin),
                                               'updown-margin' + str(margin))
            # out_file_grm_genic_margin = os.path.join(out_dir_grm_genic, 'genic-margin' + str(margin),
            #                                             'genic-margin' + str(margin))
            # out_file_grm_nongenic_margin = os.path.join(out_dir_grm_nongenic, 'nongenic-margin50',
            #                                            'nongenic-margin50')

            # Compute upstream/downstream GRM using unrelated individuals
            preprocessing.gcta_grm(in_file=in_file_pruned,
                                       par_input='--bfile',
                                       out_file=out_file_grm_updown,
                                       gcta=config.gcta,
                                       per_chr=False,
                                       ncpus=config.nbproc,
                                       other_gcta_par=["--extract", str(filesnp_updown_margin),
                                                       "--keep", config.keep_ind],
                                       sbatch=config.use_sbatch
                                       )

    # ========= Genic / 0-20k and 20-50k upstream and downstream / non-genic  =============
    if (20 in margins and 50 in margins):
        # SNPs located in the genic regions, with margin = 50
        listsnp_margin1 = os.path.join(out_dir_grm_genic, 'genic-margin' + str(50) + '.snplist')
        readsnp_genic_margin1 = pd.read_table(listsnp_margin1)
        readsnp_genic_margin1 = readsnp_genic_margin1.set_index(readsnp_genic_margin1.iloc[:, 0])

        # SNPs located in the genic regions, with margin = 20
        listsnp_margin2 = os.path.join(out_dir_grm_genic, 'genic-margin' + str(20) + '.snplist')
        readsnp_genic_margin2 = pd.read_table(listsnp_margin2)
        readsnp_genic_margin2 = readsnp_genic_margin2.set_index(readsnp_genic_margin2.iloc[:, 0])

        # difference between the two previous sets gives the SNPs in the 20k-50k regions upstream and downstream
        filesnp_updown_margin2 = os.path.join(out_dir_grm_genic, 'updown-margin' + "20-50" + '.snplist')
        readsnp_updown_margin2 = readsnp_genic_margin1.index.difference(readsnp_genic_margin2.index)
        readsnp_updown_margin2 = pd.DataFrame(data = readsnp_updown_margin2)
        readsnp_updown_margin2.to_csv(filesnp_updown_margin2, sep='\t', index=False)

        # snp in 0-20k margin
        filesnp_updown_margin1 = os.path.join(out_dir_grm_genic, 'updown-margin' + str(20) + '.snplist')
        out_file_grm_updown2 = os.path.join(out_dir_grm_genic, 'updown-margin' + "20-50",
                                   'updown-margin' + "20-50")


        # Compute upstream/downstream GRM using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_updown2,
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--extract", str(filesnp_updown_margin2),
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch
                               )


if __name__ == '__main__':
    main(sys.argv[1])
