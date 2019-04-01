#!/usr/bin/env python3

"""Download gene lists from the web"""

import os
import sys
import re
from urllib.request import urlopen
import config_dataset
# from configall import EXTERNAL_DIR, SELGENES_CNSEXPRESSION, SELGENES_NEURODEV

import pymysql
# pymysql.install_as_MySQLdb()
# import MySQLdb


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    rs_file = os.path.join(config.external_dir, 'refseq19.txt')
    genes_file = os.path.join(config.external_dir, 'hg19genes.txt')
    genes_file2 = os.path.join(config.external_dir, 'hg19genes_strand.txt')
    genic_file = os.path.join(config.external_dir, 'genic.txt')
    if not (os.path.exists(rs_file) and os.path.exists(genes_file) and os.path.exists(genic_file)):
        # get hg19 genes from UCSC server
        con = pymysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', db='hg19')
        cur = con.cursor()
        cur.execute("SELECT chrom,txStart,txEnd,strand,cdsStart,cdsEnd,exonStarts,exonEnds,exonCount,name,name2"
                    " from refGene")
        tmp = cur.fetchall()
        con.close()

        # remove data in non straightforward chromosomes
        result = []
        for row in tmp:
            if row[0][0:4] != "chrU" and '_' not in row[0]:
                result.append(row)

        # save data in refseq19.txt
        file = open(rs_file, 'w')
        for row in result:
            file.write("%s\n" % ("\t".join(map(str, row))))
        file.close()

        # extract only chr, start bp, end bp and refseq gene names, save to hg19genes.txt
        file = open(genes_file, 'w')
        for row in result:
            file.write("%s\n" % ("\t".join(map(str, [row[0], row[1], row[2], row[10]]))))
        file.close()

        # extract only chr, start bp, end bp, strand, and refseq gene names, save to hg19genes_strand.txt
        file = open(genes_file2, 'w')
        for row in result:
            file.write("%s\n" % ("\t".join(map(str, [row[0], row[1], row[2], row[3], row[4], row[5],
                                                     row[6], row[7], row[8], row[10]]))))
        file.close()

        # extract the coordinates of the region including the first exon and intron
        intron1exon1_file = os.path.join(config.external_dir, 'hg19genes_intron1exon1.txt')
        # notintron1exon1_file = os.path.join(config.external_dir, 'hg19genes_notintron1exon1.txt')
        filein1ex1 = open(intron1exon1_file, 'w')
        # filenotin1ex1 = open(notintron1exon1_file, 'w')
        for row in result:
            exon_starts = (re.split(',', row[6].decode()))
            exon_stops = (re.split(',', row[7].decode()))

            if len(exon_stops) > 2:
                strand = row[3]
                if strand == "-":
                    # exon1 = coordinates of last exons in the list
                    exon1 = (exon_starts[len(exon_starts)-2],
                             exon_stops[len(exon_stops) - 2])
                    intron1 = (int(exon_stops[len(exon_stops) - 3]) + 1,
                               int(exon_starts[len(exon_starts) - 2]) - 1)

                    exon1intron1 = (intron1[0], int(exon1[1]))
                    # notexon1intron1 = (int(exon_starts[0]),
                    #                    int(exon_stops[len(exon_stops) - 3]))

                else:
                    # exon1 = coordinates of first exons in the list
                    exon1 = (exon_starts[0],
                             exon_stops[0])
                    intron1 = (int(exon_stops[0]) + 1,
                               int(exon_starts[1]) - 1)

                    exon1intron1 = (int(exon1[0]), int(intron1[1]))
                    # notexon1intron1 = (int(exon_starts[1]),
                    #                    int(exon_stops[len(exon_stops) - 2]))

                filein1ex1.write("%s\n" % ("\t".join(map(str, [row[0], exon1intron1[0], exon1intron1[1], row[10]]))))
                # filenotin1ex1.write("%s\n" % ("\t".join(map(str, [row[0], notexon1intron1[0], notexon1intron1[1],
                # row[10]]))))

            else:
                exon1 = (exon_starts[0],
                         exon_stops[0])
                exon1intron1 = exon1

        filein1ex1.close()
        # filenotin1ex1.close()

        # save a list of all gene names to us in genic enrichment
        file = open(genic_file, 'w')
        sorted_result = []
        for row in result:
            sorted_result.append(row[10])
        sorted_result = sorted(sorted_result)
        for row in sorted_result:
            file.write("%s\n" % str(row))
        file.close()
    if not os.path.exists(config.selgenes_cnsexpression):
        # get genes used for cns expression enrichment
        response = urlopen('https://github.com/r03ert0/ENIGMA-GCTA/raw/master/lists/'
                           'cnsexpression.txt')
        file = open(config.selgenes_cnsexpression, 'w')
        file.write(response.read().decode("utf-8"))
        file.close()
    if not os.path.exists(config.selgenes_neurodev):
        # get genes used for neurodev enrichment
        response = urlopen('https://github.com/r03ert0/ENIGMA-GCTA/raw/master/lists/neurodev.txt')
        file = open(config.selgenes_neurodev, 'w')
        file.write(response.read().decode("utf-8"))
        file.close()


if __name__ == '__main__':
    main(sys.argv[1])
