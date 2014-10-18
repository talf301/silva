import os
import sys
import logging

"""
Take in a silva output and plot/analyse comparison of allele 
to silva score.
"""

class Variant:
    def __init__(self, rank, score, clss, gene, tx, chrom, pos, id, ref, alt, af, info):
        self.rank = rank
        self.score = score
        self.clss = clss
        self.gene = gene
        self.tx = tx
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.af = af
        self.info = info

    """Load a results file into a list of variants"""
    @staticmethod
    def load_res_file(file):
        variants = []
        for line in open(file):
            line = line.strip()
            if line.startswith('#'): continue
            line = line.split()
            info = line[-1]
            # Grab the Allele count and total from info to calculate AF
            ac = next(x for x in info.split(';') if x.startswith('AC')).split('=')[1]
            an = next(x for x in info.split(';') if x.startswith('AN')).split('=')[1]
            af = float(ac)/float(an)

            v = Variant(float(line[0]), float(line[1]), line[2], line[3], line[4], line[5],
                    line[6], line[7], line[8], line[9], af, info)
            variants.append(v)
        return variants


if __name__ == '__main__':
    variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/1000gp_rare_results.txt')
    for v in variants[:10]: print v.af

