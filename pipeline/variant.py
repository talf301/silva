import os
import sys
import logging

"""
Take in a silva output and plot/analyse comparison of allele 
to silva score.
"""

class Variant:
    def __init__(self, rank, score, clss, gene, tx, chrom, pos, id, ref, alt, af, info, ac_het, ac_hom, an):
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
        self.ac_het = ac_het
        self.ac_hom = ac_hom
        self.an = an
    """Load a results file into a list of variants"""
    @staticmethod
    def load_res_file(file):
        counter = 0
        variants = []
        for line in open(file):
            line = line.strip()
            if line.startswith('#'): continue
            line = line.split()
            info = line[-1]
            # Grab the Allele count and total from info to calculate AF
            # Or just grab the actual AF
            try:
                af = next(x for x in info.split(';') if x.startswith('AF')).split('=')[1]
                ac_adj = next(x for x in info.split(';') if x.startswith('AC_Adj')).split('=')[1]
                #ac = next(x for x in info.split(';') if x.startswith('AC')).split('=')[1]
                ac_het = next(x for x in info.split(';') if x.startswith('AC_Hom')).split('=')[1]
                ac_hom = next(x for x in info.split(';') if x.startswith('AC_Het')).split('=')[1]
                an = next(x for x in info.split(';') if x.startswith('AN')).split('=')[1]
            except StopIteration:
                logging.warning("Line without freq info " + info)
                continue
            af = af.split(',')
            if len(af) > 1:
                counter += 1
            af = 2 * float(ac_adj) / float(an)
            v = Variant(float(line[0]), float(line[1]), line[2], line[3], line[4], line[5],
                    line[6], line[7], line[8], line[9], af, info, ac_het, ac_hom, an)
            variants.append(v)
        logging.debug("%d multiple allele lines were ignored." % counter)
        return variants


if __name__ == '__main__':
    variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/1000gp_rare_results.txt')
    for v in variants[:10]: print v.af

