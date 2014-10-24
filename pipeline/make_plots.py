"""
This is a script for analyzing the results of running silva on a collection of synonymous
as well as the matched, randomly generated variants created from these.
"""


import os
import sys
import logging

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from variant import Variant

from argparse import ArgumentParser
NUM_INTERVALS=33

__author__ = 'Tal Friedman (talf301@gmail.com)'

"""Plot AF vs. silva score, given a list of variants"""
def plot_freq(variants, out):
    afs = [x.af for x in variants]
    scores = [x.score for x in variants]
    plt.scatter(scores, afs)
    plt.xscale('log', nonposy='clip')
    plt.xlim(0.001,1)
    plt.ylim(-0.01, 0.06)
    plt.show()
    plt.savefig(out + 'freq.png')
    plt.close()
    #f = open('/dupa-filer/talf/silva-pipeline/test.out', 'w')
    #for a,s in zip(afs, scores):
        #f.write('\t'.join([str(a),str(s)]) + '\n')
    #f.close()

def plot_dists(variants, random, out):
    #plt.hist([x.score for x in variants if x.score > 0.001], bins=10 ** np.linspace(np.log10(.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig(out + 'original_dist.png')
    plt.close()
    #plt.hist([x.score for x in random if x.score > 0.001], bins=10 ** np.linspace(np.log10(.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig(out + 'random_dist.png')
    plt.close()

def plot_thresh_dist(variants, random, out):
    path = []
    rand = []
    points = []
    # go from 10^-3 to 10^-1 by log inc of 10^0.125
    for i in range(NUM_INTERVALS):
        x = 10 ** (-3 + float(i)/((NUM_INTERVALS-1)/2))
        path_count = sum(v.score > x for v in variants)
        path.append(path_count)
        rand_count = sum(v.score > x for v in random)
        rand.append(rand_count)
        points.append(rand_count - path_count)
    plt.scatter([10 ** (-3 + float(x)/((NUM_INTERVALS-1)/2)) for x in range(NUM_INTERVALS)], points)
    plt.xscale('log')
    plt.xlim(10**-3.1, 10**-0.9)
    plt.savefig(out + 'thresh_dist.png')
    plt.close()
    # Find the max and print percentage at that point
    pts = [(10 ** (-3 + float(i)/((NUM_INTERVALS-1)/2)), points[i], i) for i in range(NUM_INTERVALS)]
    maxpt, pt, ci = max(pts, key=lambda x: x[1])
    with open(out + "statistics.txt", 'a') as file:
        file.write("The maximum threshold occurs at %.3f, with %d more random mutations (%d) than true polymorphisms (%d)." % (maxpt, pt, rand[i], path[i]))
        file.write("This suggests a %f%% rejection rate at the max threshold." % (float(pt)/path[i] * 100))
        file.write("We estimate a synonymous substitution rejection rate of %f%% (%d/%d SNPs)." % ((float(pt)/len(variants) * 100), pt, len(variants)))

def publish_t_test(variants, random, out):
    vav = sum(v.score for v in variants) / len(variants)
    rav = sum(v.score for v in random) / len(random)
    t,prob = st.ttest_ind([v.score for v in variants], [v.score for v in random], equal_var = False)
    with open(out + "statistics.txt", 'a') as file:
        file.write("In total there are %d original variants and %d matched random variants" % (len(variants), len(random)))
        file.write("Mean silva score for real variants: %.5f\n" % vav)
        file.write("Mean silva score for matched generated variants: %.5f\n" % rav)
        file.write("p-value and actual t-score based on Welch's t-test: %.8f, %.3f" % (prob, t))

def script(res1, res2, out, **kwargs):
    variants = Variant.load_res_file(res1)
    rand_variants = Variant.load_res_file(res2)
    logging.info("Finished loading files.")
    # We want to:
    # Plot frequency for real variants
    plot_freq(variants, out)
    # Plot the distribution comparison in a smart way
    plot_dists(variants, rand_variants, out)
    # Do a t-test and print the results of it (p-val, t score)
    publish_t_test(variants, rand_variants, out)
    # Do the threshold distance plot
    # Find the max point in the above plot and get estimated percentage of selection
    plot_thresh_dist(variants, rand_variants, out)
    # Anything else i'd like to add

def parse_args(args):
    parser = ArgumentParser(description=__doc__.strip())
    
    parser.add_argument('res1', metavar='ORIG_RES', 
        help='Results from original variants')
    parser.add_argument('res2', metavar='RAND_RES',
        help='Results from randomly generated variants')
    parser.add_argument('out', metavar='OUT_DIR',
        help='Directory in which to put resulting plots and analysis')
    return parser.parse_args(args)
    

def main(args = sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level='INFO')
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())

#if __name__ == '__main__':
    #variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/1000gp_rare_results.txt')
    #rand_variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/silva/new_rand_results.txt')
    #print sum(v.af < 0.05 for v in variants)
    #print sum(v.af < 0.05 for v in rand_variants)
    #print len(variants)
    #print len(rand_variants)
    #t,prob = st.ttest_ind([v.score for v in variants if v.af < 0.05], [v.score for v in rand_variants if v.af < 0.05], equal_var = False)
    #print("%.8f" % prob)
    #print t
    #plot_freq(variants)
    #plot_dists(variants, rand_variants)
    #plot_thresh_dist(variants, rand_variants)
