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
from scipy.stats import chisqprob
from variant import Variant

from argparse import ArgumentParser

__author__ = 'Tal Friedman (talf301@gmail.com)'

NUM_INTERVALS = 129
HARMFUL_THRESH = 0.05
COMMON_THRESH = 0.05

"""Test all common and harmful variants (defined by
COMMON_THRESH and HARMFUL_THRESH) to see if they
are in Hardy Weinberg equilibrium. Report the percentage
of variants which are not and some statistics about it
to the statistics file.
"""
def test_common_harmful(variants, out):
    # Get common and harmful variants, do tests
    vars = [x for x in variants if x.score > HARMFUL_THRESH and x.freq > COMMON_THRESH]
    test_stats = [(x.score, hw_test(x)) for x in vars]

    # Calculate and write the percentage with significant deviations from HW
    dev = float(sum(1 for _, _, pval in test_stats if pval < 0.05)) / len(vars)
    with open(out + 'statistics.txt', 'a') as file:
        file.write("%f%% of common harmful variants show significant deviations from\n"
                "Hardy-Weinberg equilibrium based on Pearson's chi-square test" % (100 * dev))

    # Plot chi test stat vs. score
    plt.scatter([score for score,_,_ in test_stats], [chi for _, chi, _ in test_stats])
    plt.xlabel('Silva Score')
    plt.ylabel('chi square test statistic')
    plt.savefig(open + 'chisq_vs_silva.png')
    plt.clf()
    
    # Plot silva score considered vs. % out of HW
    dens = [sum(1 for score,_,_ in test_stats if score > 0.05 * float(thresh)) for thresh in range(1,20)]
    nums = [sum(1 for score,_,pval in test_stats if pval < 0.05 and score > 0.05 * float(thresh)) for thresh in range(1,20)]
    percents = [float(num)/float(den) for num, den in zip(nums, dens)]
    plt.scatter([0.05 * float(thresh) for thresh in range(1,20)], percents)
    plt.xlabel('Silva Score cutoff')
    plt.ylabel('Percentage of common variants deviating from HW')
    plt.savefig(open + 'percent_vs_silva.png')
    plt.clf()

"""Calculate chi-square value for Pearson's Chi Square test
used on hardy weinberg equilibrium
"""
def hw_test(variant):
    p = variant.af
    q = 1 - p
    n = variant.an
    # E contains the expectations for each genotype
    E = [p*p*n, 2*p*q*n, q*q*n]
    # O contains the actual values
    O = [variant.ac_hom, variant.ac_het, variant.an - variant.ac_hom - variant.ac_het]
    chi = sum(((o - e) ** 2) / e for e,o in zip(E,O)) 
    pval = chisqprob(chi, 1)
    return chi, pval

"""Plot AF vs. silva score, given a list of variants"""
def plot_freq(variants, out):
    afs = [x.af for x in variants]
    scores = [x.score for x in variants]
    plt.scatter(scores, afs)
    plt.xscale('log', nonposy='clip')
    plt.xlim(0.001,1)
    plt.xlabel('Silva Score')
    plt.ylabel('Allele Frequency')
    plt.savefig(out + 'freq.png')
    plt.clf()
    
def plot_dists(variants, random, out):
    # Original distributions
    plot_dist(variants, 0, 1, out + 'original_dist_allfreq.png')
    plot_dist(random, 0, 1, out + 'random_dist_allfreq.png')
    # Very rare
    plot_dist(variants, 0, 0.01, out + 'original_dist_0-0.01freq.png')
    plot_dist(random, 0, 0.01, out + 'random_dist_0-0.01freq.png')
    # Medium rare
    plot_dist(variants, 0.01, 0.05, out + 'original_dist_0.01-0.05freq.png')
    plot_dist(random, 0.01, 0.05, out + 'random_dist_0.01-0.05freq.png')
    # Common
    plot_dist(variants, 0.05, 1, out + 'original_dist_0.05-1freq.png')
    plot_dist(random, 0.05, 1, out + 'random_dist_0.05-1freq.png')


"""Plot a histogram of scores, taking only variants within the given range of AF"""
def plot_dist(variants, min_freq, max_freq, out):
    plt.hist([x.score for x in variants if x.af > min_freq and x.af < max_freq], bins=10 ** np.linspace(np.log10(0.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig(out)
    plt.clf()

def plot_dists_nz(variants, random, out):
    # Original distributions
    plot_dist_nz(variants, 0, 1, out + 'original_dist_allfreq_nz.png')
    plot_dist_nz(random, 0, 1, out + 'random_dist_allfreq_nz.png')
    # Very rare
    plot_dist_nz(variants, 0, 0.01, out + 'original_dist_0-0.01freq_nz.png')
    plot_dist_nz(random, 0, 0.01, out + 'random_dist_0-0.01freq_nz.png')
    # Medium rare
    plot_dist_nz(variants, 0.01, 0.05, out + 'original_dist_0.01-0.05freq_nz.png')
    plot_dist_nz(random, 0.01, 0.05, out + 'random_dist_0.01-0.05freq_nz.png')
    # Common
    plot_dist_nz(variants, 0.05, 1, out + 'original_dist_0.05-1freq_nz.png')
    plot_dist_nz(random, 0.05, 1, out + 'random_dist_0.05-1freq_nz.png')


"""Plot a histogram of scores, taking only variants within the given range of AF, and scores above
some small threshold.
"""
def plot_dist_nz(variants, min_freq, max_freq, out):
    plt.hist([x.score for x in variants if x.af > min_freq and x.af < max_freq and x.score > 0.01], bins=10 ** np.linspace(np.log10(0.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig(out)
    plt.clf()

def plot_thresh_dist(variants, random, out):
    path = []
    rand = []
    points = []
    # go from 10^-3 to 10^-1 by log inc of 10^(2/(number of intervals-1))
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
    plt.xlabel('Silva score threshold')
    plt.ylabel('# more random than real')
    plt.savefig(out + 'thresh_dist.png')
    plt.clf()
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
    # Clear statistics file
    with open(out + 'statistics.txt', 'w') as file:
        file.write("")
    
    #Load variants
    variants = Variant.load_res_file(res1)
    rand_variants = Variant.load_res_file(res2)
    logging.info("Finished loading files.")

    # We want to:
    # Plot frequency for real variants
    plot_freq(variants, out)
    # Plot the various distribution comparisons
    plot_dists(variants, rand_variants, out)
    # Plot the distribution comparisons ignoring variants with very small scores
    plot_dists_nz(variants, rand_variants, out)
    # Do a t-test and print the results of it (p-val, t score)
    publish_t_test(variants, rand_variants, out)
    # Do the threshold distance plot
    # Find the max point in the plot and get estimated percentage of selection
    plot_thresh_dist(variants, rand_variants, out)
    # Do HW stuff
    test_common_harmful(variants, out)
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

