import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from variant import Variant

"""Plot AF vs. silva score, given a list of variants"""
def plot_freq(variants):
    afs = [x.af for x in variants]
    scores = [x.score for x in variants]
    plt.scatter(scores, afs)
    plt.xscale('log', nonposy='clip')
    plt.xlim(0.001,1)
    plt.ylim(-0.01, 0.06)
    plt.show()
    plt.savefig('test.png')
    #f = open('/dupa-filer/talf/silva-pipeline/test.out', 'w')
    #for a,s in zip(afs, scores):
        #f.write('\t'.join([str(a),str(s)]) + '\n')
    #f.close()

def plot_dists(variants, random):
    plt.hist([x.score for x in variants if x.score > 0.001], bins=10 ** np.linspace(np.log10(.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig('test3.png')
    plt.close()
    plt.hist([x.score for x in random if x.score > 0.001], bins=10 ** np.linspace(np.log10(.001), np.log10(1.0), 20), log=True)
    plt.xscale('log')
    plt.savefig('test2.png')

def plot_thresh_dist(variants, random):
    points = []
    # go from 10^-3 to 10^-1 by log inc of 10^0.125
    for i in range(17):
        x = 10 ** (-3 + float(i)/8)
        path_count = sum(v.score > x for v in variants)
        rand_count = sum(v.score > x for v in random)
        points.append(rand_count - path_count)
    plt.scatter([10 ** (-3 + float(x)/8) for x in range(17)], points)
    plt.xscale('log')
    plt.xlim(10**-3.1, 10**-0.9)
    plt.savefig('test4.png')

if __name__ == '__main__':
    variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/1000gp_rare_results.txt')
    rand_variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/silva/new_rand_results.txt')
    print sum(v.af < 0.05 for v in variants)
    print sum(v.af < 0.05 for v in rand_variants)
    print len(variants)
    print len(rand_variants)
    t,prob = st.ttest_ind([v.score for v in variants if v.af < 0.05], [v.score for v in rand_variants if v.af < 0.05], equal_var = False)
    print("%.8f" % prob)
    print t
    #plot_freq(variants)
    #plot_dists(variants, rand_variants)
    plot_thresh_dist(variants, rand_variants)
