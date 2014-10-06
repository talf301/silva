import matplotlib.pyplot as plt
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

if __name__ == '__main__':
    variants = Variant.load_res_file('/dupa-filer/talf/silva-pipeline/1000gp_rare_results.txt')
    plot_freq(variants)
