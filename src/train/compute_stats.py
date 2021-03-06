#!/usr/bin/env python

"""
Usage: $0 < CONTROL.MAT > CONTROL.STATS

Report stats for each feature.
"""
from __future__ import with_statement, division

import sys

from numpy import array

if sys.stdin.isatty():
    print __doc__
    sys.exit(1)

features = []
sums = []
sums2 = []
N = 0

for line in sys.stdin:
    line = line.strip()
    if not line: continue
    
    if line.startswith('#'):
        # Header
        features = [token.strip().strip('#') for token in line.split()]
    else:
        # Data line
        values = array([float(val) for val in line.split()])
        if N == 0:
            sums = values
            sums2 = values * values
        else:
            sums += values
            sums2 += values * values

        N += 1

for feature, sum, sum2 in zip(features, sums, sums2):
    # Don't print stats for class
    if feature == 'class': continue
    avg = sum / N
    print "%s\t%.4f\t%.4f" % (feature, avg, (N / (N - 1)) * ((1 / N) * sum2 - (avg * avg)))
