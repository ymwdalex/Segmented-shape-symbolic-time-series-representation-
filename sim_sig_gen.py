#!/usr/bin/env python

#########################################################
#                                                       #
# Segmented Shape-Symbolic Time series Representation   #
#                                                       #
# __author__ = "Zhe Sun"
# __copyright__ = "Copyright 2013, Target-holding B.V."
# __license__ = "FreeBSD"
# __version__ = "1.0.1"
# __email__ = "zhe.sun@target-holding.com"
#
#########################################################

#----------------------------------------------------------
# 
# This module provide artificial signal generator
#   This file is transplanted from MATLAB code. 
#   Due to time issue, I don't refactoring the MATLAB code
#
# History
#   2013-12-20: create the file
#   2014-01-02: refactoring the code
#
# Task:
#
# TRICK labels show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import random
import time
import numpy as np

#-----------------------------------------
# monotone_randomsignal: generate synthetic time series with the length of "length", 
#                       "nsegments" monotonous segments and with at least "minsegmentsize" samples
# return: 
#   @ yy: random monotonous signal
#   @ cc: signal type of each data point
#   @ R: the range of each segment
#
# usage: monotone_randomsignal(10,3000,200)
#------------------------------------------
def monotone_randomsignal(nsegments, length, minsegmentsize=None):
    
    if minsegmentsize is None:
        minsegmentsize=1

    random.seed(time.time())

    # 7 typical signal, see the paper esann
    signals = [\
        lambda x: np.zeros(len(x)), \
        lambda x: x, \
        lambda x: -x, \
        lambda x: np.power(x,2), \
        lambda x: -np.power(x,2), \
        lambda x: -np.power(1-x, 2), \
        lambda x: np.power(1-x, 2), \
    ]

    getsignal = lambda sigType, n: signals[sigType](np.linspace(0,1,n))

    # methods to create 'subsegments', e.g. one sine wave consists of 4 different parabola segments.
    # TODO: What is equidist??? Do we need it? the MATLAB code is as the following line.
    # equidist = @(n,labs) [reshape(repmat(labs(1:end-1), floor(n/length(labs)), 1),1,(length(labs)-1)*floor(n/length(labs))) , 
    #                        repmat(labs(end), 1, n-(length(labs)-1)*floor(n/length(labs)))]; 
    # TRICK: still, the second argument of "arange" is exclusive 
    equidistb = lambda n, nlabs: np.array(np.arange(0, nlabs)) * n//nlabs
    equidistt = lambda n, nlabs: np.append(np.array(np.arange(1, nlabs)) * n//nlabs, n-1)

    labelings = [\
        lambda n: np.zeros(1,n), \
        lambda n: np.ones(1,n), \
        lambda n: 2*np.ones(1,n), \
        lambda n: 6*np.ones(1,n), \
        lambda n: 4*np.ones(1,n), \
        lambda n: 3*np.ones(1,n), \
        lambda n: 5*np.ones(1,n),\
    ]

    nsubsegments = [1] * 7
    getlabeling = lambda i, n: (labelings[i])(n)
    negl = [0,2,1,5,6,3,4];
    negatelabels = lambda labs: negl[labs]

    spread      = lambda nsegments: np.cumsum(np.random.rand(nsegments,1))
    formatsegs  = lambda segs, length: np.append(segs.T, length) - np.insert(segs.T, 0, 0)
    randsegs_help0 = lambda nsegments, length, minsegmentsize, s: ((np.arange(1, nsegments).T)*minsegmentsize)+np.round((length-nsegments*minsegmentsize)*s[0:-1]/s[-1])   
    randsegs_help1 = lambda nsegments, length, minsegmentsize: randsegs_help0(nsegments,length,minsegmentsize,spread(nsegments))
    randlens    = lambda nsegments, length, minsegmentsize: formatsegs(randsegs_help1(nsegments,length,minsegmentsize), length)
    # TRICK: different from MATLAB code, the start index in Python is 0
    randfuncs   = lambda nsegments: np.floor(len(signals)*np.random.rand(nsegments,1))

    lens = randlens(nsegments, length, minsegmentsize).astype(int)
    labs = randfuncs(nsegments).astype(int)
    yy = np.zeros((length))
    cc = np.zeros((length))
    last = np.random.rand()*2
    lasti = 0 # TRICK: start with 0
    R = np.zeros((2, nsegments))
    Ri = 0

    for i in range(nsegments):
        l = lens[i]
        f = 0.5 + 0.5 * np.random.rand()
        if(np.random.rand()<0.5):
            f = f * -1
            cc[lasti:lasti+l-1] = negatelabels(getlabeling(labs[i][0], l));
        else:
            cc[lasti:lasti+l-1] = getlabeling(labs[i][0], l);

        part = f * getsignal(labs[i], l);

        # append the random generated signal to "yy"
        # TRICK: do not need minus 1
        yy[lasti:lasti+l] = part-part[1] + last
        
        nsub = nsubsegments[labs[i][0]]
        # TRICK: python has to use range function, the second argument is exclusive
        R[:,range(Ri,(Ri+nsub))] = np.vstack([equidistb(l, nsub), equidistt(l, nsub)]) + lasti
        lasti=lasti+l;
        Ri = Ri+nsub

        # update last index
        last = yy[lasti-1];

    return yy, cc, R
    
