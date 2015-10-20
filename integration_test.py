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
# integrationTest.py: this file mainly test the behaviour of SegmentClass
#                     by using the simulation signal generator or real dijk data
#
# History
#    2013-12-20: create the file
#    2014-01-05: add idijk data
# Task:
#
# TRICK labels show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import random
import array
import time
import logging
import csv
import numpy as np
import scipy.io

# DDSC module
from ssstsr import SSSTSR_Class
from fitting import *
from sim_sig_gen import monotone_randomsignal
import ssstsr_publib as publib
from plot_helper import *
from dist_helper import *

def monotonousSigTest():
    # generate a synthetic signal which composed by 10 basic shape and has the length of 3000 
    test_ts, _, _ = monotone_randomsignal(10,3000,200)

    a = SSSTSR_Class(ts=test_ts)
    a.build_seg(plotSwitch = False, minSize=50)
    a.seg_encode()
    
    fig = plt.figure(figsize=(24, 4))
    plot_segts_fit(ts=a.get_ts(), seg_ts=a.get_seg_ts(), seg_fit=a.get_fit_list(), imshow=True, shapelist = a.get_symbol_list)
    
    print a.get_symbol_list
    print ''.join(a.shapelist)


def hamming_dist_test():
    print "------------------"
    print hamming_match('abc', 'vbbbbbcccccv')
    print "------------------"
    print hamming_match_all('abc', 'bbbbbcccccv')
    print "------------------"
    print hamming_match_all('abc', 'abcbbbbcccccv')
    print "------------------"
    print hamming_match('abc', 'bcabcbbcccccv')
    print "------------------"
    print hamming_match('abc', 'bcabcbabc')
    print "------------------"
    print hamming_match_all('abc', 'bcabcbabc')    
    print "------------------"
    print hamming_match_best_effort('abc', 'bcabcbabc', 3)    

def fit_shape_test():
    print fit_shape(np.array([random.random() for _ in range(10)]))        

def smooth_test():
    t=np.linspace(-4,4,100)
    x=np.sin(t)
    xn=x+np.random.randn(len(t))*0.1
    y=smooth(x, window_len=5,window='flat')
    plt.plot(x, 'r')
    plt.plot(xn, 'b')
    plt.plot(y, 'g')
     
    plt.show()

def bottomUpGAPTest():
    c = SSSTSR_Class(ts=map(float, range(10) + range(10,0,-3) + range(2,20,3) + range(15,3,-1) + (np.cos(np.arange(0,6,0.2))*5).tolist()),
                     zero_thresh=0.05)
    c.build_seg(plotSwitch = True)
    pass
  
def logtest():
    c = SSSTSR_Class(ts=np.cos(np.arange(0,20,0.1)), logger_level=logging.DEBUG, start=0, end=10)
    c.build_seg(plotSwitch = True)
    pass

if __name__ == '__main__':
#     monotonousSigTest()
    # bottomUpGAPTest()
#     print distance_table
#     print sym_avg_dist
#     print edit_dist('eeba', 'abac')
#     print edit_dist('abc', 'cba')
#     print edit_dist('cbc', 'eba')
#     print edit_dist('recoginze', 'recognize')

#     print_distance_table()
#     print nw_dist('ab', 'ba')
#     print nw_dist('ac', 'ca')
#     print nw_dist('be', 'eb')
#  
 
    print leven_match('aba', "c abba c", 4)
    print edit_dist('Levenshtein', 'Lenvinsten')
    pass
