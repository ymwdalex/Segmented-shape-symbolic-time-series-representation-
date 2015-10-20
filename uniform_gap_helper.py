#!/usr/bin/env python

#########################################################
#                                                         #
# Segmented Shape-Symbolic Time series Representation    #
#                                                         #
# __author__ = "Zhe Sun"
# __copyright__ = "Copyright 2013, Target-holding B.V."
# __license__ = "FreeBSD"
# __version__ = "1.0.1"
# __email__ = "zhe.sun@target-holding.com"
#
#########################################################

#----------------------------------------------------------
# 
# uniform_gap_helper: this file contains some helper function of generating uniform
#                     random time series
# History
#    2013-01-04: create the file
#                run save_uni_gap_lib and create cache file 
#
# Task:
#
#
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import math
import logging
import json
import time
import random
import os.path
import numpy as np
import sys

import ssstsr_publib as publib
from bottomup import bottomup_seg

#---------------------------------
# Statistic the error of uniform random signal
#    @ minSize: the lenght of minimal size of segments 
#     
#---------------------------------    
def gap_uniform(lenTS, minSize = None, num_round = None):
    if minSize is None:
        minSize = 2
    # See paper: we estimate error of segmentation by taking the mean of 10 uniformly randomly time series
    if num_round is None:
        num_round = 10

    random.seed(time.time())

    # estimate error from uniform random time series
    mat = np.array([])
    for i in range(num_round):
        # generate uniform time series in the range of [0,1]
        uniformTS = np.array([random.uniform(0,1) for _ in range(lenTS)])

        # bottomUp: use initSize=2 and k=1 to iterate every possible number of segments
        _, uniformSegCost, _ = bottomup_seg(ts=uniformTS, max_err=sys.float_info.max, k=1, init_size=minSize)

        # the testing K must be in [1, len(uniformTS) // 2]
        assert len(uniformSegCost) == (len(uniformTS) // minSize - 1)

        # the length of uniformSegCost
        mat = mat.reshape(i, len(uniformSegCost))
        mat = np.vstack([mat, np.array(uniformSegCost)])
        pass

    # calculate the weighted GAP and get smoothing factor GAP
    mat = np.log(mat)
    sk = (num_round - 1) * publib.std(mat,axis=0) / num_round * np.sqrt(1 + 1/num_round)
    weightGAP = np.mean(mat, axis=0) - sk

    # return num_round uniform random experiments
    return weightGAP

#---------------------------------
# The foldername where save the cache file
#---------------------------------   
foldername = os.path.join(os.path.dirname(__file__), 'uni_gap_lib')

#---------------------------------
# save_uni_gap_lib: calculate the error of uniform random signal, and save into cache files
#    the naming convention of the files are "ugl_(start)_(end).txt"
#    The files are saved under the folder "./uni_gap_lib/"
#    @ start: the start number of the length of time series
#    @ end: the end number of the length of time series
#    @ num_round: the repetition time of calculating uniform random signal error, default value is 10
#---------------------------------   
def save_uni_gap_lib(start, end,  num_round):
    
    assert start > 2 and end > start

    logger = logging.getLogger("uniform_gap_helper.save_uni_gap_lib")
    logger.setLevel(logging.INFO)
    logger.addHandler(publib.console_handle)  
    
    logger.error("test")
    if not os.path.exists(foldername):
        try:
            os.makedirs(foldername)
        except OSError:
            logger.error("make directory failed:\t" + foldername)
            pass    

    uni_gap_lib = {}
    for n in range(start, end):
        print n
        for minSize in range(2, int(math.sqrt(n)), 1):
            gap = gap_uniform(n, minSize, num_round).tolist()
            key = str((n, minSize))
            uni_gap_lib[key] = gap

    filename = 'ugl_' + str(start) + "_" + str(end) + ".txt"
    wholename = os.path.join(foldername, filename)
    
    # TODO: what if the file has already been there?
    with open(wholename, 'w') as outfile:
        try:
            json.dump(uni_gap_lib, outfile)
        except:
            logger.error("create cache file failed:\t" + wholename)
            pass
    
def load_uni_gap_lib():
    uni_gap_lib = {}
    print os.listdir (foldername)
    for filename in os.listdir (foldername):
        print os.path.join(foldername, filename)
        with open(os.path.join(foldername, filename), 'r') as infile:
            uni_gap_lib.update(json.load(infile))
            
    return uni_gap_lib
    
# if __name__ == '__main__':
#     step = 5
#     for i in range(10, 15, step):
#         save_uni_gap_lib(i, i+step, num_round=50)
