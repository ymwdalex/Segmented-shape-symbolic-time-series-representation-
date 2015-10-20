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
# bottomup: this module provide basic functions of bottomup algorithm
#
# History
#    2014-01-07: create the file
#
# Task:
#
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import logging
import time
import matplotlib.pyplot as plt

import ssstsr_publib as publib
from fitting import *
from plot_helper import plot_segts_fit

#-------------------------------------
# Generate bottomUp initial segments: 
#    concatenate the adjacent points to form the initial segment Refactor: minimal size is an input parameter
#    The length of last part may be larger than initSize
#    input:
#        @ ts: time series
#        @ initSize: the minimal length of initial segments
#    output:
#        @ segTS: return inital segment
#--------------------------------------
def init_bottomup_segs(ts, initSize):
    segTS = []         # a list of segmentation result
    lenTS = len(ts)
    numInitSeg = lenTS // initSize # number of initial segments when using bottomUp segmentation
    for i in range(numInitSeg-1):
        segTS.append( ( i*initSize, (i+1)*initSize) )

    if lenTS < initSize:
        # extreme situation: length of ts is even less than initSize
        segTS.append( (0, lenTS) )
    else:
        segTS.append( ( (numInitSeg-1)*initSize, lenTS ) ) # append the last part

    return segTS

#--------------------------------
# merge segment:
#    In this application, the length of segment is strictly larger than 2
#    seg2 should be on the right of seg1
#--------------------------------
def merge_seg(seg1, seg2):
    logger = logging.getLogger("SegmentClass.merge_seg")
    logger.setLevel(logging.INFO)
    logger.addHandler(publib.console_handle)     

    assert seg1[0] < seg1[1] and seg2[0] < seg2[1] and seg1[1] == seg2[0]

    return (seg1[0], seg2[1])
    pass

#--------------------------------
# bottom up segmentation: generate segments with minimal number of segments k and maximum error max_err
# 
# @ input:
#        initSize: the step size of the initial segments
#        max_err: if the merge cost is greater than
#        PLOT_DEBUG: plot the bottomUp procedure
#        k: minimal number of segments
# @ output:
#        startK: 
#        endK:  
#        segTS: the list of 
#        segCostKList: 
#--------------------------------
def bottomup_seg(ts, max_err, init_size=None, k=None, PLOT_DEBUG=None):
    logger = logging.getLogger("bottomup.bottomup_seg")
    logger.setLevel(logging.INFO)
    logger.addHandler(publib.console_handle)     
    
    # assign default parameter
    if init_size == None:
        init_size = 2
    if k == None:
        k = 1
    if PLOT_DEBUG == None:
        PLOT_DEBUG = False

    assert k > 0

    # initial bottomUp segmentation
    seg_ts = init_bottomup_segs(ts, init_size)

    # computer initial segment error
    seg_cost_list, seg_fit_list = calc_segments_cost(ts, seg_ts)
    
    # Compute the costs of merging each initial segment
    merge_cost = [0,]* (len(seg_ts)-1)  # init merge_cost as 0
    for i in range(0, len(seg_ts)-1):
        merge_cost[i] = calc_seg_err(ts, merge_seg(seg_ts[i], seg_ts[i+1]))[0]

    if PLOT_DEBUG:
        fig = plt.figure(figsize=(24, 4))
        # plt.get_current_fig_manager().window.wm_geometry("400x600+20+40")
        plt.ion()

    # record the segmentation error for each K
    # TRICK: initial as numpy vector because it needs to support vector operation in GAP statistic
    seg_cost_K_list = np.array([])

    # Bottom-Up concatenation
    # two break condition: 
    #     1. merge cost is larger than max error
    #    2. number of segments is less than k
    while len(seg_ts) > k and min(merge_cost) < max_err:
        # plot segments and linear least square sum regression
        if PLOT_DEBUG:
            plot_segts_fit(ts, seg_ts, seg_fit_list)
            fig.canvas.draw()
            time.sleep(0.02)
        
        # find the minimal merge cost. If there are more than 1 minimal value, return the first one
        minInd = merge_cost.index(min(merge_cost))

        # merge and update the segment with minimal merge_cost, remove the merged segment
        seg_ts[minInd] = merge_seg(seg_ts[minInd], seg_ts[minInd+1])
        del seg_ts[minInd+1]

        # update segment cost, remove the merged cost
        residual, seg_fit = calc_seg_err(ts, seg_ts[minInd])
        seg_cost_list[minInd] = residual
        del seg_cost_list[minInd+1]

        # update the linear least square sum regression
        seg_fit_list[minInd] = seg_fit
        del seg_fit_list[minInd+1]
    
        seg_cost_K_list = np.append(seg_cost_K_list, sum(seg_cost_list))

        # update merge_cost
        del merge_cost[minInd]         
        if minInd < len(seg_ts) - 1:
            # not the last segment
            merge_cost[minInd]= calc_seg_err(ts, merge_seg(seg_ts[minInd], seg_ts[minInd+1]))[0]
        if minInd > 0:        
            # not the first segment
            merge_cost[minInd-1] = calc_seg_err(ts, merge_seg(seg_ts[minInd-1], seg_ts[minInd]))[0]

    # Plot the final result of bottomUp segmentation
    if PLOT_DEBUG:
        plot_segts_fit(ts, seg_ts, seg_fit_list)
        fig.canvas.draw()
        plt.ioff()
        # raw_input("Press Enter to continue...")

    # TRICK: reserve the seg_cost_K_list. 
    # So, seg_cost_K_list is the error from minimal number of segments to maximum number of segments
    seg_cost_K_list = seg_cost_K_list[::-1]
    return seg_ts, seg_cost_K_list, seg_fit_list
    pass

