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
# SSSTSR_Class: this class segments the time series into monotonous pieces, represent 
#         as symbolic word
# History
#    2013-12-12: create the file
#    2013-12-17: finish bottomup_seg and GAP_decide_k statics, testing...
#    2013-12-23: fix the bug of calculating the length of segment in the function calc_segments_cost
#    2014-01-02: use MonotonBase class fit_shape function to regression
#                rename the module name according to python code convention PEP 8 
#    2014-01-03: add logging
#    2014-01-07: split bottomup, plot functions into an independent file
#
# Task:
#        3. Distance: Levenshtein distance family
#        4. Match        
#        6. remove slope line shape from shape library
#
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import sys
import random
import time
import logging
import numpy as np
from matplotlib import pyplot as plt

import ssstsr_publib as publib
from fitting import *
from sim_sig_gen import monotone_randomsignal
from bottomup import bottomup_seg
from uniform_gap_helper import *
from plot_helper import *

class SSSTSR_Class:

    #------------------------------------
    # class variable
    #------------------------------------
    #__default_ts_len =  50
    #__default_ts = [random.random() for r in range(__default_ts_len)]
    __default_ts_len =  3000
    __default_ts, _, _ = monotone_randomsignal(10,3000,200)
    __default_max_err = sys.float_info.max

    #------------------------------------
    # constructors:
    #     initial the instance with time series ts[start:end]
    #     pay attention in ts[start:end], end is exclusive!!!
    #------------------------------------
    def __init__(self, ts = None, start=None, end=None, max_err=None, logger_level=None, smooth_window_len = None, zero_thresh=0.05):
        # configure logging
        if logger_level is None:
            logger_level = logging.INFO
        
        logger = logging.getLogger("SSSTSR_Class.__init__")
        logger.setLevel(logger_level)
        logger.addHandler(publib.console_handle)

        # Initial rest parameters 
        if ts is None:
            ts = SegmentClass.__default_ts
        if len(ts) == 0:
            raise ValueError("Invalid input: time series should not be empty!")

        if start is None:
            start = 0

        if end is None:
            # pay attention: a[0:end] means a[0] .. a[end-1] !!!
            end = len(ts)
        if end - start < 2:
            raise ValueError("Invalid input: minimal length of time series should be larger than 2!")

        if max_err is None:
            max_err = SSSTSR_Class.__default_max_err
            
        if smooth_window_len is None:
            smooth_window_len = 0

        # initialize instance variables
        self.__loggerLevel = logger_level    
        self.__ts = np.array(ts[start:end])        # time series
        self.__lenTS = len(self.__ts)            # length of time series
        self.__segTS = []                        # segments
        self.__seg_fit_list = []                    # the fitting result of segments
        self.__smooth_ts = publib.smooth(self.__ts, window_len=smooth_window_len) # use hamming windows by default
        self.__normTS = publib.scale(self.__smooth_ts)
        self.__maxErr = max_err                 # max err used in bottom up
        self.__zero_thresh = zero_thresh

        self.__symbollist = []
        self.__word = ""


    #---------------------------------
    # Helper functions
    #---------------------------------
    def getclass_default_ts_len(self):
        return SegmentClass.__default_ts_len

    def getclass_default_ts(self):
        return SegmentClass.__default_ts

    def getclass_default_max_err(self):
        return SegmentClass.__default_max_err

    def get_logger_level(self):
        return self.__loggerLevel

    def get_logger_name(self):
        return self.__loggerName

    def get_ts(self):
        return (self.__ts).tolist()    # time series

    def get_smooth_ts(self):
        return (self.__smooth_ts).tolist()    # time series

    def get_len_ts(self):
        return self.__lenTS     # length of time series

    def get_seg_ts(self):
        return self.__segTS  # segments

    def get_normts(self):
        return (self.__normTS).tolist()

    def get_max_err(self):
        return self.__maxErr

    def get_fit_list(self):
        return self.__seg_fit_list

    def get_symbol_list(self):
        return self.__symbollist
 
    def get_word(self):
        return self.__word 
 
 
    #---------------------------------
    # GAP_decide_k function: use GAP_decide_k 
    #    
    #     
    #---------------------------------    

    def GAP_decide_k(self, minSize):
        logger = logging.getLogger("SSSTSR_Class.GAP_decide_k")
        logger.setLevel(self.__loggerLevel)
        logger.addHandler(publib.console_handle)        

        # calculate the error of uniform random time series
        uniformSegErr = gap_uniform(len(self.__normTS), minSize)

        # get segmentation error of scaled time series
        # TRICK: input normalize time series!
        _, ts_segCost, _ = bottomup_seg(max_err=sys.float_info.max, ts=self.__normTS, k=1, init_size=minSize)

        # TRICK: sometime regression is perfect matching, so ts_segCost is very tiny
        # cutOffPoint = (ts_segCost > 1e-10).sum()
        cutOffPoint = len(ts_segCost)

        # compute gap
        gap_K = np.log(ts_segCost[range(cutOffPoint)]) - uniformSegErr[range(cutOffPoint)]
        # gap_K_minus = ts_segCost[range(cutOffPoint)] - uniformSegErr[range(cutOffPoint)]
        
        # compute gap with first derivative, second derivative
        gapDiff1 = -np.diff(gap_K)
        gapDiff2 = -np.diff(gapDiff1)

        #-----------------------------------------------------
        # TRICK: the index of gapDiff2 start from 0, but it actually means 1
        # GAP_K = np.argmax(gapDiff2) + 1

        #-----------------------------------------------------
        # TRICK: In the ideal situation, the first derivative should become zeros from the point of K
        #         So, we add a "nearly zero" region, that gap_1st_dev < T (T is hardcode, 0.05 is this case)
        #         Then, find the K after which the first derivatives are strictly "nearly zero"
        # NOTE: nearlyZeroThresh is a very important internal parameter!
        nearlyZeroThresh = self.__zero_thresh # TODO: find it automaticlly?
        temp = (gapDiff1 >= nearlyZeroThresh) 

        if ((np.nonzero(temp[::-1]))[0]).size == 0:
            # all gapDiff1 could be less than nearlyZeroThresh
            GAP_K = GAP_K = np.argmax(gapDiff2) + 1
        else:
            GAP_K = len(gapDiff1) - (np.nonzero(temp[::-1]))[0][0]
        #logger.info("GAP choose " + str(GAP_K) + " as the number of segments.")
        
        # plot GAP_decide_k, 1st diff, 2nd diff curve
        # self.plot_gap_curve(gap_K, gapDiff1, gapDiff2, ts_segCost, uniformSegErr)
        # self.plot_gap_curve(gap_K, gapDiff1, gapDiff2, ts_segCost, uniformSegErr, gap_K_minus)
        #self.plot_gap_curve(gap_K, gapDiff1, gapDiff2)

        return GAP_K

    #---------------------------------
    # build the segmentations for the time series
    #     @ plotSwitch: the switch of if plotting the animation of building the segments
    #
    # TODO: support SWAB segment methods
    #---------------------------------    
    def build_seg(self, minSize = None, plotSwitch = None):
        if plotSwitch is None:
            plotSwitch = False
        if minSize is None:
            minSize = 2
        
        GAP_k = self.GAP_decide_k(minSize)
        
        # by using the K generated by GAP_decide_k statics, segment the time series again
        # BUGFIX: bottomup_seg(ts=self.__ts...., should input self.__normTS!
        #segTS, _, seg_fit_list, seg_symbol_list = bottomup_seg(ts=self.__ts, max_err=sys.float_info.max, init_size=minSize, k=GAP_k, PLOT_DEBUG=plotSwitch)
        segTS, _, seg_fit_list = bottomup_seg(ts=self.__normTS, max_err=sys.float_info.max, init_size=minSize, k=GAP_k, PLOT_DEBUG=plotSwitch)

        if plotSwitch:
            plot_segts_fit(self.__normTS, seg_ts=segTS, seg_fit=seg_fit_list)
            plt.show(block=False)

        self.__seg_fit_list = seg_fit_list 
        self.__segTS = segTS
        pass

    #---------------------------------
    # Encode the segments into shape and code 
    #---------------------------------    
    def seg_encode(self):
        for seg in self.__segTS:
            start = seg[0]
            end = seg[1]
            y = np.array(self.__normTS[start:end])
            _, _, shape_ind, shape_dir = fit_shape(y)
            
            self.__symbollist.append(shape_convert_tbl[(shape_lib_2[shape_ind], shape_dir)])
            
        self.__word = "".join(self.__symbollist)
        pass


    #---------------------------------
    # Plotting functions
    #---------------------------------

    # plot initialized time series
    def plot_ts(self):
        plt.plot(range(self.__lenTS), self.__ts, 'bo-')    
        plt.show()

    # plot initialized time series
    def plot_smooth_ts(self):
        plt.plot(range(self.__lenTS), self.__smoothTS, 'bo-')    
        plt.show()

    # plot scaled time series
    def plot_norm_ts(self):
        plt.plot(range(self.__lenTS), self.__normTS, 'bo-')    
        plt.show()

    # plot GAP_decide_k, GAP_decide_k 1st difference and GAP_decide_k 2nd difference
    #  @ start, end indicate the interesting part
    def plot_gap_curve(self, gap_K, gapDiff1, gapDiff2, tsCost=None, uniformCost=None, gap_K_minus=None, start=None, end=None):

        if start is None:
            start = 1
        if end is None:
            end_gap = len(gap_K)+1
            end_gap1 = len(gapDiff1)+1
            end_gap2 = len(gapDiff2)+1
            end_tscost = len(gap_K)+1
            end_uniformcost = len(gap_K)+1
        else:
            end_gap = end
            end_gap1 = end
            end_gap2 = end
            end_tscost = end
            end_uniformcost = end

        _, ax  = plt.subplots()

        ax.plot(range(start,end_gap), gap_K[start-1:end_gap-1], 'r*-', label='GAP_decide_k log')
        if not (gap_K_minus is None):
            ax.plot(range(start,end_gap), gap_K_minus[start-1:end_gap-1], 'r--', label='GAP_decide_k minus')
        if not (tsCost is None):
            ax.plot(range(start,end_tscost), tsCost[start-1:end_tscost-1], 'k--', label='TS error')
        if not (uniformCost is None):
            ax.plot(range(start,end_uniformcost), uniformCost[start-1:end_uniformcost-1], 'm--', label='Uniform random error')
        ax.plot(range(start,end_gap1), gapDiff1[start-1:end_gap1-1], 'bo-', label='GAP_decide_k 1st diff')
        ax.plot(range(start,end_gap2), gapDiff2[start-1:end_gap2-1], 'gd-', label='GAP_decide_k 2nd diff')

        # add legend
        legend = ax.legend(loc='upper right', shadow=True)    
        frame = legend.get_frame()
        frame.set_facecolor('0.90')

        plt.show()    

    # plot (smoothed) time series, fitting curve and symbol
    def plot(self):
        # original time series
        fig = plt.figure(figsize=(24, 4))
        
        # TRICK: do not input self.get_fit_list(), since fit-list store the fitting curve of normalized signal
        #        we want to output the fitting of original signal here 
        plot_segts_fit(ts=self.get_smooth_ts(), seg_ts=self.get_seg_ts(), 
                       imshow=False, shape_symbol_list = self.get_symbol_list())
        plt.show(block=False) 