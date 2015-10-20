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
# ssstsr_publib: this module provide some public functions
#
# History
#    2013-12-24: create the file
#    2014-01-02: add S shape B, mean_square_error function 
#    2014-01-06: create distance table
#                add plotting functions
#    
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import numpy as np
import logging
import matplotlib.pyplot as plt
import random
import math

#---------------------------------------------------
# logging configuration
#---------------------------------------------------
# create console handler and set level
console_handle = logging.StreamHandler()
console_handle.setLevel(logging.INFO)        
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# add formatter to console handle
console_handle.setFormatter(formatter)

#---------------------------------------------------
# normalize the vector with mean 0 and std 1
#    @ ts: np vector
#---------------------------------------------------
def normalize(ts):
    return (ts - np.mean(ts)) / np.std(ts, ddof=1)
    pass

#---------------------------------------------------
# Encapsulate np.std. This std is unbiased standard deviation
#---------------------------------------------------
def std(x, axis=None):
    if axis is None:
        return np.std(x, ddof=1)
    else:
        return np.std(x, ddof=1, axis=axis)

#---------------------------------------------------
# scale the vector within the range [0,1]
#    @ ts: np vector 
#---------------------------------------------------
def scale(ts):
    minV = min(ts)
    maxV = max(ts)

    # TRICK: np.divide is vector element-wise float division
    return np.divide(ts-minV, maxV-minV)
    pass

#---------------------------------------------------
# 
# smooth the data using a window with requested size.
# source: http://wiki.scipy.org/Cookbook/SignalSmooth
# 
# This method is based on the convolution of a scaled window with the signal.
# The signal is prepared by introducing reflected copies of the signal 
# (with the window size) in both ends so that transient parts are minimized
# in the begining and end part of the output signal.
# 
# input:
#     x: the input signal 
#     window_len: the dimension of the smoothing window; should be an odd integer
#     window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
#         flat window will produce a moving average smoothing.
# 
# output:
#     the smoothed signal
#     
# example:
# 
# t=linspace(-2,2,0.1)
# x=sin(t)+randn(len(t))*0.1
# y=smooth(x)
# 
# see also: 
# 
# np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
# scipy.signal.lfilter
# 
# TODO: the window parameter could be the window itself if an array instead of a string
# NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y
# Zhe's NOTE: it's not right!!! when window_len is a odd number, the result is not correct, see the last sentence of this function
#---------------------------------------------------
def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    # TRICK: I think there is a bug of the original code
    # s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-(window_len+1):-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    # TRICK: the length of smoothed time series is len(x)+window_size
    y=np.convolve(w/w.sum(),s,mode='valid')
   
    # TRICK: original code is not right
    # return y[(window_len/2-1):-(window_len/2)]
    return y[math.floor(window_len/2.0-1) : math.floor(-(window_len/2.0))]
    
#---------------------------------------------------
# get mean square error of two time series
#    @ x, y: np vectors 
#---------------------------------------------------
def mean_square_err(x, y):
    return np.sum(np.power(x - y, 2))
    pass

#---------------------------------------------------
# build shape library
#---------------------------------------------------
def shape_flat(x): return np.zeros(len(x))                    # flat
def shape_linear(x): return normalize(x)                        # linear
def shape_leftparab(x): return normalize(np.power(x, 2))           # left parab
def shape_rightparab(x): return normalize(-1 * np.power((1-x), 2))  # right parab
def shape_s_shape1(x): return normalize(-1 * np.cos(np.pi*x))     # S shape A
# TODO: what's the correct name of this kind of shape?
def shape_s_shape2(x): return normalize(-1*np.arccos(2*x-1))     # S shape B, diagonal symmetric with S shape A

  


