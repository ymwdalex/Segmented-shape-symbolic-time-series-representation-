#!/usr/bin/env python

#########################################################
#                                                         #
# Segmented Shape-Symbolic Time series Representation    #
#                                                         #
# __author__ = "Zhe Sun"
# __copyright__ = "Copyright 2013, Target-holding B.V."
# __license__ = "GPL"
# __version__ = "1.0.1"
# __email__ = "zhe.sun@target-holding.com"
#
#########################################################

#===============================================================================
# montonbase: this module provides provide a suite of helper functions to do the
#             maximum likelihood shape estimate related function
#
# History
#    2013-12-24: create the file
#    2014-01-02: finish basic function
#    2014-
#
# Task:
#    1.  
#
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#===============================================================================

import random
import numpy as np
import logging

import ssstsr_publib as publib

# a hash table of shape series with length of n, n is the key
shape_cache = {} 

# shape_lib = [shape_0, shape_linear, shape_leftparab, shape_rightparab, shape_s_shape1]
# we do not use flat shape in the shape library. See the paper
# shape_lib = [publib.shape_linear, publib.shape_leftparab, publib.shape_rightparab, publib.shape_s_shape1, publib.shape_s_shape2]
shape_lib = [publib.shape_linear, publib.shape_leftparab, publib.shape_rightparab]
shape_lib_len = len(shape_lib)

# keep the same coding as MATLAB code
# RC = [[0,-1,-1]; [-1,1,2]; [-1,6,4]; [-1,3,5]; [-1,7,8]];
shape_convert_tbl = {(publib.shape_flat, 'flat'): 'a',
                     (publib.shape_linear, 'inc'): 'b',
                     (publib.shape_linear, 'dec'): 'e',
                     (publib.shape_leftparab, 'inc'): 'c',
                     (publib.shape_leftparab, 'dec'): 'f',
                     (publib.shape_rightparab, 'inc'): 'd',
                     (publib.shape_rightparab, 'dec'): 'g',
                     (publib.shape_s_shape1, 'inc'): 'h',
                     (publib.shape_s_shape1, 'dec'): 'j',
                     (publib.shape_s_shape2, 'inc'): 'i',
                     (publib.shape_s_shape2, 'dec'): 'k'}
# shape_convert_tbl = {(publib.shape_flat, "flat"): 0,
#                      (publib.shape_linear, "inc"): 1,
#                      (publib.shape_linear, "dec"): 2,
#                      (publib.shape_leftparab, "inc"): 6,
#                      (publib.shape_leftparab, "dec"): 4,
#                      (publib.shape_rightparab, "inc"): 3,
#                      (publib.shape_rightparab, "dec"): 5,
#                      (publib.shape_s_shape1, "inc"): 7,
#                      (publib.shape_s_shape1, "dec"): 8,
#                      (publib.shape_s_shape2, "inc"): 9,
#                      (publib.shape_s_shape2, "dec"): 10}


# TRICK: this list just appends shape_flat at the end of the list "shape_lib" so that the shape encoding can be done easier
shape_lib_2 = [publib.shape_linear, publib.shape_leftparab, publib.shape_rightparab, publib.shape_flat]

# thresh: the thresh is used for detecting flat shape. Please refer the paper
thresh_flat = 0.26

#-----------------------------------------------
# init_shape: build the shape series with length of n, and store in the shape_cache
def init_shape(n):
    if n in shape_cache:
        return shape_cache[n]
    else:
        x = np.linspace(0, 1, n)
        shapes = np.empty([shape_lib_len, n])

        # iterate the shape library
        for i in range(shape_lib_len): 
            shapeFun = shape_lib[i]
            shapes[i] = shapeFun(x)

        # store 
        shape_cache[n] = shapes
        return shapes

#-------------------------------------
# fit_shape: find the best fit shape
#    input:
#        @ y: a piece of time series (segment), the data type is np.array
#    return: 
#        @ fit_y: the fitted time series
#        @ likelihood: the likelihood of estimation
#        @ shape_ind: the shape index of the list shape_lib_2
#        @ shape_dir: "flat", "inc" or "dec"
#-------------------------------------
def fit_shape(y):

    mean_y = np.mean(y)
     
    if len(y) < 5:
        return np.ones(len(y)) * np.mean(y), 0, shape_lib_len, "flat"

    # load or build shape series with the same length as y
    n = len(y)
    if n in shape_cache:
        shapes = shape_cache[n]
    else:
        x = np.linspace(0, 1, n)
        shapes = np.empty([shape_lib_len, n])

        # iterate the shape library
        for i in range(shape_lib_len): 
            shapeFun = shape_lib[i]
            shapes[i] = shapeFun(x)

        # store 
        shape_cache[n] = shapes
    
    # normalize y, pay attention that standard deviation could be very small
    std_y = publib.std(y)
    if std_y < 1e-10:
        std_y = 1
    normal_y = (y - mean_y) / std_y

    # Linear least square estimation
    numerator = np.dot(shapes, normal_y) 
    # TRICK: we do not calculate norm(shape.^2, 2), save computation time
    # denominator = np.sum(shapes * shapes, axis=1)
    denominator = np.ones(shape_lib_len) * len(y)
    theta = numerator / denominator
    
    # maximum likelihood estimation:
    # Calculate the error with standard shapes, The matlab code is below:
    # err   = sum((MB.shapes.*repmat(theta,1,length(yy))- repmat(yyn,size(MB.shapes,1), 1)).^2,2);        
    err = np.sum(np.power((shapes.T * theta).T - normal_y, 2), axis=1)

    # find the minimum one        
    likelihood, shape_ind = -1 * err.min(), err.argmin()

    # TRICK: flat shape is judged by |\theta|, see the paper         
    if abs(theta[shape_ind]) < thresh_flat:
        shape_ind = shape_lib_len
        shape_dir = 'flat'
        fit_y = np.zeros(len(y))
        likelihood = -1 * np.var(y)
    elif theta[shape_ind] < 0:
        shape_dir = 'dec'
        fit_y = theta[shape_ind] * shapes[shape_ind]
    else:
        shape_dir = 'inc'
        fit_y = theta[shape_ind] * shapes[shape_ind]
    
    fit_y = fit_y * std_y + mean_y

    return fit_y, likelihood, shape_ind, shape_dir


#===============================================================================
#
# TRICK: we use LLSE to fit the shapes (see shape libraries in the module ssstsr_publib), but not ax+b
#        Comment this code, and the code following is what we use
#===============================================================================
#     #--------------------------------
#     # create segment:
#     #     A function which takes in time series and returns a linear segment approximation of it and the the approximation error
#     #    @ seg: input the segment, seg[0] is the start index, seg[1] is the end index (inclusive!!!)
#     #     @ (theta, alpha): return the linear approximation coefficient
#     #    @ residuals: error
#     #--------------------------------
#     def calc_seg_err(self, ts, seg):
# 
#         self.logger.info("calc_seg_err")
# 
#         startPoint = seg[0]
#         endPoint = seg[1]
#         assert startPoint >= 0 and startPoint < len(ts) and endPoint >= 0 and endPoint < len(ts) and endPoint > startPoint        
# 
#         #-------------------
#         # pay attention for Python beginner like me:
#         # a = [1,2,3,4,5]
#         # a[1:3] --> [2,3]
#         #-------------------
# 
#         # TRICK: if the segment is (4,9), we input (0,1,2,3,4,5) into linalg.lstsq
#         x = np.array(range(endPoint-startPoint+1))
#         y = np.array(ts[startPoint:endPoint+1])
#         A = np.vstack([x, np.ones(len(x))]).T
# 
#         [theta, alpha], residuals, _, _ = np.linalg.lstsq(A, y)
#         
#         self.logger.debug(theta, alpha, residuals)
#         return (residuals, theta, alpha)
#===============================================================================

#--------------------------------
# calc_seg_err: A function which takes in segment and returns a linear segment 
#               approximation of it and the the approximation error
#    input:
#        @ ts: time series
#        @ seg: input the segment, seg[0] is the start index, seg[1] is the end index (exclusive!!!)
#    output:
#        @ residuals: mean square error between the segment and fitted curve
#        @ fit_y: fitted curve
#--------------------------------
def calc_seg_err(ts, seg):
    logger = logging.getLogger("fitting.calc_seg_err")
    logger.setLevel(logging.INFO)  
    logger.addHandler(publib.console_handle)   

    startpoint = seg[0]
    endpoint = seg[1]
    assert startpoint >= 0 and startpoint < len(ts) and endpoint >= 0 and endpoint <= len(ts) and endpoint - 1 > startpoint        

    y = np.array(ts[startpoint:endpoint])

    fit_y, _, shape_ind, shape_dir = fit_shape(y)
    residuals = publib.mean_square_err(fit_y, y)

    return residuals, fit_y
    

#-------------------------------------
# Calculate residuals, theta and alpha for every segments
#    input:
#        @ts: time series, list
#        @ seg_ts: list of tuple, segments of the time series
#    output:
#        @ seg_cost_list: list of residuals of each segments
#        @ seg_fit_list: list of fitting curve of each segments
#--------------------------------------
def calc_segments_cost(ts, seg_ts):
    seg_cost_list = []
    seg_fit_list = []

    for seg in seg_ts:
        residual, seg_fit = calc_seg_err(ts, seg)
        seg_cost_list.append(residual)
        seg_fit_list.append(seg_fit)

    return seg_cost_list, seg_fit_list

