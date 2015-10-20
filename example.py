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
# example.py: this module provide some testing function and examples of using this package
#
# History
#    2014-01-07: create the file
#    
# TRICK labels: show some parts which are not easy to understand
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
import threading

# DDSC module
from ssstsr import SSSTSR_Class
from fitting import *
from sim_sig_gen import monotone_randomsignal
import ssstsr_publib as publib
from plot_helper import *
from dist_helper import *

# def thread_plot(ssstrs_obj):
#     fig = plt.figure(figsize=(24, 4))
#     plot_segts_fit(ts=ssstrs_obj.get_smooth_ts(), seg_ts=ssstrs_obj.get_seg_ts(), seg_fit=ssstrs_obj.get_fit_list(), 
#                    imshow=False, shape_symbol_list = ssstrs_obj.get_symbol_list())
#     plt.show() 
#     detach_display()
# 
# def run_thread_plot(ssstrs_obj):
#     t = threading.Thread(target=thread_plot, args = (ssstrs_obj,))
#     t.daemon = True
#     t.start()

def idijktest():
    min_size = 5
    smooth_window_len = 5
    zero_thresh = 0.001
    
    # load idijk data
    # mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_east-Humidity.mat')
    # mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_west-Rainfall.mat')
    # mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_east-Radiation.mat')
    # mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_east-Temperature')
    # mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_east-Winddirection')
    mat = scipy.io.loadmat('D:\Dropbox\TH\DDSC\Matlab\data\weatherstation_east-Windspeed')
    
    # fetch data and time tags
    data = mat['values'][0]
    times = mat['times'][0]
    
    # adjust smooth_window_len to get different level of smoothing
    a = SSSTSR_Class(ts=data, smooth_window_len = smooth_window_len, start=0, end=5000, zero_thresh=zero_thresh)
    a.build_seg(plotSwitch = False, minSize=min_size)
    a.seg_encode()
    a.plot()

    print "Word of the database:\t" + a.get_word()
    print "The number of symbols:\t" + str(len(a.get_word()))
    print "-----------------------------------"           
    
    piece_range = range(5200, 5400)
#     piece_range = range(1600, 1800)
#     piece_range = range(5000, 5200)

    piece = [data[i] for i in piece_range]
    piece_noise = piece + np.random.randn(len(piece))*2
    b = SSSTSR_Class(ts=piece, smooth_window_len = smooth_window_len, zero_thresh=zero_thresh)
    b.build_seg(plotSwitch = True, minSize=min_size)
    b.seg_encode()
    b.plot()
    print "The keyword is \"" + b.get_word() + "\", contains " + str(len(b.get_word())) + " symbols" 
    print "-----------------------------------"   

    # search the best n matching
    n = 3
    print "Searching the best " + str(n) + " matching time series pieces..."    

# Hamming distance matching    
    match_pos_list, match_dist_list = hamming_match_best_effort(b.get_word(), a.get_word(), n)
    for i in range(n):
        pos = match_pos_list[i]
        dist = match_dist_list[i]
        match_str_startpos = (a.get_seg_ts()[pos])[0]
        match_str_endpos = (a.get_seg_ts()[pos+len(b.get_word())-1])[1]
        print "Matching " + str(i+1) + ":\t" + a.get_word()[pos:pos+len(b.get_word())-1] + \
                "[" + str(match_str_startpos) + ":" + str(match_str_endpos) + "] dist=" + str(dist) 
  
    # plot matching results   
    fig = plt.figure()
    ax = plt.subplot(n+1,1,1)
    ax.plot(piece_range, b.get_smooth_ts(), linewidth=1, color='r')
    ax.set_title("keyword time series: " + b.get_word())
    ts = a.get_smooth_ts()
    segs = a.get_seg_ts()
    title_text = []
    for i in range(n):
        pos = match_pos_list[i]
        dist = match_dist_list[i]
        match_str_startpos = (segs[pos])[0]
        match_str_endpos = (segs[pos+len(b.get_word())])[1]
          
        ax = plt.subplot(n+1,1,i+2)
        ax.plot(range(match_str_startpos, match_str_endpos), ts[match_str_startpos:match_str_endpos], linewidth=1, color='b')
        ax.set_title("Matching " + str(i) + ": " + a.get_word()[pos:pos+len(b.get_word())] + ", dist=" + str(dist))
              
    fig.subplots_adjust(hspace=1)
    plt.show(block=True) 

    # Hamming distance matching    
#     match_pos_list, match_dist_list = hamming_match_best_effort(b.get_word(), a.get_word(), n)

    # Levenshtein distance matching   
#     match_pos_list, matching_words, match_dist_list, _ = leven_match(b.get_word(), a.get_word(), n)
#  
#     for i in range(n):
#         len_matched_word = len(matching_words[i])
# #         len_matched_word = len(b.get_word())
#         pos = match_pos_list[i]
#         dist = match_dist_list[i]
#         match_str_startpos = (a.get_seg_ts()[pos])[0]
#         match_str_endpos = (a.get_seg_ts()[pos+len_matched_word-1])[1]
#         print "Matching " + str(i+1) + ":\t" + a.get_word()[pos:pos+len_matched_word-1] + \
#                 "[" + str(match_str_startpos) + ":" + str(match_str_endpos) + "] dist=" + str(dist) 
#   
#     # plot matching results   
#     fig = plt.figure()
#     ax = plt.subplot(n+1,1,1)
#     ax.plot(piece_range, b.get_smooth_ts(), linewidth=1, color='r')
#     ax.set_title("keyword time series: " + b.get_word())
#     ts = a.get_smooth_ts()
#     segs = a.get_seg_ts()
#     title_text = []
#     for i in range(n):
#         pos = match_pos_list[i]
#         dist = match_dist_list[i]
#         match_str_startpos = (segs[pos])[0]
#         match_str_endpos = (segs[pos+len_matched_word])[1]
#           
#         ax = plt.subplot(n+1,1,i+2)
#         ax.plot(range(match_str_startpos, match_str_endpos), ts[match_str_startpos:match_str_endpos], linewidth=1, color='b')
#         ax.set_title("Matching " + str(i) + ": " + a.get_word()[pos:pos+len_matched_word] + ", dist=" + str(dist))
#               
#     fig.subplots_adjust(hspace=1)
#     plt.show(block=True)     

if __name__ == '__main__':
    idijktest()
    