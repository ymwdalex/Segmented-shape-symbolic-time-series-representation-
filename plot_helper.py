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
# plot_helper: this module provide some plotting functions
#
# History
#    2014-01-07: create the file (refactoring from segmentclass.py)
#    
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

from matplotlib import pyplot as plt
from fitting import fit_shape

#---------------------------------------------
# plot segments and least square sum regression segment
#    input:
#    @ ts: time series
#    @ seg_ts: segments of time series, list of tuple, i.e, [(0,5),(5,20)...]
#    @ seg_fit: the fitting shape of segment, list of list
#    @ shape_symbol_list: the symbols of each segment, list of string
#    @ imshow: if the figure show immediately. When plot the animation of bottomup,
#              imshow=False, otherwise imshow=True
#----------------------------------------------

# TRICK: declare fig. matplotlib.pyplot is a state
def plot_segts_fit(ts, seg_ts, seg_fit = None, imshow=False, shape_symbol_list=None):
    plt.clf()

    # iterate each segment and plot
    for i in range(len(seg_ts)):
        seg = seg_ts[i]
        startpoint = seg[0]
        endpoint = seg[1]

        plt.plot(range(startpoint, endpoint), ts[startpoint:endpoint], linewidth=1, color='b')
        if seg_fit is None:
            plt.plot(range(startpoint, endpoint), fit_shape(ts[startpoint:endpoint])[0], linewidth=1, color='r')
        else:
            plt.plot(range(startpoint, endpoint), seg_fit[i], linewidth=1, color='r')

    # plot vertical line to divide segments, just draw the line in the middle of two segments
    for seg in seg_ts[:-1]:
        # seg[1] is the exclusive(!!!) end point 
        seg_boundary = seg[1] - 0.5 
        plt.axvline(seg_boundary, linewidth=2, color='k')

    titlestr = "#Segments={k}".format(k=str(len(seg_ts)))

    # add the text of shape symbol as well
    if shape_symbol_list:
        ymin, ymax = plt.ylim()
        ypos = ymin + (ymax-ymin) * 0.9
        for i in range(len(seg_ts)):
            seg = seg_ts[i]
            # text_ypos  = (seg[0] + seg[1])/2.005
            text_ypos = seg[0] + 1
            plt.text(text_ypos, ypos, shape_symbol_list[i], fontsize=14, color='green')

        # display the word in the figure
        titlestr = titlestr + ", symbolic representation = \"" + "".join(shape_symbol_list) + "\""

    plt.title(titlestr, fontsize=15)
    
    # TRICK: control if show the figure immediately
    if imshow:
        plt.show()