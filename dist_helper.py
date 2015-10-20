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
# dist_helper: this module provide some distance measurement function and matching function
#
# History
#    2014-01-07: create the file
#    
# TRICK labels: show some parts which are not easy to understand
#
# The code follows the code style PEP 8
#--------------------------------------------------------

import operator
import numpy as np
import sys
from ssstsr_publib import *

# lookup table: from number to symbol
num_to_sym_tbl = {0:'a', 1:'b', 2:'e', 3:'c', 4:'f', 5:'d', 6:'g'}
# inverse lookup table: from symbol to number
sym_to_num_tbl = {'a':0, 'b':1, 'e':2, 'c':3, 'f':4, 'd':5, 'g':6}

#---------------------------------------------------
# build distance table: the distance between different symbols
# use dictionary since its faster ( O(1) when lookups)
#--------------------------------------------------- 

def get_distance_table(n):
    x=np.linspace(0,1,n)
     
    shape_abs_val = [shape_flat, shape_linear, shape_linear, shape_leftparab, shape_leftparab, shape_rightparab, shape_rightparab]
    shape_sgn = np.array([1,1,-1,1,-1,1,-1])
     
    shapes = []
    for i in range(len(shape_abs_val)):
        shapes.append((shape_abs_val[i])(x) * shape_sgn[i])
             
    distMat = {}
    for i in range(len(shape_abs_val)):
        sym1 = num_to_sym_tbl[i]
        for j in range(len(shape_abs_val)):
            sym2 = num_to_sym_tbl[j]
            if i == j:
                distMat[(sym1,sym2)] = float(0)
            else:
                # Euclidean distance of n points
                # distMat[(sym1,sym2)] = np.sqrt(np.sum(np.power((shapes[i] - shapes[j]), 2))) / n
                # Manhattan distance
                distMat[(sym1,sym2)] = np.sum(np.abs(shapes[i] - shapes[j])) / n
    
    return distMat
    pass

#---------------------------------------------------
# distance_table version A: calculate the real 4 level distances between standard shapes
#    0:         distance between exactly the same shape
#    0.0079728: different shape, but with the same monotonous trend; both increase or decrease
#    0.0316069: different shape, one shape is flat
#    0.0612027: different shape, opposite monotonous trend; one increase while another decrease
#--------------------------------------------------- 
distance_table = get_distance_table(50000)

def print_distance_table():
    for i in range(len(num_to_sym_tbl)):
        sym1 = num_to_sym_tbl[i]
        for j in range(len(num_to_sym_tbl)):
            sym2 = num_to_sym_tbl[j]
            sys.stdout.write("(" + sym1 + "," + sym2 + "):" +  str(distance_table[(sym1,sym2)]) + "\t")
        print ""
        #"{0:.2f}".format(13.949999999999999)

#---------------------------------------------------
# distance_table version B: assign four level distances (0,1,2,3) to different shape
#    0: distance between exactly the same shape
#    1: different shape, but with the same monotonous trend; both increase or decrease
#    2: different shape, one shape is flat
#    3: different shape, opposite monotonous trend; one increase while another decrease
#--------------------------------------------------- 

# use number as symbolic representation
# distance_table = { (0,0):0, (0,1):2, (0,2):2, (0,3):2, (0,4):2, (0,5):2, (0,6):2,
#                    (1,0):2, (1,1):0, (1,2):3, (1,3):1, (1,4):3, (1,5):1, (1,6):3,
#                    (2,0):2, (2,1):3, (2,2):0, (2,3):3, (2,4):1, (2,5):3, (2,6):1,
#                    (3,0):2, (3,1):1, (3,2):3, (3,3):0, (3,4):3, (3,5):1, (3,6):3,
#                    (4,0):2, (4,1):3, (4,2):1, (4,3):3, (4,4):0, (4,5):3, (4,6):1,
#                    (5,0):2, (5,1):1, (5,2):3, (5,3):1, (5,4):3, (5,5):0, (5,6):3,
#                    (6,0):2, (6,1):3, (6,2):1, (6,3):3, (6,4):1, (6,5):3, (6,6):0}

# use string as symbolic representation
# distance_table = { ('a','a'):0, ('a','b'):2, ('a','e'):2, ('a','c'):2, ('a','f'):2, ('a','d'):2, ('a','g'):2,
#                    ('b','a'):2, ('b','b'):0, ('b','e'):3, ('b','c'):1, ('b','f'):3, ('b','d'):1, ('b','g'):3,
#                    ('e','a'):2, ('e','b'):3, ('e','e'):0, ('e','c'):3, ('e','f'):1, ('e','d'):3, ('e','g'):1,
#                    ('c','a'):2, ('c','b'):1, ('c','e'):3, ('c','c'):0, ('c','f'):3, ('c','d'):1, ('c','g'):3,
#                    ('f','a'):2, ('f','b'):3, ('f','e'):1, ('f','c'):3, ('f','f'):0, ('f','d'):3, ('f','g'):1,
#                    ('d','a'):2, ('d','b'):1, ('d','e'):3, ('d','c'):1, ('d','f'):3, ('d','d'):0, ('d','g'):3,
#                    ('g','a'):2, ('g','b'):3, ('g','e'):1, ('g','c'):3, ('g','f'):1, ('g','d'):3, ('g','g'):0}

sym_avg_dist = sum(v for (_, v) in distance_table.items()) / float(len(distance_table))


#---------------------------------
# calculating the levenstein distance of two string
#     input:
#        @ two strings
#     output:
#        @ levenstein distance with float type
#     example:
#         >>> edit_dist('eeba', 'abac')
#         3
#         >>> edit_dist('abc', 'cba')
#         2
#         >>> edit_dist('cbc', 'eba')
#         2
#         >>> edit_dist('recoginze', 'recognize')
#         2
#         >>> edit_dist('sailn', 'failing')
#         3
#---------------------------------    
def edit_dist(s1, s2):
    len1 = len(s1)
    len2 = len(s2) 
     
    # for all i and j, d[i,j] will hold the Levenshtein distance between
    # the first i characters of s and the first j characters of t;
    # note that d has (m+1)*(n+1) values
    matrix = [[i+j for j in range(len2 + 1)] for i in range(len1 + 1)]
     
    for row in range(len1):
        for col in range(len2):
                substitute_cost = (s1[row] != s2[col])
                matrix[row+1][col+1] = min(matrix[row+1][col]+1, # delete
                                           matrix[row][col+1]+1, # insert
                                           matrix[row][col]+substitute_cost)   # substitution
                         
    return matrix[len1][len2]


#---------------------------------
# Needleman-Wunch distance: the cost of substitute are arbitrary (cost table in this function)
#                           the insert and delete cost are the same, also called "gap cost"
# this function calculating the Needleman-Wunch distance of two string
#     gap cost is sym_avg_dist, and substitution cost is from distance_table
#     input:
#        @ normalized: normalized by the length of the string if normalized is true
#     output:
#        @ Needleman-Wunch distance with float type
#---------------------------------    
def nw_dist(s1, s2):
    len1 = len(s1)
    len2 = len(s2) 
     
    matrix = [[(i+j)*sym_avg_dist for j in range(len2 + 1)] for i in range(len1 + 1)]
     
    for row in range(len1):
        for col in range(len2):
            if s1[row] == s2[col]:
                # if new characters are the same
                matrix[row+1][col+1] = min(matrix[row+1][col]+sym_avg_dist, # delete
                                           matrix[row][col+1]+sym_avg_dist, # insert
                                           matrix[row][col])   # substitution
            else:
                matrix[row+1][col+1] = min(matrix[row+1][col]+sym_avg_dist, # delete
                                           matrix[row][col+1]+sym_avg_dist, # insert
                                           matrix[row][col]+distance_table[(s1[row], s2[col])])   # substitution
                         
    return matrix[len1][len2]


#---------------------------------
# leven_match: approximate substring matching
#     input:
#        @ 
#     output:
#        @ the position of first matching substring, if no matching, return -1
#---------------------------------    
def leven_match(keyword, word, best_n=1):
    len_keyword = len(keyword)
    len_word = len(word) 

    if best_n <= 0 or best_n > len_word-len_keyword+1:
        raise ValueError, "best_n must be greater or equal than 1 and less equal than maximum number of possible substrings!"
    
    # initial value of edit distance 
    matrix = [[0 for j in range(len_word + 1)] for i in range(len_keyword + 1)]
    for j in range(len_keyword + 1):
        matrix[j][0] = j
    
    # the direction of backward searching
    dir = [[0 for j in range(len_word + 1)] for i in range(len_keyword + 1)]
     
    for row in range(len_keyword):
        for col in range(len_word):
#             if keyword[row] == word[col]:
#                 # if new characters are the same
#                 complist = np.array([
#                             matrix[row][col],  # substitution                                     
#                             matrix[row+1][col]+1,  # delete
#                             matrix[row][col+1]+1  # insert
#                             ])  
#             else:
#                 complist = np.array([
#                             matrix[row][col]+1,  # substitution, go diagonal                                   
#                             matrix[row+1][col]+1,  # delete, go right
#                             matrix[row][col+1]+1  # insert, go down
#                             ])  
                
            # TRICK: the cost of substitution is on the first place of the list, because
            #        when backward search, we prefer the substring has equal length
            if keyword[row] == word[col]:
                # if new characters are the same
                complist = np.array([
                            matrix[row][col],  # substitution                                     
                            matrix[row+1][col]+sym_avg_dist,  # delete
                            matrix[row][col+1]+sym_avg_dist  # insert
                            ])  
            else:                
                complist = np.array([
                            matrix[row][col]+distance_table[(keyword[row], word[col])],  # substitution                                     
                            matrix[row+1][col]+sym_avg_dist,  # delete
                            matrix[row][col+1]+sym_avg_dist  # insert
                            ])                  
     
                          
            matrix[row+1][col+1] = np.min(complist)
            dir[row+1][col+1] = np.argmin(complist) + 1
                
    # get the best_n minimal distances
    lenvdist_list = matrix[len_keyword]
    matching_words = []
    min_dist_list = np.argsort(lenvdist_list)[0:best_n]
    row_end = len_keyword
    action = []
    matching_pos = []
    
    for col_end in min_dist_list:
        search_record = []
        # backward search three direction
        col = col_end
        row = row_end
        while row != 0:
            search_record.append(dir[row][col])
            # reverse operation of substitute, delete and insert 
            if dir[row][col] == 1:
                # substitution, diagonal
                col = col - 1
                row = row - 1
            elif dir[row][col] == 2:
                # delete, go left
                col = col - 1
            elif dir[row][col] == 3:
                # insert, go up 
                row = row - 1
            else:
                raise ValueError, "it is odd: the backward search direction must be in [1,2,3]"
                
        # TRICK: pay attention to the index of the word
        matching_words.append(word[col:col_end])
        matching_pos.append(col)
        
        # editing action
        action_dict = {1:'s', 2:'i', 3:'d'}
        action.append("".join([action_dict[i] for i in search_record[::-1]]))
        
    return matching_pos, matching_words, [lenvdist_list[i] for i in min_dist_list], action

#---------------------------------
# calculating the hamming distance of two string
#     input:
#        @ normalized: normalized by the length of the string if normalized is true
#     output:
#        @ hamming distance with float type
#---------------------------------    
def hamming_dist(str1, str2, normalized=True):
    len_str1 = len(str1)
    len_str2 = len(str2) 
    
    if len_str1 != len_str2:
        raise ValueError, "the length of two strings must be the same!"
    
    dist = 0
    for i in range(len_str1):
        try:
            str_dist = distance_table[(str1[i], str2[i])]
        except KeyError:
            #TODO: use logging
            # print "invalid shape symbol (must be in a-g):\t" + str1[i] + "\t" + str2[i]
            return sys.float_info.max
            
        dist = dist + str_dist
    
    if normalized:
        dist = dist / float(len_str1)
    return dist

#---------------------------------
# match: looking up the first matching (distance less equal than the threshold) substring based on the keyword
#     input:
#        @ dist_thresh: default value is 0
#     output:
#        @ the position of first matching substring, if no matching, return -1
#---------------------------------    
def hamming_match(keyword, word, dist_thresh=0.0):
    len_keyword = len(keyword)
    len_str = len(word) 
    
    if len_keyword > len_str:
        raise ValueError, "the length of the keyword should be less than the string!"
    
    for i in range(len_str-len_keyword+1):
        dist = hamming_dist(keyword, word[i:i+len_keyword])
        if hamming_dist(keyword, word[i:i+len_keyword]) <= dist_thresh:
            return i
        # elif dist == sys.float_info.max:
        #    print "invalid shape symbol (must be in a-g) in the position:\t" + str(i)
            
    # no matching
    return -1

#---------------------------------
# match: looking up all matching (distance less equal than the threshold) substrings based on the keyword
#    output: a list of substring index
#---------------------------------
def hamming_match_all(keyword, word, dist_thresh=0.0):
    len_keyword = len(keyword)
    len_str = len(word) 
    
    if len_keyword > len_str:
        raise ValueError, "the length of the keyword should be less than the string!"
    
    matching_pos = []
    for i in range(len_str-len_keyword+1):
        dist = hamming_dist(keyword, word[i:i+len_keyword])
        if dist <= dist_thresh:
            matching_pos.append(i) 
        # elif dist == sys.float_info.max:
        #    print "invalid shape symbol (must be in a-g) in the position:\t" + str(i)             
    
    return matching_pos

#---------------------------------
# match_best_effort: looking up the best n matching substrings based on the keyword
#    output: a list of substring index
#---------------------------------
def hamming_match_best_effort(keyword, word, best_n=1):
    
    len_keyword = len(keyword)
    len_str = len(word) 
    
    if len_keyword > len_str:
        raise ValueError, "the length of the keyword must be less than the string!"

    if best_n <= 0 or best_n > len_str-len_keyword+1:
        raise ValueError, "best_n must be greater or equal than 1 and less equal than maximum number of possible substrings!"
    
    # store all distance
    distlist = []
    for i in range(len_str-len_keyword+1):
        distlist.append(hamming_dist(keyword, word[i:i+len_keyword]))
    
    # get sorted index
    vals = np.array(distlist)    
    matching_pos = np.argsort(vals)[0:best_n]

    # return best_n matching result: positions and distances
    return matching_pos, [distlist[i] for i in matching_pos]

