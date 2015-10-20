#########################################################
#                                                         
# Segmented Shape-Symbolic Time series Representation    
#                                                         #
# __author__ = "Zhe Sun"
# __copyright__ = "Copyright 2013, Target-holding B.V."
# __license__ = "FreeBSD"
# __email__ = "zhe.sun@target-holding.com"
#
#########################################################

This package implements Segmented Shape-Symbolic Time series Representation method 
(different from the original algorithm in future).

File list:
	|
	|---- ssstsr.py: 		the class of Segmented Shape-Symbolic Time series Representation
	|---- ssstsr_public.py:	provide some public functions, normalize, smooth, basic shape, etc 
	|---- bottomup.py: 		provide bottom-up segments merging functions
	|---- uniform_gap_helper.py:	provide functions which calculate gap statistic error of 
	|								the uniformsignals 
	|---- dist_helper.py:	provide different distance metrics and matching 
	|						functions. Support self-define distance table, hamming distance so far 
	|---- fitting.py:		provide fitting and symbolic functions, which fit the time series to 7 basic shapes. 
	|						They can calculate residuals and estimate the maximum likelihood shape
	|---- plot_helper.py:	provide some plotting helper functions
	|---- sim_sig_gen.py:	generate random signal which composed by monotonous pieces
	|---- integration_test.py:	integration test functions 	
	|---- example.py:		a example of using livedijk data, segment, symbolic representation and matching
	

TODO:
	1. Lenvenstein distance --- Done
	2. uniform error cache store and retrieve
	3. minSize and point of inflection
	4. unit test and integrate test?
	