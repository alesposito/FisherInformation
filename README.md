# FisherInformation

# Numerical evaluation of Fisher Information Matrix

FPT files are related to our publication (MATLAB, codes working since Matlab 2013, last tested with Matlab 2018a)
https://doi.org/10.1371/journal.pone.0077392
"Maximizing the Biochemical Resolving Power of Fluorescence Microscopy"
Alessandro Esposito , Marina Popleteeva, Ashok R. Venkitaraman
2013 PLOS ONE

the main file is fpt_optimize_time_gates.m
which handles the execution of code and the visualization of the results
fpt_fvalue.m computes the Fisher Information (I) and the figure of merit (F) for the evaluation of a lifetime value on a borad ranges of lifetimes
fpt_tg_bu.m partitions the photon arrival statistics onto an histogram of higher resolutions to identify the optimal partition (bottom-up method). Each bin is plit into two symmetric bins until there is a significant gain of Fisher information. 
fpt_tg_fcost.m cost function used to optimize the edges of histrograms of even bins

fpt_optimize_time_gates.m computes F and I for a range of fluorescence lifetimes using:
1) Use fpt_tg_bu.m to identify an optimal partition (optim_par and FO). The bins are of uneven width but they are all multiples of the smallest bin width 
2) Use the in-built function fminsearch and the bespoke fpt_tg_fcost.m to optimize the boundaries of optim_par, to obtain a refine_par with Fisher information FF.
3) Evaluate the Fisher information (FE) for an a histrogam with bins of even size (even_par). The number of bins is the same identified with step #1
4) Evaluate the Fisher information (FR) for a dense reference histrogram of equal bins (dense_par, e.g. n=256) as for TCSPC

fpt_optimize_time_gates.m then visualize the F-value curves and the partitions resulting from these optimization processes
