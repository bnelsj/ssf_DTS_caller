import sys, os, tables, pysam, glob
from optparse import OptionParser
from collections import defaultdict

import numpy as np
import Bio.Statistics.lowess as biostats
import ml_get_cpn as ml_methods

from wssd_pw_common import *
from kitz_wssd.wssd_common_v2 import *
from math import *
from sys import stderr

from wnd_cp_data import wnd_cp_indiv, dCGH
from ssf_caller import ssf_caller, callset
from plot import line_plot, plot

from null_distribution import null_distribution

from get_windowed_variance import get_windowed_variance
import pdb


if __name__=='__main__':

    opts = OptionParser()
    opts.add_option('','--in_DTSs',dest='fn_in_DTS')
    opts.add_option('','--ref_DTS',dest='fn_ref_DTS')
    opts.add_option('','--in_prefix',dest='in_prefix')
    opts.add_option('','--in_chrs',dest='in_chrs')
    opts.add_option('','--contigs',dest='fn_contigs')
    opts.add_option('','--GC_contents_dir',dest='GC_contents_dir')
    opts.add_option('','--gene_tabix',dest='fn_gene_tabix')
    opts.add_option('','--GC_tabix',dest='fn_GC_tabix')
    opts.add_option('','--limit_to_wnds',dest='limit_to_wnds',)
    opts.add_option('','--out_viz_dir',dest='out_viz_dir')
    opts.add_option('','--out_call_dir',dest='out_call_dir')
    #opts.add_option('','--out_DNA_copy_dir',dest='out_dna_copy_dir')
    opts.add_option('','--output',dest='fn_output')
    opts.add_option('','--chunk_len',dest='chunk_len',type=int,default=1000000)
    opts.add_option('','--out_append',dest='out_append',default='')
    opts.add_option('','--binomial_smooth_k_times',dest='bin_smooth_k',
                    default=0,
                    type=int)
    opts.add_option('','--color_file',dest='color_file',
                    default="/net/gs/vol1/home/psudmant/EEE_Lab/1000G/\
                            analysis_files/colors/colors")
    #opts.add_option('','--k_mult_transitions',dest='k_mult_transitions',
    #                default=5.0,type=float)
    opts.add_option('','--cutoff_scale',dest='cutoff_scale',type=float)
    opts.add_option('','--max_merge_dif',dest='max_merge_dif',type=float)
    opts.add_option('','--merge_type',dest='merge_type',default='hierarchical')
    opts.add_option('','--plot_lims',dest='plot_lims',default=None)
    opts.add_option('','--assess_segmentation_boundaries',dest='assess_boundaries',
                    default=None)
    opts.add_option('','--window_size',dest='window_size',default=None)
    opts.add_option('','--do_plot',dest='do_plot',default=False,
                    action='store_true')
    opts.add_option('','--GC_lambda',dest='GC_lambda',
                    default="lambda gc_by_seg:np.array(gc_by_seg['max'])>1.0")
    opts.add_option('','--mask_track',dest='mask_track')
    opts.add_option('','--mask_contigs',dest='mask_contigs')
    opts.add_option('','--do_smooth_GC_biased_loci',dest='fn_smooth_GC_biased_loci',
                    default=None)
    opts.add_option('','--gaps',dest='fn_gaps',default=None)
    opts.add_option('','--segdups',dest='fn_dups',default=None)
    opts.add_option('','--superdup_track',dest='fn_superdups',default=None)
    opts.add_option('','--fdr',dest='fdr',default=0.00,type=float)
    opts.add_option('','--auto_detect_SDs',dest='auto_detect_SDs',default=False, action='store_true')

    (o, args) = opts.parse_args()
        
    GC_filt_lambda=eval(o.GC_lambda)
     
    chunk_len=o.chunk_len
    g_to_group = lambda  g: g.split("-")[0]
    
    bp_cutoff_scale = o.cutoff_scale
    rand_indiv = None
     
    tbx_gaps = pysam.Tabixfile(o.fn_gaps)
    tbx_dups = pysam.Tabixfile(o.fn_dups)
    tbx_gc   = pysam.Tabixfile(o.fn_GC_tabix)
    
    chrs = o.in_chrs.rstrip(":").split(":")
    #used to be only_wnds
    
    fn_contigs = o.fn_contigs
    wnd =  int(o.window_size)

    ##
    cutoff_scale=float(bp_cutoff_scale/wnd)
    max_merge = o.max_merge_dif 
    cp_data = dCGH(o.fn_in_DTS,o.fn_ref_DTS,o.fn_contigs,wnd)    
    
    null_dist = null_distribution()

    segment_callset = callset()
    
    caller_by_chr = {}
    
    for chr in chrs:
        print >>stderr,"%s..."%chr
        magnitude_vect = cp_data.get_cps_by_chr(chr)
        starts_vect,ends_vect = cp_data.get_wnds_by_chr(chr) 
        
        gapped_wnds = cp_data.get_overlapping_wnds(chr,tbx_gaps) 
        segdup_wnds = cp_data.get_overlapping_wnds(chr,tbx_dups)
        null_filter_wnds = [gapped_wnds, segdup_wnds]
        
        if o.auto_detect_SDs:
            cp_dup_wnds = cp_data.get_cp_dup_loci(chr)
            null_filter_wnds+=cp_dup_wnds

        null_dist.add(magnitude_vect, null_filter_wnds)
        
        print>>stderr, "cutoff_scale:%f"%( cutoff_scale )
        caller = ssf_caller(chr,magnitude_vect,starts_vect,ends_vect,cutoff_scale,use_means=True,max_merge=max_merge,scale_width=1)
        caller_by_chr[chr] = caller
        
        #plotter = line_plot(chr,caller,o.fn_gene_tabix,o.fn_dup_tabix,o.plot_lims,o.out_viz_dir) 
        #segment_callset_chr[chr] = 
        segment_callset += caller.get_callset([tbx_gaps]) 
    
    segment_callset.get_p_values(null_dist)
   
    if o.do_plot: 
        for chr, caller in caller_by_chr.iteritems():
            plotter = plot(chr,caller,o.fn_gene_tabix, o.fn_superdups ,o.plot_lims,o.out_viz_dir,segment_callset) 
            #plotter.plot_all(chunk_len=1000000,bp_start=150000000,add_heatmap=False)
            plotter.plot_all(chunk_len=1000000,bp_start=0,add_heatmap=False)

    fn_out = "%s"%(o.fn_output)
    print fn_out
    segment_callset.output(fn_out,t_stats=False)

