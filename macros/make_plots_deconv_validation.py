#! /usr/bin/env python
from __future__ import print_function
import os
import sys
import re
#import pdb
import numpy as np
from ROOT import TFile
from eophelper.plotting import *
from eophelper.util import *

def main(argv):

    # Config - modify based on local setup
    # ---------------------------------
    tag = "deconv_inputs_20171005" # file tag
    eop_dir = "/Users/joakim/eoverp/"
    input_dir = os.path.join(eop_dir, "results/eoverp_cross_checks/eop_files/")
    output_dir = os.path.join(input_dir, tag+"_validation_plots/")
    # Options
    do_lcw = True
    do_cells = False
    # ---------------------------------

    ensure_dir(output_dir)

    input_file = TFile.Open(input_dir+"/"+tag+".root", "READ")
    folders = [input_file.Get("ValidationHistograms_EM")]
    if do_lcw:
        folders.append(input_file.Get("ValidationHistograms_LCW"))
    if do_cells:
        folders.append(input_file.Get("ValidationHistograms_Cells"))

    #_ Default setup (low-mu samples)
    file_tag = "overlay_lowmu_data_mc"
    leg_labels = ["Data 2015, Low-#LT#font[50]{#mu}#GT", "Pythia MinBias"]
    lumi_label = "#sqrt{s} = 13 TeV, 1.6 nb^{-1}"
    atlas_label = "Internal"
    #atlas_label = "Work-in-progress"
    options = "deconv_"
    do_data = True # plot data as solid black markers
    do_data_mc = True # plot data_mc_ratio
    invert_data_mc_ratio = True # invert data mc ratio order
    do_ratio = True # add a ratio plot under main plot
    do_norm = False # normalize all plots to unity
    do_norm_mc_to_data = False #
    do_auto_scale = False # scale y-axis appropriately
    draw_style = "hist][" # pass args to h.Draw("args"), NB: not applicable for do_data
    draw_avg = False
    plot_tags = []
    marker_type = "data_mc" #sets the markers and colors

    # plot deconvolution inputs
    hlist_deconv = ["DeconvolutionInputs_EM/h2_EM_diff_zfData_zfMC",
                    "DeconvolutionInputs_EM/h2_EM_avgeop_cor_data_over_mc"]
    if do_lcw:
        hlist_deconv.append("DeconvolutionInputs_LCW/h2_LCW_diff_zfData_zfMC")
        hlist_deconv.append("DeconvolutionInputs_LCW/h2_LCW_avgeop_cor_data_over_mc")

    for hname in hlist_deconv:
        h = input_file.Get(hname)
        plot_name = re.sub("DeconvolutionInputs_.+/", "", hname)
        plot_opts = PlotOptions(plot_name,
                                options,
                                draw_style='E',
                                )

        plot_output_path = output_dir+"/"+plot_name
        save_hist(h, plot_opts, plot_output_path)

    # plot validation histograms
    for folder in folders:
        for key in folder.GetListOfKeys():

            if "_zf_ratio" in key.GetName():
                plot_name = key.GetName()
                plot_opts = PlotOptions(plot_name,
                                        options=options,
                                        draw_style='E',
                                        )

                plot_output_path = output_dir+"/"+file_tag+"_"+plot_name
                h = folder.Get(key.GetName())
                save_hist(h, plot_opts, plot_output_path)

            classname = key.GetClassName()
            cl = gROOT.GetClass(classname)
            if not cl:
                continue
            if not "TH1" in str(cl):
                continue
            if "data" in str(key):
                hdata_name = key.GetName()
                hdata = folder.Get(hdata_name)
                hmc_name = re.sub("data", "mc", hdata_name)
                hmc = folder.Get(hmc_name)

                plot_name = re.sub("_data", "", hdata_name)

                plot_opts = PlotOptions(plot_name+"Run1paper",
                                        options,
                                        do_data,
                                        do_data_mc,
                                        invert_data_mc_ratio,
                                        do_ratio,
                                        do_norm,
                                        do_norm_mc_to_data,
                                        do_auto_scale,
                                        draw_style,
                                        marker_type,
                                        leg_labels,
                                        plot_tags,
                                        atlas_label,
                                        lumi_label,
                                        draw_avg)

                plot_output_path = output_dir+"/"+file_tag+"_"+plot_name
                save_hists_overlay([hdata, hmc], plot_opts, plot_output_path)

if __name__ == "__main__":
    main(sys.argv)
