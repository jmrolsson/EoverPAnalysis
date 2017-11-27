#! /usr/bin/env python
from __future__ import print_function
import os
import sys
import re
#import pdb
import numpy as np
from array import array
from ROOT import TFile, TH1D, TH2F
from eophelper.util import *

def main(argv):

    # Config - modify based on local setup
    # ---------------------------------
    tag = "deconv_inputs_20171005" # file tag
    eop_dir = "/Users/joakim/eoverp/"
    input_dir = os.path.join(eop_dir, "results/eoverp_cross_checks/eop_files/EoverPAnalysis_outputs")
    output_dir = os.path.join(eop_dir, "results/eoverp_cross_checks/eop_files/")
    # Input data and MC files
    data = os.path.join(input_dir, "condor_all_eop_lowmu_runII_general_data_20170530_deconvinputs.root")
    mc = os.path.join(input_dir, "condor_all_eop_lowmu_runII_general_mc_merged_20170530_deconvinputs.root")
    # Options
    write_raw_hists = True
    do_lcw = True
    do_cells = False
    # ---------------------------------

    print("Input data-file:", data)
    print("Input MC-file:", mc)

    datafile = TFile.Open(data, "READ")
    mcfile = TFile.Open(mc, "READ")

    ensure_dir(output_dir)
    output = os.path.join(output_dir, tag+".root")
    outputfile = TFile.Open(output, "RECREATE")
    outputfile.cd()

    p_bins = [.5, .8, 1.2, 1.8, 2.2, 2.8, 3.4, 4.2, 5., 6., 7., 9., 12., 15., 20., 30.]
    eta_bins = [0., .6, 1.1, 1.4, 1.5, 1.8, 1.9, 2.3]

    selection_tags = ["EM"]
    energy_calib = ["ClusterEnergy"]
    if do_lcw:
        selection_tags.append("LCW")
        energy_calib.append("ClusterEnergyLCW")
    if do_cells:
        selection_tags.append("Cells")
        energy_calib.append("CellEnergy")

    # Save numbers for cross-checks
    n_trks_data_numbers = np.zeros( (len(selection_tags), len(p_bins),len(eta_bins)) )
    n_trks_mc_numbers = np.zeros( (len(selection_tags), len(p_bins),len(eta_bins)) )
    n_trks_leq0_data_numbers = np.zeros( (len(selection_tags), len(p_bins),len(eta_bins)) )
    n_trks_leq0_mc_numbers = np.zeros( (len(selection_tags), len(p_bins),len(eta_bins)) )
    zf_numbers = np.zeros( (len(selection_tags), len(p_bins),len(eta_bins)) )
    eop_cor_numbers = np.zeros( (len(selection_tags), len(p_bins), len(eta_bins)) )
    err_eop_cor_numbers = np.zeros( (len(selection_tags), len(p_bins), len(eta_bins)) )

    for s in xrange(len(energy_calib)):

        selection = "EoverP_LoosePrimaryTrks_"+energy_calib[s]+"_Run1paper"

        raw_dir = outputfile.mkdir("RawHistograms_"+selection_tags[s])
        valid_dir = outputfile.mkdir("ValidationHistograms_"+selection_tags[s])
        deconv_dir = outputfile.mkdir("DeconvolutionInputs_"+selection_tags[s])

        # Zero fractions

        # Input histos for zero fraction calculations
        h2_N_data = datafile.Get(selection+"/trk_n_E_200")
        h2_N_l0_data = datafile.Get(selection+"/trk_n_E_200_l0")
        h2_N_eq0_data = datafile.Get(selection+"/trk_n_E_200_eq0")
        h2_N_leq0_data = datafile.Get(selection+"/trk_n_E_200_leq0")

        if write_raw_hists:
            raw_dir.cd()
            append_name(h2_N_data, "_data")
            append_name(h2_N_l0_data, "_data")
            append_name(h2_N_eq0_data, "_data")
            append_name(h2_N_leq0_data, "_data")
            h2_N_data.Write()
            h2_N_l0_data.Write()
            h2_N_eq0_data.Write()
            h2_N_leq0_data.Write()

        h2_N_mc = mcfile.Get(selection+"/trk_n_E_200")
        h2_N_l0_mc = mcfile.Get(selection+"/trk_n_E_200_l0")
        h2_N_eq0_mc = mcfile.Get(selection+"/trk_n_E_200_eq0")
        h2_N_leq0_mc = mcfile.Get(selection+"/trk_n_E_200_leq0")

        if write_raw_hists:
            raw_dir.cd()
            append_name(h2_N_mc, "_mc")
            append_name(h2_N_l0_mc, "_mc")
            append_name(h2_N_eq0_mc, "_mc")
            append_name(h2_N_leq0_mc, "_mc")
            h2_N_mc.Write()
            h2_N_l0_mc.Write()
            h2_N_eq0_mc.Write()
            h2_N_leq0_mc.Write()

        # Output histos for zero fraction calculations
        h2_N_leq0_over_N_data = TH2F("h2_"+selection_tags[s]+"_N_leq0_over_N_data", "h2_N_leq0_over_N_data", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_N_leq0_over_N_data.GetXaxis().SetTitle("p [GeV]")
        h2_N_leq0_over_N_data.GetYaxis().SetTitle("|#eta|")
        h2_N_leq0_over_N_data.GetZaxis().SetTitle("N_{Data}(E#leq0)/N_{Data}")
        h2_N_leq0_over_N_mc = TH2F("h2_"+selection_tags[s]+"_N_leq0_over_N_mc", "h2_N_leq0_over_N_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_N_leq0_over_N_mc.GetXaxis().SetTitle("p [GeV]")
        h2_N_leq0_over_N_mc.GetYaxis().SetTitle("|#eta|")
        h2_N_leq0_over_N_mc.GetZaxis().SetTitle("N_{MC}(E#leq0)/N_{MC}")
        h2_diff_zfData_zfMC = TH2F("h2_"+selection_tags[s]+"_diff_zfData_zfMC", "h2_diff_zfData_zfMC", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_diff_zfData_zfMC.GetXaxis().SetTitle("p [GeV]")
        h2_diff_zfData_zfMC.GetYaxis().SetTitle("|#eta|")
        h2_diff_zfData_zfMC.GetZaxis().SetTitle("(N(E#leq0)/N)_{Data} - (N(E#leq0)/N)_{MC}")

        h2_N_l0_over_N_eq0_data = TH2F("h2_"+selection_tags[s]+"_N_l0_over_N_eq0_data", "h2_N_l0_over_N_eq0_data", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_N_l0_over_N_eq0_data.GetXaxis().SetTitle("p [GeV]")
        h2_N_l0_over_N_eq0_data.GetYaxis().SetTitle("|#eta|")
        h2_N_l0_over_N_eq0_data.GetZaxis().SetTitle("N_{Data}(E#leq0)/N_{Data}")
        h2_N_l0_over_N_eq0_mc = TH2F("h2_"+selection_tags[s]+"_N_l0_over_N_eq0_mc", "h2_N_l0_over_N_eq0_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_N_l0_over_N_eq0_mc.GetXaxis().SetTitle("p [GeV]")
        h2_N_l0_over_N_eq0_mc.GetYaxis().SetTitle("|#eta|")
        h2_N_l0_over_N_eq0_mc.GetZaxis().SetTitle("N_{MC}(E#leq0)/N_{MC}")

        # Background corrected average E/p

        # Output histos for background corrected average E/p
        h2_eop_raw_data = TH2F("h2_"+selection_tags[s]+"_avgeop_raw_data", "h2_avgeop_raw_data", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_raw_data.GetXaxis().SetTitle("p [GeV]")
        h2_eop_raw_data.GetYaxis().SetTitle("|#eta|")
        h2_eop_raw_data.GetZaxis().SetTitle("<E/p>^{Data}_{RAW}")
        h2_eop_raw_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_raw_mc", "h2_avgeop_raw_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_raw_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_raw_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_raw_mc.GetZaxis().SetTitle("<E/p>^{MC}_{RAW}")
        h2_eop_raw_data_over_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_raw_data_over_mc", "h2_avgeop_raw_data_over_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_raw_data_over_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_raw_data_over_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_raw_data_over_mc.GetZaxis().SetTitle("<E/p>^{Data}_{RAW}/<E/p>^{MC}_{RAW}")

        h2_eop_bg_data = TH2F("h2_"+selection_tags[s]+"_avgeop_bg_data", "h2_avgeop_bg_data", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_bg_data.GetXaxis().SetTitle("p [GeV]")
        h2_eop_bg_data.GetYaxis().SetTitle("|#eta|")
        h2_eop_bg_data.GetZaxis().SetTitle("<E/p>^{Data}_{BG}")
        h2_eop_bg_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_bg_mc", "h2_avgeop_bg_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_bg_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_bg_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_bg_mc.GetZaxis().SetTitle("<E/p>^{MC}_{BG}")
        h2_eop_bg_data_over_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_bg_data_over_mc", "h2_avgeop_bg_data_over_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_bg_data_over_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_bg_data_over_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_bg_data_over_mc.GetZaxis().SetTitle("<E/p>^{Data}_{BG}/<E/p>^{MC}_{BG}")

        h2_eop_cor_data = TH2F("h2_"+selection_tags[s]+"_avgeop_cor_data", "h2_avgeop_cor_data", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_cor_data.GetXaxis().SetTitle("p [GeV]")
        h2_eop_cor_data.GetYaxis().SetTitle("|#eta|")
        h2_eop_cor_data.GetZaxis().SetTitle("<E/p>^{Data}_{COR}")
        h2_eop_cor_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_cor_mc", "h2_avgeop_cor_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_cor_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_cor_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_cor_mc.GetZaxis().SetTitle("<E/p>^{MC}_{COR}")
        h2_eop_cor_data_over_mc = TH2F("h2_"+selection_tags[s]+"_avgeop_cor_data_over_mc", "h2_avgeop_cor_data_over_mc", len(p_bins)-1, array('d', p_bins), len(eta_bins)-1, array('d', eta_bins))
        h2_eop_cor_data_over_mc.GetXaxis().SetTitle("p [GeV]")
        h2_eop_cor_data_over_mc.GetYaxis().SetTitle("|#eta|")
        h2_eop_cor_data_over_mc.GetZaxis().SetTitle("<E/p>^{Data}_{COR}/<E/p>^{MC}_{COR}")

        for i in xrange(1, len(p_bins)):

            p_str = 'pG{:d}L{:d}'.format((int)(p_bins[i-1]*1000), (int)(p_bins[i]*1000))

            # Raw data (TProfile)
            h_eop_raw_data = datafile.Get(selection+"/"+p_str+"_eop_Total_"+energy_calib[s]+"_0_200")
            h_eop_bg_data = datafile.Get(selection+"/"+p_str+"_eop_EM_BG_"+energy_calib[s]+"")

            # Raw MC (TProfile)
            h_eop_raw_mc = mcfile.Get(selection+"/"+p_str+"_eop_Total_"+energy_calib[s]+"_0_200")
            h_eop_bg_mc = mcfile.Get(selection+"/"+p_str+"_eop_EM_BG_"+energy_calib[s]+"")

            if write_raw_hists:
                raw_dir.cd()
                append_name(h_eop_raw_data, "_data")
                append_name(h_eop_bg_data, "_data")
                append_name(h_eop_raw_mc, "_mc")
                append_name(h_eop_bg_mc, "_mc")
                h_eop_raw_data.Write()
                h_eop_bg_data.Write()
                h_eop_raw_mc.Write()
                h_eop_bg_mc.Write()

            # histograms of ZF vs. eta in p-bins
            h_N_leq0_over_N_vs_eta_data = TH1D(selection_tags[s]+"_"+p_str+"_N_leq0_over_N_vs_eta_data", p_str+"_N_leq0_over_N_vs_eta_data", len(eta_bins)-1, array('d', eta_bins))
            h_N_leq0_over_N_vs_eta_mc = TH1D(selection_tags[s]+"_"+p_str+"_N_leq0_over_N_vs_eta_mc", p_str+"_N_leq0_over_N_vs_eta_mc", len(eta_bins)-1, array('d', eta_bins))

            h_zf_vs_eta = TH1D(selection_tags[s]+"_"+p_str+"_zf_vs_eta", p_str+"_zf_vs_eta", len(eta_bins)-1, array('d', eta_bins))
            h_zf_ratio_vs_eta = TH1D(selection_tags[s]+"_"+p_str+"_zf_ratio_vs_eta", p_str+"_zf_ratio_vs_eta", len(eta_bins)-1, array('d', eta_bins))

            # histograms of <E/p> vs. eta in p-bins
            h_eop_raw_vs_eta_data = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_raw_vs_eta_data", p_str+"_avgeop_raw_vs_eta_data", len(eta_bins)-1, array('d', eta_bins))
            h_eop_bg_vs_eta_data = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_bg_vs_eta_data", p_str+"_avgeop_bg_vs_eta_data", len(eta_bins)-1, array('d', eta_bins))
            h_eop_cor_vs_eta_data = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_cor_vs_eta_data", p_str+"_avgeop_cor_vs_eta_data", len(eta_bins)-1, array('d', eta_bins))
            h_eop_raw_vs_eta_mc = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_raw_vs_eta_mc", p_str+"_avgeop_raw_vs_eta_mc", len(eta_bins)-1, array('d', eta_bins))
            h_eop_bg_vs_eta_mc = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_bg_vs_eta_mc", p_str+"_avgeop_bg_vs_eta_mc", len(eta_bins)-1, array('d', eta_bins))
            h_eop_cor_vs_eta_mc = TH1D(selection_tags[s]+"_"+p_str+"_avgeop_cor_vs_eta_mc", p_str+"_avgeop_cor_vs_eta_mc", len(eta_bins)-1, array('d', eta_bins))

            for j in xrange(1, len(eta_bins)):
                p_eta_str = 'pG{:d}L{:d}_etaG{:d}L{:d}'.format(
                    (int)(p_bins[i-1]*1000), (int)(p_bins[i]*1000),
                    (int)(eta_bins[j-1]*10), (int)(eta_bins[j]*10))

                # 1D E/p distributions

                #Data
                h1_eop_raw_data = datafile.Get(selection+"/"+p_eta_str+"_eop_Total_"+energy_calib[s]+"_0_200")
                h1_eop_bg_data = datafile.Get(selection+"/"+p_eta_str+"_eop_EM_BG_"+energy_calib[s]+"")

                #MC
                h1_eop_raw_mc = mcfile.Get(selection+"/"+p_eta_str+"_eop_Total_"+energy_calib[s]+"_0_200")
                h1_eop_bg_mc = mcfile.Get(selection+"/"+p_eta_str+"_eop_EM_BG_"+energy_calib[s]+"")

                if write_raw_hists:
                    raw_dir.cd()
                    append_name(h1_eop_raw_data, "_data")
                    append_name(h1_eop_bg_data, "_data")
                    append_name(h1_eop_raw_mc, "_mc")
                    append_name(h1_eop_bg_mc, "_mc")
                    h1_eop_raw_data.Write()
                    h1_eop_bg_data.Write()
                    h1_eop_raw_mc.Write()
                    h1_eop_bg_mc.Write()

                # Zero fractions

                # Data
                N_data = h2_N_data.GetBinContent(i, j)
                N_l0_data = h2_N_l0_data.GetBinContent(i, j)
                N_eq0_data = h2_N_eq0_data.GetBinContent(i, j)
                N_leq0_data = h2_N_leq0_data.GetBinContent(i, j)
                n_trks_data_numbers[s, i, j] = N_data
                n_trks_leq0_data_numbers[s, i, j] = N_leq0_data

                N_leq0_over_N_data = 0
                err_N_leq0_over_N_data = 0
                if N_data > 0 and N_leq0_data > 0:
                    N_leq0_over_N_data = N_leq0_data / N_data
                    err_N_data = h2_N_data.GetBinError(i, j)
                    err_N_leq0_data = h2_N_leq0_data.GetBinError(i, j)
                    err_N_leq0_over_N_data = frac_err(N_leq0_data, N_data, err_N_leq0_data, err_N_data)
                    h2_N_leq0_over_N_data.SetBinContent(i, j, N_leq0_over_N_data)
                    h2_N_leq0_over_N_data.SetBinError(i, j, err_N_leq0_over_N_data)

                    valid_dir.cd()
                    h_N_leq0_over_N_vs_eta_data.SetBinContent(j, N_leq0_over_N_data)
                    h_N_leq0_over_N_vs_eta_data.SetBinError(j, err_N_leq0_over_N_data)

                N_l0_over_N_eq0_data = 0
                err_N_l0_over_N_eq0_data = 0
                if N_l0_data > 0 and N_eq0_data > 0:
                    N_l0_over_N_eq0_data = N_l0_data / N_eq0_data
                    err_N_l0_data = h2_N_l0_data.GetBinError(i, j)
                    err_N_eq0_data = h2_N_eq0_data.GetBinError(i, j)
                    err_N_l0_over_N_eq0_data = frac_err(N_l0_data, N_eq0_data, err_N_l0_data, err_N_eq0_data)
                    h2_N_l0_over_N_eq0_data.SetBinContent(i, j, N_l0_over_N_eq0_data)
                    h2_N_l0_over_N_eq0_data.SetBinError(i, j, err_N_l0_over_N_eq0_data)

                # MC
                N_mc = h2_N_mc.GetBinContent(i, j)
                N_l0_mc = h2_N_l0_mc.GetBinContent(i, j)
                N_eq0_mc = h2_N_eq0_mc.GetBinContent(i, j)
                N_leq0_mc = h2_N_leq0_mc.GetBinContent(i, j)
                n_trks_mc_numbers[s, i, j] = N_data
                n_trks_leq0_mc_numbers[s, i, j] = N_leq0_data

                N_leq0_over_N_mc = 0
                err_N_leq0_over_N_mc = 0
                if N_mc > 0 and N_leq0_mc > 0:
                    N_leq0_over_N_mc = N_leq0_mc / N_mc
                    err_N_mc = h2_N_mc.GetBinError(i, j)
                    err_N_leq0_mc = h2_N_leq0_mc.GetBinError(i, j)
                    err_N_leq0_over_N_mc = frac_err(N_leq0_mc, N_mc, err_N_leq0_mc, err_N_mc)
                    h2_N_leq0_over_N_mc.SetBinContent(i, j, N_leq0_over_N_mc)
                    h2_N_leq0_over_N_mc.SetBinError(i, j, err_N_leq0_over_N_mc)

                    valid_dir.cd()
                    h_N_leq0_over_N_vs_eta_mc.SetBinContent(j, N_leq0_over_N_mc)
                    h_N_leq0_over_N_vs_eta_mc.SetBinError(j, err_N_leq0_over_N_mc)

                N_l0_over_N_eq0_mc = 0
                err_N_l0_over_N_eq0_mc = 0
                if N_l0_mc > 0 and N_eq0_mc > 0:
                    N_l0_over_N_eq0_mc = N_l0_mc / N_eq0_mc
                    err_N_l0_mc = h2_N_l0_mc.GetBinError(i, j)
                    err_N_eq0_mc = h2_N_eq0_mc.GetBinError(i, j)
                    err_N_l0_over_N_eq0_mc = frac_err(N_l0_mc, N_eq0_mc, err_N_l0_mc, err_N_eq0_mc)
                    h2_N_l0_over_N_eq0_mc.SetBinContent(i, j, N_l0_over_N_eq0_mc)
                    h2_N_l0_over_N_eq0_mc.SetBinError(i, j, err_N_l0_over_N_eq0_mc)

                if N_leq0_over_N_mc > 0:
                    # deconv_zf = np.abs( (N_leq0_over_N_data - N_leq0_over_N_mc) / N_leq0_over_N_mc )
                    deconv_zf = N_leq0_over_N_data - N_leq0_over_N_mc
                    h2_diff_zfData_zfMC.SetBinContent(i, j, 0)
                    h2_diff_zfData_zfMC.SetBinError(i, j, deconv_zf)

                    # Save <E/p> vs. eta histograms
                    valid_dir.cd()
                    h_zf_vs_eta.SetBinContent(j, deconv_zf)
                    h_zf_vs_eta.SetBinError(j, 0.)
                    h_zf_ratio_vs_eta.SetBinContent(j, N_leq0_over_N_data / N_leq0_over_N_mc)
                    h_zf_ratio_vs_eta.SetBinError(j, frac_err(N_leq0_over_N_data, N_leq0_over_N_mc, err_N_leq0_over_N_data, err_N_leq0_over_N_mc))

                    # save numbers for cross-checks
                    zf_numbers[s, i, j] = deconv_zf

                # Background corrected average E/p

                # Data

                avg_eop_raw_data = h_eop_raw_data.GetBinContent(j)
                err_eop_raw_data = np.abs(h_eop_raw_data.GetBinError(j))
                avg_eop_bg_data = h_eop_bg_data.GetBinContent(j)
                err_eop_bg_data = np.abs(h_eop_bg_data.GetBinError(j))
                avg_eop_cor_data = avg_eop_raw_data - 4./3.*avg_eop_bg_data
                err_eop_cor_data = np.sqrt(err_eop_raw_data**2 + (4./3.*err_eop_bg_data)**2)

                h2_eop_raw_data.SetBinContent(i, j, avg_eop_raw_data)
                h2_eop_raw_data.SetBinError(i, j, err_eop_raw_data)
                h2_eop_bg_data.SetBinContent(i, j, avg_eop_bg_data)
                h2_eop_bg_data.SetBinError(i, j, err_eop_bg_data)
                h2_eop_cor_data.SetBinContent(i, j, avg_eop_cor_data)
                h2_eop_cor_data.SetBinError(i, j, err_eop_cor_data)

                h_eop_raw_vs_eta_data.SetBinContent(j, avg_eop_raw_data)
                h_eop_raw_vs_eta_data.SetBinError(j, err_eop_raw_data)
                h_eop_bg_vs_eta_data.SetBinContent(j, avg_eop_bg_data)
                h_eop_bg_vs_eta_data.SetBinError(j, err_eop_bg_data)
                h_eop_cor_vs_eta_data.SetBinContent(j, avg_eop_cor_data)
                h_eop_cor_vs_eta_data.SetBinError(j, err_eop_cor_data)

                # MC

                avg_eop_raw_mc = h_eop_raw_mc.GetBinContent(j)
                err_eop_raw_mc = np.abs(h_eop_raw_mc.GetBinError(j))
                avg_eop_bg_mc = h_eop_bg_mc.GetBinContent(j)
                err_eop_bg_mc = np.abs(h_eop_bg_mc.GetBinError(j))
                avg_eop_cor_mc = avg_eop_raw_mc - 4./3.*avg_eop_bg_mc
                err_eop_cor_mc = np.sqrt(err_eop_raw_mc**2 + (4./3.*err_eop_bg_mc)**2)

                h2_eop_raw_mc.SetBinContent(i, j, avg_eop_raw_mc)
                h2_eop_raw_mc.SetBinError(i, j, err_eop_raw_mc)
                h2_eop_bg_mc.SetBinContent(i, j, avg_eop_bg_mc)
                h2_eop_bg_mc.SetBinError(i, j, err_eop_bg_mc)
                h2_eop_cor_mc.SetBinContent(i, j, avg_eop_cor_mc)
                h2_eop_cor_mc.SetBinError(i, j, err_eop_cor_mc)

                h_eop_raw_vs_eta_mc.SetBinContent(j, avg_eop_raw_mc)
                h_eop_raw_vs_eta_mc.SetBinError(j, err_eop_raw_mc)
                h_eop_bg_vs_eta_mc.SetBinContent(j, avg_eop_bg_mc)
                h_eop_bg_vs_eta_mc.SetBinError(j, err_eop_bg_mc)
                h_eop_cor_vs_eta_mc.SetBinContent(j, avg_eop_cor_mc)
                h_eop_cor_vs_eta_mc.SetBinError(j, err_eop_cor_mc)

                # Data over MC
                if avg_eop_raw_data > 0 and avg_eop_raw_mc > 0:
                    ratio_raw = avg_eop_raw_data/avg_eop_raw_mc
                    err_ratio_raw = frac_err(avg_eop_raw_data, avg_eop_raw_mc, err_eop_raw_data, err_eop_raw_mc)
                    h2_eop_raw_data_over_mc.SetBinContent(i, j, ratio_raw)
                    h2_eop_raw_data_over_mc.SetBinError(i, j, err_ratio_raw)

                if avg_eop_bg_data > 0 and avg_eop_bg_mc > 0:
                    ratio_bg = avg_eop_bg_data/avg_eop_bg_mc
                    err_ratio_bg = frac_err(avg_eop_bg_data, avg_eop_bg_mc, err_eop_bg_data, err_eop_bg_mc)
                    h2_eop_bg_data_over_mc.SetBinContent(i, j, ratio_bg)
                    h2_eop_bg_data_over_mc.SetBinError(i, j, err_ratio_bg)

                if avg_eop_cor_data > 0 and avg_eop_cor_mc > 0:
                    ratio_cor = avg_eop_cor_data/avg_eop_cor_mc
                    err_ratio_cor = frac_err(avg_eop_cor_data, avg_eop_cor_mc, err_eop_cor_data, err_eop_cor_mc)
                    h2_eop_cor_data_over_mc.SetBinContent(i, j, ratio_cor)
                    h2_eop_cor_data_over_mc.SetBinError(i, j, err_ratio_cor)
                    # save numbers for cross-checks
                    eop_cor_numbers[s][i][j] = ratio_cor
                    err_eop_cor_numbers[s][i][j] = err_ratio_cor

                # print (p_eta_str)
                # print ('<E/p>_RAW = {:f} +/- {:f}'.format(avg_eop_raw, err_eop_raw))
                # print ('4./3.*<E/p>_BG = {:f} +/- {:f}'.format(avg_eop_bg, err_eop_bg))
                # print ('<E/p>_COR = {:f} +/- {:f}'.format(avg_eop_cor, err_eop_cor))

            # Save <E/p> vs. eta histograms
            valid_dir.cd()
            h_N_leq0_over_N_vs_eta_data.Write()
            h_N_leq0_over_N_vs_eta_mc.Write()
            h_zf_vs_eta.Write()
            h_zf_ratio_vs_eta.Write()
            h_eop_raw_vs_eta_data.Write()
            h_eop_bg_vs_eta_data.Write()
            h_eop_cor_vs_eta_data.Write()
            h_eop_raw_vs_eta_mc.Write()
            h_eop_bg_vs_eta_mc.Write()
            h_eop_cor_vs_eta_mc.Write()

        # loop over again to create <E/p> vs. p in eta-bins
        for i in xrange(1, len(eta_bins)):

            eta_str = 'etaG{:d}L{:d}'.format((int)(eta_bins[i-1]*10), (int)(eta_bins[i]*10))

            # Raw data (TProfile)
            h_eop_raw_data = datafile.Get(selection+"/"+eta_str+"_eop_Total_"+energy_calib[s]+"_0_200")
            h_eop_bg_data = datafile.Get(selection+"/"+eta_str+"_eop_EM_BG_"+energy_calib[s]+"")

            # MC
            h_eop_raw_mc = mcfile.Get(selection+"/"+eta_str+"_eop_Total_"+energy_calib[s]+"_0_200")
            h_eop_bg_mc = mcfile.Get(selection+"/"+eta_str+"_eop_EM_BG_"+energy_calib[s]+"")

            if write_raw_hists:
                raw_dir.cd()
                append_name(h_eop_raw_data, "_data")
                append_name(h_eop_bg_data, "_data")
                append_name(h_eop_raw_mc, "_mc")
                append_name(h_eop_bg_mc, "_mc")
                h_eop_raw_data.Write()
                h_eop_bg_data.Write()
                h_eop_raw_mc.Write()
                h_eop_bg_mc.Write()

            # histograms of ZF vs. eta in p-bins
            h_N_leq0_over_N_vs_p_data = TH1D(selection_tags[s]+"_"+eta_str+"_N_leq0_over_N_vs_p_data", eta_str+"_N_leq0_over_N_vs_p_data", len(p_bins)-1, array('d', p_bins))
            h_N_leq0_over_N_vs_p_mc = TH1D(selection_tags[s]+"_"+eta_str+"_N_leq0_over_N_vs_p_mc", eta_str+"_N_leq0_over_N_vs_p_mc", len(p_bins)-1, array('d', p_bins))

            h_zf_vs_p = TH1D(selection_tags[s]+"_"+eta_str+"_zf_vs_p", eta_str+"_zf_vs_p", len(p_bins)-1, array('d', p_bins))
            h_zf_ratio_vs_p = TH1D(selection_tags[s]+"_"+eta_str+"_zf_ratio_vs_p", eta_str+"_zf_ratio_vs_p", len(p_bins)-1, array('d', p_bins))

            # histograms of <E/p> vs. p in eta-bins
            h_eop_raw_vs_p_data = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_raw_vs_p_data", eta_str+"_avgeop_raw_vs_p_data", len(p_bins)-1, array('d', p_bins))
            h_eop_bg_vs_p_data = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_bg_vs_p_data", eta_str+"_avgeop_bg_vs_p_data", len(p_bins)-1, array('d', p_bins))
            h_eop_cor_vs_p_data = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_cor_vs_p_data", eta_str+"_avgeop_cor_vs_p_data", len(p_bins)-1, array('d', p_bins))
            h_eop_raw_vs_p_mc = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_raw_vs_p_mc", eta_str+"_avgeop_raw_vs_p_mc", len(p_bins)-1, array('d', p_bins))
            h_eop_bg_vs_p_mc = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_bg_vs_p_mc", eta_str+"_avgeop_bg_vs_p_mc", len(p_bins)-1, array('d', p_bins))
            h_eop_cor_vs_p_mc = TH1D(selection_tags[s]+"_"+eta_str+"_avgeop_cor_vs_p_mc", eta_str+"_avgeop_cor_vs_p_mc", len(p_bins)-1, array('d', p_bins))

            for j in xrange(1, len(p_bins)):
                p_eta_str = 'pG{:d}L{:d}_etaG{:d}L{:d}'.format(
                    (int)(p_bins[j-1]*1000), (int)(p_bins[j]*1000),
                    (int)(eta_bins[i-1]*10), (int)(eta_bins[i]*10))

                # Zero fractions

                # Data
                N_data = h2_N_data.GetBinContent(j, i)
                N_leq0_data = h2_N_leq0_data.GetBinContent(j, i)

                N_leq0_over_N_data = 0
                err_N_leq0_over_N_data = 0
                if N_data > 0 and N_leq0_data > 0:
                    N_leq0_over_N_data = N_leq0_data / N_data
                    err_N_data = h2_N_data.GetBinError(j, i)
                    err_N_leq0_data = h2_N_leq0_data.GetBinError(j, i)
                    err_N_leq0_over_N_data = frac_err(N_leq0_data, N_data, err_N_leq0_data, err_N_data)

                    valid_dir.cd()
                    h_N_leq0_over_N_vs_p_data.SetBinContent(j, N_leq0_over_N_data)
                    h_N_leq0_over_N_vs_p_data.SetBinError(j, err_N_leq0_over_N_data)

                # MC

                N_mc = h2_N_mc.GetBinContent(j, i)
                N_leq0_mc = h2_N_leq0_mc.GetBinContent(j, i)

                N_leq0_over_N_mc = 0
                err_N_leq0_over_N_mc = 0
                if N_mc > 0 and N_leq0_mc > 0:
                    N_leq0_over_N_mc = N_leq0_mc / N_mc
                    err_N_mc = h2_N_mc.GetBinError(j, i)
                    err_N_leq0_mc = h2_N_leq0_mc.GetBinError(j, i)
                    err_N_leq0_over_N_mc = frac_err(N_leq0_mc, N_mc, err_N_leq0_mc, err_N_mc)

                    valid_dir.cd()
                    h_N_leq0_over_N_vs_p_mc.SetBinContent(j, N_leq0_over_N_mc)
                    h_N_leq0_over_N_vs_p_mc.SetBinError(j, err_N_leq0_over_N_mc)

                if N_leq0_over_N_mc > 0:
                    deconv_zf = N_leq0_over_N_data - N_leq0_over_N_mc

                    # Save <E/p> vs. eta histograms
                    valid_dir.cd()
                    h_zf_vs_p.SetBinContent(j, deconv_zf)
                    h_zf_vs_p.SetBinError(j, 0.)
                    h_zf_ratio_vs_p.SetBinContent(j, N_leq0_over_N_data / N_leq0_over_N_mc)
                    h_zf_ratio_vs_p.SetBinError(j, frac_err(N_leq0_over_N_data, N_leq0_over_N_mc, err_N_leq0_over_N_data, err_N_leq0_over_N_mc))

                # Background corrected average E/p

                # Data

                avg_eop_raw_data = h_eop_raw_data.GetBinContent(j)
                err_eop_raw_data = np.abs(h_eop_raw_data.GetBinError(j))
                avg_eop_bg_data = h_eop_bg_data.GetBinContent(j)
                err_eop_bg_data = np.abs(h_eop_bg_data.GetBinError(j))
                avg_eop_cor_data = avg_eop_raw_data - (4./3.)*avg_eop_bg_data
                err_eop_cor_data = np.sqrt(err_eop_raw_data**2 + (4./3.*err_eop_bg_data)**2)

                h_eop_raw_vs_p_data.SetBinContent(j, avg_eop_raw_data)
                h_eop_raw_vs_p_data.SetBinError(j, err_eop_raw_data)
                h_eop_bg_vs_p_data.SetBinContent(j, avg_eop_bg_data)
                h_eop_bg_vs_p_data.SetBinError(j, err_eop_bg_data)
                h_eop_cor_vs_p_data.SetBinContent(j, avg_eop_cor_data)
                h_eop_cor_vs_p_data.SetBinError(j, err_eop_cor_data)

                # MC

                avg_eop_raw_mc = h_eop_raw_mc.GetBinContent(j)
                err_eop_raw_mc = np.abs(h_eop_raw_mc.GetBinError(j))
                avg_eop_bg_mc = h_eop_bg_mc.GetBinContent(j)
                err_eop_bg_mc = np.abs(h_eop_bg_mc.GetBinError(j))
                avg_eop_cor_mc = avg_eop_raw_mc - 4./3.*avg_eop_bg_mc
                err_eop_cor_mc = np.sqrt(err_eop_raw_mc**2 + (4./3.*err_eop_bg_mc)**2)

                h_eop_raw_vs_p_mc.SetBinContent(j, avg_eop_raw_mc)
                h_eop_raw_vs_p_mc.SetBinError(j, err_eop_raw_mc)
                h_eop_bg_vs_p_mc.SetBinContent(j, avg_eop_bg_mc)
                h_eop_bg_vs_p_mc.SetBinError(j, err_eop_bg_mc)
                h_eop_cor_vs_p_mc.SetBinContent(j, avg_eop_cor_mc)
                h_eop_cor_vs_p_mc.SetBinError(j, err_eop_cor_mc)

            # Save <E/p> vs. p histograms
            valid_dir.cd()
            h_N_leq0_over_N_vs_p_data.Write()
            h_N_leq0_over_N_vs_p_mc.Write()
            h_zf_vs_p.Write()
            h_zf_ratio_vs_p.Write()
            h_eop_raw_vs_p_data.Write()
            h_eop_bg_vs_p_data.Write()
            h_eop_cor_vs_p_data.Write()
            h_eop_raw_vs_p_mc.Write()
            h_eop_bg_vs_p_mc.Write()
            h_eop_cor_vs_p_mc.Write()

        # Save 2D histograms for deconvolution input
        deconv_dir.cd()
        h2_N_leq0_over_N_data.Write()
        h2_N_leq0_over_N_mc.Write()
        # h2_N_leq0_over_N_data_over_mc.Write()
        h2_diff_zfData_zfMC.Write()

        h2_N_l0_over_N_eq0_data.Write()
        h2_N_l0_over_N_eq0_mc.Write()
        # h2_N_l0_over_N_eq0_data_over_mc.Write()

        h2_eop_raw_data.Write()
        h2_eop_raw_mc.Write()
        h2_eop_raw_data_over_mc.Write()

        h2_eop_bg_data.Write()
        h2_eop_bg_mc.Write()
        h2_eop_bg_data_over_mc.Write()

        h2_eop_cor_data.Write()
        h2_eop_cor_mc.Write()
        h2_eop_cor_data_over_mc.Write()

    with open(output_dir+"/"+tag+".log", "w") as f:
        f.write("# p, eta, eopcor, eopcor_err, zf, N_data, N_leq0_data, N_mc, N_leq0_mc\n")
        for s in xrange(len(selection_tags)):
            f.write("\n# "+selection_tags[s]+"\n")
            for i in xrange(1, len(p_bins)):
                for j in xrange(1, len(eta_bins)):
                    f.write("\n{}, {}, {}, {}, {}, {}, {}, {}, {}".format(p_bins[i], eta_bins[j], eop_cor_numbers[s][i][j],
                                                                  err_eop_cor_numbers[s][i][j], zf_numbers[s][i][j],
                                                                  n_trks_data_numbers[s][i][j], n_trks_leq0_data_numbers[s][i][j],
                                                                  n_trks_mc_numbers[s][i][j], n_trks_leq0_mc_numbers[s][i][j],
                                                                  ))

    outputfile.Close()

def frac_err(a, b, err_a, err_b):
    return np.abs(a/b)*np.sqrt( (err_a/a)**2 + (err_b/b)**2 )

def append_name(TObj, a = ''):
    TObj.SetName(TObj.GetName()+a)

if __name__ == "__main__":
    main(sys.argv)
