#! /usr/bin/env python
from __future__ import print_function
import sys
import os
import re
import argparse
import numpy as np

from ROOT import TFile, TCanvas, TH2D

def main(argv):

    path = '/Users/joakim/eoverp/results/eoverp_cross_checks/'
    tag = 'September15'

    ##### ---- Combined ND,SD,DD ----

    # Millie
    m_d_lcw = TFile.Open(os.path.join(path, 'eop_files/millie/data_LCWscale_13092017.root'))
    m_mc_lcw = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCWscale_13092017.root'))
    # m_d_em = TFile.Open(os.path.join(path, 'eop_files/millie/data_EMscale_13092017.root'))
    # m_mc_em = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EMscale_13092017.root.root'))

    # Joakim
    j_deconv = TFile.Open(os.path.join(path, 'eop_files/deconv_inputs_20171005.root'))

    # Output files
    f_diff_d_eop_lcw = open(os.path.join(path,'diffs/diff_data_avgeop_lcw_'+tag+'.txt'), 'w')
    f_diff_mc_eop_lcw = open(os.path.join(path,'diffs/diff_mc_avgeop_lcw_'+tag+'.txt'), 'w')
    # f_diff_d_eop_em = open(os.path.join(path,'diffs/diff_data_avgeop_em_'+tag+'.txt'), 'w')
    # f_diff_mc_eop_em = open(os.path.join(path,'diffs/diff_mc_avgeop_em_'+tag+'.txt'), 'w')

    # Run diff
    print_comparison(m_d_lcw, j_deconv, f_diff_d_eop_lcw)
    print_comparison(m_mc_lcw, j_deconv, f_diff_mc_eop_lcw, isMC=True)
    # print_comparison(m_d_em, j_deconv, f_diff_d_eop_em, isEM=True)
    # print_comparison(m_mc_em, j_deconv, f_diff_mc_eop_em, isEM=True, isMC=True)

    ## ----------------------------------------------------------------------------------------------

    #### ---- ND ----
    m_mc_lcw_nd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_ND.root'))
    m_mc_lcw_nd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_ND.root'))
    m_mc_em_nd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_ND.root'))
    m_mc_em_nd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_ND.root'))
    j_deconv_nd = TFile.Open(os.path.join(path, 'eop_files/deconv_inputs_20170529_joakim_ND.root'))

    f_diff_mc_eop_lcw_nd = open(os.path.join(path,'diffs/diff_mc_avgeop_lcw_nd_'+tag+'.txt'), 'w')
    f_diff_mc_eop_em_nd = open(os.path.join(path,'diffs/diff_mc_avgeop_em_nd_'+tag+'.txt'), 'w')

    print_comparison(m_mc_lcw_nd, j_deconv_nd, f_diff_mc_eop_lcw_nd, isMC=True, oldEtaBins=True)
    print_comparison(m_mc_em_nd, j_deconv_nd, f_diff_mc_eop_em_nd, isEM=True, isMC=True, oldEtaBins=True)

    #### ---- SD ----
    m_mc_lcw_sd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_SD.root'))
    m_mc_lcw_sd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_SD.root'))
    m_mc_em_sd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_SD.root'))
    m_mc_em_sd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_SD.root'))
    j_deconv_sd = TFile.Open(os.path.join(path, 'eop_files/deconv_inputs_20170529_joakim_SD.root'))

    f_diff_mc_eop_lcw_sd = open(os.path.join(path,'diffs/diff_mc_avgeop_lcw_sd_'+tag+'.txt'), 'w')
    f_diff_mc_eop_em_sd = open(os.path.join(path,'diffs/diff_mc_avgeop_em_sd_'+tag+'.txt'), 'w')

    print_comparison(m_mc_lcw_sd, j_deconv_sd, f_diff_mc_eop_lcw_sd, isMC=True, oldEtaBins=True)
    print_comparison(m_mc_em_sd, j_deconv_sd, f_diff_mc_eop_em_sd, isEM=True, isMC=True, oldEtaBins=True)

    #### ---- DD ----
    m_mc_lcw_dd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_DD.root'))
    m_mc_lcw_dd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_LCW_milliem_01062017_DD.root'))
    m_mc_em_dd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_DD.root'))
    m_mc_em_dd = TFile.Open(os.path.join(path, 'eop_files/millie/mc_EM_milliem_01062017_DD.root'))
    j_deconv_dd = TFile.Open(os.path.join(path, 'eop_files/deconv_inputs_20170529_joakim_DD.root'))

    f_diff_mc_eop_lcw_dd = open(os.path.join(path,'diffs/diff_mc_avgeop_lcw_dd_'+tag+'.txt'), 'w')
    f_diff_mc_eop_em_dd = open(os.path.join(path,'diffs/diff_mc_avgeop_em_dd_'+tag+'.txt'), 'w')

    print_comparison(m_mc_lcw_dd, j_deconv_dd, f_diff_mc_eop_lcw_dd, isMC=True, oldEtaBins=True)
    print_comparison(m_mc_em_dd, j_deconv_dd, f_diff_mc_eop_em_dd, isEM=True, isMC=True, oldEtaBins=True)

def print_comparison(m_file, j_file, f, isEM=False, isMC=False, oldEtaBins=False):

    # eta,p binning
    m_eta = ['0p0_eta_0p6', '0p6_eta_1p1', '1p1_eta_1p4', '1p4_eta_1p5', '1p5_eta_1p8', '1p8_eta_1p9', '1p9_eta_2p3']
    if oldEtaBins:
        m_eta = ['|#eta|<0.6', '0.6<|#eta|<1.1', '1.1<|#eta|<1.4', '1.4<|#eta|<1.5', '1.5<|#eta|<1.8', '1.8<|#eta|<1.9', '1.9<|#eta|<2.3']
    j_eta = ['etaG0L6', 'etaG6L11', 'etaG11L14', 'etaG14L15', 'etaG15L18', 'etaG18L19', 'etaG19L23']
    p_bins = [.5, .8, 1.2, 1.8, 2.2, 2.8, 3.4, 4.2, 5., 6., 7., 9., 12., 15., 20., 30.]
    eta_bins = [0., .6, 1.1, 1.4, 1.5, 1.8, 1.9, 2.3]

    print ("{: <20}, {: <20}, {: <20}, {: <10}, {: <10}, {: <10}, {: <10}".format('selection', 'eta', 'p', 'millie', 'joakim', 'diff', 'frac [%]'), file=f)
    for i in xrange(0, len(m_eta)):

        h_m = [m_file.Get(os.path.join(m_eta[i], 'hp_EoP_mean_raw')),
               m_file.Get(os.path.join(m_eta[i], 'hp_EoP_mean_bkd')),
               m_file.Get(os.path.join(m_eta[i], 'h1_EoP_mean_bkdcorr')),
               m_file.Get(os.path.join(m_eta[i], 'h1_zeroFraction'))]

        if isEM:
            e_tag = 'ClusterEnergy'
            folder_tag = 'EM'
        else:
            e_tag = 'ClusterEnergyLCW'
            folder_tag = 'LCW'

        if isMC:
            d_mc_tag = 'mc'
        else:
            d_mc_tag = 'data'

        # print('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_avgeop_cor_vs_p_'+d_mc_tag)
        # exit()

        # h_j = [j_file.Get('RawHistograms_'+folder_tag+'/'+j_eta[i]+'_eop_Total_'+e_tag+'_0_200_'+d_mc_tag),
        #        j_file.Get('RawHistograms_'+folder_tag+'/'+j_eta[i]+'_eop_EM_BG_'+e_tag+'_'+d_mc_tag),
        #        j_file.Get('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_avgeop_cor_vs_p_'+d_mc_tag)]
        h_j = [#j_file.Get('RawHistograms_'+folder_tag+'/'+j_eta[i]+'_eop_Total_'+e_tag+'_0_200_'+d_mc_tag),
               # j_file.Get('RawHistograms_'+folder_tag+'/'+j_eta[i]+'_eop_EM_BG_'+e_tag+'_'+d_mc_tag),
               j_file.Get('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_avgeop_raw_vs_p_'+d_mc_tag),
               j_file.Get('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_avgeop_bg_vs_p_'+d_mc_tag),
               j_file.Get('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_avgeop_cor_vs_p_'+d_mc_tag),
               j_file.Get('ValidationHistograms_'+folder_tag+'/'+folder_tag+'_'+j_eta[i]+'_N_leq0_over_N_vs_p_'+d_mc_tag)]

        labels = ["<E/p>_RAW", "<E/p>_BG", "<E/p>_COR", "ZeroFrac"]

        eta_label = m_eta[i]

        for j in xrange(1, len(p_bins)):

            for hindex in xrange(0, len(h_j)):

                m_bin = h_m[hindex].GetBinContent(j)
                j_bin = h_j[hindex].GetBinContent(j)
                diff_bin = m_bin - j_bin
                if (j_bin > 0):
                    frac_bin = np.abs(m_bin / j_bin)
                else:
                    frac_bin = 0

                p_label = "{:.1f}<p/[GeV]<{:.1f}".format(p_bins[j-1], p_bins[j])

                print("{: <20}, {: <20}, {: <20}, {:10.8f}, {:10.8f}, {:10.4f}, {:10.0f}".format(labels[hindex], eta_label, p_label, m_bin, j_bin, diff_bin, frac_bin*100), file=f)
            print("", file=f)

if __name__ == "__main__":
    main(sys.argv)
