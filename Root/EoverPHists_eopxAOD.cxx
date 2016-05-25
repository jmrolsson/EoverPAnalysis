#include "EoverP/EoverPHists_eopxAOD.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <regex>
#include <math.h>

EoverPHists_eopxAOD :: EoverPHists_eopxAOD (std::string name, std::string detailStr, std::string trkExtrapol, bool doBgSubtr, bool doEMcalib, bool doLCWcalib, bool doCells, bool doCaloEM, bool doCaloHAD, bool doCaloTotal, float trkIsoDRmax, float trkIsoPfrac, float LarEmax, float TileEfrac, std::string Ebins, bool doEbinsArray, std::string EbinsArray, std::string EtaAbsbins, bool doEtaAbsbinsArray, std::string EtaAbsbinsArray, bool doTrkPcut, float trkPmin, float trkPmax, bool doTrkEtacut, float trkEtamin, float trkEtamax) :
  HistogramManager(name, detailStr)
{
  m_trkExtrapol = trkExtrapol; 
  m_doBgSubtr = doBgSubtr;
  m_doEMcalib = doEMcalib;
  m_doLCWcalib = doLCWcalib;
  m_doCells = doCells;
  m_doCaloEM = doCaloEM;
  m_doCaloHAD = doCaloHAD;
  m_doCaloTotal = doCaloTotal;
  m_trkIsoDRmax = trkIsoDRmax;
  m_trkIsoPfrac = trkIsoPfrac;
  m_LarEmax = LarEmax;
  m_TileEfrac = TileEfrac;
  m_Ebins = Ebins;
  m_doEbinsArray = doEbinsArray;
  m_EbinsArray = EbinsArray; 
  m_EtaAbsbins = EtaAbsbins; 
  m_doEtaAbsbinsArray = doEtaAbsbinsArray;
  m_EtaAbsbinsArray = EtaAbsbinsArray; 
  m_doTrkPcut = doTrkPcut;
  m_trkPmin = trkPmin;
  m_trkPmax = trkPmax; 
  m_doTrkEtacut = doTrkEtacut;
  m_trkEtamin = trkEtamin;
  m_trkEtamax = trkEtamax; 
}

EoverPHists_eopxAOD :: ~EoverPHists_eopxAOD () {}

StatusCode EoverPHists_eopxAOD::initialize()
{

  // number of bins and ranges for histograms
  int nBinsMu = 50;        float minMu = -0.5;            float maxMu = 49.5;
  int nBinsNPV = 50;       float minNPV = -0.5;           float maxNPV = 49.5;
  int nBinsDR = 60;        float minDR = 0;               float maxDR = 3;
  int nBinsPhi = 128;      float minPhi = -3.2;           float maxPhi = 3.2; 

  int nBinsTrkN = 200;     float minTrkN = -0.5;          float maxTrkN = 199.5;

  int nBinsEop = 250;      float minEop = -5;             float maxEop = 20;

  int nBinsE = 300;        float minE = 0;                float maxE = 30;
  std::vector<double> Ebins = str2vec(m_Ebins);
  if (Ebins.size() > 2) {
    nBinsE = (int) Ebins[0];
    minE = (float) Ebins[1];
    maxE = (float) Ebins[2];
  }

  int nBinsEta = 100;      float minEta = -2.5;           float maxEta = 2.5;
  int nBinsEtaAbs = 50;    float minEtaAbs = 0;           float maxEtaAbs = 2.5;
  std::vector<double> EtaAbsbins = str2vec(m_EtaAbsbins);
  if (EtaAbsbins.size() > 2) {
    nBinsEtaAbs = (int) Ebins[0];
    minEtaAbs = (float) Ebins[1];
    maxEtaAbs = (float) Ebins[2];
    nBinsEta = (int) 2*Ebins[0];
    maxEta = (float) Ebins[2];
    minEta = -maxEta;
  }

  std::vector<double> EbinsArray;
  int nEbinsArray = 0;
  if (m_doEbinsArray){
    EbinsArray = str2vec(m_EbinsArray);
    nEbinsArray = EbinsArray.size()-1;
  }

  std::vector<double> EtaAbsbinsArray;
  int nEtaAbsbinsArray = 0;
  if (m_doEtaAbsbinsArray){
    EtaAbsbinsArray = str2vec(m_EtaAbsbinsArray);
    nEtaAbsbinsArray = EtaAbsbinsArray.size()-1;
  }

  //// Book histograms

  // event level plots
  m_mu = book(m_name, "mu", "#mu", nBinsMu, minMu, maxMu); 
  m_mu_avg = book(m_name, "mu_avg", "<#mu>", nBinsMu, minMu, maxMu); 
  m_npv = book(m_name, "npv", "NPV", nBinsNPV, minNPV, maxNPV); 

  // track plots
  m_trk_n_nocut= book(m_name, "trk_n_nocut", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n= book(m_name, "trk_n", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_p = book(m_name, "trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE); 
  m_trk_eta = book(m_name, "trk_eta", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_eta_abs = book(m_name, "trk_eta_abs", "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs); 
  m_trk_phi = book(m_name, "trk_phi", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_DR_CALO_ID = book(m_name, std::string("trk_DR_"+m_trkExtrapol+"_ID"), std::string("#Delta R_{trk}("+m_trkExtrapol+", ID)"), nBinsDR, minDR, maxDR); 
  m_trk_DEta_CALO_ID = book(m_name, std::string("trk_DEta_"+m_trkExtrapol+"_ID"), std::string("#Delta #eta_{trk}("+m_trkExtrapol+", ID)"), nBinsEtaAbs, minEtaAbs, maxEtaAbs); 
  m_trk_DPhi_CALO_ID = book(m_name, std::string("trk_DPhi_"+m_trkExtrapol+"_ID"), std::string("#Delta #phi_{trk}("+m_trkExtrapol+", ID)"), nBinsPhi, minPhi, maxPhi); 

  m_trk_ntrks_maxDR01 = book(m_name, "trk_ntrks_maxDR01", "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02 = book(m_name, "trk_ntrks_maxDR02", "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03 = book(m_name, "trk_ntrks_maxDR03", "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04 = book(m_name, "trk_ntrks_maxDR04", "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05 = book(m_name, "trk_ntrks_maxDR05", "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06 = book(m_name, "trk_ntrks_maxDR06", "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07 = book(m_name, "trk_ntrks_maxDR07", "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08 = book(m_name, "trk_ntrks_maxDR08", "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09 = book(m_name, "trk_ntrks_maxDR09", "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10 = book(m_name, "trk_ntrks_maxDR10", "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_maxDR01_vs_trk_p = book(m_name, "trk_ntrks_maxDR01_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02_vs_trk_p = book(m_name, "trk_ntrks_maxDR02_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03_vs_trk_p = book(m_name, "trk_ntrks_maxDR03_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04_vs_trk_p = book(m_name, "trk_ntrks_maxDR04_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05_vs_trk_p = book(m_name, "trk_ntrks_maxDR05_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06_vs_trk_p = book(m_name, "trk_ntrks_maxDR06_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07_vs_trk_p = book(m_name, "trk_ntrks_maxDR07_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08_vs_trk_p = book(m_name, "trk_ntrks_maxDR08_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09_vs_trk_p = book(m_name, "trk_ntrks_maxDR09_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10_vs_trk_p = book(m_name, "trk_ntrks_maxDR10_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p", "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p", "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsE, minE, maxE); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsE, minE, maxE); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsE, minE, maxE); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{avg} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{leading} / p_{trk}", 100, 0, 5.0); 

  // sum of the energy of clusters matched to tracks vs. calo layer
  m_ClusterEnergy_100_vs_layer = book(m_name, "ClusterEnergy_100_vs_layer", "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 
  m_ClusterEnergy_200_vs_layer = book(m_name, "ClusterEnergy_200_vs_layer", "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 
  m_ClusterEnergy_100_vs_layer_passTrkIso = book(m_name, "ClusterEnergy_100_vs_layer_passTrkIso", "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 
  m_ClusterEnergy_200_vs_layer_passTrkIso = book(m_name, "ClusterEnergy_200_vs_layer_passTrkIso", "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 

  // calo layer with the highest cluster energy
  m_ClusterEnergy_100_highestEnergyLayer = book (m_name, "ClusterEnergy_100_highestEnergyLayer", "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_200_highestEnergyLayer = book (m_name, "ClusterEnergy_200_highestEnergyLayer", "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_100_highestEnergyLayer_vs_E = book (m_name, "ClusterEnergy_100_highestEnergyLayer_vs_E", "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_200_highestEnergyLayer_vs_E = book (m_name, "ClusterEnergy_200_highestEnergyLayer_vs_E", "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_100_highestEnergyLayer_passTrkIso = book (m_name, "ClusterEnergy_100_highestEnergyLayer_passTrkIso", "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_200_highestEnergyLayer_passTrkIso = book (m_name, "ClusterEnergy_200_highestEnergyLayer_passTrkIso", "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_100_highestEnergyLayer_passTrkIso_vs_E = book (m_name, "ClusterEnergy_100_highestEnergyLayer_passTrkIso_vs_E", "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 
  m_ClusterEnergy_200_highestEnergyLayer_passTrkIso_vs_E = book (m_name, "ClusterEnergy_200_highestEnergyLayer_passTrkIso_vs_E", "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 

  m_trk_TileEfrac_100 = book (m_name, "trk_TileEfrac_100", "Tile Energy Fraction", 100, 0, 5.0); 
  m_trk_TileEfrac_100_vs_trk_p = book (m_name, "trk_TileEfrac_100_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "Tile Energy Fraction", 100, 0, 5.0); 
  m_trk_TileEfrac_200 = book (m_name, "trk_TileEfrac_200", "Tile Energy Fraction", 100, 0, 5.0); 
  m_trk_TileEfrac_200_vs_trk_p = book (m_name, "trk_TileEfrac_200_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "Tile Energy Fraction", 100, 0, 5.0); 

  m_trk_SumTileLayers_over_HAD_100 = book (m_name, "trk_SumTileLayers_over_HAD_100", "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_100 = book (m_name, "trk_SumLarLayers_over_EM_100", "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumTileLarHECLayers_over_HAD_100 = book (m_name, "trk_SumTileLarHECLayers_over_HAD_100", "#Sigma Tile+LarEM layers / EM", 100, 0, 5.0); 
  m_trk_SumLarEMLayers_over_EM_100 = book (m_name, "trk_SumLarEMLayers_over_EM_100", "#Sigma LarEM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_100 = book (m_name, "trk_EMandHAD_over_Total_100", "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_100 = book (m_name, "trk_SumAllLayers_over_Total_100", "#Sigma all layers / Total", 100, 0, 5.0); 
  m_trk_SumTileLayers_over_HAD_200 = book (m_name, "trk_SumTileLayers_over_HAD_200", "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_200 = book (m_name, "trk_SumLarLayers_over_EM_200", "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumTileLarHECLayers_over_HAD_200 = book (m_name, "trk_SumTileLarHECLayers_over_HAD_200", "#Sigma Tile+LarEM layers / EM", 100, 0, 5.0); 
  m_trk_SumLarEMLayers_over_EM_200 = book (m_name, "trk_SumLarEMLayers_over_EM_200", "#Sigma LarEM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_200 = book (m_name, "trk_EMandHAD_over_Total_200", "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_200 = book (m_name, "trk_SumAllLayers_over_Total_200", "#Sigma all layers / Total", 100, 0, 5.0); 

  m_trk_DR_CALO_ID_vs_trk_p = book(m_name, std::string("trk_DR_"+m_trkExtrapol+"_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta R_{trk}("+m_trkExtrapol+", ID)"), nBinsDR, minDR, maxDR); 
  m_trk_DEta_CALO_ID_vs_trk_p = book(m_name, std::string("trk_DEta_"+m_trkExtrapol+"_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta #eta_{trk}("+m_trkExtrapol+", ID)"), nBinsEtaAbs, minEtaAbs, maxEtaAbs); 
  m_trk_DPhi_CALO_ID_vs_trk_p = book(m_name, std::string("trk_DPhi_"+m_trkExtrapol+"_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta #phi_{trk}("+m_trkExtrapol+", ID)"), nBinsPhi, minPhi, maxPhi); 

  // clusters calibrated at the EM-scale
  if (m_doEMcalib) { 
    // all calo clusters associated with a selected track (dR < trkIsoDRmax)
    if (m_doCaloTotal) {
      // dR_matched < 0.1
      m_trk_matched_Total_ClusterEnergy_100 = book(m_name, "trk_Total_ClusterEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_Total_ClusterEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_Total_ClusterEnergy_200 = book(m_name, "trk_Total_ClusterEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_Total_ClusterEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloEM) {
      // dR_matched < 0.1
      m_trk_matched_EM_ClusterEnergy_100 = book(m_name, "trk_EM_ClusterEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_EM_ClusterEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_EM_ClusterEnergy_200 = book(m_name, "trk_EM_ClusterEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_EM_ClusterEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloHAD) {
      // dR_matched < 0.1
      m_trk_matched_HAD_ClusterEnergy_100 = book(m_name, "trk_HAD_ClusterEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_HAD_ClusterEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_HAD_ClusterEnergy_200 = book(m_name, "trk_HAD_ClusterEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_HAD_ClusterEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doBgSubtr) {
      m_eop_EM_BG_ClusterEnergy = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy"), "E/p Background", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_BG_ClusterEnergy_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_ClusterEnergy_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_BG_ClusterEnergy_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_ClusterEnergy_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergy_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergy_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergy_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergy_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergy_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
  } // END doEMCalib

  // clusters calibrated at the LCW-scale
  if (m_doLCWcalib) { 
    // all calo clusters associated with a selected track (dR < trkIsoDRmax)
    if (m_doCaloTotal) {
      // dR_matched < 0.1
      m_trk_matched_Total_ClusterEnergyLCW_100 = book(m_name, "trk_Total_ClusterEnergyLCW_0_100", "E", nBinsE, minE, maxE);
      m_eop_Total_ClusterEnergyLCW_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_Total_ClusterEnergyLCW_200 = book(m_name, "trk_Total_ClusterEnergyLCW_0_200", "E", nBinsE, minE, maxE);
      m_eop_Total_ClusterEnergyLCW_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_ClusterEnergyLCW_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_ClusterEnergyLCW_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloEM) {
      // dR_matched < 0.1
      m_trk_matched_EM_ClusterEnergyLCW_100 = book(m_name, "trk_EM_ClusterEnergyLCW_0_100", "E", nBinsE, minE, maxE);
      m_eop_EM_ClusterEnergyLCW_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_EM_ClusterEnergyLCW_200 = book(m_name, "trk_EM_ClusterEnergyLCW_0_200", "E", nBinsE, minE, maxE);
      m_eop_EM_ClusterEnergyLCW_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_ClusterEnergyLCW_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_ClusterEnergyLCW_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloHAD) {
      // dR_matched < 0.1
      m_trk_matched_HAD_ClusterEnergyLCW_100 = book(m_name, "trk_HAD_ClusterEnergyLCW_0_100", "E", nBinsE, minE, maxE);
      m_eop_HAD_ClusterEnergyLCW_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergyLCW_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergyLCW_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_HAD_ClusterEnergyLCW_200 = book(m_name, "trk_HAD_ClusterEnergyLCW_0_200", "E", nBinsE, minE, maxE);
      m_eop_HAD_ClusterEnergyLCW_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergyLCW_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_ClusterEnergyLCW_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_ClusterEnergyLCW_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_ClusterEnergyLCW_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doBgSubtr) {
      m_eop_EM_BG_ClusterEnergyLCW = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW"), "E/p Background", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_BG_ClusterEnergyLCW_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_ClusterEnergyLCW_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_BG_ClusterEnergyLCW_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_ClusterEnergyLCW_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergyLCW_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergyLCW_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergyLCW_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_ClusterEnergyLCW_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_ClusterEnergyLCW_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
  }

  // cells 
  if (m_doCells) { 
    // all calo clusters associated with a selected track (dR < trkIsoDRmax)
    if (m_doCaloTotal) {
      // dR_matched < 0.1
      m_trk_matched_Total_CellEnergy_100 = book(m_name, "trk_Total_CellEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_Total_CellEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_Total_CellEnergy_200 = book(m_name, "trk_Total_CellEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_Total_CellEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_Total_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_Total_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_Total_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_Total_CellEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_Total_CellEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloEM) {
      // dR_matched < 0.1
      m_trk_matched_EM_CellEnergy_100 = book(m_name, "trk_EM_CellEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_EM_CellEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_EM_CellEnergy_200 = book(m_name, "trk_EM_CellEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_EM_CellEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_CellEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_CellEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doCaloHAD) {
      // dR_matched < 0.1
      m_trk_matched_HAD_CellEnergy_100 = book(m_name, "trk_HAD_CellEnergy_0_100", "E", nBinsE, minE, maxE);
      m_eop_HAD_CellEnergy_100 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_CellEnergy_100_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_CellEnergy_100_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_100_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_100_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_100_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_100_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
      // dR_matched < 0.2
      m_trk_matched_HAD_CellEnergy_200 = book(m_name, "trk_HAD_CellEnergy_0_200", "E", nBinsE, minE, maxE);
      m_eop_HAD_CellEnergy_200 = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200"), "E/p", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_HAD_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_CellEnergy_200_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_HAD_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_HAD_CellEnergy_200_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_200_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_200_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_200_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_HAD_CellEnergy_200_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_HAD_CellEnergy_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
    if (m_doBgSubtr) {
      m_eop_EM_BG_CellEnergy = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy"), "E/p Background", nBinsEop, minEop, maxEop);
      if (m_doEbinsArray) m_eop_EM_BG_CellEnergy_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_CellEnergy_vs_trkP = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
      if(m_doEtaAbsbinsArray) m_eop_EM_BG_CellEnergy_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_trkEta"), "|#eta_{trk}|", nEtaAbsbinsArray, &EtaAbsbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
      else m_eop_EM_BG_CellEnergy_vs_trkEta = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_trkEta"), "|#eta_{trk}|", nBinsEtaAbs, minEtaAbs, maxEtaAbs, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_CellEnergy_vs_trkPhi = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_CellEnergy_vs_mu = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_CellEnergy_vs_mu_avg = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
      m_eop_EM_BG_CellEnergy_vs_npv = book(m_name, std::string("eop_trkEtaPhi_"+m_trkExtrapol+"_EM_BG_CellEnergy_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    }
  }

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHists_eopxAOD::execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, float eventWeight )
{

  // get pileup
  float mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  float mu_avg(-1.);
  if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) ) {
    mu_avg = eventInfo->averageInteractionsPerCrossing();
  }
  m_mu -> Fill(mu, eventWeight);
  m_mu_avg -> Fill(mu_avg, eventWeight);

  // get number of primary vtxs 
  float npv = HelperFunctions::countPrimaryVertices(vtxs, 2);
  m_npv -> Fill(npv, eventWeight);

  // number of tracks
  m_trk_n_nocut-> Fill(trks->size());

  // track-cluster plots
  float trk_n = 0.;
  float trk_p = -1.;
  float trk_pt = 0;
  float trk_etaID = 0; 
  float trk_phiID = 0; 
  float trk_etaCALO = 0; 
  float trk_phiCALO = 0; 
  float dR_CALO_ID = 0; 
  float dEta_CALO_ID = 0; 
  float dPhi_CALO_ID = 0; 
  float trk2_p = -1.;
  float trk2_etaCALO = 0; 
  float trk2_phiCALO = 0; 
  float trk_trk2_dR = 0; 
  float trk_trk2_dEta = 0; 
  float trk_trk2_dPhi = 0; 
  float surr_trk_sum_p = 0.;
  float surr_trk_avg_p = 0.;
  float surr_trk_leading_p = 0.;

  float trk_matched_Total_ClusterEnergy_100 = 0; 
  float trk_matched_EM_ClusterEnergy_100 = 0; 
  float trk_matched_HAD_ClusterEnergy_100 = 0; 
  float trk_matched_Total_ClusterEnergy_200 = 0; 
  float trk_matched_EM_ClusterEnergy_200 = 0; 
  float trk_matched_HAD_ClusterEnergy_200 = 0; 
  float trk_matched_Total_ClusterEnergyLCW_100 = 0; 
  float trk_matched_EM_ClusterEnergyLCW_100 = 0; 
  float trk_matched_HAD_ClusterEnergyLCW_100 = 0; 
  float trk_matched_Total_ClusterEnergyLCW_200 = 0; 
  float trk_matched_EM_ClusterEnergyLCW_200 = 0; 
  float trk_matched_HAD_ClusterEnergyLCW_200 = 0; 
  float trk_matched_Total_CellEnergy_100 = 0; 
  float trk_matched_EM_CellEnergy_100 = 0; 
  float trk_matched_HAD_CellEnergy_100 = 0; 
  float trk_matched_Total_CellEnergy_200 = 0; 
  float trk_matched_EM_CellEnergy_200 = 0; 
  float trk_matched_HAD_CellEnergy_200 = 0; 
  float trk_matched_Tile_ClusterEnergy_100 = 0; 
  float trk_matched_Tile_ClusterEnergy_200 = 0; 
  float trk_matched_Lar_ClusterEnergy_100 = 0; 
  float trk_matched_Lar_ClusterEnergy_200 = 0; 
  float trk_matched_TileLarHEC_ClusterEnergy_100 = 0; 
  float trk_matched_TileLarHEC_ClusterEnergy_200 = 0; 
  float trk_matched_LarEM_ClusterEnergy_100 = 0; 
  float trk_matched_LarEM_ClusterEnergy_200 = 0; 

  // loop over tracks and clusters around a selected track, plot how many tracks (clusters) fall within a certain dR
  float maxDR_ranges[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();

  // track-track and track-cluster (outer loop over all tracks)
  for( ; trk_itr != trk_end; ++trk_itr ) {

    // count the number of tracks passing all selections
    trk_n += 1;
    const xAOD::TrackParticle* trk = (*trk_itr);
    if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
    trk_pt  = trk->pt()/1e3;
    trk_etaID = trk->eta();
    trk_phiID = trk->phi();
    trk_etaCALO = trk->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
    trk_phiCALO = trk->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));

    if (m_doTrkPcut) {
      if (trk_p < m_trkPmin) continue;
      if (trk_p > m_trkPmax) continue;
    }

    if (m_doTrkEtacut) {
      if (TMath::Abs(trk_etaCALO) < m_trkEtamin) continue;
      if (TMath::Abs(trk_etaCALO) > m_trkEtamax) continue;
    }

    // apply TileCal cuts
    trk_matched_Total_ClusterEnergy_100 = trk->auxdata<float>("CALO_Total_ClusterEnergy_0_100")/1e3; 
    trk_matched_HAD_ClusterEnergy_100 = trk->auxdata<float>("CALO_HAD_ClusterEnergy_0_100")/1e3; 
    trk_matched_EM_ClusterEnergy_100 = trk->auxdata<float>("CALO_EM_ClusterEnergy_0_100")/1e3; 
    trk_matched_Total_ClusterEnergy_200 = trk->auxdata<float>("CALO_Total_ClusterEnergy_0_200")/1e3; 
    trk_matched_HAD_ClusterEnergy_200 = trk->auxdata<float>("CALO_HAD_ClusterEnergy_0_200")/1e3; 
    trk_matched_EM_ClusterEnergy_200 = trk->auxdata<float>("CALO_EM_ClusterEnergy_0_200")/1e3; 
    trk_matched_Tile_ClusterEnergy_100 = 0; 
    trk_matched_Tile_ClusterEnergy_200 = 0; 
    trk_matched_Lar_ClusterEnergy_100 = 0; 
    trk_matched_Lar_ClusterEnergy_200 = 0; 
    trk_matched_TileLarHEC_ClusterEnergy_100 = 0; 
    trk_matched_TileLarHEC_ClusterEnergy_200 = 0; 
    trk_matched_LarEM_ClusterEnergy_100 = 0; 
    trk_matched_LarEM_ClusterEnergy_200 = 0; 
    for (unsigned int i=0; i<m_layer_tile.size(); i++) {
      trk_matched_Tile_ClusterEnergy_100 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_tile[i]+"_100"))/1e3; 
      trk_matched_Tile_ClusterEnergy_200 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_tile[i]+"_200"))/1e3; 
    }
    for (unsigned int i=0; i<m_layer_lar.size(); i++) {
      trk_matched_Lar_ClusterEnergy_100 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_lar[i]+"_100"))/1e3; 
      trk_matched_Lar_ClusterEnergy_200 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_lar[i]+"_200"))/1e3; 
    }
    for (unsigned int i=0; i<m_layer_tile_larhec.size(); i++) {
      trk_matched_TileLarHEC_ClusterEnergy_100 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_tile_larhec[i]+"_100"))/1e3; 
      trk_matched_TileLarHEC_ClusterEnergy_200 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_tile_larhec[i]+"_200"))/1e3; 
    }
    for (unsigned int i=0; i<m_layer_larem.size(); i++) {
      trk_matched_LarEM_ClusterEnergy_100 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_larem[i]+"_100"))/1e3; 
      trk_matched_LarEM_ClusterEnergy_200 += trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer_larem[i]+"_200"))/1e3; 
    }
    float trk_TileEfrac_100 = trk_matched_Tile_ClusterEnergy_100/trk_matched_Total_ClusterEnergy_100;    
    float trk_TileEfrac_200 = trk_matched_Tile_ClusterEnergy_200/trk_matched_Total_ClusterEnergy_200;    
    if (trk_matched_EM_ClusterEnergy_200 > m_LarEmax) continue; // energy loss in LAr 
    if (trk_TileEfrac_200 < m_TileEfrac) continue;

    dEta_CALO_ID = TMath::Abs(trk_etaCALO - trk_etaID);
    if (TMath::Abs(trk_phiCALO - trk_phiID) < TMath::Pi())
      dPhi_CALO_ID = TMath::Abs(trk_phiCALO - trk_phiID);
    else
      dPhi_CALO_ID = 2*TMath::Pi() - TMath::Abs(trk_phiCALO - trk_phiID);

    dR_CALO_ID = sqrt( pow(dEta_CALO_ID, 2) + pow(dPhi_CALO_ID, 2) );

    m_trk_TileEfrac_100 -> Fill(trk_TileEfrac_100, eventWeight);
    m_trk_TileEfrac_100_vs_trk_p -> Fill(trk_p, trk_TileEfrac_100, eventWeight);
    m_trk_TileEfrac_200 -> Fill(trk_TileEfrac_200, eventWeight);
    m_trk_TileEfrac_200_vs_trk_p -> Fill(trk_p, trk_TileEfrac_200, eventWeight);
    m_trk_SumTileLayers_over_HAD_100 -> Fill(trk_matched_Tile_ClusterEnergy_100/trk_matched_HAD_ClusterEnergy_100, eventWeight);
    m_trk_SumTileLayers_over_HAD_200 -> Fill(trk_matched_Tile_ClusterEnergy_200/trk_matched_HAD_ClusterEnergy_200, eventWeight);
    m_trk_SumLarLayers_over_EM_100 -> Fill(trk_matched_Lar_ClusterEnergy_100/trk_matched_EM_ClusterEnergy_100, eventWeight);
    m_trk_SumLarLayers_over_EM_200 -> Fill(trk_matched_Lar_ClusterEnergy_200/trk_matched_EM_ClusterEnergy_200, eventWeight);
    m_trk_SumTileLarHECLayers_over_HAD_100 -> Fill(trk_matched_TileLarHEC_ClusterEnergy_100/trk_matched_HAD_ClusterEnergy_100, eventWeight);
    m_trk_SumTileLarHECLayers_over_HAD_200 -> Fill(trk_matched_TileLarHEC_ClusterEnergy_200/trk_matched_HAD_ClusterEnergy_200, eventWeight);
    m_trk_SumLarEMLayers_over_EM_100 -> Fill(trk_matched_LarEM_ClusterEnergy_100/trk_matched_EM_ClusterEnergy_100, eventWeight);
    m_trk_SumLarEMLayers_over_EM_200 -> Fill(trk_matched_LarEM_ClusterEnergy_200/trk_matched_EM_ClusterEnergy_200, eventWeight);
    m_trk_EMandHAD_over_Total_100 -> Fill((trk_matched_EM_ClusterEnergy_100+trk_matched_HAD_ClusterEnergy_100)/trk_matched_Total_ClusterEnergy_100, eventWeight);
    m_trk_EMandHAD_over_Total_200 -> Fill((trk_matched_EM_ClusterEnergy_200+trk_matched_HAD_ClusterEnergy_200)/trk_matched_Total_ClusterEnergy_200, eventWeight);

    // basic track properties
    m_trk_p -> Fill(trk_p, eventWeight); 
    m_trk_eta -> Fill(trk_etaCALO, eventWeight); 
    m_trk_eta_abs -> Fill(TMath::Abs(trk_etaCALO), eventWeight); 
    m_trk_phi -> Fill(trk_phiCALO, eventWeight); 
    m_trk_DR_CALO_ID -> Fill(dR_CALO_ID, eventWeight); 
    m_trk_DR_CALO_ID_vs_trk_p -> Fill(trk_p, dR_CALO_ID, eventWeight);
    m_trk_DEta_CALO_ID -> Fill(dEta_CALO_ID, eventWeight); 
    m_trk_DEta_CALO_ID_vs_trk_p -> Fill(trk_p, dEta_CALO_ID, eventWeight);
    m_trk_DPhi_CALO_ID -> Fill(dPhi_CALO_ID, eventWeight); 
    m_trk_DPhi_CALO_ID_vs_trk_p -> Fill(trk_p, dPhi_CALO_ID, eventWeight);

    // keep track of the number of tracks within a certain dR from the selected track
    int trk_ntrks_maxDR[10] = {0};
    trk2_p = 0.;
    surr_trk_sum_p = 0.;
    surr_trk_leading_p = 0.;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      const xAOD::TrackParticle* trk2 = (*trk2_itr);
      trk2_etaCALO = trk2->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
      trk2_phiCALO = trk2->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));
      trk_trk2_dEta= TMath::Abs(trk2_etaCALO - trk_etaCALO);
      if (TMath::Abs(trk2_phiCALO - trk_phiCALO) < TMath::Pi())
        trk_trk2_dPhi = TMath::Abs(trk2_phiCALO - trk_phiCALO);
      else
        trk_trk2_dPhi = 2*TMath::Pi() - TMath::Abs(trk2_phiCALO - trk_phiCALO);
      trk_trk2_dR = sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );

      if (trk_itr != trk2_itr) { // do not double count the selected track 

        // count the number of track wihin a given radius
        for (int i = 0; i < 10; ++i) {
          if (trk_trk2_dR < maxDR_ranges[i]) {
            trk_ntrks_maxDR[i]++;
          }
        }

        // track isolation
        if (trk_trk2_dR < m_trkIsoDRmax) { // check if trk2 falls within DRmax of trk 

          // calculate the leading and avg p of the surrounding tracks
          if (TMath::Abs(trk2->qOverP())>0.) trk2_p = (1./TMath::Abs(trk2->qOverP()))/1e3; 
          surr_trk_sum_p += trk2_p;
          if (trk2_p > surr_trk_leading_p) {
            surr_trk_leading_p = trk2_p;
          }
        }
      }
    } // END looping trk2

    // plots of the number of close-by tracks within a certain radius
    m_trk_ntrks_maxDR01 -> Fill(trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02 -> Fill(trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03 -> Fill(trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04 -> Fill(trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05 -> Fill(trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06 -> Fill(trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07 -> Fill(trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08 -> Fill(trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09 -> Fill(trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10 -> Fill(trk_ntrks_maxDR[9], eventWeight);

    // plots of the number of close-by tracks vs p of the selected track
    m_trk_ntrks_maxDR01_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[9], eventWeight);

    // plots of leading p and avg p of the tracks surrounding the selected track 
    surr_trk_avg_p = surr_trk_sum_p/(trks->size()-1); // calculate avg p of the surrounding tracks 
    m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p -> Fill(surr_trk_sum_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p -> Fill(surr_trk_avg_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p -> Fill(surr_trk_leading_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_sum_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_avg_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_leading_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p -> Fill(surr_trk_sum_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p -> Fill(surr_trk_avg_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p -> Fill(surr_trk_leading_p/trk_p, eventWeight);

    int trk_highClusE_100_layer = 0;
    int trk_highClusE_200_layer = 0;
    float trk_highClusE_100 = -1e8;
    float trk_highClusE_200 = -1e8;
    float trk_clusE_100_tmp = -1e8;
    float trk_clusE_200_tmp = -1e8;
    float trk_clusE_100_sum = 0.;
    float trk_clusE_200_sum = 0.;
    for (int i=0; i<m_layer.size(); i++) { 
      trk_clusE_100_tmp = trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer[i]+"_100"))/1e3; 
      trk_clusE_200_tmp = trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer[i]+"_200"))/1e3; 
      trk_clusE_100_sum += trk_clusE_100_tmp; 
      trk_clusE_200_sum += trk_clusE_200_tmp; 
      m_ClusterEnergy_100_vs_layer -> Fill(i, trk_clusE_100_tmp, eventWeight); 
      m_ClusterEnergy_200_vs_layer -> Fill(i, trk_clusE_200_tmp, eventWeight); 
      
      if (trk_clusE_100_tmp > trk_highClusE_100) {
        trk_highClusE_100 = trk_clusE_100_tmp;
        trk_highClusE_100_layer = i; 
      }
      if (trk_clusE_200_tmp > trk_highClusE_200) {
        trk_highClusE_200 = trk_clusE_200_tmp;
        trk_highClusE_200_layer = i; 
      }
    }
    m_ClusterEnergy_100_highestEnergyLayer -> Fill(trk_highClusE_100_layer, eventWeight);
    m_ClusterEnergy_200_highestEnergyLayer -> Fill(trk_highClusE_200_layer, eventWeight);
    m_ClusterEnergy_100_highestEnergyLayer_vs_E -> Fill(trk_highClusE_100, trk_highClusE_100_layer, eventWeight);
    m_ClusterEnergy_200_highestEnergyLayer_vs_E -> Fill(trk_highClusE_200, trk_highClusE_200_layer, eventWeight);

    m_trk_SumAllLayers_over_Total_100 -> Fill(trk_clusE_100_sum/trk_matched_Total_ClusterEnergy_100, eventWeight);
    m_trk_SumAllLayers_over_Total_200 -> Fill(trk_clusE_200_sum/trk_matched_Total_ClusterEnergy_200, eventWeight);

    // check track isolation requirement
    if (TMath::Abs(surr_trk_sum_p/trk_p) > m_trkIsoPfrac) return StatusCode::SUCCESS;

    trk_highClusE_100_layer = 0;
    trk_highClusE_200_layer = 0;
    trk_highClusE_100 = -1e6;
    trk_highClusE_200 = -1e6;
    for (int i=0; i<m_layer.size(); i++) { 
      trk_clusE_100_tmp = trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer[i]+"_100"))/1e3; 
      trk_clusE_200_tmp = trk->auxdata<float>(std::string("CALO_ClusterEnergy_"+m_layer[i]+"_200"))/1e3; 
      m_ClusterEnergy_100_vs_layer_passTrkIso -> Fill(i, trk_clusE_100_tmp, eventWeight); 
      m_ClusterEnergy_200_vs_layer_passTrkIso -> Fill(i, trk_clusE_200_tmp, eventWeight); 

      if (trk_clusE_100_tmp > trk_highClusE_100) {
        trk_highClusE_100 = trk_clusE_100_tmp;
        trk_highClusE_100_layer = i; 
      }
      if (trk_clusE_200_tmp > trk_highClusE_200) {
        trk_highClusE_200 = trk_clusE_200_tmp;
        trk_highClusE_200_layer = i; 
      }
    }
    m_ClusterEnergy_100_highestEnergyLayer_passTrkIso -> Fill(trk_highClusE_100_layer, eventWeight);
    m_ClusterEnergy_200_highestEnergyLayer_passTrkIso -> Fill(trk_highClusE_200_layer, eventWeight);
    m_ClusterEnergy_100_highestEnergyLayer_passTrkIso_vs_E -> Fill(trk_highClusE_100, trk_highClusE_100_layer, eventWeight);
    m_ClusterEnergy_200_highestEnergyLayer_passTrkIso_vs_E -> Fill(trk_highClusE_200, trk_highClusE_200_layer, eventWeight);

    // fill eoverp
    if (m_doEMcalib) { 
      if (m_doCaloTotal) {
        // dR_matched < 0.1
        trk_matched_Total_ClusterEnergy_100 = trk->auxdata<float>("CALO_Total_ClusterEnergy_0_100")/1e3; 
        m_trk_matched_Total_ClusterEnergy_100 -> Fill(trk_matched_Total_ClusterEnergy_100, eventWeight); 
        m_eop_Total_ClusterEnergy_100 -> Fill(trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_mu -> Fill(mu, trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_100_vs_npv -> Fill(npv, trk_matched_Total_ClusterEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_Total_ClusterEnergy_200 = trk->auxdata<float>("CALO_Total_ClusterEnergy_0_200")/1e3; 
        m_trk_matched_Total_ClusterEnergy_200 -> Fill(trk_matched_Total_ClusterEnergy_200, eventWeight); 
        m_eop_Total_ClusterEnergy_200 -> Fill(trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_mu -> Fill(mu, trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergy_200_vs_npv -> Fill(npv, trk_matched_Total_ClusterEnergy_200/trk_p, eventWeight); 
      }
      if (m_doCaloEM) {
        // dR_matched < 0.1
        trk_matched_EM_ClusterEnergy_100 = trk->auxdata<float>("CALO_EM_ClusterEnergy_0_100")/1e3; 
        m_trk_matched_EM_ClusterEnergy_100 -> Fill(trk_matched_EM_ClusterEnergy_100, eventWeight); 
        m_eop_EM_ClusterEnergy_100 -> Fill(trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_mu -> Fill(mu, trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_100_vs_npv -> Fill(npv, trk_matched_EM_ClusterEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_EM_ClusterEnergy_200 = trk->auxdata<float>("CALO_EM_ClusterEnergy_0_200")/1e3; 
        m_trk_matched_EM_ClusterEnergy_200 -> Fill(trk_matched_EM_ClusterEnergy_200, eventWeight); 
        m_eop_EM_ClusterEnergy_200 -> Fill(trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_mu -> Fill(mu, trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergy_200_vs_npv -> Fill(npv, trk_matched_EM_ClusterEnergy_200/trk_p, eventWeight); 
      }
      if (m_doCaloHAD) {
        // dR_matched < 0.1
        trk_matched_HAD_ClusterEnergy_100 = trk->auxdata<float>("CALO_HAD_ClusterEnergy_0_100")/1e3; 
        m_trk_matched_HAD_ClusterEnergy_100 -> Fill(trk_matched_HAD_ClusterEnergy_100, eventWeight); 
        m_eop_HAD_ClusterEnergy_100 -> Fill(trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_mu -> Fill(mu, trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_100_vs_npv -> Fill(npv, trk_matched_HAD_ClusterEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_HAD_ClusterEnergy_200 = trk->auxdata<float>("CALO_HAD_ClusterEnergy_0_200")/1e3; 
        m_trk_matched_HAD_ClusterEnergy_200 -> Fill(trk_matched_HAD_ClusterEnergy_200, eventWeight); 
        m_eop_HAD_ClusterEnergy_200 -> Fill(trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_mu -> Fill(mu, trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergy_200_vs_npv -> Fill(npv, trk_matched_HAD_ClusterEnergy_200/trk_p, eventWeight); 
      }
      if (m_doBgSubtr) {
        if (trk_matched_EM_ClusterEnergy_100 < 1.1 && trk_matched_HAD_ClusterEnergy_100/trk_p > 0.4) {
          m_eop_EM_BG_ClusterEnergy -> Fill( (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_trkP -> Fill(trk_p, (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_trkPhi -> Fill(trk_phiCALO, (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_mu -> Fill(mu, (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_mu_avg -> Fill(mu_avg, (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergy_vs_npv -> Fill(npv, (trk_matched_EM_ClusterEnergy_200 - trk_matched_EM_ClusterEnergy_100)/trk_p, eventWeight); 
        }
      }
    } // END doEMcalib

    if (m_doLCWcalib) { 
      if (m_doCaloTotal) {
        // dR_matched < 0.1
        trk_matched_Total_ClusterEnergyLCW_100 = trk->auxdata<float>("CALO_Total_ClusterEnergyLCW_0_100")/1e3; 
        m_trk_matched_Total_ClusterEnergyLCW_100 -> Fill(trk_matched_Total_ClusterEnergyLCW_100, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100 -> Fill(trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_trkP -> Fill(trk_p, trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_mu -> Fill(mu, trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_100_vs_npv -> Fill(npv, trk_matched_Total_ClusterEnergyLCW_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_Total_ClusterEnergyLCW_200 = trk->auxdata<float>("CALO_Total_ClusterEnergyLCW_0_200")/1e3; 
        m_trk_matched_Total_ClusterEnergyLCW_200 -> Fill(trk_matched_Total_ClusterEnergyLCW_200, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200 -> Fill(trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_trkP -> Fill(trk_p, trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_mu -> Fill(mu, trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_Total_ClusterEnergyLCW_200_vs_npv -> Fill(npv, trk_matched_Total_ClusterEnergyLCW_200/trk_p, eventWeight); 
      }
      if (m_doCaloEM) {
        // dR_matched < 0.1
        trk_matched_EM_ClusterEnergyLCW_100 = trk->auxdata<float>("CALO_EM_ClusterEnergyLCW_0_100")/1e3; 
        m_trk_matched_EM_ClusterEnergyLCW_100 -> Fill(trk_matched_EM_ClusterEnergyLCW_100, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100 -> Fill(trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_trkP -> Fill(trk_p, trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_mu -> Fill(mu, trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_100_vs_npv -> Fill(npv, trk_matched_EM_ClusterEnergyLCW_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_EM_ClusterEnergyLCW_200 = trk->auxdata<float>("CALO_EM_ClusterEnergyLCW_0_200")/1e3; 
        m_trk_matched_EM_ClusterEnergyLCW_200 -> Fill(trk_matched_EM_ClusterEnergyLCW_200, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200 -> Fill(trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_trkP -> Fill(trk_p, trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_mu -> Fill(mu, trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_EM_ClusterEnergyLCW_200_vs_npv -> Fill(npv, trk_matched_EM_ClusterEnergyLCW_200/trk_p, eventWeight); 
      }
      if (m_doCaloHAD) {
        // dR_matched < 0.1
        trk_matched_HAD_ClusterEnergyLCW_100 = trk->auxdata<float>("CALO_HAD_ClusterEnergyLCW_0_100")/1e3; 
        m_trk_matched_HAD_ClusterEnergyLCW_100 -> Fill(trk_matched_HAD_ClusterEnergyLCW_100, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100 -> Fill(trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_trkP -> Fill(trk_p, trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_mu -> Fill(mu, trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_100_vs_npv -> Fill(npv, trk_matched_HAD_ClusterEnergyLCW_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_HAD_ClusterEnergyLCW_200 = trk->auxdata<float>("CALO_HAD_ClusterEnergyLCW_0_200")/1e3; 
        m_trk_matched_HAD_ClusterEnergyLCW_200 -> Fill(trk_matched_HAD_ClusterEnergyLCW_200, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200 -> Fill(trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_trkP -> Fill(trk_p, trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_mu -> Fill(mu, trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
        m_eop_HAD_ClusterEnergyLCW_200_vs_npv -> Fill(npv, trk_matched_HAD_ClusterEnergyLCW_200/trk_p, eventWeight); 
      }
      if (m_doBgSubtr) {
        if (trk_matched_EM_ClusterEnergyLCW_100 < 1.1 && trk_matched_HAD_ClusterEnergyLCW_100/trk_p > 0.4) {
          m_eop_EM_BG_ClusterEnergyLCW -> Fill( (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW -> Fill( (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_trkP -> Fill(trk_p, (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_trkPhi -> Fill(trk_phiCALO, (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_mu -> Fill(mu, (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_mu_avg -> Fill(mu_avg, (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
          m_eop_EM_BG_ClusterEnergyLCW_vs_npv -> Fill(npv, (trk_matched_EM_ClusterEnergyLCW_200 - trk_matched_EM_ClusterEnergyLCW_100)/trk_p, eventWeight); 
        }
      }
    } // END doLCWcalib

    if (m_doCells) { 
      if (m_doCaloTotal) {
        // dR_matched < 0.1
        trk_matched_Total_CellEnergy_100 = trk->auxdata<float>("CALO_Total_CellEnergy_0_100")/1e3; 
        m_trk_matched_Total_CellEnergy_100 -> Fill(trk_matched_Total_CellEnergy_100, eventWeight); 
        m_eop_Total_CellEnergy_100 -> Fill(trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_mu -> Fill(mu, trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_100_vs_npv -> Fill(npv, trk_matched_Total_CellEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_Total_CellEnergy_200 = trk->auxdata<float>("CALO_Total_CellEnergy_0_200")/1e3; 
        m_trk_matched_Total_CellEnergy_200 -> Fill(trk_matched_Total_CellEnergy_200, eventWeight); 
        m_eop_Total_CellEnergy_200 -> Fill(trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_mu -> Fill(mu, trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
        m_eop_Total_CellEnergy_200_vs_npv -> Fill(npv, trk_matched_Total_CellEnergy_200/trk_p, eventWeight); 
      }
      if (m_doCaloEM) {
        // dR_matched < 0.1
        trk_matched_EM_CellEnergy_100 = trk->auxdata<float>("CALO_EM_CellEnergy_0_100")/1e3; 
        m_trk_matched_EM_CellEnergy_100 -> Fill(trk_matched_EM_CellEnergy_100, eventWeight); 
        m_eop_EM_CellEnergy_100 -> Fill(trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_mu -> Fill(mu, trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_100_vs_npv -> Fill(npv, trk_matched_EM_CellEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_EM_CellEnergy_200 = trk->auxdata<float>("CALO_EM_CellEnergy_0_200")/1e3; 
        m_trk_matched_EM_CellEnergy_200 -> Fill(trk_matched_EM_CellEnergy_200, eventWeight); 
        m_eop_EM_CellEnergy_200 -> Fill(trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_mu -> Fill(mu, trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
        m_eop_EM_CellEnergy_200_vs_npv -> Fill(npv, trk_matched_EM_CellEnergy_200/trk_p, eventWeight); 
      }
      if (m_doCaloHAD) {
        // dR_matched < 0.1
        trk_matched_HAD_CellEnergy_100 = trk->auxdata<float>("CALO_HAD_CellEnergy_0_100")/1e3; 
        m_trk_matched_HAD_CellEnergy_100 -> Fill(trk_matched_HAD_CellEnergy_100, eventWeight); 
        m_eop_HAD_CellEnergy_100 -> Fill(trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_trkP -> Fill(trk_p, trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_mu -> Fill(mu, trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_100_vs_npv -> Fill(npv, trk_matched_HAD_CellEnergy_100/trk_p, eventWeight); 
        // dR_matched < 0.2
        trk_matched_HAD_CellEnergy_200 = trk->auxdata<float>("CALO_HAD_CellEnergy_0_200")/1e3; 
        m_trk_matched_HAD_CellEnergy_200 -> Fill(trk_matched_HAD_CellEnergy_200, eventWeight); 
        m_eop_HAD_CellEnergy_200 -> Fill(trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_trkP -> Fill(trk_p, trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_trkPhi -> Fill(trk_phiCALO, trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_mu -> Fill(mu, trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_mu_avg -> Fill(mu_avg, trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
        m_eop_HAD_CellEnergy_200_vs_npv -> Fill(npv, trk_matched_HAD_CellEnergy_200/trk_p, eventWeight); 
      }
    } // END doCells
    if (m_doBgSubtr) {
      if (trk_matched_EM_CellEnergy_100 < 1.1 && trk_matched_HAD_CellEnergy_100/trk_p > 0.4) {
        m_eop_EM_BG_CellEnergy -> Fill( (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_trkP -> Fill(trk_p, (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_trkEta -> Fill(TMath::Abs(trk_etaCALO), (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_trkPhi -> Fill(trk_phiCALO, (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_mu -> Fill(mu, (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_mu_avg -> Fill(mu_avg, (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
        m_eop_EM_BG_CellEnergy_vs_npv -> Fill(npv, (trk_matched_EM_CellEnergy_200 - trk_matched_EM_CellEnergy_100)/trk_p, eventWeight); 
      }
    }

  } // END loop tracks 

  m_trk_n -> Fill(trk_n);

  return StatusCode::SUCCESS;
}

Double_t* EoverPHists_eopxAOD::linspace(float a, float b, unsigned int n) {
  Double_t *vec = new Double_t[n];
  vec[0] = a;
  if ((a == b) || (n == 1)) return vec;
  Double_t diff = b-a;
  Double_t step=diff/(n-1);
  for(unsigned int i=1; i<n; i++)
    vec[i] = vec[i-1] + step;
  return vec;
}

Double_t* EoverPHists_eopxAOD::logspace(float a, float b, unsigned int n) {
  Double_t *vec = linspace(TMath::Log10(a), TMath::Log10(b), n);
  for(unsigned int i=0; i<n; i++)
    vec[i] = float(TMath::Nint(TMath::Power(10, vec[i])*10))/10; 
  return vec;
}

std::vector<double> EoverPHists_eopxAOD::str2vec(std::string str)
{
  std::vector<double> vec;
  str.erase (std::remove (str.begin(), str.end(), ' '), str.end()); // remove whitespace
  std::stringstream ss(str);
  std::string token; // split string at ','
  while ( std::getline(ss, token, ','))
    vec.push_back(std::stof(token));
  return vec;
}

