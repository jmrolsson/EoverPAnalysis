// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#include "EoverPAnalysis/EoverPHists.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <regex>
#include <sstream>
#include <cmath>

EoverPHists :: EoverPHists (std::string name, std::string detailStr, std::string energyCalib, bool doCaloTotal, bool doCaloEM, bool doCaloHAD, bool doBgSubtr, bool doTileLayer, std::string Ebins, bool doEbinsArray, std::string EbinsArray, std::string Etabins, bool doEtabinsArray, std::string EtabinsArray, bool doExtraEtaEnergyBinHists) : HistogramManager(name, detailStr)
{
  m_energyCalib = energyCalib;
  m_doCaloTotal = doCaloTotal;
  m_doCaloEM = doCaloEM;
  m_doCaloHAD = doCaloHAD;
  m_doBgSubtr = doBgSubtr;
  m_doTileLayer = doTileLayer;
  m_Ebins = Ebins;
  m_doEbinsArray = doEbinsArray;
  m_EbinsArray = EbinsArray; 
  m_Etabins = Etabins; 
  m_doEtabinsArray = doEtabinsArray;
  m_EtabinsArray = EtabinsArray; 
  m_doExtraEtaEnergyBinHists = doExtraEtaEnergyBinHists;
}

EoverPHists :: ~EoverPHists () {}

StatusCode EoverPHists::initialize()
{
  // number of bins and ranges for histograms
  unsigned int nBinsMu = 50;          float minMu = 0;                  float maxMu = 50;
  unsigned int nBinsNPV = 50;         float minNPV = -0.5;              float maxNPV = 49.5;
  unsigned int nBinsDR = 60;          float minDR = 0;                  float maxDR = 3;
  unsigned int nBinsPhi = 32;         float minPhi = -TMath::Pi();      float maxPhi = TMath::Pi(); 
  unsigned int nBinsPhiExtra = 1024;  float minPhiExtra = -TMath::Pi(); float maxPhiExtra = TMath::Pi(); 
  unsigned int nBinsPhiExtra2 = 800;  float minPhiExtra2 = -4.0;        float maxPhiExtra2 = 4.0; 
  unsigned int nBinsEop = 300;        float minEop = -4;                float maxEop = 20;
  unsigned int nBinsEop_l = 255;      float minEop_l = -100;            float maxEop_l = 5000;

  unsigned int nBinsE = 500;          float minE = 0;                   float maxE = 50;
  std::vector<double> Ebins = str2vec(m_Ebins);
  if (Ebins.size() > 2) {
    nBinsE = (int) Ebins[0];
    minE = (float) Ebins[1];
    maxE = (float) Ebins[2];
  }

  nEbinsArray = 0;
  if (m_doEbinsArray){
    EbinsArray = str2vec(m_EbinsArray);
    nEbinsArray = EbinsArray.size()-1;
  }
  if (nEbinsArray <= 1)
    m_doEbinsArray = false;

  unsigned int nBinsEta = 100;      float minEta = -2.5;           float maxEta = 2.5;
  std::vector<double> Etabins = str2vec(m_Etabins);
  if (Etabins.size() > 2) {
    nBinsEta = (int) Etabins[0];
    minEta = (float) Etabins[1];
    maxEta = (float) Etabins[2];
  }
  if (minEta < 0) m_doEtaAbs = false;
  else m_doEtaAbs = true;

  nEtabinsArray = 0;
  if (m_doEtabinsArray){
    EtabinsArray = str2vec(m_EtabinsArray);
    nEtabinsArray = EtabinsArray.size()-1;
  }
  if (nEtabinsArray <= 1)
    m_doEtabinsArray = false;

  // track kinematics
  m_trk_p = book(m_name, "trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE); 
  if (m_doEbinsArray) m_trk_p_array = book(m_name, "trk_p_EbinsArray", "p_{trk} [GeV]", nEbinsArray, &EbinsArray[0]); 
  m_trk_pt = book(m_name, "trk_pt", "p_{T,trk} [GeV]", nBinsE, minE, maxE); 
  if (m_doEbinsArray) m_trk_pt_array = book(m_name, "trk_pt_EbinsArray", "p_{T,trk} [GeV]", nEbinsArray, &EbinsArray[0]); 
  m_trk_etaID = book(m_name, "trk_etaID", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEMB2 = book(m_name, "trk_etaEMB2", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEME2 = book(m_name, "trk_etaEME2", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  if (m_doEtabinsArray) m_trk_eta_array = book(m_name, "trk_eta_EtabinsArray", "#eta_{trk}", nEtabinsArray, &EtabinsArray[0]); 
  m_trk_phiID = book(m_name, "trk_phiID", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEMB2 = book(m_name, "trk_phiEMB2", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEME2 = book(m_name, "trk_phiEME2", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phi_extra = book(m_name, "trk_phi_extra", "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra); 
  m_trk_phi_extra2 = book(m_name, "trk_phi_extra2", "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2); 
  m_trk_DR_EMB2_ID = book(m_name, std::string("trk_DR_EMB2_ID"), std::string("#Delta R_{trk}(EMB2, ID)"), nBinsDR, minDR, maxDR); 
  m_trk_DEta_EMB2_ID = book(m_name, std::string("trk_DEta_EMB2_ID"), std::string("#Delta #eta_{trk}(EMB2, ID)"), nBinsEta, minEta, maxEta); 
  m_trk_DPhi_EMB2_ID = book(m_name, std::string("trk_DPhi_EMB2_ID"), std::string("#Delta #phi_{trk}(EMB2, ID)"), nBinsPhi, minPhi, maxPhi); 
  m_trk_DR_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DR_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta R_{trk}(EMB2, ID)"), nBinsDR, minDR, maxDR); 
  m_trk_DEta_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DEta_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta #eta_{trk}(EMB2, ID)"), nBinsEta, minEta, maxEta); 
  m_trk_DPhi_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DPhi_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, std::string("#Delta #phi_{trk}(EMB2, ID)"), nBinsPhi, minPhi, maxPhi); 

  // Tile energy fractions
  m_trk_TileEfrac_100 = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_100"), "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_100_vs_trk_p = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_100_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_200 = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_200"), "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_200_vs_trk_p = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_200_vs_trk_p"), "p_{trk} [GeV]", nBinsE, minE, maxE, "E(tile)/E(total)", 100, 0, 5.0); 

  // basic "sanity" checks of the E/p xAOD derivation
  m_trk_sumE_Tile_100 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Tile_100"), "E(Tile)", nBinsE, minE, maxE); 
  m_trk_sumE_Tile_200 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Tile_200"), "E(Tile)", nBinsE, minE, maxE); 
  m_trk_sumE_Lar_100 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Lar_100"), "E(Lar)", nBinsE, minE, maxE); 
  m_trk_sumE_Lar_200 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Lar_200"), "E(Lar)", nBinsE, minE, maxE); 
  m_trk_SumTileLayers_over_HAD_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumTileLayers_over_HAD_100"), "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumLarLayers_over_EM_100"), "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumHADLayers_over_HAD_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumHADLayers_over_HAD_100"), "#Sigma Tile+EM layers / EM", 100, 0, 5.0); 
  m_trk_SumEMLayers_over_EM_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumEMLayers_over_EM_100"), "#Sigma EM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_100 = book (m_name, std::string("trk_"+m_energyCalib+"_EMandHAD_over_Total_100"), "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumAllLayers_over_Total_100"), "#Sigma all layers / Total", 100, 0, 5.0); 
  m_trk_SumTileLayers_over_HAD_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumTileLayers_over_HAD_200"), "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumLarLayers_over_EM_200"), "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumHADLayers_over_HAD_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumHADLayers_over_HAD_200"), "#Sigma Tile+EM layers / EM", 100, 0, 5.0); 
  m_trk_SumEMLayers_over_EM_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumEMLayers_over_EM_200"), "#Sigma EM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_200 = book (m_name, std::string("trk_"+m_energyCalib+"_EMandHAD_over_Total_200"), "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumAllLayers_over_Total_200"), "#Sigma all layers / Total", 100, 0, 5.0); 

  // calo layer with the highest cluster energy
  m_trk_E_100_highElayer = book (m_name, std::string("trk_"+m_energyCalib+"_E_100_highElayer"), "Calorimeter layer with highest E", 21, 0, 21); 
  m_trk_E_200_highElayer = book (m_name, std::string("trk_"+m_energyCalib+"_E_200_highElayer"), "Calorimeter layer with highest E", 21, 0, 21); 
  m_trk_E_100_highElayer_vs_E = book (m_name, std::string("trk_"+m_energyCalib+"_E_100_highElayer_vs_E"), "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 
  m_trk_E_200_highElayer_vs_E = book (m_name, std::string("trk_"+m_energyCalib+"_E_200_highElayer_vs_E"), "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 

  // sum of the energy of clusters matched to tracks vs. calo layer
  m_trk_E_100_vs_layer = book(m_name, std::string("trk_"+m_energyCalib+"_E_100_vs_layer"), "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 
  m_trk_E_200_vs_layer = book(m_name, std::string("trk_"+m_energyCalib+"_E_200_vs_layer"), "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 

  // Zero fraction histograms
  if (m_doEbinsArray && m_doEtabinsArray) {
    m_trk_n_E_200 = book(m_name, "trk_n_E_200", "p_{trk} [GeV]", nEbinsArray, &EbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_l0 = book(m_name, "trk_n_E_200_l0", "p_{trk} [GeV]", nEbinsArray, &EbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_eq0 = book(m_name, "trk_n_E_200_eq0", "p_{trk} [GeV]", nEbinsArray, &EbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_leq0 = book(m_name, "trk_n_E_200_leq0", "p_{trk} [GeV]", nEbinsArray, &EbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
  }

  // create E/p histograms

  // total calorimeter energy
  if (m_doCaloTotal) {
    // dR(trk,cluster) < 0.1
    m_trk_E_Total_100 = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_trk_E_Total_100_vs_mu_avg = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E", nBinsE, minE, maxE);

    m_eop_Total_100 = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_Total_100_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_Total_100_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_Total_100_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_Total_100_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_trkPhi = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_trkPhiID = book(m_name, std::string("eop_ID_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_trkPhi_extra = book(m_name, std::string("eop_extra_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_trkPhi_extra2 = book(m_name, std::string("eop_extra2_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_mu = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_mu_avg = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_npv = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_100_vs_highE_layer = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_Total_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_Total_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // dR(trk,cluster) < 0.2
    m_trk_E_Total_200 = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_trk_E_Total_200_vs_mu_avg = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E", nBinsE, minE, maxE);
    m_eop_Total_200 = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_l = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_l"), "E/p", nBinsEop_l, minEop_l, maxEop_l);
    if (m_doEbinsArray) m_eop_Total_200_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_Total_200_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_Total_200_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_Total_200_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_trkPhi = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_trkPhiID = book(m_name, std::string("eop_ID_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_trkPhi_extra = book(m_name, std::string("eop_extra_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_trkPhi_extra2 = book(m_name, std::string("eop_extra2_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_mu = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_mu_avg = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_npv = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_vs_highE_layer = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_Total_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_Total_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
  } // doCaloTotal

  // EM calorimeter (EMB+EMEC)
  if (m_doCaloEM) {
    // dR(trk,cluster) < 0.1
    m_trk_E_EM_100 = book(m_name, std::string("trk_E_EM_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_EM_100 = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_EM_100_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_100_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_EM_100_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_100_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_100_vs_trkPhi = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_100_vs_mu = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_100_vs_mu_avg = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_100_vs_npv = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_100_vs_highE_layer = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_EM_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_EM_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // dR(trk,cluster) < 0.2
    m_trk_E_EM_200 = book(m_name, std::string("trk_E_EM_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_EM_200 = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_EM_200_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_200_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_EM_200_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_200_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_vs_trkPhi = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_vs_mu = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_vs_mu_avg = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_vs_npv = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_vs_highE_layer = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_EM_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_EM_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }

    // No clusters in HAD matched to trk 
    m_eop_Total_200_noHAD = book(m_name, std::string("eop_Total_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_noHAD = book(m_name, std::string("eop_EM_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_noHAD = book(m_name, std::string("eop_HAD_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);

  } // doCaloEM

  // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
  if (m_doCaloHAD) {
    // dR(trk,cluster) < 0.1
    m_trk_E_HAD_100 = book(m_name, std::string("trk_E_HAD_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_HAD_100 = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_HAD_100_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_HAD_100_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_HAD_100_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_HAD_100_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_100_vs_trkPhi = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_100_vs_mu = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_100_vs_mu_avg = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_100_vs_npv = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_100_vs_highE_layer = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_HAD_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_HAD_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // dR(trk,cluster) < 0.2
    m_trk_E_HAD_200 = book(m_name, std::string("trk_E_HAD_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_HAD_200 = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_HAD_200_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_HAD_200_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_HAD_200_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_HAD_200_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_vs_trkPhi = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_vs_mu = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_vs_mu_avg = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_vs_npv = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_vs_highE_layer = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_HAD_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_HAD_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }

    // MIP requirement for HAD calorimeter
    m_eop_Total_200_MIP = book(m_name, std::string("eop_Total_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_MIP = book(m_name, std::string("eop_EM_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_MIP = book(m_name, std::string("eop_HAD_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);

  } // doCaloHAD

  // background subtraction, the way it was done in run 1
  if (m_doBgSubtr) {
    m_eop_EM_BG = book(m_name, std::string("eop_EM_BG_"+m_energyCalib), "E/p Background", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_EM_BG_vs_trkP = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_BG_vs_trkP = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_EM_BG_vs_trkEta = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_EM_BG_vs_trkEta = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_BG_vs_trkPhi = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_BG_vs_mu = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_BG_vs_mu_avg = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_BG_vs_npv = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_EM_BG_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_BG_"+m_energyCalib), "E/p", nBinsEop, minEop, maxEop);
          m_eop_EM_BG_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
  } // END doBgSubr

  // E/p where E is maximum in either Tile layer A, BC, or D
  if (m_doTileLayer) {
    // dR(trk,cluster) < 0.1
    // max in Tile Layer A
    m_trk_E_highTileA_100 = book(m_name, std::string("trk_E_highTileA_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileA_100 = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileA_100_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileA_100_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileA_100_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileA_100_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_100_vs_trkPhi = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_100_vs_mu = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_100_vs_mu_avg = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_100_vs_npv = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_100_vs_TileA_E = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_TileA_E"), "E(Tile, layer A)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileA_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileA_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // max in Tile Layer BC
    m_trk_E_highTileB_100 = book(m_name, std::string("trk_E_highTileB_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileB_100 = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileB_100_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileB_100_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileB_100_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileB_100_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_100_vs_trkPhi = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_100_vs_mu = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_100_vs_mu_avg = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_100_vs_npv = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_100_vs_TileB_E = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_TileB_E"), "E(Tile, layer B)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileB_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileB_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // max in Tile Layer D
    m_trk_E_highTileD_100 = book(m_name, std::string("trk_E_highTileD_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileD_100 = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileD_100_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileD_100_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileD_100_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileD_100_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_100_vs_trkPhi = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_100_vs_mu = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_100_vs_mu_avg = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_100_vs_npv = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_100_vs_TileD_E = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_TileD_E"), "E(Tile, layer D)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileD_100_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray); 
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileD_100_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // dR(trk,cluster) < 0.2
    // max in Tile Layer A
    m_trk_E_highTileA_200 = book(m_name, std::string("trk_E_highTileA_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileA_200 = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileA_200_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileA_200_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileA_200_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileA_200_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_200_vs_trkPhi = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_200_vs_mu = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_200_vs_mu_avg = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_200_vs_npv = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileA_200_vs_TileA_E = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_TileA_E"), "E(Tile, layer A)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileA_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray);
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileA_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // max in Tile Layer B
    m_trk_E_highTileB_200 = book(m_name, std::string("trk_E_highTileB_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileB_200 = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileB_200_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileB_200_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileB_200_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileB_200_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_200_vs_trkPhi = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_200_vs_mu = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_200_vs_mu_avg = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_200_vs_npv = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileB_200_vs_TileB_E = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_TileB_E"), "E(Tile, layer B)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileB_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray);
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileB_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
    // max in Tile Layer D
    m_trk_E_highTileD_200 = book(m_name, std::string("trk_E_highTileD_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileD_200 = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doEbinsArray) m_eop_highTileD_200_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nEbinsArray, &EbinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileD_200_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if(m_doEtabinsArray) m_eop_highTileD_200_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "E/p", nBinsEop, minEop, maxEop);
    else m_eop_highTileD_200_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_200_vs_trkPhi = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_200_vs_mu = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_200_vs_mu_avg = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_200_vs_npv = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "E/p", nBinsEop, minEop, maxEop);
    m_eop_highTileD_200_vs_TileD_E = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_TileD_E"), "E(Tile, layer D)", nBinsE, minE, maxE, "E/p", nBinsEop, minEop, maxEop);
    if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
      m_eop_highTileD_200_EtaEnergyRanges = std::vector<std::vector<TH1F*> >(nEbinsArray);
      char buffer [200];
      for (unsigned int i=0; i<nEbinsArray; i++) {
        for (unsigned int j=0; j<nEtabinsArray; j++) {
          std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((EbinsArray[i]*1000)), (int)std::round((EbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
          TH1F* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
          m_eop_highTileD_200_EtaEnergyRanges[i].push_back(tmp_hist);
        }
      }
    }
  } // END doTileLayer

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHists::execute( const xAOD::TrackParticle* trk, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, float eventWeight )
{
  // get pileup
  float mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  float mu_avg(-1.);
  if( eventInfo->isAvailable< float >( "corrected_averageInteractionsPerCrossing" ) &&
      !eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )) 
    mu_avg = eventInfo->auxdata< float >( "corrected_averageInteractionsPerCrossing" );
  else if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) )
    mu_avg = eventInfo->averageInteractionsPerCrossing();

  // get number of primary vtxs 
  float npv = HelperFunctions::countPrimaryVertices(vtxs, 2);

  float trk_p = 0;
  if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
  float trk_pt = trk->pt()/1e3;
  // coordinates of the track in the ID
  float trk_etaID = trk->eta();
  float trk_phiID = trk->phi();
  // coordinates of the track extrapolated to the specified calorimeter layer
  float trk_etaEMB2 = trk->auxdata<float>("CALO_trkEta_EMB2");
  float trk_phiEMB2 = trk->auxdata<float>("CALO_trkPhi_EMB2");
  float trk_etaEME2 = trk->auxdata<float>("CALO_trkEta_EME2");
  float trk_phiEME2 = trk->auxdata<float>("CALO_trkPhi_EME2");
  float dEta_EMB2_ID = TMath::Abs(trk_etaEMB2 - trk_etaID);
  float dPhi_EMB2_ID = TMath::Abs(trk_phiEMB2 - trk_phiID);
  if (dPhi_EMB2_ID > TMath::Pi())
    dPhi_EMB2_ID = 2*TMath::Pi() - dPhi_EMB2_ID;
  float dR_EMB2_ID = sqrt( pow(dEta_EMB2_ID, 2) + pow(dPhi_EMB2_ID, 2) );

  if (m_doEtaAbs) {
    trk_etaID = TMath::Abs(trk_etaID);
    trk_etaEMB2 = TMath::Abs(trk_etaEMB2);
  }

  // fill histos with basic track properties
  m_trk_p -> Fill(trk_p, eventWeight); 
  if (m_doEbinsArray) m_trk_p_array -> Fill(trk_p, eventWeight); 
  m_trk_pt -> Fill(trk_pt, eventWeight); 
  if (m_doEbinsArray) m_trk_pt_array -> Fill(trk_pt, eventWeight); 
  m_trk_etaID -> Fill(trk_etaID, eventWeight); 
  m_trk_etaEMB2 -> Fill(trk_etaEMB2, eventWeight); 
  m_trk_etaEME2 -> Fill(trk_etaEME2, eventWeight); 
  if (m_doEtabinsArray) m_trk_eta_array -> Fill(trk_etaID, eventWeight); 
  m_trk_phiID -> Fill(trk_phiID, eventWeight); 
  m_trk_phiEMB2 -> Fill(trk_phiEMB2, eventWeight); 
  m_trk_phiEME2 -> Fill(trk_phiEME2, eventWeight); 
  m_trk_phi_extra -> Fill(trk_phiID, eventWeight); 
  m_trk_phi_extra2 -> Fill(trk_phiID, eventWeight); 
  m_trk_DR_EMB2_ID -> Fill(dR_EMB2_ID, eventWeight); 
  m_trk_DEta_EMB2_ID -> Fill(dEta_EMB2_ID, eventWeight); 
  m_trk_DPhi_EMB2_ID -> Fill(dPhi_EMB2_ID, eventWeight); 
  m_trk_DR_EMB2_ID_vs_trk_p -> Fill(trk_p, dR_EMB2_ID, eventWeight);
  m_trk_DEta_EMB2_ID_vs_trk_p -> Fill(trk_p, dEta_EMB2_ID, eventWeight);
  m_trk_DPhi_EMB2_ID_vs_trk_p -> Fill(trk_p, dPhi_EMB2_ID, eventWeight);

  // cluster energy associated with the track
  float trk_E_Total_100 = trk->auxdata<float>(std::string("CALO_Total_"+m_energyCalib+"_0_100"))/1e3; 
  float trk_E_Total_200 = trk->auxdata<float>(std::string("CALO_Total_"+m_energyCalib+"_0_200"))/1e3; 
  float trk_E_EM_100 = trk->auxdata<float>(std::string("CALO_EM_"+m_energyCalib+"_0_100"))/1e3; 
  float trk_E_EM_200 = trk->auxdata<float>(std::string("CALO_EM_"+m_energyCalib+"_0_200"))/1e3; 
  float trk_E_HAD_100 = trk->auxdata<float>(std::string("CALO_HAD_"+m_energyCalib+"_0_100"))/1e3; 
  float trk_E_HAD_200 = trk->auxdata<float>(std::string("CALO_HAD_"+m_energyCalib+"_0_200"))/1e3; 

  // calculate energy fractions in different layers
  float trk_sumE_Tile_100 = 0; 
  float trk_sumE_Tile_200 = 0; 
  float trk_sumE_Lar_100 = 0; 
  float trk_sumE_Lar_200 = 0; 
  float trk_sumE_HAD_100 = 0; 
  float trk_sumE_HAD_200 = 0; 
  float trk_sumE_EM_100 = 0; 
  float trk_sumE_EM_200 = 0; 
  for (unsigned int i=0; i<m_layer_tile.size(); i++) {
    float trk_E_100_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_100"))/1e3; 
    if (trk_E_100_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
      trk_sumE_Tile_100 += trk_E_100_tmp; 
    float trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_200"))/1e3; 
    if (trk_E_200_tmp > 0.)
      trk_sumE_Tile_200 += trk_E_200_tmp; 
  }
  for (unsigned int i=0; i<m_layer_lar.size(); i++) {
    trk_sumE_Lar_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_100"))/1e3; 
    trk_sumE_Lar_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_200"))/1e3; 
  }
  for (unsigned int i=0; i<m_layer_had.size(); i++) {
    trk_sumE_HAD_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_had[i]+"_100"))/1e3; 
    trk_sumE_HAD_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_had[i]+"_200"))/1e3; 
  }
  for (unsigned int i=0; i<m_layer_em.size(); i++) {
    trk_sumE_EM_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_em[i]+"_100"))/1e3; 
    trk_sumE_EM_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_em[i]+"_200"))/1e3; 
  }

  // basic "sanity" checks of the E/p xAOD derivation
  m_trk_sumE_Tile_100 -> Fill(trk_sumE_Tile_100, eventWeight);
  m_trk_sumE_Tile_200 -> Fill(trk_sumE_Tile_200, eventWeight);
  m_trk_sumE_Lar_100 -> Fill(trk_sumE_Lar_100, eventWeight);
  m_trk_sumE_Lar_200 -> Fill(trk_sumE_Lar_200, eventWeight);
  m_trk_SumTileLayers_over_HAD_100 -> Fill(trk_sumE_Tile_100/trk_E_HAD_100, eventWeight);
  m_trk_SumTileLayers_over_HAD_200 -> Fill(trk_sumE_Tile_200/trk_E_HAD_200, eventWeight);
  m_trk_SumLarLayers_over_EM_100 -> Fill(trk_sumE_Lar_100/trk_E_EM_100, eventWeight);
  m_trk_SumLarLayers_over_EM_200 -> Fill(trk_sumE_Lar_200/trk_E_EM_200, eventWeight);
  m_trk_SumHADLayers_over_HAD_100 -> Fill(trk_sumE_HAD_100/trk_E_HAD_100, eventWeight);
  m_trk_SumHADLayers_over_HAD_200 -> Fill(trk_sumE_HAD_200/trk_E_HAD_200, eventWeight);
  m_trk_SumEMLayers_over_EM_100 -> Fill(trk_sumE_EM_100/trk_E_EM_100, eventWeight);
  m_trk_SumEMLayers_over_EM_200 -> Fill(trk_sumE_EM_200/trk_E_EM_200, eventWeight);
  m_trk_EMandHAD_over_Total_100 -> Fill((trk_E_EM_100+trk_E_HAD_100)/trk_E_Total_100, eventWeight);
  m_trk_EMandHAD_over_Total_200 -> Fill((trk_E_EM_200+trk_E_HAD_200)/trk_E_Total_200, eventWeight);

  // determine which layer has the highest energy deposit 
  int trk_highE_100_layer = 0;
  int trk_highE_200_layer = 0;
  float trk_highE_100 = -1e8;
  float trk_highE_200 = -1e8;
  float trk_E_100_tmp = -1e8;
  float trk_E_200_tmp = -1e8;
  float trk_sumE_Total_100 = 0.;
  float trk_sumE_Total_200 = 0.;
  float trk_sumEg0_Total_100 = 0.;
  float trk_sumEg0_Total_200 = 0.;
  for (unsigned int i=0; i<m_layer.size(); i++) { 
    trk_E_100_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_100"))/1e3; 
    trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_200"))/1e3; 
    trk_sumE_Total_100 += trk_E_100_tmp; 
    trk_sumE_Total_200 += trk_E_200_tmp; 
    if (trk_E_100_tmp > 0.) 
      trk_sumEg0_Total_100 += trk_E_100_tmp; 
    if (trk_E_200_tmp > 0.) 
      trk_sumEg0_Total_200 += trk_E_200_tmp; 
    m_trk_E_100_vs_layer -> Fill(i, trk_E_100_tmp, eventWeight); 
    m_trk_E_200_vs_layer -> Fill(i, trk_E_200_tmp, eventWeight); 

    if (trk_E_100_tmp > trk_highE_100) {
      trk_highE_100 = trk_E_100_tmp;
      trk_highE_100_layer = i; 
    }
    if (trk_E_200_tmp > trk_highE_200) {
      trk_highE_200 = trk_E_200_tmp;
      trk_highE_200_layer = i; 
    }
  }

  float trk_TileEfrac_100 = trk_sumE_Tile_100/trk_sumEg0_Total_100;    
  float trk_TileEfrac_200 = trk_sumE_Tile_200/trk_sumEg0_Total_200;    

  int trk_p_i = -1;
  int trk_eta_i = -1;
  if (m_doEtabinsArray && m_doEbinsArray && m_doExtraEtaEnergyBinHists) {
    for (unsigned int i=0; i<nEbinsArray; i++) {
      if (trk_p > EbinsArray[i] && trk_p < EbinsArray[i+1]) 
        trk_p_i = i;
    }
    for (unsigned int i=0; i<nEtabinsArray; i++) {
      if (trk_etaID > EtabinsArray[i] && trk_etaID < EtabinsArray[i+1]) 
        trk_eta_i = i; 
    }
  }

  // fill histos with Tile energy fractions
  m_trk_TileEfrac_100 -> Fill(trk_TileEfrac_100, eventWeight);
  m_trk_TileEfrac_100_vs_trk_p -> Fill(trk_p, trk_TileEfrac_100, eventWeight);
  m_trk_TileEfrac_200 -> Fill(trk_TileEfrac_200, eventWeight);
  m_trk_TileEfrac_200_vs_trk_p -> Fill(trk_p, trk_TileEfrac_200, eventWeight);
  // "sanity checks" check so that the energy in all clusters sum to the Total
  m_trk_SumAllLayers_over_Total_100 -> Fill(trk_sumE_Total_100/trk_E_Total_100, eventWeight);
  m_trk_SumAllLayers_over_Total_200 -> Fill(trk_sumE_Total_200/trk_E_Total_200, eventWeight);

  // fill histos of highest energy layer
  m_trk_E_100_highElayer -> Fill(trk_highE_100_layer, eventWeight);
  m_trk_E_200_highElayer -> Fill(trk_highE_200_layer, eventWeight);
  m_trk_E_100_highElayer_vs_E -> Fill(trk_highE_100, trk_highE_100_layer, eventWeight);
  m_trk_E_200_highElayer_vs_E -> Fill(trk_highE_200, trk_highE_200_layer, eventWeight);

  // fill E/p histograms
  // Total calorimeter energy
  if (m_doCaloTotal) {
    // dR(trk,cluster) < 0.1
    m_trk_E_Total_100 -> Fill(trk_E_Total_100, eventWeight); 
    m_trk_E_Total_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100, eventWeight); 
    m_eop_Total_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhiID -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi_extra -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi_extra2 -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_Total_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_Total_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    // dR(trk,cluster) < 0.2
 
    // fill zero fraction histograms
    m_trk_n_E_200 -> Fill(trk_p, trk_etaID, eventWeight);
    if (trk_E_Total_200 < 0) m_trk_n_E_200_l0 -> Fill(trk_p, trk_etaID, eventWeight);
    if (trk_E_Total_200 == 0) m_trk_n_E_200_eq0 -> Fill(trk_p, trk_etaID, eventWeight);
    if (trk_E_Total_200 <= 0) m_trk_n_E_200_leq0 -> Fill(trk_p, trk_etaID, eventWeight);

    m_trk_E_Total_200 -> Fill(trk_E_Total_200, eventWeight); 
    m_trk_E_Total_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200, eventWeight); 
    m_eop_Total_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_l -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhiID -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi_extra -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi_extra2 -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_Total_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_Total_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
  } // END doCaloTotal

  // EM calorimeter (EMB+EMEC)
  if (m_doCaloEM) {
    // dR(trk,cluster) < 0.1
    m_trk_E_EM_100 -> Fill(trk_E_EM_100, eventWeight); 
    m_eop_EM_100 -> Fill(trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkP -> Fill(trk_p, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkEta -> Fill(trk_etaID, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkPhi -> Fill(trk_phiID, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_mu -> Fill(mu, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_mu_avg -> Fill(mu_avg, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_npv -> Fill(npv, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_EM_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_EM_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_EM_100/trk_p, eventWeight); 
    // dR(trk,cluster) < 0.2
    m_trk_E_EM_200 -> Fill(trk_E_EM_200, eventWeight); 
    m_eop_EM_200 -> Fill(trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkP -> Fill(trk_p, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkEta -> Fill(trk_etaID, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkPhi -> Fill(trk_phiID, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_mu -> Fill(mu, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_mu_avg -> Fill(mu_avg, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_npv -> Fill(npv, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_EM_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_EM_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_EM_200/trk_p, eventWeight); 

    // No clusters in HAD matched to trk 
    if (trk_E_HAD_200 <=0) {
      m_eop_Total_200_noHAD -> Fill(trk_E_Total_200/trk_p, eventWeight);
      m_eop_EM_200_noHAD -> Fill(trk_E_EM_200/trk_p, eventWeight);
      m_eop_HAD_200_noHAD -> Fill(trk_E_HAD_200/trk_p, eventWeight);
    }

  } // END doCaloEM

  // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
  if (m_doCaloHAD) {
    // dR(trk,cluster) < 0.1
    m_trk_E_HAD_100 -> Fill(trk_E_HAD_100, eventWeight); 
    m_eop_HAD_100 -> Fill(trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkP -> Fill(trk_p, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkEta -> Fill(trk_etaID, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkPhi -> Fill(trk_phiID, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_mu -> Fill(mu, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_mu_avg -> Fill(mu_avg, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_npv -> Fill(npv, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_HAD_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_HAD_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_HAD_100/trk_p, eventWeight); 
    // dR(trk,cluster) < 0.2
    m_trk_E_HAD_200 -> Fill(trk_E_HAD_200, eventWeight); 
    m_eop_HAD_200 -> Fill(trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkP -> Fill(trk_p, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkEta -> Fill(trk_etaID, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkPhi -> Fill(trk_phiID, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_mu -> Fill(mu, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_mu_avg -> Fill(mu_avg, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_npv -> Fill(npv, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_HAD_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_HAD_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_HAD_200/trk_p, eventWeight); 
    // MIP requirement for HAD calorimeter
    if (trk_E_EM_100 < 1.1 && trk_E_HAD_100/trk_p > 0.4 && trk_E_HAD_100/trk_p < 0.9) {
      m_eop_Total_200_MIP -> Fill(trk_E_Total_200/trk_p, eventWeight);
      m_eop_EM_200_MIP -> Fill(trk_E_EM_200/trk_p, eventWeight);
      m_eop_HAD_200_MIP -> Fill(trk_E_HAD_200/trk_p, eventWeight);
    }

  } // END doCaloHAD

  // background subtraction, the way it was done in run 1
  if (m_doBgSubtr) {
    if (trk_E_EM_100 < 1.1 && trk_E_HAD_100/trk_p > 0.4 && trk_E_HAD_100/trk_p < 0.9) {
      m_eop_EM_BG -> Fill( (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkP -> Fill(trk_p, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkEta -> Fill(trk_etaID, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkPhi -> Fill(trk_phiID, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_mu -> Fill(mu, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_mu_avg -> Fill(mu_avg, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_npv -> Fill(npv, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_EM_BG_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill((trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
    }
  } // END doBgSubtr


  // E/p where E is maximum in either Tile layer A, BC, or D
  if (m_doTileLayer) {
    float highTileLayer_100 = -1e8;
    float highTileLayer_200 = -1e8;
    int highTileLayer_i_100 = -1; 
    int highTileLayer_i_200 = -1; 
    if (trk_highE_100_layer == 12 || trk_highE_100_layer == 18) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 1;                         
    } if (trk_highE_100_layer == 13 || trk_highE_100_layer == 19) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 2;                         
    } if (trk_highE_100_layer == 14 || trk_highE_100_layer == 20) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 3;                         
    } if (trk_highE_200_layer == 12 || trk_highE_200_layer == 18) {
      highTileLayer_200 = trk_highE_200;             
      highTileLayer_i_200 = 1;                         
    } if (trk_highE_200_layer == 13 || trk_highE_200_layer == 19) {
      highTileLayer_200 = trk_highE_200;             
      highTileLayer_i_200 = 2;                         
    } if (trk_highE_200_layer == 14 || trk_highE_200_layer == 20) {
      highTileLayer_200 = trk_highE_200;
      highTileLayer_i_200 = 3;
    }
    // dR(trk,cluster) < 0.1
    // max in Tile Layer A
    if (highTileLayer_i_100 == 1) { 
      m_trk_E_highTileA_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileA_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_TileA_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileA_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // max in Tile Layer BC
    if (highTileLayer_i_100 == 2) { 
      m_trk_E_highTileB_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileB_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_TileB_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileB_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // max in Tile Layer D
    if (highTileLayer_i_100 == 3) { 
      m_trk_E_highTileD_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileD_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_TileD_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileD_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // dR(trk,cluster) < 0.2
    // max in Tile Layer A
    if (highTileLayer_i_200 == 1) { 
      m_trk_E_highTileA_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileA_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_TileA_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileA_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
    // max in Tile Layer BC
    if (highTileLayer_i_200 == 2) { 
      m_trk_E_highTileB_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileB_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_TileB_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileB_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
    // max in Tile Layer D
    if (highTileLayer_i_200 == 3) { 
      m_trk_E_highTileD_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileD_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_TileD_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileD_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
  } // END doTileLayer

  return StatusCode::SUCCESS;
}

// Double_t* EoverPHists::linspace(float a, float b, unsigned int n) {
//   Double_t *vec = new Double_t[n];
//   vec[0] = a;
//   if ((a == b) || (n == 1)) return vec;
//   Double_t diff = b-a;
//   Double_t step=diff/(n-1);
//   for(unsigned int i=1; i<n; i++)
//     vec[i] = vec[i-1] + step;
//   return vec;
// }
//
// Double_t* EoverPHists::logspace(float a, float b, unsigned int n) {
//   Double_t *vec = linspace(TMath::Log10(a), TMath::Log10(b), n);
//   for(unsigned int i=0; i<n; i++)
//     vec[i] = float(TMath::Nint(TMath::Power(10, vec[i])*10))/10; 
//   return vec;
// }

std::vector<double> EoverPHists::str2vec(std::string str)
{
  std::vector<double> vec;
  if (str.size() > 0) {
    str.erase (std::remove (str.begin(), str.end(), ' '), str.end()); // remove whitespace
    std::stringstream ss(str);
    std::string token; // split string at ','
    while ( std::getline(ss, token, ','))
      vec.push_back(std::stof(token));
  }
  return vec;
}
