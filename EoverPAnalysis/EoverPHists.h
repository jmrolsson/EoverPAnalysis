// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#ifndef EoverPAnalysis_EoverPHists_H
#define EoverPAnalysis_EoverPHists_H

#include "EoverPAnalysis/HistogramManagerD.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

class EoverPHists : public HistogramManagerD
{
  public:
    EoverPHists(std::string name, std::string detailStr, std::string energyCalib = "ClusterEnergy", bool doCaloTotal = true, bool doCaloEM = false, bool doCaloHAD = false, bool doBgSubtr = true, bool doTileLayer = false, std::string Pbins = "", bool doPbinsArray = false, std::string PbinsArray = "", std::string Etabins = "", bool doEtabinsArray = false, std::string EtabinsArray = "", bool doExtraEtaEnergyBinHists = false);
    ~EoverPHists();

    StatusCode initialize();

    StatusCode execute( const xAOD::TrackParticle* trk, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, double eventWeight );
    using HistogramManagerD::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManagerD::execute; // overload

    // Double_t* linspace(double a, double b, unsigned int n);
    // Double_t* logspace(double a, double b, unsigned int n);
    std::vector<double> str2vec(std::string str);

  protected:

    // bools to control which histograms are filled
    bool m_doCaloTotal; //!
    bool m_doCaloEM; //!
    bool m_doCaloHAD; //!
    bool m_doBgSubtr; //!
    bool m_doTileLayer; //!

    std::string m_energyCalib; //! what type of energy calibration to use (EM, LCW, cells)
    std::string m_Pbins; //!
    bool m_doPbinsArray; //!
    std::string m_PbinsArray; //!
    std::string m_Etabins; //!
    bool m_doEtaAbs; //!
    bool m_doEtabinsArray; //!
    std::string m_EtabinsArray; //!
    std::vector<double> PbinsArray = {0}; //!
    unsigned int nPbinsArray; //!
    std::vector<double> EtabinsArray = {0}; //!
    unsigned int nEtabinsArray; //!

    // save separate 1D histograms for the specified track eta and energy ranges
    bool m_doExtraEtaEnergyBinHists; //!

  private:

    const std::vector<std::string> m_layer = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3", "TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of all the calo layers
    const std::vector<std::string> m_layer_lar = {"EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of lar layers only
    const std::vector<std::string> m_layer_tile = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of tile layers only
    const std::vector<std::string> m_layer_had = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of HAD (Tile+HEC) layers only
    const std::vector<std::string> m_layer_em = {"EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3"}; //! array of EM layers only

    // track kinematics
    TH1D* m_trk_p; //!
    TH1D* m_trk_p_array; //!
    TH1D* m_trk_pt; //!
    TH1D* m_trk_pt_array; //!
    TH1D* m_trk_etaID; //!
    TH1D* m_trk_etaEMB2; //!
    TH1D* m_trk_etaEME2; //!
    TH1D* m_trk_eta_array; //!
    TH1D* m_trk_phiID; //!
    TH1D* m_trk_phiEMB2; //!
    TH1D* m_trk_phiEME2; //!
    TH1D* m_trk_phi_extra; //!
    TH1D* m_trk_phi_extra2; //!
    TH2D* m_trk_p_vs_eta; //!
    TH2D* m_trk_p_vs_etaEMB2; //!
    TH2D* m_trk_p_vs_etaEME2; //!
    TH2D* m_trk_p_vs_eta_array; //!

    TH1D* m_trk_d0;
    TH1D* m_trk_d0_s;
    TH1D* m_trk_z0;
    TH1D* m_trk_z0_s;
    TH1D* m_trk_z0sinT;
    TH1D* m_trk_chi2Prob;
    TH1D* m_trk_charge;
    TH1D* m_trk_nSi; //!
    TH1D* m_trk_nSiAndDead; //!
    TH1D* m_trk_nSiDead; //!
    TH1D* m_trk_nSCT; //!
    TH1D* m_trk_nPix; //!
    TH1D* m_trk_nPixHoles; //!
    TH1D* m_trk_nPixDead; //!
    TH1D* m_trk_nBL; //!
    TH1D* m_trk_nTRT; //!
    TH2D* m_trk_nTRT_vs_p; //!
    TH2D* m_trk_nTRT_vs_eta; //!

    TH1D* m_trk_DR_EMB2_ID; //!
    TH1D* m_trk_DEta_EMB2_ID; //!
    TH1D* m_trk_DPhi_EMB2_ID; //!
    TH2D* m_trk_DR_EMB2_ID_vs_trk_p; //!
    TH2D* m_trk_DEta_EMB2_ID_vs_trk_p; //!
    TH2D* m_trk_DPhi_EMB2_ID_vs_trk_p; //!

    // Tile energy fractions
    TH1D* m_trk_TileEfrac_100; //! 
    TH2D* m_trk_TileEfrac_100_vs_trk_p; //! 
    TH1D* m_trk_TileEfrac_200; //! 
    TH2D* m_trk_TileEfrac_200_vs_trk_p; //! 

    // basic "sanity" checks of the E/p xAOD derivation
    TH1D* m_trk_sumE_Tile_100; //!
    TH1D* m_trk_sumE_Tile_200; //!
    TH1D* m_trk_sumE_Lar_100; //!
    TH1D* m_trk_sumE_Lar_200; //!
    TH1D* m_trk_SumTileLayers_over_HAD_100; //! 
    TH1D* m_trk_SumLarLayers_over_EM_100; //!
    TH1D* m_trk_SumHADLayers_over_HAD_100; //!
    TH1D* m_trk_SumEMLayers_over_EM_100; //!
    TH1D* m_trk_EMandHAD_over_Total_100; //! 
    TH1D* m_trk_SumAllLayers_over_Total_100; //! 
    TH1D* m_trk_SumTileLayers_over_HAD_200; //! 
    TH1D* m_trk_SumLarLayers_over_EM_200; //!
    TH1D* m_trk_SumHADLayers_over_HAD_200; //!
    TH1D* m_trk_SumEMLayers_over_EM_200; //!
    TH1D* m_trk_EMandHAD_over_Total_200; //! 
    TH1D* m_trk_SumAllLayers_over_Total_200; //! 

    // calo layer with the highest cluster energy
    TH1D* m_trk_E_100_highElayer; //!
    TH1D* m_trk_E_200_highElayer; //!
    TH2D* m_trk_E_100_highElayer_vs_E; //!
    TH2D* m_trk_E_200_highElayer_vs_E; //!

    // sum of the energy of clusters matched to tracks vs. calo layer
    TH2D* m_trk_E_100_vs_layer; //! 
    TH2D* m_trk_E_200_vs_layer; //! 

    // zero fraction histograms
    TH2D* m_trk_n_E_200; //!
    TH2D* m_trk_n_E_200_l0; //!
    TH2D* m_trk_n_E_200_eq0; //!
    TH2D* m_trk_n_E_200_leq0; //!

    // E/p histos

    // total calorimeter energy
    TH1D* m_trk_E_Total_100; //!
    TH2D* m_trk_E_Total_100_vs_mu_avg; //!
    TH1D* m_eop_Total_100; //!
    TH2D* m_eop_Total_100_vs_trkP; //!
    TH2D* m_eop_Total_100_vs_trkEta; //!
    TH2D* m_eop_Total_100_vs_trkPhi; //!
    TH2D* m_eop_Total_100_vs_trkPhiID; //!
    TH2D* m_eop_Total_100_vs_trkPhi_extra; //!
    TH2D* m_eop_Total_100_vs_trkPhi_extra2; //!
    TH2D* m_eop_Total_100_vs_mu; //!
    TH2D* m_eop_Total_100_vs_mu_avg; //!
    TH2D* m_eop_Total_100_vs_npv; //!
    TH2D* m_eop_Total_100_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_Total_100_EtaEnergyRanges; //!

    TH1D* m_trk_E_Total_200; //!
    TH2D* m_trk_E_Total_200_vs_mu_avg; //!
    TH1D* m_eop_Total_200; //!
    TH1D* m_eop_Total_200_l; //!
    TH2D* m_eop_Total_200_vs_trkP; //!
    TH2D* m_eop_Total_200_vs_trkEta; //!
    TH2D* m_eop_Total_200_vs_trkPhi; //!
    TH2D* m_eop_Total_200_vs_trkPhiID; //!
    TH2D* m_eop_Total_200_vs_trkPhi_extra; //!
    TH2D* m_eop_Total_200_vs_trkPhi_extra2; //!
    TH2D* m_eop_Total_200_vs_mu; //!
    TH2D* m_eop_Total_200_vs_mu_avg; //!
    TH2D* m_eop_Total_200_vs_npv; //!
    TH2D* m_eop_Total_200_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_Total_200_EtaEnergyRanges; //!

    // EM calorimeter (EMB+EMEC)
    TH1D* m_trk_E_EM_100; //!
    TH1D* m_eop_EM_100; //!
    TH2D* m_eop_EM_100_vs_trkP; //!
    TH2D* m_eop_EM_100_vs_trkEta; //!
    TH2D* m_eop_EM_100_vs_trkPhi; //!
    TH2D* m_eop_EM_100_vs_mu; //!
    TH2D* m_eop_EM_100_vs_mu_avg; //!
    TH2D* m_eop_EM_100_vs_npv; //!
    TH2D* m_eop_EM_100_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_EM_100_EtaEnergyRanges; //!

    TH1D* m_trk_E_EM_200; //!
    TH1D* m_eop_EM_200; //!
    TH2D* m_eop_EM_200_vs_trkP; //!
    TH2D* m_eop_EM_200_vs_trkEta; //!
    TH2D* m_eop_EM_200_vs_trkPhi; //!
    TH2D* m_eop_EM_200_vs_mu; //!
    TH2D* m_eop_EM_200_vs_mu_avg; //!
    TH2D* m_eop_EM_200_vs_npv; //!
    TH2D* m_eop_EM_200_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_EM_200_EtaEnergyRanges; //!

    // No clusters in HAD matched to trk 
    TH1D* m_eop_Total_200_noHAD; //!
    TH1D* m_eop_EM_200_noHAD; //!
    TH1D* m_eop_HAD_200_noHAD; //!

    // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
    TH1D* m_trk_E_HAD_100; //!
    TH1D* m_eop_HAD_100; //!
    TH2D* m_eop_HAD_100_vs_trkP; //!
    TH2D* m_eop_HAD_100_vs_trkEta; //!
    TH2D* m_eop_HAD_100_vs_trkPhi; //!
    TH2D* m_eop_HAD_100_vs_mu; //!
    TH2D* m_eop_HAD_100_vs_mu_avg; //!
    TH2D* m_eop_HAD_100_vs_npv; //!
    TH2D* m_eop_HAD_100_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_HAD_100_EtaEnergyRanges; //!

    TH1D* m_trk_E_HAD_200; //!
    TH1D* m_eop_HAD_200; //!
    TH2D* m_eop_HAD_200_vs_trkP; //!
    TH2D* m_eop_HAD_200_vs_trkEta; //!
    TH2D* m_eop_HAD_200_vs_trkPhi; //!
    TH2D* m_eop_HAD_200_vs_mu; //!
    TH2D* m_eop_HAD_200_vs_mu_avg; //!
    TH2D* m_eop_HAD_200_vs_npv; //!
    TH2D* m_eop_HAD_200_vs_highE_layer; //!
    std::vector<std::vector<TH1D*> > m_eop_HAD_200_EtaEnergyRanges; //!

    // MIP requirement for HAD calorimeter
    TH1D* m_eop_Total_200_MIP; //!
    TH1D* m_eop_EM_200_MIP; //!
    TH1D* m_eop_HAD_200_MIP; //!

    // background subtraction, the way it was done in run 1
    TH1D* m_eop_EM_BG; //!
    TH2D* m_eop_EM_BG_vs_trkP; //!
    TH2D* m_eop_EM_BG_vs_trkEta; //!
    TH2D* m_eop_EM_BG_vs_trkPhi; //!
    TH2D* m_eop_EM_BG_vs_mu; //!
    TH2D* m_eop_EM_BG_vs_mu_avg; //!
    TH2D* m_eop_EM_BG_vs_npv; //!
    std::vector<std::vector<TH1D*> > m_eop_EM_BG_EtaEnergyRanges; //!

    // E/p where E is maximum in either Tile layer A, BC, or D
    // max in Tile Layer A
    TH1D* m_trk_E_highTileA_100; //!
    TH1D* m_eop_highTileA_100; //!
    TH2D* m_eop_highTileA_100_vs_trkP; //!
    TH2D* m_eop_highTileA_100_vs_trkEta; //!
    TH2D* m_eop_highTileA_100_vs_trkPhi; //!
    TH2D* m_eop_highTileA_100_vs_mu; //!
    TH2D* m_eop_highTileA_100_vs_mu_avg; //!
    TH2D* m_eop_highTileA_100_vs_npv; //!
    TH2D* m_eop_highTileA_100_vs_TileA_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileA_100_EtaEnergyRanges; //!
    TH1D* m_trk_E_highTileA_200; //!
    TH1D* m_eop_highTileA_200; //!
    TH2D* m_eop_highTileA_200_vs_trkP; //!
    TH2D* m_eop_highTileA_200_vs_trkEta; //!
    TH2D* m_eop_highTileA_200_vs_trkPhi; //!
    TH2D* m_eop_highTileA_200_vs_mu; //!
    TH2D* m_eop_highTileA_200_vs_mu_avg; //!
    TH2D* m_eop_highTileA_200_vs_npv; //!
    TH2D* m_eop_highTileA_200_vs_TileA_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileA_200_EtaEnergyRanges; //!
    // max in Tile Layer BC
    TH1D* m_trk_E_highTileB_100; //!
    TH1D* m_eop_highTileB_100; //!
    TH2D* m_eop_highTileB_100_vs_trkP; //!
    TH2D* m_eop_highTileB_100_vs_trkEta; //!
    TH2D* m_eop_highTileB_100_vs_trkPhi; //!
    TH2D* m_eop_highTileB_100_vs_mu; //!
    TH2D* m_eop_highTileB_100_vs_mu_avg; //!
    TH2D* m_eop_highTileB_100_vs_npv; //!
    TH2D* m_eop_highTileB_100_vs_TileB_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileB_100_EtaEnergyRanges; //!
    TH1D* m_trk_E_highTileB_200; //!
    TH1D* m_eop_highTileB_200; //!
    TH2D* m_eop_highTileB_200_vs_trkP; //!
    TH2D* m_eop_highTileB_200_vs_trkEta; //!
    TH2D* m_eop_highTileB_200_vs_trkPhi; //!
    TH2D* m_eop_highTileB_200_vs_mu; //!
    TH2D* m_eop_highTileB_200_vs_mu_avg; //!
    TH2D* m_eop_highTileB_200_vs_npv; //!
    TH2D* m_eop_highTileB_200_vs_TileB_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileB_200_EtaEnergyRanges; //!
    // max in Tile Layer D
    TH1D* m_trk_E_highTileD_100; //!
    TH1D* m_eop_highTileD_100; //!
    TH2D* m_eop_highTileD_100_vs_trkP; //!
    TH2D* m_eop_highTileD_100_vs_trkEta; //!
    TH2D* m_eop_highTileD_100_vs_trkPhi; //!
    TH2D* m_eop_highTileD_100_vs_mu; //!
    TH2D* m_eop_highTileD_100_vs_mu_avg; //!
    TH2D* m_eop_highTileD_100_vs_npv; //!
    TH2D* m_eop_highTileD_100_vs_TileD_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileD_100_EtaEnergyRanges; //!
    TH1D* m_trk_E_highTileD_200; //!
    TH1D* m_eop_highTileD_200; //!
    TH2D* m_eop_highTileD_200_vs_trkP; //!
    TH2D* m_eop_highTileD_200_vs_trkEta; //!
    TH2D* m_eop_highTileD_200_vs_trkPhi; //!
    TH2D* m_eop_highTileD_200_vs_mu; //!
    TH2D* m_eop_highTileD_200_vs_mu_avg; //!
    TH2D* m_eop_highTileD_200_vs_npv; //!
    TH2D* m_eop_highTileD_200_vs_TileD_E; //!
    std::vector<std::vector<TH1D*> > m_eop_highTileD_200_EtaEnergyRanges; //!

};

#endif
