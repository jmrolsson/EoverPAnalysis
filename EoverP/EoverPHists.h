// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#ifndef EoverP_EoverPHists_H
#define EoverP_EoverPHists_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

class EoverPHists : public HistogramManager
{
  public:
    EoverPHists(std::string name, std::string detailStr, std::string energyCalib = "ClusterEnergy", std::string trkExtrapol = "EMB2", bool doCaloTotal = true, bool doCaloEM = false, bool doCaloHAD = false, bool doBgSubtr = true, bool doTileLayer = false, std::string Ebins = "", bool doEbinsArray = false, std::string EbinsArray = "", std::string Etabins = "", bool doEtabinsArray = false, std::string EtabinsArray = "", bool doExtraEtaEnergyBinHists = false);
    ~EoverPHists();

    StatusCode initialize();

    StatusCode execute( const xAOD::TrackParticle* trk, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload

    // Double_t* linspace(float a, float b, unsigned int n);
    // Double_t* logspace(float a, float b, unsigned int n);
    std::vector<double> str2vec(std::string str);

  protected:

    // bools to control which histograms are filled
    bool m_doCaloTotal; //!
    bool m_doCaloEM; //!
    bool m_doCaloHAD; //!
    bool m_doBgSubtr; //!
    bool m_doTileLayer; //!

    std::string m_trkExtrapol; //! layer where tracks are extrapolated
    std::string m_energyCalib; //! what type of energy calibration to use (EM, LCW, cells)
    std::string m_Ebins; //!
    bool m_doEbinsArray; //!
    std::string m_EbinsArray; //!
    std::string m_Etabins; //!
    bool m_doEtaAbs; //!
    bool m_doEtabinsArray; //!
    std::string m_EtabinsArray; //!
    std::vector<double> EbinsArray = {0}; //!
    unsigned int nEbinsArray; //!
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
    TH1F* m_trk_p; //!
    TH1F* m_trk_p_array; //!
    TH1F* m_trk_eta; //!
    TH1F* m_trk_eta_array; //!
    TH1F* m_trk_phi; //!

    TH1F* m_trk_DR_CALO_ID; //!
    TH1F* m_trk_DEta_CALO_ID; //!
    TH1F* m_trk_DPhi_CALO_ID; //!
    TH2F* m_trk_DR_CALO_ID_vs_trk_p; //!
    TH2F* m_trk_DEta_CALO_ID_vs_trk_p; //!
    TH2F* m_trk_DPhi_CALO_ID_vs_trk_p; //!

    // Tile energy fractions
    TH1F* m_trk_TileEfrac_100; //! 
    TH2F* m_trk_TileEfrac_100_vs_trk_p; //! 
    TH1F* m_trk_TileEfrac_200; //! 
    TH2F* m_trk_TileEfrac_200_vs_trk_p; //! 

    // basic "sanity" checks of the E/p xAOD derivation
    TH1F* m_trk_sumE_Tile_100; //!
    TH1F* m_trk_sumE_Tile_200; //!
    TH1F* m_trk_sumE_Lar_100; //!
    TH1F* m_trk_sumE_Lar_200; //!
    TH1F* m_trk_SumTileLayers_over_HAD_100; //! 
    TH1F* m_trk_SumLarLayers_over_EM_100; //!
    TH1F* m_trk_SumHADLayers_over_HAD_100; //!
    TH1F* m_trk_SumEMLayers_over_EM_100; //!
    TH1F* m_trk_EMandHAD_over_Total_100; //! 
    TH1F* m_trk_SumAllLayers_over_Total_100; //! 
    TH1F* m_trk_SumTileLayers_over_HAD_200; //! 
    TH1F* m_trk_SumLarLayers_over_EM_200; //!
    TH1F* m_trk_SumHADLayers_over_HAD_200; //!
    TH1F* m_trk_SumEMLayers_over_EM_200; //!
    TH1F* m_trk_EMandHAD_over_Total_200; //! 
    TH1F* m_trk_SumAllLayers_over_Total_200; //! 

    // calo layer with the highest cluster energy
    TH1F* m_trk_E_100_highElayer; //!
    TH1F* m_trk_E_200_highElayer; //!
    TH2F* m_trk_E_100_highElayer_vs_E; //!
    TH2F* m_trk_E_200_highElayer_vs_E; //!

    // sum of the energy of clusters matched to tracks vs. calo layer
    TH2F* m_trk_E_100_vs_layer; //! 
    TH2F* m_trk_E_200_vs_layer; //! 

    // E/p histos

    // total calorimeter energy
    TH1F* m_trk_E_Total_100; //!
    TH1F* m_eop_Total_100; //!
    TH2F* m_eop_Total_100_vs_trkP; //!
    TH2F* m_eop_Total_100_vs_trkEta; //!
    TH2F* m_eop_Total_100_vs_trkPhi; //!
    TH2F* m_eop_Total_100_vs_mu; //!
    TH2F* m_eop_Total_100_vs_mu_avg; //!
    TH2F* m_eop_Total_100_vs_npv; //!
    TH2F* m_eop_Total_100_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_Total_100_EtaEnergyRanges; //!

    TH1F* m_trk_E_Total_200; //!
    TH1F* m_eop_Total_200; //!
    TH2F* m_eop_Total_200_vs_trkP; //!
    TH2F* m_eop_Total_200_vs_trkEta; //!
    TH2F* m_eop_Total_200_vs_trkPhi; //!
    TH2F* m_eop_Total_200_vs_mu; //!
    TH2F* m_eop_Total_200_vs_mu_avg; //!
    TH2F* m_eop_Total_200_vs_npv; //!
    TH2F* m_eop_Total_200_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_Total_200_EtaEnergyRanges; //!

    // EM calorimeter (EMB+EMEC)
    TH1F* m_trk_E_EM_100; //!
    TH1F* m_eop_EM_100; //!
    TH2F* m_eop_EM_100_vs_trkP; //!
    TH2F* m_eop_EM_100_vs_trkEta; //!
    TH2F* m_eop_EM_100_vs_trkPhi; //!
    TH2F* m_eop_EM_100_vs_mu; //!
    TH2F* m_eop_EM_100_vs_mu_avg; //!
    TH2F* m_eop_EM_100_vs_npv; //!
    TH2F* m_eop_EM_100_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_EM_100_EtaEnergyRanges; //!

    TH1F* m_trk_E_EM_200; //!
    TH1F* m_eop_EM_200; //!
    TH2F* m_eop_EM_200_vs_trkP; //!
    TH2F* m_eop_EM_200_vs_trkEta; //!
    TH2F* m_eop_EM_200_vs_trkPhi; //!
    TH2F* m_eop_EM_200_vs_mu; //!
    TH2F* m_eop_EM_200_vs_mu_avg; //!
    TH2F* m_eop_EM_200_vs_npv; //!
    TH2F* m_eop_EM_200_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_EM_200_EtaEnergyRanges; //!

    // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
    TH1F* m_trk_E_HAD_100; //!
    TH1F* m_eop_HAD_100; //!
    TH2F* m_eop_HAD_100_vs_trkP; //!
    TH2F* m_eop_HAD_100_vs_trkEta; //!
    TH2F* m_eop_HAD_100_vs_trkPhi; //!
    TH2F* m_eop_HAD_100_vs_mu; //!
    TH2F* m_eop_HAD_100_vs_mu_avg; //!
    TH2F* m_eop_HAD_100_vs_npv; //!
    TH2F* m_eop_HAD_100_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_HAD_100_EtaEnergyRanges; //!

    TH1F* m_trk_E_HAD_200; //!
    TH1F* m_eop_HAD_200; //!
    TH2F* m_eop_HAD_200_vs_trkP; //!
    TH2F* m_eop_HAD_200_vs_trkEta; //!
    TH2F* m_eop_HAD_200_vs_trkPhi; //!
    TH2F* m_eop_HAD_200_vs_mu; //!
    TH2F* m_eop_HAD_200_vs_mu_avg; //!
    TH2F* m_eop_HAD_200_vs_npv; //!
    TH2F* m_eop_HAD_200_vs_highE_layer; //!
    std::vector<std::vector<TH1F*> > m_eop_HAD_200_EtaEnergyRanges; //!

    // background subtraction, the way it was done in run 1
    TH1F* m_eop_EM_BG; //!
    TH2F* m_eop_EM_BG_vs_trkP; //!
    TH2F* m_eop_EM_BG_vs_trkEta; //!
    TH2F* m_eop_EM_BG_vs_trkPhi; //!
    TH2F* m_eop_EM_BG_vs_mu; //!
    TH2F* m_eop_EM_BG_vs_mu_avg; //!
    TH2F* m_eop_EM_BG_vs_npv; //!
    std::vector<std::vector<TH1F*> > m_eop_EM_BG_EtaEnergyRanges; //!

    // E/p where E is maximum in either Tile layer A, BC, or D
    // max in Tile Layer A
    TH1F* m_trk_E_highTileA_100; //!
    TH1F* m_eop_highTileA_100; //!
    TH2F* m_eop_highTileA_100_vs_trkP; //!
    TH2F* m_eop_highTileA_100_vs_trkEta; //!
    TH2F* m_eop_highTileA_100_vs_trkPhi; //!
    TH2F* m_eop_highTileA_100_vs_mu; //!
    TH2F* m_eop_highTileA_100_vs_mu_avg; //!
    TH2F* m_eop_highTileA_100_vs_npv; //!
    TH2F* m_eop_highTileA_100_vs_TileA_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileA_100_EtaEnergyRanges; //!
    TH1F* m_trk_E_highTileA_200; //!
    TH1F* m_eop_highTileA_200; //!
    TH2F* m_eop_highTileA_200_vs_trkP; //!
    TH2F* m_eop_highTileA_200_vs_trkEta; //!
    TH2F* m_eop_highTileA_200_vs_trkPhi; //!
    TH2F* m_eop_highTileA_200_vs_mu; //!
    TH2F* m_eop_highTileA_200_vs_mu_avg; //!
    TH2F* m_eop_highTileA_200_vs_npv; //!
    TH2F* m_eop_highTileA_200_vs_TileA_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileA_200_EtaEnergyRanges; //!
    // max in Tile Layer BC
    TH1F* m_trk_E_highTileB_100; //!
    TH1F* m_eop_highTileB_100; //!
    TH2F* m_eop_highTileB_100_vs_trkP; //!
    TH2F* m_eop_highTileB_100_vs_trkEta; //!
    TH2F* m_eop_highTileB_100_vs_trkPhi; //!
    TH2F* m_eop_highTileB_100_vs_mu; //!
    TH2F* m_eop_highTileB_100_vs_mu_avg; //!
    TH2F* m_eop_highTileB_100_vs_npv; //!
    TH2F* m_eop_highTileB_100_vs_TileB_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileB_100_EtaEnergyRanges; //!
    TH1F* m_trk_E_highTileB_200; //!
    TH1F* m_eop_highTileB_200; //!
    TH2F* m_eop_highTileB_200_vs_trkP; //!
    TH2F* m_eop_highTileB_200_vs_trkEta; //!
    TH2F* m_eop_highTileB_200_vs_trkPhi; //!
    TH2F* m_eop_highTileB_200_vs_mu; //!
    TH2F* m_eop_highTileB_200_vs_mu_avg; //!
    TH2F* m_eop_highTileB_200_vs_npv; //!
    TH2F* m_eop_highTileB_200_vs_TileB_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileB_200_EtaEnergyRanges; //!
    // max in Tile Layer D
    TH1F* m_trk_E_highTileD_100; //!
    TH1F* m_eop_highTileD_100; //!
    TH2F* m_eop_highTileD_100_vs_trkP; //!
    TH2F* m_eop_highTileD_100_vs_trkEta; //!
    TH2F* m_eop_highTileD_100_vs_trkPhi; //!
    TH2F* m_eop_highTileD_100_vs_mu; //!
    TH2F* m_eop_highTileD_100_vs_mu_avg; //!
    TH2F* m_eop_highTileD_100_vs_npv; //!
    TH2F* m_eop_highTileD_100_vs_TileD_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileD_100_EtaEnergyRanges; //!
    TH1F* m_trk_E_highTileD_200; //!
    TH1F* m_eop_highTileD_200; //!
    TH2F* m_eop_highTileD_200_vs_trkP; //!
    TH2F* m_eop_highTileD_200_vs_trkEta; //!
    TH2F* m_eop_highTileD_200_vs_trkPhi; //!
    TH2F* m_eop_highTileD_200_vs_mu; //!
    TH2F* m_eop_highTileD_200_vs_mu_avg; //!
    TH2F* m_eop_highTileD_200_vs_npv; //!
    TH2F* m_eop_highTileD_200_vs_TileD_E; //!
    std::vector<std::vector<TH1F*> > m_eop_highTileD_200_EtaEnergyRanges; //!

};

#endif
