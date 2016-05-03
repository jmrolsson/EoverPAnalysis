#ifndef EoverP_EoverPHists_eopxAOD_H
#define EoverP_EoverPHists_eopxAOD_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

class EoverPHists_eopxAOD : public HistogramManager
{
  public:
    EoverPHists_eopxAOD(std::string name, std::string detailStr, std::string trkExtrapol, bool doBgSubtr, bool doEMcalib, bool doLCWcalib, bool doCells, bool doCaloEM, bool doCaloHAD, bool doCaloTotal, float trkIsoDRmax, float trkIsoPfrac, float LarEmax, float TileEfrac, bool doTrkPcut, float trkPmin, float trkPmax, bool doTrkEtacut, float trkEtamin, float trkEtamax);
    ~EoverPHists_eopxAOD();

    StatusCode initialize();

    StatusCode execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload

    Double_t* linspace(float a, float b, unsigned int n);
    Double_t* logspace(float a, float b, unsigned int n);

  protected:
    // bools to control which histograms are filled
    bool m_fillDebugging; //!

  private:

    //LAYER = PreSamplerB/E, EMB1/2/3, EME1/2/3, HEC0/1/2/3, TileBar0/1/2, TileGap1/2/3 and TileExt0/1/2
    const std::vector<std::string> m_layer = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3", "TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of all the calo layers
    const std::vector<std::string> m_layer_lar = {"EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of lar layers only
    const std::vector<std::string> m_layer_tile = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of tile layers only
    std::string m_trkExtrapol; //! layer where tracks are extrapolated
    bool m_doBgSubtr; //!
    bool m_doEMcalib; //!
    bool m_doLCWcalib; //!
    bool m_doCells; //!
    bool m_doCaloEM; //!
    bool m_doCaloHAD; //!
    bool m_doCaloTotal; //!
    float m_trkIsoDRmax; //! track isolation max DR
    float m_trkIsoPfrac; //! track isolation max p fraction
    float m_LarEmax; //! maximum cluster energy in LAr
    float m_TileEfrac; //! minimum energy fraction in Tile calorimeter
    bool m_doTrkPcut; //! apply track momentum cut
    float m_trkPmin; //!
    float m_trkPmax; //! 
    bool m_doTrkEtacut; //! apply track eta cut
    float m_trkEtamin; //! 
    float m_trkEtamax; //! 

    // event level plots
    TH1F* m_mu; //!
    TH1F* m_mu_avg; //!
    TH1F* m_npv; //!

    // track plots
    TH1F* m_trk_n_nocut; //!
    TH1F* m_trk_n; //!
    TH1F* m_trk_p; //!
    TH1F* m_trk_eta; //!
    TH1F* m_trk_eta_abs; //!
    TH1F* m_trk_phi; //!

    TH1F* m_trk_DR_CALO_ID; //!
    TH1F* m_trk_DEta_CALO_ID; //!
    TH1F* m_trk_DPhi_CALO_ID; //!

    TH1F* m_trk_ntrks_maxDR01; //!
    TH1F* m_trk_ntrks_maxDR02; //!
    TH1F* m_trk_ntrks_maxDR03; //!
    TH1F* m_trk_ntrks_maxDR04; //!
    TH1F* m_trk_ntrks_maxDR05; //!
    TH1F* m_trk_ntrks_maxDR06; //!
    TH1F* m_trk_ntrks_maxDR07; //!
    TH1F* m_trk_ntrks_maxDR08; //!
    TH1F* m_trk_ntrks_maxDR09; //!
    TH1F* m_trk_ntrks_maxDR10; //!

    TH2F* m_trk_ntrks_maxDR01_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR02_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR03_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR04_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR05_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR06_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR07_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR08_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR09_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR10_vs_trk_p; //!

    TH2F* m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p; //!

    TH2F* m_ClusterEnergy_100_vs_layer; //! 
    TH2F* m_ClusterEnergy_200_vs_layer; //! 
    TH2F* m_ClusterEnergy_100_vs_layer_passTrkIso; //! 
    TH2F* m_ClusterEnergy_200_vs_layer_passTrkIso; //! 
    TH1F* m_ClusterEnergy_100_highestEnergyLayer; //!
    TH1F* m_ClusterEnergy_200_highestEnergyLayer; //!
    TH1F* m_ClusterEnergy_100_highestEnergyLayer_passTrkIso; //!
    TH1F* m_ClusterEnergy_200_highestEnergyLayer_passTrkIso; //!

    TH1F* m_trk_TileEfrac_100; //! 
    TH2F* m_trk_TileEfrac_100_vs_trk_p; //! 
    TH1F* m_trk_TileEfrac_200; //! 
    TH2F* m_trk_TileEfrac_200_vs_trk_p; //! 

    TH1F* m_trk_SumTileLayers_over_HAD_100; //! 
    TH1F* m_trk_SumLarLayers_over_EM_100; //!
    TH1F* m_trk_EMandHAD_over_Total_100; //! 
    TH1F* m_trk_SumAllLayers_over_Total_100; //! 
    TH1F* m_trk_SumTileLayers_over_HAD_200; //! 
    TH1F* m_trk_SumLarLayers_over_EM_200; //!
    TH1F* m_trk_EMandHAD_over_Total_200; //! 
    TH1F* m_trk_SumAllLayers_over_Total_200; //! 

    TH2F* m_trk_DR_CALO_ID_vs_trk_p; //!
    TH2F* m_trk_DEta_CALO_ID_vs_trk_p; //!
    TH2F* m_trk_DPhi_CALO_ID_vs_trk_p; //!

    // clusters calibrated at the EM-scale
    TH1F* m_trk_matched_Total_ClusterEnergy_200; //!
    TH1F* m_eop_Total_ClusterEnergy_200; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_trkP; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_trkEta; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_mu; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_Total_ClusterEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_EM_ClusterEnergy_200; //!
    TH1F* m_eop_EM_ClusterEnergy_200; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_trkP; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_trkEta; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_mu; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_EM_ClusterEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_HAD_ClusterEnergy_200; //!
    TH1F* m_eop_HAD_ClusterEnergy_200; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_trkP; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_trkEta; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_mu; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_HAD_ClusterEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_Total_ClusterEnergy_100; //!
    TH1F* m_eop_Total_ClusterEnergy_100; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_trkP; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_trkEta; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_mu; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_Total_ClusterEnergy_100_vs_npv; //!

    TH1F* m_trk_matched_EM_ClusterEnergy_100; //!
    TH1F* m_eop_EM_ClusterEnergy_100; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_trkP; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_trkEta; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_mu; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_EM_ClusterEnergy_100_vs_npv; //!

    TH1F* m_trk_matched_HAD_ClusterEnergy_100; //!
    TH1F* m_eop_HAD_ClusterEnergy_100; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_trkP; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_trkEta; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_mu; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_HAD_ClusterEnergy_100_vs_npv; //!

    // clusters calibrated at the LCW-scale
    TH1F* m_trk_matched_Total_ClusterEnergyLCW_200; //!
    TH1F* m_eop_Total_ClusterEnergyLCW_200; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_trkP; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_trkEta; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_trkPhi; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_mu; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_mu_avg; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_200_vs_npv; //!

    TH1F* m_trk_matched_EM_ClusterEnergyLCW_200; //!
    TH1F* m_eop_EM_ClusterEnergyLCW_200; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_trkP; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_trkEta; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_trkPhi; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_mu; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_mu_avg; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_200_vs_npv; //!

    TH1F* m_trk_matched_HAD_ClusterEnergyLCW_200; //!
    TH1F* m_eop_HAD_ClusterEnergyLCW_200; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_trkP; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_trkEta; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_trkPhi; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_mu; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_mu_avg; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_200_vs_npv; //!

    TH1F* m_trk_matched_Total_ClusterEnergyLCW_100; //!
    TH1F* m_eop_Total_ClusterEnergyLCW_100; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_trkP; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_trkEta; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_trkPhi; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_mu; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_mu_avg; //!
    TH2F* m_eop_Total_ClusterEnergyLCW_100_vs_npv; //!

    TH1F* m_trk_matched_EM_ClusterEnergyLCW_100; //!
    TH1F* m_eop_EM_ClusterEnergyLCW_100; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_trkP; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_trkEta; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_trkPhi; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_mu; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_mu_avg; //!
    TH2F* m_eop_EM_ClusterEnergyLCW_100_vs_npv; //!

    TH1F* m_trk_matched_HAD_ClusterEnergyLCW_100; //!
    TH1F* m_eop_HAD_ClusterEnergyLCW_100; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_trkP; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_trkEta; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_trkPhi; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_mu; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_mu_avg; //!
    TH2F* m_eop_HAD_ClusterEnergyLCW_100_vs_npv; //!

    // cells 
    TH1F* m_trk_matched_Total_CellEnergy_200; //!
    TH1F* m_eop_Total_CellEnergy_200; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_trkP; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_trkEta; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_mu; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_Total_CellEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_EM_CellEnergy_200; //!
    TH1F* m_eop_EM_CellEnergy_200; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_trkP; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_trkEta; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_mu; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_EM_CellEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_HAD_CellEnergy_200; //!
    TH1F* m_eop_HAD_CellEnergy_200; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_trkP; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_trkEta; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_trkPhi; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_mu; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_mu_avg; //!
    TH2F* m_eop_HAD_CellEnergy_200_vs_npv; //!

    TH1F* m_trk_matched_Total_CellEnergy_100; //!
    TH1F* m_eop_Total_CellEnergy_100; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_trkP; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_trkEta; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_mu; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_Total_CellEnergy_100_vs_npv; //!

    TH1F* m_trk_matched_EM_CellEnergy_100; //!
    TH1F* m_eop_EM_CellEnergy_100; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_trkP; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_trkEta; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_mu; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_EM_CellEnergy_100_vs_npv; //!

    TH1F* m_trk_matched_HAD_CellEnergy_100; //!
    TH1F* m_eop_HAD_CellEnergy_100; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_trkP; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_trkEta; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_trkPhi; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_mu; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_mu_avg; //!
    TH2F* m_eop_HAD_CellEnergy_100_vs_npv; //!

    // background subtraction
    TH1F* m_eop_EM_BG_ClusterEnergy; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_trkP; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_trkEta; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_trkPhi; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_mu; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_mu_avg; //!
    TH2F* m_eop_EM_BG_ClusterEnergy_vs_npv; //!

    TH1F* m_eop_EM_BG_ClusterEnergyLCW; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_trkP; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_trkEta; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_trkPhi; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_mu; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_mu_avg; //!
    TH2F* m_eop_EM_BG_ClusterEnergyLCW_vs_npv; //!

    TH1F* m_eop_EM_BG_CellEnergy; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_trkP; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_trkEta; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_trkPhi; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_mu; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_mu_avg; //!
    TH2F* m_eop_EM_BG_CellEnergy_vs_npv; //!
};

#endif
