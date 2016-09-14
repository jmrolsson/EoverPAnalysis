// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#ifndef EoverPAnalysis_EoverPAnalysis_H
#define EoverPAnalysis_EoverPAnalysis_H

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// Histograms
#include "EoverPAnalysis/EoverPHists.h"
#include "EoverPAnalysis/EoverPHistsTrks.h"

class EoverPAnalysis : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    std::string m_inTrackContainerName;

    // configuration variables
    std::string m_detailStr;

    // calorimeter layer to which tracks are extrapolated (for track eta, phi)
    std::string m_trkExtrapol ="EMB2";

    // track isolation settings
    float m_trkIsoDRmax = .4;
    float m_trkIsoPfrac = 0.;

    // what plots to make (eop from clusters in entire calo, EM calo, HAD calo, with background subtraction, for maximum energy in tile)
    bool m_doCaloTotal= true;
    bool m_doCaloEM= false;
    bool m_doCaloHAD= false;
    bool m_doBgSubtr = false;
    bool m_doTileLayer = false;

    // global track p cuts
    bool m_doTrkPcut = false;
    float m_trkPmin = 0.;
    float m_trkPmax = 1e8; 

    // global track p cuts
    bool m_doTrkEtacut = false;
    float m_trkEtamin = 0.;
    float m_trkEtamax = 1e8; 

    // make plots with Tile specific cuts
    bool m_doTileCuts = false;
    float m_LarEmax = 1e8;
    float m_TileEfracmin = -1;

    // user defined energy (p) bins
    std::string m_Ebins = "";
    bool m_doEbinsArray = false;
    std::string m_EbinsArray = "";

    // user defined eta bins
    std::string m_Etabins = "";
    bool m_doEtabinsArray = false;
    std::string m_EtabinsArray = "";

    // turn on extra E/p histograms for each of the track eta and p bins
    bool m_doExtraEtaEnergyBinHists = false;

    // make all plots for different TileEfrac cuts
    bool m_doGlobalTileEfracRanges = false;

    // make all plots for different track p and eta selections 
    bool m_doGlobalEnergyRanges = false; 
    bool m_doGlobalEtaRanges = false; 

    // extra ranges for comparisons with Run 1 studies
    bool m_doGlobalExtraRanges = false; 

    // extra histograms for track isolation testing
    bool m_doTrkIsoHists = true;

    // turn on cutflows
    bool m_useCutFlow = false; 

    // energy calibration, either "ClusterEnergy", "ClusterEnergyLCW", or "CellEnergy"
    std::string m_energyCalib = "ClusterEnergy";

    // pileup reweighting
    bool m_doCustomPUreweighting = false;

  private:

    // cutflow
    TH1D* m_cutflowHist; //!
    TH1D* m_cutflowHistW; //!
    int   m_cutflow_bin; //!

    // pileup reweighting
    TH1D* m_puwHist; //! 

    /* object-level cutflow */

    TH1D* m_trk_cutflowHist_1;  //!
    TH1D* m_trk_cutflowHist_eop;  //!

    int m_trk_cutflow_eop_all_bin; //!
    int m_trk_cutflow_eop_pass_iso_bin; //!
    int m_trk_cutflow_eop_pass_p_bin; //!
    int m_trk_cutflow_eop_pass_eta_bin; //!
    int m_trk_cutflow_eop_pass_larEmax_bin; //!
    int m_trk_cutflow_eop_pass_tileEfrac_bin; //!
    int m_trk_cutflow_eop_all; //!
    int m_trk_cutflow_eop_pass_p; //!
    int m_trk_cutflow_eop_pass_eta; //!
    int m_trk_cutflow_eop_pass_iso; //!
    int m_trk_cutflow_eop_pass_larEmax; //!
    int m_trk_cutflow_eop_pass_tileEfrac; //!

    int m_numEvent;         //!
    int m_numEventPass;     //!
    int m_weightNumEventPass; //!

    // number of tracks per event, after each selection
    TH1D* m_trk_n_all; //!
    TH1D* m_trk_n_pass_p; //!
    TH1D* m_trk_n_pass_eta; //!
    TH1D* m_trk_n_pass_iso; //!
    TH1D* m_trk_n_pass_larEmax; //!
    TH1D* m_trk_n_pass_tileEfrac; //!
    int m_trk_n_all_tmp; //!
    int m_trk_n_pass_p_tmp; //!
    int m_trk_n_pass_eta_tmp; //!
    int m_trk_n_pass_iso_tmp; //!
    int m_trk_n_pass_larEmax_tmp; //!
    int m_trk_n_pass_tileEfrac_tmp; //!

    // list of calo layers
    const std::vector<std::string> m_layer = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3", "TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of all the calo layers
    const std::vector<std::string> m_layer_lar = {"EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of lar layers only
    const std::vector<std::string> m_layer_tile = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of tile layers only


    // histogram making classes 
    EoverPHists* m_plots_eop; //!
    EoverPHistsTrks* m_plots_eop_trks; //!

    EoverPHists* m_plots_eop_TileEfrac000; //!
    EoverPHists* m_plots_eop_TileEfrac010; //!
    EoverPHists* m_plots_eop_TileEfrac030; //!
    EoverPHists* m_plots_eop_TileEfrac050; //!
    EoverPHists* m_plots_eop_TileEfrac060; //!
    EoverPHists* m_plots_eop_TileEfrac070; //!
    EoverPHists* m_plots_eop_TileEfrac075; //!
    EoverPHists* m_plots_eop_TileEfrac080; //!

    EoverPHists* m_plots_eop_pL4000; //!
    EoverPHists* m_plots_eop_pG4000L8000; //!
    EoverPHists* m_plots_eop_pG8000L12000; //!
    EoverPHists* m_plots_eop_pG12000; //!

    EoverPHists* m_plots_eop_etaL05; //!
    EoverPHists* m_plots_eop_etaG05L07; //!
    EoverPHists* m_plots_eop_etaG07; //!

    // extra ranges for comparisons with Run 1 studies
    EoverPHists* m_plots_eop_pG1200L1800; //!
    EoverPHists* m_plots_eop_pG1800L2200; //!
    EoverPHists* m_plots_eop_pG2200L2800; //!
    EoverPHists* m_plots_eop_pG2800L3600; //!
    EoverPHists* m_plots_eop_pG3600L4600; //!
    EoverPHists* m_plots_eop_pG4600L5600; //!
    EoverPHists* m_plots_eop_etaL06_pG2200L4600; //!
    EoverPHists* m_plots_eop_etaL06_pG4600L50000; //!
    EoverPHists* m_plots_eop_etaL06; //!
    EoverPHists* m_plots_eop_etaG06L11; //!
    EoverPHists* m_plots_eop_etaG11L14; //!
    EoverPHists* m_plots_eop_etaG14L15; //!
    EoverPHists* m_plots_eop_etaG15L18; //!

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
  public:
    // Tree *myTree; //!
    // TH1 *myHist; //!

    // this is a standard constructor
    EoverPAnalysis (std::string className = "EoverPAnalysis");

    // these are the functions inherited from Algorithm
    virtual EL::StatusCode setupJob (EL::Job& job);
    virtual EL::StatusCode fileExecute ();
    virtual EL::StatusCode histInitialize ();
    virtual EL::StatusCode changeInput (bool firstFile);
    virtual EL::StatusCode initialize ();
    virtual EL::StatusCode execute ();
    virtual EL::StatusCode postExecute ();
    virtual EL::StatusCode finalize ();
    virtual EL::StatusCode histFinalize ();

    /// @cond
    // this is needed to distribute the algorithm to the workers
    ClassDef(EoverPAnalysis, 1);
    /// @endcond
};

#endif
