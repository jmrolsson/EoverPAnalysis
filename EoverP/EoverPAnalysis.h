#ifndef EoverP_EoverPAnalysis_H
#define EoverP_EoverPAnalysis_H

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// Histograms
#include "EoverP/EoverPHists.h"

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

    // what type of clusters/cells to use for the energy sum
    bool m_doEMcalib = true;
    bool m_doLCWcalib = false;
    bool m_doCells = false;

    // track isolation settings
    float m_trkIsoDRmax = .4;
    float m_trkIsoPfrac = 0.;

    // what plots to make (eop from clusters in entire calo, EM calo, HAD calo, with background subtraction, for maximum energy in tile)
    bool m_doCaloTotal= true;
    bool m_doCaloEM= false;
    bool m_doCaloHAD= false;
    bool m_doBgSubtr = true;
    bool m_doTileLayer = true;

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

    // make plots for the energy (p) and eta bins as defined above
    // (if no bins are specified by the user, then the default bins are used)
    bool m_doEtaEnergyRanges = false;

  private:
    const std::vector<std::string> m_layer = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3", "TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of all the calo layers
    const std::vector<std::string> m_layer_lar = {"EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of lar layers only
    const std::vector<std::string> m_layer_tile = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of tile layers only

    std::string m_energyCalib = "ClusterEnergy";
    EoverPHists* m_plots_eop; //!

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
