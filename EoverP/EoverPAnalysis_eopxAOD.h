#ifndef EoverP_EoverPAnalysis_eopxAOD_H
#define EoverP_EoverPAnalysis_eopxAOD_H

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// Histograms
#include "EoverP/EoverPHists_eopxAOD.h"

class EoverPAnalysis_eopxAOD : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    std::string m_inTrackContainerName;

    // configuration variables
    std::string m_detailStr;

    std::string m_trkExtrapol ="EMB2";
    bool m_doBgSubtr = true;
    bool m_doEMcalib = true;
    bool m_doLCWcalib = false;
    bool m_doCells = false;
    bool m_doCaloEM= false;
    bool m_doCaloHAD= false;
    bool m_doCaloTotal= true;
    bool m_doEtaPranges = false;

    float m_trkIsoDRmax = .4;
    float m_trkIsoPfrac = .0;

    bool m_doTrkPcut = false;
    float m_trkPmin = 0.;
    float m_trkPmax = 1e8; 

    bool m_doTrkEtacut = false;
    float m_trkEtamin = 0.;
    float m_trkEtamax = 1e8; 

    bool m_doTileCuts = false;
    float m_LarEmax = 1e8;
    float m_TileEfrac = -1;

    std::string m_Ebins = "";
    bool m_doEbinsArray = false;
    std::string m_EbinsArray = "";

    std::string m_EtaAbsbins = "";
    bool m_doEtaAbsbinsArray = false;
    std::string m_EtaAbsbinsArray = "";

  private:
    EoverPHists_eopxAOD* m_plots_eop; //!
    EoverPHists_eopxAOD* m_plots_eop_etaL06; //!
    EoverPHists_eopxAOD* m_plots_eop_etaG06L11; //!
    EoverPHists_eopxAOD* m_plots_eop_etaG18L19; //!
    EoverPHists_eopxAOD* m_plots_eop_etaG19L23; //!
    EoverPHists_eopxAOD* m_plots_eop_pG1200L1800; //!
    EoverPHists_eopxAOD* m_plots_eop_pG1800L2200; //!
    EoverPHists_eopxAOD* m_plots_eop_pG2200L2800; //!
    EoverPHists_eopxAOD* m_plots_eop_pG2800L3600; //!
    EoverPHists_eopxAOD* m_plots_eop_pG2000; //!
    EoverPHists_eopxAOD* m_plots_eop_etaL06_pG2000; //!
    EoverPHists_eopxAOD* m_plots_eop_etaL20_pG2000_LarL1000; //!
    EoverPHists_eopxAOD* m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075; //!

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
  public:
    // Tree *myTree; //!
    // TH1 *myHist; //!

    // this is a standard constructor
    EoverPAnalysis_eopxAOD (std::string className = "EoverPAnalysis_eopxAOD");

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
    ClassDef(EoverPAnalysis_eopxAOD, 1);
    /// @endcond
};

#endif
