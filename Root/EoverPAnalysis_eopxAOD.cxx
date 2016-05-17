#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"

#include <EoverP/EoverPAnalysis_eopxAOD.h>

#include <xAODAnaHelpers/tools/ReturnCheck.h>

// this is needed to distribute the algorithm to the workers
ClassImp(EoverPAnalysis_eopxAOD)

EoverPAnalysis_eopxAOD :: EoverPAnalysis_eopxAOD (std::string className) :
    Algorithm(className),
    m_plots_eop(nullptr),
    m_plots_eop_etaL06(nullptr),
    m_plots_eop_etaG06L11(nullptr),
    m_plots_eop_etaG18L19(nullptr),
    m_plots_eop_etaG19L23(nullptr),
    m_plots_eop_pG1200L1800(nullptr),
    m_plots_eop_pG1800L2200(nullptr),
    m_plots_eop_pG2200L2800(nullptr),
    m_plots_eop_pG2800L3600(nullptr),
    m_plots_eop_pG2000(nullptr),
    m_plots_eop_etaL06_pG2000(nullptr),
    m_plots_eop_etaL20_pG2000_LarL1000(nullptr),
    m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075(nullptr)
{
  m_inTrackContainerName    = "";
  m_detailStr               = "";
  m_debug                   = false;
}

EL::StatusCode EoverPAnalysis_eopxAOD :: setupJob (EL::Job& job)
{
  job.useXAOD();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("EoverPAnalysis_eopxAOD").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis_eopxAOD :: histInitialize ()
{

  Info("histInitialize()", "%s", m_name.c_str() );
  RETURN_CHECK("xAH::Algorithm::algInitialize()", xAH::Algorithm::algInitialize(), "");
  // needed here and not in initalize since this is called first
  if( m_inTrackContainerName.empty() || m_detailStr.empty() ){
    Error("configure()", "One or more required configuration values are empty");
    return EL::StatusCode::FAILURE;
  }

  // declare class and add histograms to output
  m_plots_eop = new EoverPHists_eopxAOD(m_name, m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, m_doTrkPcut, m_trkPmin, m_trkPmax, m_doTrkEtacut, m_trkEtamin, m_trkEtamax);
  if (m_doEtaPranges) {
    m_plots_eop_etaL06 = new EoverPHists_eopxAOD(m_name+"_etaL06", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, false, 0.0, 1e8, true, 0.0, 0.6);
    m_plots_eop_etaG06L11 = new EoverPHists_eopxAOD(m_name+"_etaG06L11", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, false, 0.0, 1e8, true, 0.6, 1.1);
    m_plots_eop_etaG18L19 = new EoverPHists_eopxAOD(m_name+"_etaG18L19", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, false, 0.0, 1e8, true, 1.8, 1.9);
    m_plots_eop_etaG19L23 = new EoverPHists_eopxAOD(m_name+"_etaG19L23", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, false, 0.0, 1e8, true, 1.9, 2.3);
    m_plots_eop_pG1200L1800 = new EoverPHists_eopxAOD(m_name+"_pG1200L1800", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 1.2, 1.8, false, 0.0, 1e8);
    m_plots_eop_pG1800L2200 = new EoverPHists_eopxAOD(m_name+"_pG1800L2200", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 1.8, 2.2, false, 0.0, 1e8);
    m_plots_eop_pG2200L2800 = new EoverPHists_eopxAOD(m_name+"_pG2200L2800", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 2.2, 2.8, false, 0.0, 1e8);
    m_plots_eop_pG2800L3600 = new EoverPHists_eopxAOD(m_name+"_pG2800L3600", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 2.8, 3.6, false, 0.0, 1e8);
    m_plots_eop_pG2000 = new EoverPHists_eopxAOD(m_name+"_pG2000", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 2.0, 1e8, false, 0.0, 1e8);
    m_plots_eop_etaL06_pG2000 = new EoverPHists_eopxAOD(m_name+"_etaL06_pG2000", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, m_LarEmax, m_TileEfrac, true, 2.0, 1e8, true, 0.0, 0.6);
  }
  if (m_doTileCuts) {
    m_plots_eop_etaL20_pG2000_LarL1000 = new EoverPHists_eopxAOD(m_name+"_etaL20_pG2000_LarL1000", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, 1.0, -1, true, 2.0, 1e8, true, 0.0, 2.0);
    m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075 = new EoverPHists_eopxAOD(m_name+"_etaL20_pG2000_LarL1000_TileEfrac075", m_detailStr, m_trkExtrapol, m_doBgSubtr, m_doEMcalib, m_doLCWcalib, m_doCells, m_doCaloEM, m_doCaloHAD, m_doCaloTotal, m_trkIsoDRmax, m_trkIsoPfrac, 1.0, 0.75, true, 2.0, 1e8, true, 0.0, 2.0);
  }

  RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop -> initialize(), "");
  m_plots_eop -> record( wk() );

  if (m_doEtaPranges) {
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL06 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG06L11 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG18L19 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG19L23 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG1200L1800 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG1800L2200 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG2200L2800 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG2800L3600 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG2000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL06_pG2000 -> initialize(), "");
    m_plots_eop_etaL06 -> record( wk() );
    m_plots_eop_etaG06L11 -> record( wk() );
    m_plots_eop_etaG18L19 -> record( wk() );
    m_plots_eop_etaG19L23 -> record( wk() );
    m_plots_eop_pG1200L1800 -> record( wk() );
    m_plots_eop_pG1800L2200 -> record( wk() );
    m_plots_eop_pG2200L2800 -> record( wk() );
    m_plots_eop_pG2800L3600 -> record( wk() );
    m_plots_eop_pG2000 -> record( wk() );
    m_plots_eop_etaL06_pG2000 -> record( wk() );
  }
  if (m_doTileCuts) {
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL20_pG2000_LarL1000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075 -> initialize(), "");
    m_plots_eop_etaL20_pG2000_LarL1000 -> record( wk() );
    m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075 -> record( wk() );
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis_eopxAOD :: fileExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis_eopxAOD :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode EoverPAnalysis_eopxAOD :: initialize ()
{
  Info("initialize()", "EoverPAnalysis_eopxAOD");
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis_eopxAOD :: execute ()
{
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");

  float eventWeight(1);
  if( eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
    eventWeight = eventInfo->auxdecor< float >( "mcEventWeight" );
  }

  const xAOD::VertexContainer *vtxs(nullptr);
  RETURN_CHECK("TrackVertexSelection::execute()", HelperFunctions::retrieve(vtxs, "PrimaryVertices", m_event, m_store, m_verbose) ,"");

  const xAOD::TrackParticleContainer* trks(nullptr);
  RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", HelperFunctions::retrieve(trks, m_inTrackContainerName, m_event, m_store, m_verbose) ,"");

  // make eop plots for the eopxAOD samples
  RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop->execute(trks, vtxs, eventInfo, eventWeight), "");
  if (m_doEtaPranges) {
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaL06->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaG06L11->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaG18L19->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaG19L23->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_pG1200L1800->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_pG1800L2200->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_pG2200L2800->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_pG2800L3600->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_pG2000->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaL06_pG2000->execute(trks, vtxs, eventInfo, eventWeight), "");
  }
  if (m_doTileCuts) {
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaL20_pG2000_LarL1000->execute(trks, vtxs, eventInfo, eventWeight), "");
    RETURN_CHECK("EoverPAnalysis_eopxAOD::execute()", m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075->execute(trks, vtxs, eventInfo, eventWeight), "");
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis_eopxAOD :: postExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis_eopxAOD :: finalize () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis_eopxAOD :: histFinalize ()
{
  // clean up memory
  if(m_plots_eop) delete m_plots_eop;
  if(m_plots_eop_etaL06) delete m_plots_eop_etaL06;
  if(m_plots_eop_etaG06L11) delete m_plots_eop_etaG06L11;
  if(m_plots_eop_etaG18L19) delete m_plots_eop_etaG18L19;
  if(m_plots_eop_etaG19L23) delete m_plots_eop_etaG19L23;
  if(m_plots_eop_pG1200L1800) delete m_plots_eop_pG1200L1800;
  if(m_plots_eop_pG1800L2200) delete m_plots_eop_pG1800L2200;
  if(m_plots_eop_pG2200L2800) delete m_plots_eop_pG2200L2800;
  if(m_plots_eop_pG2800L3600) delete m_plots_eop_pG2800L3600;
  if(m_plots_eop_pG2000) delete m_plots_eop_pG2000;
  if(m_plots_eop_etaL06_pG2000) delete m_plots_eop_etaL06_pG2000;
  if(m_plots_eop_etaL20_pG2000_LarL1000) delete m_plots_eop_etaL20_pG2000_LarL1000;
  if(m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075) delete m_plots_eop_etaL20_pG2000_LarL1000_TileEfracG075;

  RETURN_CHECK("xAH::Algorithm::algFinalize()", xAH::Algorithm::algFinalize(), "");
  return EL::StatusCode::SUCCESS;
}
