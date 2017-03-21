// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <EoverPAnalysis/EoverPAnalysis.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

// this is needed to distribute the algorithm to the workers
ClassImp(EoverPAnalysis)

  EoverPAnalysis :: EoverPAnalysis (std::string className) :
    Algorithm(className),
    // cutflows
    m_cutflowHist(nullptr),
    m_cutflowHistW(nullptr),
    m_trk_cutflowHist_1(nullptr),
    m_trk_cutflowHist_eop(nullptr),
    // pileup reweighting
    m_puwHist(nullptr),
    m_ptHist(nullptr),
    // number of tracks per event, after each selection
    m_trk_n_all(nullptr),
    m_trk_n_pass_extrapol(nullptr),
    m_trk_n_pass_trk1etaphi(nullptr),
    m_trk_n_pass_iso(nullptr),
    m_trk_n_pass_p(nullptr),
    m_trk_n_pass_pG500(nullptr),
    m_trk_n_pass_pG800(nullptr),
    m_trk_n_pass_pG1200(nullptr),
    m_trk_n_pass_pG2200(nullptr),
    m_trk_n_pass_pG3400(nullptr),
    m_trk_n_pass_pG5000(nullptr),
    m_trk_n_pass_eta(nullptr),
    m_trk_n_pass_etaL06(nullptr),
    m_trk_n_pass_etaG06L15(nullptr),
    m_trk_n_pass_etaG15L23(nullptr),
    m_trk_n_pass_larEmax(nullptr),
    m_trk_n_pass_tileEfrac(nullptr),
    // main histograms
    m_plots_eop(nullptr),
    // extra histograms for track isolation testing, turn on with m_doTrkIsoHists
    m_plots_eop_trks(nullptr),
    // extra plots, turn on with m_doGlobalTileEfracRanges
    m_plots_eop_TileEfrac000(nullptr),
    m_plots_eop_TileEfrac010(nullptr),
    m_plots_eop_TileEfrac030(nullptr),
    m_plots_eop_TileEfrac050(nullptr),
    m_plots_eop_TileEfrac060(nullptr),
    m_plots_eop_TileEfrac070(nullptr),
    m_plots_eop_TileEfrac075(nullptr),
    m_plots_eop_TileEfrac080(nullptr),
    // extra plots with different trk_p cuts, turn on with m_doGlobalEnergyRanges
    m_plots_eop_pL4000(nullptr),
    m_plots_eop_pG4000L8000(nullptr),
    m_plots_eop_pG8000L12000(nullptr),
    m_plots_eop_pG12000(nullptr),
    // extra plots with different trk_eta cuts, turn on with m_doGlobalEtaRanges
    m_plots_eop_etaL05(nullptr),
    m_plots_eop_etaG05L07(nullptr),
    m_plots_eop_etaG07(nullptr),
    // extra ranges for comparisons with Run 1 studies
    m_plots_eop_pG1200L1800(nullptr),
    m_plots_eop_pG1800L2200(nullptr),
    m_plots_eop_pG2200L2800(nullptr),
    m_plots_eop_pG2800L3400(nullptr),
    m_plots_eop_pG3400L4200(nullptr),
    m_plots_eop_pG4200L5000(nullptr),
    m_plots_eop_etaL06_pG2200L4200(nullptr),
    m_plots_eop_etaL06_pG4200L50000(nullptr),
    m_plots_eop_etaL06(nullptr),
    m_plots_eop_etaG06L11(nullptr),
    m_plots_eop_etaG11L14(nullptr),
    m_plots_eop_etaG14L15(nullptr),
    m_plots_eop_etaG15L18(nullptr),
    m_plots_eop_etaG18L19(nullptr),
    m_plots_eop_etaG19L23(nullptr)
{
  m_inTrackContainerName    = "";
  m_detailStr               = "";
  m_debug                   = false;
}

EL::StatusCode EoverPAnalysis :: setupJob (EL::Job& job)
{
  job.useXAOD();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("EoverPAnalysis").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: histInitialize ()
{

  Info("histInitialize()", "%s", m_name.c_str() );
  RETURN_CHECK("xAH::Algorithm::algInitialize()", xAH::Algorithm::algInitialize(), "");
  // needed here and not in initalize since this is called first
  if( m_inTrackContainerName.empty() || m_detailStr.empty() ){
    Error("configure()", "One or more required configuration values are empty");
    return EL::StatusCode::FAILURE;
  }

  // declare class and add histograms to output
  m_plots_eop = new EoverPHists(m_name, m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
  RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop -> initialize(), "");
  m_plots_eop -> record( wk() );

  if (m_doTrkIsoHists) {
    // extra histograms for track isolation testing
    m_plots_eop_trks = new EoverPHistsTrks(m_name, m_detailStr, m_trkIsoDRmax, m_trkIsoPfrac, m_doTrkPcut, m_trkPmin, m_trkPmax, m_doTrkEtacut, m_trkEtamin, m_trkEtamax, false); 
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_trks -> initialize(), "");
    m_plots_eop_trks -> record( wk() );
  }

  if (m_doGlobalTileEfracRanges) {

    m_plots_eop_TileEfrac000 = new EoverPHists(m_name+"_TileEfrac000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac010 = new EoverPHists(m_name+"_TileEfrac010", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac030 = new EoverPHists(m_name+"_TileEfrac030", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac050 = new EoverPHists(m_name+"_TileEfrac050", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac060 = new EoverPHists(m_name+"_TileEfrac060", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac070 = new EoverPHists(m_name+"_TileEfrac070", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac075 = new EoverPHists(m_name+"_TileEfrac075", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac080 = new EoverPHists(m_name+"_TileEfrac080", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac010 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac030 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac050 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac060 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac070 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac075 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_TileEfrac080 -> initialize(), "");

    m_plots_eop_TileEfrac000 -> record( wk() );
    m_plots_eop_TileEfrac010 -> record( wk() );
    m_plots_eop_TileEfrac030 -> record( wk() );
    m_plots_eop_TileEfrac050 -> record( wk() );
    m_plots_eop_TileEfrac060 -> record( wk() );
    m_plots_eop_TileEfrac070 -> record( wk() );
    m_plots_eop_TileEfrac075 -> record( wk() );
    m_plots_eop_TileEfrac080 -> record( wk() );
  }

  if (m_doGlobalEnergyRanges) {

    m_plots_eop_pL4000 = new EoverPHists(m_name+"_pL4000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG4000L8000 = new EoverPHists(m_name+"_pG4000L8000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG8000L12000 = new EoverPHists(m_name+"_pG8000L12000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG12000 = new EoverPHists(m_name+"_pG12000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pL4000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG4000L8000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG8000L12000 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG12000 -> initialize(), "");

    m_plots_eop_pL4000 -> record( wk() );
    m_plots_eop_pG4000L8000 -> record( wk() );
    m_plots_eop_pG8000L12000 -> record( wk() );
    m_plots_eop_pG12000 -> record( wk() );

  }

  if (m_doGlobalEtaRanges) {

    m_plots_eop_etaL05 = new EoverPHists(m_name+"_etaL05", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_etaG05L07 = new EoverPHists(m_name+"_etaG05L07", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_etaG07 = new EoverPHists(m_name+"_etaG07", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL05 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG05L07 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG07 -> initialize(), "");

    m_plots_eop_etaL05 -> record( wk() );
    m_plots_eop_etaG05L07 -> record( wk() );
    m_plots_eop_etaG07 -> record( wk() );

  }

  if (m_doGlobalExtraRanges) {
    m_plots_eop_pG1200L1800 = new EoverPHists(m_name+"_pG1200L1800", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_pG1800L2200 = new EoverPHists(m_name+"_pG1800L2200", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_pG2200L2800 = new EoverPHists(m_name+"_pG2200L2800", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_pG2800L3400 = new EoverPHists(m_name+"_pG2800L3400", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_pG3400L4200 = new EoverPHists(m_name+"_pG3400L4200", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_pG4200L5000 = new EoverPHists(m_name+"_pG4200L5000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaL06_pG2200L4200 = new EoverPHists(m_name+"_etaL06_pG2200L4200", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaL06_pG4200L50000 = new EoverPHists(m_name+"_etaL06_pG4200L50000", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaL06 = new EoverPHists(m_name+"_etaL06", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG06L11 = new EoverPHists(m_name+"_etaG06L11", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG11L14 = new EoverPHists(m_name+"_etaG11L14", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG14L15 = new EoverPHists(m_name+"_etaG14L15", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG15L18 = new EoverPHists(m_name+"_etaG15L18", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG18L19 = new EoverPHists(m_name+"_etaG18L19", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);
    m_plots_eop_etaG19L23 = new EoverPHists(m_name+"_etaG19L23", m_detailStr, m_energyCalib, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Pbins, m_doPbinsArray, m_PbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, false);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG1200L1800-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG1800L2200-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG2200L2800-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG2800L3400-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG3400L4200-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG4200L5000-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL06_pG2200L4200-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL06_pG4200L50000-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL06-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG06L11-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG11L14-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG14L15-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG15L18-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG18L19-> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG19L23-> initialize(), "");

    m_plots_eop_pG1200L1800 -> record( wk() );
    m_plots_eop_pG1800L2200 -> record( wk() );
    m_plots_eop_pG2200L2800 -> record( wk() );
    m_plots_eop_pG2800L3400 -> record( wk() );
    m_plots_eop_pG3400L4200 -> record( wk() );
    m_plots_eop_pG4200L5000 -> record( wk() );
    m_plots_eop_etaL06_pG2200L4200 -> record( wk() );
    m_plots_eop_etaL06_pG4200L50000 -> record( wk() );
    m_plots_eop_etaL06 -> record( wk() );
    m_plots_eop_etaG06L11 -> record( wk() );
    m_plots_eop_etaG11L14 -> record( wk() );
    m_plots_eop_etaG14L15 -> record( wk() );
    m_plots_eop_etaG15L18 -> record( wk() );
    m_plots_eop_etaG18L19 -> record( wk() );
    m_plots_eop_etaG19L23 -> record( wk() );

  }

  m_trk_cutflowHist_eop = new TH1D("cutflow_trks_eop", "cutflow_trks_1", 1, 1, 2);
  m_trk_cutflowHist_eop->SetCanExtend(TH1::kAllAxes);
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop all");
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop pass trk p cuts");
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop pass trk eta cuts");
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop pass trk iso");
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop pass larEmax");
  m_trk_cutflowHist_eop->GetXaxis()->FindBin("eop pass tileEfrac");

  wk()->addOutput(m_trk_cutflowHist_eop);

  int nBinsTrkN = 200; double minTrkN = -0.5; double maxTrkN = 199.5;
  m_trk_n_all = new TH1D((std::string(m_name+"/trk_n_all")).c_str(), "trk_n_all", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_extrapol = new TH1D((std::string(m_name+"/trk_n_pass_extrapol")).c_str(), "trk_n_pass_extrapol", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_trk1etaphi = new TH1D((std::string(m_name+"/trk_n_pass_trk1etaphi")).c_str(), "trk_n_pass_trk1etaphi", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_iso = new TH1D((std::string(m_name+"/trk_n_pass_iso")).c_str(), "trk_n_pass_iso", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_p = new TH1D((std::string(m_name+"/trk_n_pass_p")).c_str(), "trk_n_pass_p", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG500 = new TH1D((std::string(m_name+"/trk_n_pass_pG500")).c_str(), "trk_n_pass_pG500", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG800 = new TH1D((std::string(m_name+"/trk_n_pass_pG800")).c_str(), "trk_n_pass_pG800", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG1200 = new TH1D((std::string(m_name+"/trk_n_pass_pG1200")).c_str(), "trk_n_pass_pG1200", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG2200 = new TH1D((std::string(m_name+"/trk_n_pass_pG2200")).c_str(), "trk_n_pass_pG2200", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG3400 = new TH1D((std::string(m_name+"/trk_n_pass_pG3400")).c_str(), "trk_n_pass_pG3400", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_pG5000 = new TH1D((std::string(m_name+"/trk_n_pass_pG5000")).c_str(), "trk_n_pass_pG5000", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_eta = new TH1D((std::string(m_name+"/trk_n_pass_eta")).c_str(), "trk_n_pass_eta", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_etaL06 = new TH1D((std::string(m_name+"/trk_n_pass_etaL06")).c_str(), "trk_n_pass_etaL06", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_etaG06L15 = new TH1D((std::string(m_name+"/trk_n_pass_etaG06L15")).c_str(), "trk_n_pass_etaG06L15", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_etaG15L23 = new TH1D((std::string(m_name+"/trk_n_pass_etaG15L23")).c_str(), "trk_n_pass_etaG15L23", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_larEmax = new TH1D((std::string(m_name+"/trk_n_pass_larEmax")).c_str(), "trk_n_pass_larEmax", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_pass_tileEfrac = new TH1D((std::string(m_name+"/trk_n_pass_tileEfrac")).c_str(), "trk_n_pass_tileEfrac", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n_all->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_extrapol->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_trk1etaphi->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_iso->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_p->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG500->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG800->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG1200->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG2200->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG3400->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_pG5000->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_eta->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_etaL06->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_etaG06L15->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_etaG15L23->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_larEmax->GetXaxis()->SetTitle("N_trks");
  m_trk_n_pass_tileEfrac->GetXaxis()->SetTitle("N_trks");

  wk()->addOutput(m_trk_n_all);
  wk()->addOutput(m_trk_n_pass_extrapol);
  wk()->addOutput(m_trk_n_pass_trk1etaphi);
  wk()->addOutput(m_trk_n_pass_iso);
  wk()->addOutput(m_trk_n_pass_p);
  wk()->addOutput(m_trk_n_pass_pG500);
  wk()->addOutput(m_trk_n_pass_pG800);
  wk()->addOutput(m_trk_n_pass_pG1200);
  wk()->addOutput(m_trk_n_pass_pG2200);
  wk()->addOutput(m_trk_n_pass_pG3400);
  wk()->addOutput(m_trk_n_pass_pG5000);
  wk()->addOutput(m_trk_n_pass_eta);
  wk()->addOutput(m_trk_n_pass_etaL06);
  wk()->addOutput(m_trk_n_pass_etaG06L15);
  wk()->addOutput(m_trk_n_pass_etaG15L23);
  wk()->addOutput(m_trk_n_pass_larEmax);
  wk()->addOutput(m_trk_n_pass_tileEfrac);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: fileExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode EoverPAnalysis :: initialize ()
{
  Info("initialize()", "EoverPAnalysis");
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  if(m_useCutFlow) {
    TFile *file = wk()->getOutputFile ("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");
    m_cutflow_bin  = m_cutflowHist->GetXaxis()->FindBin("eop pass all");
    m_cutflowHistW->GetXaxis()->FindBin("eop pass all");

    // retrieve the object cutflow
    m_trk_cutflowHist_1  = (TH1D*)file->Get("cutflow_trks_1");
    m_trk_cutflow_eop_all_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop all");
    m_trk_cutflow_eop_extrapol_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk extrapol");
    m_trk_cutflow_eop_trk1etaphi_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk1 eta,phi<1000");
    m_trk_cutflow_eop_pass_iso_bin  = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk iso");
    m_trk_cutflow_eop_pass_p_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk p cuts");
    m_trk_cutflow_eop_pass_eta_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk eta cuts");
    m_trk_cutflow_eop_pass_larEmax_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass larEmax");
    m_trk_cutflow_eop_pass_tileEfrac_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass tileEfrac");

  }

  // pileup reweighting
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");
  if (m_doCustomPUreweighting && eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )) {
    // TFile *f_prw = wk()->getOutputFile ("pileup");
    // m_puwHist = (TH1D*)f_prw->Get("pileup_weights");
    TFile *f_pileup = new TFile(std::string("$ROOTCOREBIN/data/EoverPAnalysis/"+m_pileupReweightingFile).c_str(), "READ");
    m_puwHist = (TH1D*)f_pileup->Get("h_pileupweight");
  }
  if (m_doTrkPtReweighting && eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )) {
    TFile *f_pt = new TFile(std::string("$ROOTCOREBIN/data/EoverPAnalysis/"+m_trkPtReweightingFile).c_str(), "READ");
    m_ptHist = (TH1D*)f_pt->Get("h_ptweight");
  }

  m_numEvent = 0;
  m_numEventPass = 0;
  m_weightNumEventPass = 0;

  m_trk_cutflow_eop_all = 0;
  m_trk_cutflow_eop_extrapol = 0;
  m_trk_cutflow_eop_trk1etaphi = 0;
  m_trk_cutflow_eop_pass_iso = 0;
  m_trk_cutflow_eop_pass_p = 0;
  m_trk_cutflow_eop_pass_eta = 0;
  m_trk_cutflow_eop_pass_larEmax = 0;
  m_trk_cutflow_eop_pass_tileEfrac = 0;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: execute ()
{
  // retrieve event
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");

  //  1.) the PU weight ("PileupWeight")
  //  2.) the corrected mu ("corrected_averageInteractionsPerCrossing")
  double eventWeight(1.);
  int mu_avg(1e8); // initialize with a that won't pass the selection
  if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) )
    mu_avg = eventInfo->averageInteractionsPerCrossing();
  if( eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
    eventWeight = eventInfo->auxdecor< float >( "mcEventWeight" );
    // std::cout << "eventWeight, before PRW: " << eventWeight << std::endl;
    // if( eventInfo->isAvailable< float >( "corrected_averageInteractionsPerCrossing" ) ) 
    //   mu_avg = eventInfo->auxdata< float >( "corrected_averageInteractionsPerCrossing" );
    double pileupWeight(0.);
    // std::cout << "before getting PileupWeight" << std::endl;
    if (m_doCustomPUreweighting) {
      if (mu_avg > 0. && mu_avg <= m_puwHist->GetNbinsX()) {
        pileupWeight = m_puwHist->GetBinContent(mu_avg+1);
      }
      eventWeight *= pileupWeight;
    }
  }

  m_numEvent++;

  if (mu_avg < m_mu_avg_min || mu_avg > m_mu_avg_max) return EL::StatusCode::SUCCESS;

  m_trk_n_all_tmp = 0;
  m_trk_n_pass_extrapol_tmp = 0;
  m_trk_n_pass_trk1etaphi_tmp = 0;
  m_trk_n_pass_iso_tmp = 0;
  m_trk_n_pass_p_tmp = 0;
  m_trk_n_pass_pG500_tmp = 0;
  m_trk_n_pass_pG800_tmp = 0;
  m_trk_n_pass_pG1200_tmp = 0;
  m_trk_n_pass_pG2200_tmp = 0;
  m_trk_n_pass_pG3400_tmp = 0;
  m_trk_n_pass_pG5000_tmp = 0;
  m_trk_n_pass_eta_tmp = 0;
  m_trk_n_pass_etaL06_tmp = 0;
  m_trk_n_pass_etaG06L15_tmp = 0;
  m_trk_n_pass_etaG15L23_tmp = 0;
  m_trk_n_pass_larEmax_tmp = 0;
  m_trk_n_pass_tileEfrac_tmp = 0;

  const xAOD::VertexContainer *vtxs(nullptr);
  RETURN_CHECK("TrackVertexSelection::execute()", HelperFunctions::retrieve(vtxs, "PrimaryVertices", m_event, m_store, m_verbose) ,"");

  const xAOD::TrackParticleContainer* trks(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(trks, m_inTrackContainerName, m_event, m_store, m_verbose) ,"");

  // extra histograms for track isolation testing
  if (m_doTrkIsoHists) {
    RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_trks->execute(trks, vtxs, eventInfo, eventWeight), "");
  }

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();
  // loop over all tracks only once
  for( ; trk_itr != trk_end; ++trk_itr ) {

    m_trk_cutflow_eop_all++;
    m_trk_n_all_tmp++;

    const xAOD::TrackParticle* trk = (*trk_itr);

    if (trk->auxdata<int>("CALO_extrapolation") == 0) continue;

    m_trk_cutflow_eop_extrapol++;
    m_trk_n_pass_extrapol_tmp++;

    double trk_pt = trk->pt()/1e3;
    double trk_p = 0;
    if (fabs(trk->qOverP())>0.) trk_p = (1./fabs(trk->qOverP()))/1e3; 
    // coordinates of the track in the ID
    double trk_etaID = trk->eta();
    double trk_phiID = trk->phi();
    // coordinates of the track extrapolated to the calorimeter
    // EMB2
    double trk_etaEMB2 = trk->auxdata<float>("CALO_trkEta_EMB2");
    double trk_phiEMB2 = trk->auxdata<float>("CALO_trkPhi_EMB2");
    //EME2
    double trk_etaEME2 = trk->auxdata<float>("CALO_trkEta_EME2");
    double trk_phiEME2 = trk->auxdata<float>("CALO_trkPhi_EME2");

    // check that the track is extrapolated to either EMB2 or EME2
    // (if not then trk_eta = trk_phi = -999999999)
    if ( (fabs(trk_etaEMB2) > 1000.0 || fabs(trk_phiEMB2) > 1000.0) && 
         (fabs(trk_etaEME2) > 1000.0 || fabs(trk_phiEME2) > 1000.0) )
      continue;

    m_trk_cutflow_eop_trk1etaphi++;
    m_trk_n_pass_trk1etaphi_tmp++;

    bool trk_not_isolated_EMB2 = false;
    bool trk_not_isolated_EME2 = false;

    // track isolation: p(cone of DR='m_trkIsoDRmax')/p(track) < 'm_trkIsoPfrac' 
    // float surr_trk_sum_p = 0.;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      if (trk_itr != trk2_itr) { // do not double count the selected track 
        const xAOD::TrackParticle* trk2 = (*trk2_itr);

        //EMB2
        if (fabs(trk_etaEMB2) < (double)1000.0 && fabs(trk_phiEMB2) < (double)1000.0) {
          double trk2_etaEMB2 = trk2->auxdata<float>("CALO_trkEta_EMB2");
          double trk2_phiEMB2 = trk2->auxdata<float>("CALO_trkPhi_EMB2");
          if (fabs(trk2_etaEMB2) < (double)1000.0 && fabs(trk2_phiEMB2) < (double)1000.0) {
            double trk_trk2_dR_EMB2 = deltaR(trk_etaEMB2, trk_phiEMB2, trk2_etaEMB2, trk2_phiEMB2);
            if (trk_trk2_dR_EMB2 <= m_trkIsoDRmax) trk_not_isolated_EMB2 = true;
          }
        }

        //EME2
        if (fabs(trk_etaEME2) < (double)1000.0 && fabs(trk_phiEME2) < (double)1000.0) {
          double trk2_etaEME2 = trk2->auxdata<float>("CALO_trkEta_EME2");
          double trk2_phiEME2 = trk2->auxdata<float>("CALO_trkPhi_EME2");

          if (fabs(trk2_etaEME2) < (double)1000.0 && fabs(trk2_phiEME2) < (double)1000.0) {
            double trk_trk2_dR_EME2 = deltaR(trk_etaEME2, trk_phiEME2, trk2_etaEME2, trk2_phiEME2);
            if (trk_trk2_dR_EME2 <= m_trkIsoDRmax) trk_not_isolated_EME2 = true;
          }
        }

        // //track isolation - check if trk2 falls within DRmax of trk
        // if (trk_trk2_dR_min < m_trkIsoDRmax) {  
        //   // calculate the leading and avg p of the surrounding tracks, 
        //   // used for TileCal comparisons with Run1
        //   if (fabs(trk2->qOverP())>0.) surr_trk_sum_p += (1./fabs(trk2->qOverP()))/1e3; 
        // }
      }
    } // END looping trk2

    // check track isolation requirement
    // if (fabs(surr_trk_sum_p/trk_p) > m_trkIsoPfrac) continue;
    if (trk_not_isolated_EMB2) continue;
    if (trk_not_isolated_EME2) continue;

    m_trk_cutflow_eop_pass_iso++;
    m_trk_n_pass_iso_tmp++;

    if (trk_p > 0.5)
      m_trk_n_pass_pG500_tmp++;
    if (trk_p > 0.8)
      m_trk_n_pass_pG800_tmp++;
    if (trk_p > 1.2)
      m_trk_n_pass_pG1200_tmp++;
    if (trk_p > 2.2)
      m_trk_n_pass_pG2200_tmp++;
    if (trk_p > 3.4)
      m_trk_n_pass_pG3400_tmp++;
    if (trk_p > 5.0)
      m_trk_n_pass_pG5000_tmp++;

    if (fabs(trk_etaID) < 0.6)
      m_trk_n_pass_etaL06_tmp++;
    if (fabs(trk_etaID) >= 0.6 && fabs(trk_etaID) < 1.5)
      m_trk_n_pass_etaG06L15_tmp++;
    if (fabs(trk_etaID) >= 1.5 && fabs(trk_etaID) < 2.3)
      m_trk_n_pass_etaG15L23_tmp++;

    // check track p requirement
    // if (m_doTrkPcut) {
    //   if (trk_p < m_trkPmin) continue;
    //   if (trk_p >= m_trkPmax) continue;
    // }
    m_trk_cutflow_eop_pass_p++;
    m_trk_n_pass_p_tmp++;

    // check track eta requirement
    if (m_doTrkEtacut) {
      // if (fabs(trk_etaID) < m_trkEtamin) continue;
      // if (fabs(trk_etaID) >= m_trkEtamax) continue;
      if ((double)fabs(trk_etaID) >= (double)2.3) continue;
    }
    m_trk_cutflow_eop_pass_eta++;
    m_trk_n_pass_eta_tmp++;

    double trkWeight = eventWeight;
    if (m_doTrkPtReweighting && eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
      if (trk_pt > 0. && trk_pt < 30.) {
        trkWeight *= m_ptHist->GetBinContent(m_ptHist->FindBin(trk_pt));
      }
    }

    if (m_doTileCuts) {

      // check LAr energy loss requirement
      double trk_sumE_Lar_200 = 0.; 
      for (unsigned int i=0; i<m_layer_lar.size(); i++) {
        double trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_200"))/1e3; 
        if (trk_E_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Lar_200 += trk_E_tmp;
      }
      if (trk_sumE_Lar_200 > m_LarEmax) continue;

      m_trk_cutflow_eop_pass_larEmax++;
      m_trk_n_pass_larEmax_tmp++;

      // check E(tile)/E(total) requirement
      double trk_sumE_Tile_200 = 0.; 
      for (unsigned int i=0; i<m_layer_tile.size(); i++) {
        double trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_200"))/1e3; 
        if (trk_E_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Tile_200 += trk_E_tmp;
      }
      double trk_sumE_Total_200 = 0.;
      for (unsigned int i=0; i<m_layer.size(); i++) { 
        double trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_200"))/1e3; 
        if (trk_E_200_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Total_200 += trk_E_200_tmp; 
      }
      double trk_TileEfrac_200 = 0.;
      if (trk_sumE_Total_200 > 0.)  
        trk_TileEfrac_200 = trk_sumE_Tile_200/trk_sumE_Total_200;

      // Make all the histograms for different TileEfrac selections
      if (m_doGlobalTileEfracRanges && trk_TileEfrac_200 >= 0. && trk_TileEfrac_200 < 1.) {
        if (m_doGlobalTileEfracRanges && trk_TileEfrac_200 == 0.)
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac000->execute(trk, vtxs, eventInfo, trkWeight), "");
        if (trk_TileEfrac_200 >= .1) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac010->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .3) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac030->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .5) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac050->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .6) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac060->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .7) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac070->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .75) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac075->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
        if (trk_TileEfrac_200 >= .8) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac080->execute(trk, vtxs, eventInfo, trkWeight), "");
        }
      }

      // Default TileEfrac cut
      if (trk_TileEfrac_200 < 0. || trk_TileEfrac_200 >= 1.0 || trk_TileEfrac_200 < m_TileEfracmin) continue;
      m_trk_cutflow_eop_pass_tileEfrac++;
      m_trk_n_pass_tileEfrac_tmp++;
    }

    // fill eop histograms
    RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop -> execute(trk, vtxs, eventInfo, trkWeight), "");

    // fill eop histograms for different trk p ranges
    if (m_doGlobalEnergyRanges) {
      if (trk_p < 4)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pL4000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 4 && trk_p < 8)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG4000L8000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 8 && trk_p < 12)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG8000L12000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 12)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG12000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
    }
    // fill eop histograms for different trk eta ranges
    if (m_doGlobalEtaRanges) {
      if (fabs(trk_etaID) < .5) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaL05 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= .5 && fabs(trk_etaID) < .7) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG05L07 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= .7) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG07 -> execute(trk, vtxs, eventInfo, trkWeight), "");
    }

    if (m_doGlobalExtraRanges) {
      if (trk_p >= 1.2 && trk_p < 1.8)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG1200L1800 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 1.8 && trk_p < 2.2)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG1800L2200 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 2.2 && trk_p < 2.8)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG2200L2800 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 2.8 && trk_p < 3.4)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG2800L3400 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 3.4 && trk_p < 4.2)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG3400L4200 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (trk_p >= 4.2 && trk_p < 5.)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG4200L5000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) < .6){
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaL06 -> execute(trk, vtxs, eventInfo, trkWeight), "");
        if (trk_p >= 2.2 && trk_p < 4.6)
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaL06_pG2200L4200 -> execute(trk, vtxs, eventInfo, trkWeight), "");
        if (trk_p >= 4.6 && trk_p < 50.)
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaL06_pG4200L50000 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      }
      if (fabs(trk_etaID) >= .6 && fabs(trk_etaID) < 1.1) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG06L11 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= 1.1 && fabs(trk_etaID) < 1.4) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG11L14 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= 1.4 && fabs(trk_etaID) < 1.5) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG14L15 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= 1.5 && fabs(trk_etaID) < 1.8) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG15L18 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= 1.8 && fabs(trk_etaID) < 1.9) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG18L19 -> execute(trk, vtxs, eventInfo, trkWeight), "");
      if (fabs(trk_etaID) >= 1.9 && fabs(trk_etaID) < 2.3) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG19L23 -> execute(trk, vtxs, eventInfo, trkWeight), "");
    }

  } // END looping trk

  m_numEventPass++;
  m_weightNumEventPass += eventWeight;

  m_trk_n_all->Fill(m_trk_n_all_tmp, eventWeight);
  m_trk_n_pass_extrapol->Fill(m_trk_n_pass_extrapol_tmp, eventWeight);
  m_trk_n_pass_trk1etaphi->Fill(m_trk_n_pass_trk1etaphi_tmp, eventWeight);
  m_trk_n_pass_iso->Fill(m_trk_n_pass_iso_tmp, eventWeight);
  m_trk_n_pass_p->Fill(m_trk_n_pass_p_tmp, eventWeight);
  m_trk_n_pass_pG500->Fill(m_trk_n_pass_pG500_tmp, eventWeight);
  m_trk_n_pass_pG800->Fill(m_trk_n_pass_pG800_tmp, eventWeight);
  m_trk_n_pass_pG1200->Fill(m_trk_n_pass_pG1200_tmp, eventWeight);
  m_trk_n_pass_pG2200->Fill(m_trk_n_pass_pG2200_tmp, eventWeight);
  m_trk_n_pass_pG3400->Fill(m_trk_n_pass_pG3400_tmp, eventWeight);
  m_trk_n_pass_pG5000->Fill(m_trk_n_pass_pG5000_tmp, eventWeight);
  m_trk_n_pass_eta->Fill(m_trk_n_pass_eta_tmp, eventWeight);
  m_trk_n_pass_etaL06->Fill(m_trk_n_pass_etaL06_tmp, eventWeight);
  m_trk_n_pass_etaG06L15->Fill(m_trk_n_pass_etaG06L15_tmp, eventWeight);
  m_trk_n_pass_etaG15L23->Fill(m_trk_n_pass_etaG15L23_tmp, eventWeight);
  m_trk_n_pass_larEmax->Fill(m_trk_n_pass_larEmax_tmp, eventWeight);
  m_trk_n_pass_tileEfrac->Fill(m_trk_n_pass_tileEfrac_tmp, eventWeight);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: postExecute () { return EL::StatusCode::SUCCESS; }

EL::StatusCode EoverPAnalysis :: finalize () { 

  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );

    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_all_bin, m_trk_cutflow_eop_all );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_extrapol_bin, m_trk_cutflow_eop_extrapol );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_trk1etaphi_bin, m_trk_cutflow_eop_trk1etaphi );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_iso_bin, m_trk_cutflow_eop_pass_iso );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_p_bin, m_trk_cutflow_eop_pass_p );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_eta_bin, m_trk_cutflow_eop_pass_eta );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_larEmax_bin, m_trk_cutflow_eop_pass_larEmax );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_tileEfrac_bin, m_trk_cutflow_eop_pass_tileEfrac );

  }

  return EL::StatusCode::SUCCESS; 
}

EL::StatusCode EoverPAnalysis :: histFinalize ()
{
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_all_bin, m_trk_cutflow_eop_all );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_extrapol_bin, m_trk_cutflow_eop_extrapol );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_trk1etaphi_bin, m_trk_cutflow_eop_trk1etaphi );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_iso_bin, m_trk_cutflow_eop_pass_iso );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_p_bin, m_trk_cutflow_eop_pass_p );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_eta_bin, m_trk_cutflow_eop_pass_eta );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_larEmax_bin, m_trk_cutflow_eop_pass_larEmax );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_tileEfrac_bin, m_trk_cutflow_eop_pass_tileEfrac );

  // clean up memory
  if(m_plots_eop) delete m_plots_eop;
  if(m_plots_eop_trks) delete m_plots_eop_trks;

  if(m_plots_eop_TileEfrac000) delete m_plots_eop_TileEfrac000;
  if(m_plots_eop_TileEfrac010) delete m_plots_eop_TileEfrac010;
  if(m_plots_eop_TileEfrac030) delete m_plots_eop_TileEfrac030;
  if(m_plots_eop_TileEfrac050) delete m_plots_eop_TileEfrac050;
  if(m_plots_eop_TileEfrac060) delete m_plots_eop_TileEfrac060;
  if(m_plots_eop_TileEfrac070) delete m_plots_eop_TileEfrac070;
  if(m_plots_eop_TileEfrac075) delete m_plots_eop_TileEfrac075;
  if(m_plots_eop_TileEfrac080) delete m_plots_eop_TileEfrac080;

  if(m_plots_eop_pL4000) delete m_plots_eop_pL4000;
  if(m_plots_eop_pG4000L8000) delete m_plots_eop_pG4000L8000;
  if(m_plots_eop_pG8000L12000) delete m_plots_eop_pG8000L12000;
  if(m_plots_eop_pG12000) delete m_plots_eop_pG12000;

  if(m_plots_eop_etaL05) delete m_plots_eop_etaL05;
  if(m_plots_eop_etaG05L07) delete m_plots_eop_etaG05L07;
  if(m_plots_eop_etaG07) delete m_plots_eop_etaG07;

  if(m_plots_eop_pG1200L1800) delete m_plots_eop_pG1200L1800;
  if(m_plots_eop_pG1800L2200) delete m_plots_eop_pG1800L2200;
  if(m_plots_eop_pG2200L2800) delete m_plots_eop_pG2200L2800;
  if(m_plots_eop_pG2800L3400) delete m_plots_eop_pG2800L3400;
  if(m_plots_eop_pG3400L4200) delete m_plots_eop_pG3400L4200;
  if(m_plots_eop_pG4200L5000) delete m_plots_eop_pG4200L5000;
  if(m_plots_eop_etaL06) delete m_plots_eop_etaL06;
  if(m_plots_eop_etaL06_pG2200L4200) delete m_plots_eop_etaL06_pG2200L4200;
  if(m_plots_eop_etaL06_pG4200L50000) delete m_plots_eop_etaL06_pG4200L50000;
  if(m_plots_eop_etaG06L11) delete m_plots_eop_etaG06L11;
  if(m_plots_eop_etaG11L14) delete m_plots_eop_etaG11L14;
  if(m_plots_eop_etaG14L15) delete m_plots_eop_etaG14L15;
  if(m_plots_eop_etaG15L18) delete m_plots_eop_etaG15L18;
  if(m_plots_eop_etaG18L19) delete m_plots_eop_etaG18L19;
  if(m_plots_eop_etaG19L23) delete m_plots_eop_etaG19L23;

  RETURN_CHECK("xAH::Algorithm::algFinalize()", xAH::Algorithm::algFinalize(), "");
  return EL::StatusCode::SUCCESS;
}

double EoverPAnalysis :: deltaR (double trk_eta, double trk_phi, double trk2_eta, double trk2_phi)
{
  double trk_trk2_dEta = fabs(trk2_eta - trk_eta);
  double trk_trk2_dPhi = fabs(trk2_phi - trk_phi);

  if (trk_trk2_dPhi > TMath::Pi())
    trk_trk2_dPhi = 2*TMath::Pi() - trk_trk2_dPhi;

  return sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );
}
