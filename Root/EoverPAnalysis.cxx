// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <EoverP/EoverPAnalysis.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

// this is needed to distribute the algorithm to the workers
ClassImp(EoverPAnalysis)

  EoverPAnalysis :: EoverPAnalysis (std::string className) :
    Algorithm(className),
    m_cutflowHist(nullptr),
    m_cutflowHistW(nullptr),
    m_trk_cutflowHist_1(nullptr),
    m_trk_cutflowHist_eop(nullptr),
    m_plots_eop(nullptr),
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
    m_plots_eop_pL4(nullptr),
    m_plots_eop_pG4L8(nullptr),
    m_plots_eop_pG8L12(nullptr),
    m_plots_eop_pG12(nullptr),
    // extra plots with different trk_eta cuts, turn on with m_doGlobalEtaRanges
    m_plots_eop_etaL05(nullptr),
    m_plots_eop_etaG05L07(nullptr),
    m_plots_eop_etaG07(nullptr)
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
  m_plots_eop = new EoverPHists(m_name, m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

  RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop -> initialize(), "");

  m_plots_eop -> record( wk() );

  if (m_doGlobalTileEfracRanges) {

    m_plots_eop_TileEfrac000 = new EoverPHists(m_name+"_TileEfrac000", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac010 = new EoverPHists(m_name+"_TileEfrac010", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac030 = new EoverPHists(m_name+"_TileEfrac030", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac050 = new EoverPHists(m_name+"_TileEfrac050", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac060 = new EoverPHists(m_name+"_TileEfrac060", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac070 = new EoverPHists(m_name+"_TileEfrac070", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac075 = new EoverPHists(m_name+"_TileEfrac075", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_TileEfrac080 = new EoverPHists(m_name+"_TileEfrac080", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

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

    m_plots_eop_pL4 = new EoverPHists(m_name+"_pL4", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG4L8 = new EoverPHists(m_name+"_pG4L8", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG8L12 = new EoverPHists(m_name+"_pG8L12", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_pG12 = new EoverPHists(m_name+"_pG12", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pL4 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG4L8 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG8L12 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_pG12 -> initialize(), "");

    m_plots_eop_pL4 -> record( wk() );
    m_plots_eop_pG4L8 -> record( wk() );
    m_plots_eop_pG8L12 -> record( wk() );
    m_plots_eop_pG12 -> record( wk() );

  }

  if (m_doGlobalEtaRanges) {

    m_plots_eop_etaL05 = new EoverPHists(m_name+"_etaL05", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_etaG05L07 = new EoverPHists(m_name+"_etaG05L07", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);
    m_plots_eop_etaG07 = new EoverPHists(m_name+"_etaG07", m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doExtraEtaEnergyBinHists);

    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaL05 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG05L07 -> initialize(), "");
    RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop_etaG07 -> initialize(), "");

    m_plots_eop_etaL05 -> record( wk() );
    m_plots_eop_etaG05L07 -> record( wk() );
    m_plots_eop_etaG07 -> record( wk() );

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
    m_trk_cutflow_eop_pass_p_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk p cuts");
    m_trk_cutflow_eop_pass_eta_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk eta cuts");
    m_trk_cutflow_eop_pass_iso_bin  = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass trk iso");
    m_trk_cutflow_eop_pass_larEmax_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass larEmax");
    m_trk_cutflow_eop_pass_tileEfrac_bin = m_trk_cutflowHist_1->GetXaxis()->FindBin("eop pass tileEfrac");

  }

  m_numEvent = 0;
  m_numEventPass = 0;
  m_weightNumEventPass = 0;

  m_trk_cutflow_eop_all = 0;
  m_trk_cutflow_eop_pass_p = 0;
  m_trk_cutflow_eop_pass_eta = 0;
  m_trk_cutflow_eop_pass_iso = 0;
  m_trk_cutflow_eop_pass_larEmax = 0;
  m_trk_cutflow_eop_pass_tileEfrac = 0;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: execute ()
{
  // retrieve event
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");

  float eventWeight(1);
  if( eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
    eventWeight = eventInfo->auxdecor< float >( "mcEventWeight" );
  }

  m_numEvent++;

  const xAOD::VertexContainer *vtxs(nullptr);
  RETURN_CHECK("TrackVertexSelection::execute()", HelperFunctions::retrieve(vtxs, "PrimaryVertices", m_event, m_store, m_verbose) ,"");

  const xAOD::TrackParticleContainer* trks(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(trks, m_inTrackContainerName, m_event, m_store, m_verbose) ,"");

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();
  // loop over all tracks only once
  for( ; trk_itr != trk_end; ++trk_itr ) {

    m_trk_cutflow_eop_all++;

    const xAOD::TrackParticle* trk = (*trk_itr);
    float trk_p = 0;
    if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
    float trk_etaCALO = trk->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
    float trk_phiCALO = trk->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));

    // check track p requirement
    if (m_doTrkPcut) {
      if (trk_p < m_trkPmin) continue;
      if (trk_p > m_trkPmax) continue;
    }
    m_trk_cutflow_eop_pass_p++;

    // check track eta requirement
    if (m_doTrkEtacut) {
      if (TMath::Abs(trk_etaCALO) < m_trkEtamin) continue;
      if (TMath::Abs(trk_etaCALO) > m_trkEtamax) continue;
    }
    m_trk_cutflow_eop_pass_eta++;

    // track isolation: p(cone of DR='m_trkIsoDRmax')/p(track) < 'm_trkIsoPfrac' 
    float surr_trk_sum_p = 0.;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      if (trk_itr != trk2_itr) { // do not double count the selected track 
        const xAOD::TrackParticle* trk2 = (*trk2_itr);
        float trk2_etaCALO = trk2->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
        float trk2_phiCALO = trk2->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));
        float trk_trk2_dEta= TMath::Abs(trk2_etaCALO - trk_etaCALO);
        float trk_trk2_dPhi = TMath::Abs(trk2_phiCALO - trk_phiCALO);
        if (trk_trk2_dPhi > TMath::Pi())
          trk_trk2_dPhi = 2*TMath::Pi() - trk_trk2_dPhi;
        float trk_trk2_dR = sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );

        // track isolation
        if (trk_trk2_dR < m_trkIsoDRmax) { // check if trk2 falls within DRmax of trk 
          // calculate the leading and avg p of the surrounding tracks
          if (TMath::Abs(trk2->qOverP())>0.) surr_trk_sum_p += (1./TMath::Abs(trk2->qOverP()))/1e3; 
        }
      }
    } // END looping trk2

    // check track isolation requirement
    if (TMath::Abs(surr_trk_sum_p/trk_p) > m_trkIsoPfrac) continue;

    m_trk_cutflow_eop_pass_iso++;

    if (m_doTileCuts) {

      // check LAr energy loss requirement
      float trk_sumE_Lar_200 = 0.; 
      for (unsigned int i=0; i<m_layer_lar.size(); i++) {
        float trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_200"))/1e3; 
        if (trk_E_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Lar_200 += trk_E_tmp;
      }
      if (trk_sumE_Lar_200 > m_LarEmax) continue;

      m_trk_cutflow_eop_pass_larEmax++;

      // check E(tile)/E(total) requirement
      float trk_sumE_Tile_200 = 0.; 
      for (unsigned int i=0; i<m_layer_tile.size(); i++) {
        float trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_200"))/1e3; 
        if (trk_E_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Tile_200 += trk_E_tmp;
      }
      float trk_sumE_Total_200 = 0.;
      for (unsigned int i=0; i<m_layer.size(); i++) { 
        float trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_200"))/1e3; 
        if (trk_E_200_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Total_200 += trk_E_200_tmp; 
      }
      float trk_TileEfrac_200 = 0.;
      if (trk_sumE_Total_200 > 0.)  
        trk_TileEfrac_200 = trk_sumE_Tile_200/trk_sumE_Total_200;
      
      // Make all the histograms for different TileEfrac selections
      if (m_doGlobalTileEfracRanges && trk_TileEfrac_200 >= 0. && trk_TileEfrac_200 < 1.) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac000->execute(trk, vtxs, eventInfo, eventWeight), "");
        if (trk_TileEfrac_200 >= .1) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac010->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .3) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac030->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .5) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac050->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .6) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac060->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .7) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac070->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .75) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac075->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
        if (trk_TileEfrac_200 >= .8) {
          RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_TileEfrac080->execute(trk, vtxs, eventInfo, eventWeight), "");
        }
      }

      // Default TileEfrac cut
      if (trk_TileEfrac_200 < 0. || trk_TileEfrac_200 >= 1.0 || trk_TileEfrac_200 < m_TileEfracmin) continue;
      m_trk_cutflow_eop_pass_tileEfrac++;
    }

    // fill eop histograms
    RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop -> execute(trk, vtxs, eventInfo, eventWeight), "");

    // fill eop histograms for different trk p ranges
    if (m_doGlobalEnergyRanges) {
      if (trk_p < 4)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pL4 -> execute(trk, vtxs, eventInfo, eventWeight), "");
      if (trk_p >= 4 && trk_p < 8)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG4L8 -> execute(trk, vtxs, eventInfo, eventWeight), "");
      if (trk_p >= 8 && trk_p < 12)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG8L12 -> execute(trk, vtxs, eventInfo, eventWeight), "");
      if (trk_p >= 12)
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_pG12 -> execute(trk, vtxs, eventInfo, eventWeight), "");
    }
    // fill eop histograms for different trk eta ranges
    if (m_doGlobalEtaRanges) {
      if (TMath::Abs(trk_etaCALO) < .5) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaL05 -> execute(trk, vtxs, eventInfo, eventWeight), "");
      if (TMath::Abs(trk_etaCALO) >= .5 && TMath::Abs(trk_etaCALO) < .7) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG05L07 -> execute(trk, vtxs, eventInfo, eventWeight), "");
      if (TMath::Abs(trk_etaCALO) >= .7) 
        RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop_etaG07 -> execute(trk, vtxs, eventInfo, eventWeight), "");
    }

  } // END looping trk

  m_numEventPass++;
  m_weightNumEventPass += eventWeight;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: postExecute () { return EL::StatusCode::SUCCESS; }

EL::StatusCode EoverPAnalysis :: finalize () { 

  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );

    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_all_bin, m_trk_cutflow_eop_all );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_p_bin, m_trk_cutflow_eop_pass_p );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_eta_bin, m_trk_cutflow_eop_pass_eta );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_iso_bin, m_trk_cutflow_eop_pass_iso );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_larEmax_bin, m_trk_cutflow_eop_pass_larEmax );
    m_trk_cutflowHist_1->SetBinContent( m_trk_cutflow_eop_pass_tileEfrac_bin, m_trk_cutflow_eop_pass_tileEfrac );

  }
  
  return EL::StatusCode::SUCCESS; 
}

EL::StatusCode EoverPAnalysis :: histFinalize ()
{
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_all_bin, m_trk_cutflow_eop_all );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_p_bin, m_trk_cutflow_eop_pass_p );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_eta_bin, m_trk_cutflow_eop_pass_eta );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_iso_bin, m_trk_cutflow_eop_pass_iso );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_larEmax_bin, m_trk_cutflow_eop_pass_larEmax );
  m_trk_cutflowHist_eop->SetBinContent( m_trk_cutflow_eop_pass_tileEfrac_bin, m_trk_cutflow_eop_pass_tileEfrac );

  // clean up memory
  if(m_plots_eop) delete m_plots_eop;
  
  if(m_plots_eop_TileEfrac000) delete m_plots_eop_TileEfrac000;
  if(m_plots_eop_TileEfrac010) delete m_plots_eop_TileEfrac010;
  if(m_plots_eop_TileEfrac030) delete m_plots_eop_TileEfrac030;
  if(m_plots_eop_TileEfrac050) delete m_plots_eop_TileEfrac050;
  if(m_plots_eop_TileEfrac060) delete m_plots_eop_TileEfrac060;
  if(m_plots_eop_TileEfrac070) delete m_plots_eop_TileEfrac070;
  if(m_plots_eop_TileEfrac075) delete m_plots_eop_TileEfrac075;
  if(m_plots_eop_TileEfrac080) delete m_plots_eop_TileEfrac080;

  if(m_plots_eop_pL4) delete m_plots_eop_pL4;
  if(m_plots_eop_pG4L8) delete m_plots_eop_pG4L8;
  if(m_plots_eop_pG8L12) delete m_plots_eop_pG8L12;
  if(m_plots_eop_pG12) delete m_plots_eop_pG12;

  if(m_plots_eop_etaL05) delete m_plots_eop_etaL05;
  if(m_plots_eop_etaG05L07) delete m_plots_eop_etaG05L07;
  if(m_plots_eop_etaG07) delete m_plots_eop_etaG07;

  RETURN_CHECK("xAH::Algorithm::algFinalize()", xAH::Algorithm::algFinalize(), "");
  return EL::StatusCode::SUCCESS;
}
