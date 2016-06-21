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
    m_plots_eop(nullptr)
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
  m_plots_eop = new EoverPHists(m_name, m_detailStr, m_energyCalib, m_trkExtrapol, m_doCaloTotal, m_doCaloEM, m_doCaloHAD, m_doBgSubtr, m_doTileLayer, m_Ebins, m_doEbinsArray, m_EbinsArray, m_Etabins, m_doEtabinsArray, m_EtabinsArray, m_doEtaEnergyRanges);

  RETURN_CHECK("TrackHistsAlgo::histInitialize()", m_plots_eop -> initialize(), "");

  m_plots_eop -> record( wk() );

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: fileExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode EoverPAnalysis :: initialize ()
{
  Info("initialize()", "EoverPAnalysis");
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: execute ()
{
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");

  float eventWeight(1);
  if( eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
    eventWeight = eventInfo->auxdecor< float >( "mcEventWeight" );
  }

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

    // check track eta requirement
    if (m_doTrkEtacut) {
      if (TMath::Abs(trk_etaCALO) < m_trkEtamin) continue;
      if (TMath::Abs(trk_etaCALO) > m_trkEtamax) continue;
    }

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
    // std::cout << "----> pass trk isolation" << std::endl;

    if (m_doTileCuts) {

      // check LAr energy loss requirement
      float trk_sumE_Lar_200 = 0; 
      for (unsigned int i=0; i<m_layer_lar.size(); i++) {
        float trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_200"))/1e3; 
        if (trk_E_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
          trk_sumE_Lar_200 += trk_E_tmp;
      }
      if (trk_sumE_Lar_200 > m_LarEmax) continue;

      // check E(tile)/E(total) requirement
      float trk_sumE_Tile_200 = 0; 
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
      float trk_TileEfrac_200 = trk_sumE_Tile_200/trk_sumE_Total_200; 
      if (trk_TileEfrac_200 < m_TileEfracmin) continue;

    }

    // fill eop histograms
    RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop->execute(trk, vtxs, eventInfo, eventWeight), "");

  } // END looping trk

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode EoverPAnalysis :: postExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis :: finalize () { return EL::StatusCode::SUCCESS; }
EL::StatusCode EoverPAnalysis :: histFinalize ()
{
  // clean up memory
  if(m_plots_eop) delete m_plots_eop;

  RETURN_CHECK("xAH::Algorithm::algFinalize()", xAH::Algorithm::algFinalize(), "");
  return EL::StatusCode::SUCCESS;
}
