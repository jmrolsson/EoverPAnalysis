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
  m_inClusterContainerName  = "";
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
  if( m_inTrackContainerName.empty() || m_inClusterContainerName.empty() || m_detailStr.empty() ){
    Error("configure()", "One or more required configuration values are empty");
    return EL::StatusCode::FAILURE;
  }


  // declare class and add histograms to output
  m_plots_eop = new EoverPHists(m_name, m_detailStr);
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

  const xAOD::CaloClusterContainer* ccls(nullptr);
  RETURN_CHECK("EoverPAnalysis::execute()", HelperFunctions::retrieve(ccls, m_inClusterContainerName, m_event, m_store, m_verbose) ,"");

  // make some diagnostic plots to test the track-cluster matching
  RETURN_CHECK("EoverPAnalysis::execute()", m_plots_eop->execute(trks, ccls, vtxs, eventInfo, eventWeight), "");

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  for( ; trk_itr != trk_end; ++trk_itr ) {
    RETURN_CHECK("EoverPAnalysis::execute()", this->trackClusterMatching( (*trk_itr), ccls, eventWeight ), "");
  }

  return EL::StatusCode::SUCCESS;
}

StatusCode EoverPAnalysis :: trackClusterMatching(const xAOD::TrackParticle* trk, const xAOD::CaloClusterContainer* ccls, float eventWeight)
{
  // FIXME Function for track cluster matching, at the moment this is implemented in EoverPHists
  // xAOD::CaloClusterContainer::const_iterator ccl_itr = ccls->begin();
  // xAOD::CaloClusterContainer::const_iterator ccl_end = ccls->end();
  // for( ; ccl_itr != ccl_end; ++ccl_itr ) {
  //   const xAOD::CaloCluster* ccl = (*ccl_itr);
  //   float trk_ccl_dR = trk->p4().DeltaR( ccl->p4() );
  //   // std::cout << trk_ccl_dR << std::endl;
  // }

  return StatusCode::SUCCESS;
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
