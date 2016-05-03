#ifndef EoverP_EoverPAnalysis_H
#define EoverP_EoverPAnalysis_H

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"

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
  std::string m_inClusterContainerName;

  // configuration variables
  std::string m_detailStr;

private:
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

  StatusCode trackClusterMatching (const xAOD::TrackParticle* trk, const xAOD::CaloClusterContainer* ccls, float eventWeight);
  
  /// @cond
  // this is needed to distribute the algorithm to the workers
  ClassDef(EoverPAnalysis, 1);
  /// @endcond

};

#endif
