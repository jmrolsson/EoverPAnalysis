#include "EoverP/TCMatchHists.h"

#include <math.h>

#include "xAODAnaHelpers/tools/ReturnCheck.h"

TCMatchHists :: TCMatchHists (std::string name, std::string detailStr) :
  HistogramManager(name, detailStr)
{
}

TCMatchHists :: ~TCMatchHists () {}

StatusCode TCMatchHists::initialize()
{

  // These plots are always made
  int nBinsDR = 100;       float minDR = 0;                float maxDR = 2;
  int nBinsEta = 80;       float minEta = -4.0;            float maxEta = 4.0;
  int nBinsPhi = 120;      float minPhi = -TMath::Pi();    float maxPhi = TMath::Pi();
  int nBinsTrkPt = 200;    float minTrkPt = 0;             float maxTrkPt = 20;
  int nBinsEoP = 220;      float minEoP = -10;             float maxEoP = 100;

  // 1D
  m_trk_neccl_dR = book(m_name, "trk_neccl_dR", "dR between a track and the nearest cluster", nBinsDR, minDR, maxDR);
  m_trk_allccl_dRmax = book(m_name, "trk_allccl_dRmax ", "number of clusters within dR from a track", 9, 0.1, 1);
  m_eop_all = book(m_name, "eop_all", "E/p", nBinsEoP, minEoP, maxEoP);

  // 2D
  m_trk_neccl_dR_vs_neccl_e = book(m_name, "trk_neccl_dR_vs_neccl_e", "cluster e", 100, -5, 15, "dR between a track and the nearest cluster", nBinsDR, minDR, maxDR);
  m_trk_neccl_trk_eta_vs_neccl_eta = book(m_name, "trk_neccl_trk_eta_vs_neccl_eta", "cluster #eta", nBinsEta, minEta, maxEta, "track #eta", nBinsEta, minEta, maxEta);
  m_trk_neccl_trk_phi_vs_neccl_phi = book(m_name, "trk_neccl_trk_phi_vs_neccl_phi", "cluster #phi", nBinsPhi, minPhi, maxPhi, "track #phi", nBinsPhi, minPhi, maxPhi);
  m_eop_vs_trk_pt = book(m_name, "eop_vs_trk_pt", "track p_{T}", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p", nBinsEoP, minEoP, maxEoP);
  m_eop_vs_trk_eta = book(m_name, "eop_vs_trk_eta", "track p_{T}", nBinsEta, minEta, maxEta, "E/p", nBinsEoP, minEoP, maxEoP);
  m_eop_vs_trk_phi = book(m_name, "eop_vs_trk_phi", "track p_{T}", nBinsPhi, minPhi, maxPhi, "E/p", nBinsEoP, minEoP, maxEoP);

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode TCMatchHists::execute( const xAOD::TrackParticleContainer* trks, const xAOD::CaloClusterContainer* ccls, float eventWeight )
{
  float minDR = 1e8;
  float trk_ccl_dR = 0; 
  float trk_pt = 0;
  float trk_eta = 0; 
  float trk_phi = 0; 
  float neccl_e = 0; 
  float neccl_eta = 0; 
  float neccl_phi = 0; 
  float ccl_e_sum = 0; 

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::CaloClusterContainer::const_iterator ccl_itr = ccls->begin();
  xAOD::CaloClusterContainer::const_iterator ccl_end = ccls->end();

  for( ; trk_itr != trk_end; ++trk_itr ) {
    const xAOD::TrackParticle* trk = (*trk_itr);
    trk_pt  = trk->pt()/1e3;
    trk_eta = trk->eta();
    trk_phi = trk->phi();

    for( ; ccl_itr != ccl_end; ++ccl_itr ) {
      const xAOD::CaloCluster* ccl = (*ccl_itr);
      trk_ccl_dR = trk->p4().DeltaR( ccl->p4() );
      if (trk_ccl_dR < minDR) {
        neccl_e = ccl->e()/1e3;
        neccl_eta = ccl->eta();
        neccl_phi = ccl->phi();
        minDR = trk_ccl_dR;   
      } 
      // plot number of clusters within a certain dR from a track
      // FIXME
      if (trk_ccl_dR < m_cut_trk_ccl_dRmax) {
        ccl_e_sum += ccl->e()/1e3; 
      }
    }

    // track-cluster matching histograms
    m_trk_neccl_dR -> Fill(minDR, eventWeight); 
    m_trk_neccl_dR_vs_neccl_e -> Fill(neccl_e, minDR, eventWeight); 
    m_trk_neccl_trk_eta_vs_neccl_eta -> Fill(neccl_eta, trk_eta, eventWeight);
    m_trk_neccl_trk_phi_vs_neccl_phi -> Fill(neccl_phi, trk_phi, eventWeight);

    // E/p histograms
    float eop = ccl_e_sum/trk_pt;
    m_eop_all -> Fill(eop, eventWeight); 
    m_eop_vs_trk_pt  -> Fill(trk_pt,  eop, eventWeight); 
    m_eop_vs_trk_eta -> Fill(trk_eta, eop, eventWeight); 
    m_eop_vs_trk_phi -> Fill(trk_phi, eop, eventWeight); 

    if (abs(trk_eta) < 0.6) {
    }
    if (abs(trk_eta) > 0.6 && abs(trk_eta) < 1.2) {
    }

  }

  return StatusCode::SUCCESS;
}
