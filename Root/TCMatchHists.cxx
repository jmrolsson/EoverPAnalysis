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
  // number of bins and ranges for histograms
  int nBinsDR = 100;       float minDR = 0;                float maxDR = 1;
  int nBinsEta = 80;       float minEta = -4.0;            float maxEta = 4.0;
  int nBinsPhi = 120;      float minPhi = -TMath::Pi();    float maxPhi = TMath::Pi();

  int nBinsTrkPt = 200;    float minTrkPt = 0;             float maxTrkPt = 20;
  int nBinsTrkN = 10;      float minTrkN = -0.5;           float maxTrkN = 9.5;
  int nBinsCclN = 100;     float minCclN = -0.5;           float maxCclN = 99.5;

  int nBinsEoP = 220;      float minEoP = -10;             float maxEoP = 100;

  //// 1D histograms

  //track-track plots
  m_trk_ntrks_maxDR01 = book(m_name, "trk_ntrks_maxDR01", "number of tracks within dR=0.1 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02 = book(m_name, "trk_ntrks_maxDR02", "number of tracks within dR=0.2 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03 = book(m_name, "trk_ntrks_maxDR03", "number of tracks within dR=0.3 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04 = book(m_name, "trk_ntrks_maxDR04", "number of tracks within dR=0.4 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05 = book(m_name, "trk_ntrks_maxDR05", "number of tracks within dR=0.5 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06 = book(m_name, "trk_ntrks_maxDR06", "number of tracks within dR=0.6 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07 = book(m_name, "trk_ntrks_maxDR07", "number of tracks within dR=0.7 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08 = book(m_name, "trk_ntrks_maxDR08", "number of tracks within dR=0.8 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09 = book(m_name, "trk_ntrks_maxDR09", "number of tracks within dR=0.9 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10 = book(m_name, "trk_ntrks_maxDR10", "number of tracks within dR=1.0 from a selected track", nBinsTrkN, minTrkN, maxTrkN); 
 
  // track-cluster plots
  m_trk_neccl_dR = book(m_name, "trk_neccl_dR", "dR between a track and the nearest cluster", nBinsDR, minDR, maxDR);

  m_trk_nccls_maxDR01 = book(m_name, "trk_nccls_maxDR01", "number of clusters within dR=0.1 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR02 = book(m_name, "trk_nccls_maxDR02", "number of clusters within dR=0.2 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR03 = book(m_name, "trk_nccls_maxDR03", "number of clusters within dR=0.3 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR04 = book(m_name, "trk_nccls_maxDR04", "number of clusters within dR=0.4 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR05 = book(m_name, "trk_nccls_maxDR05", "number of clusters within dR=0.5 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR06 = book(m_name, "trk_nccls_maxDR06", "number of clusters within dR=0.6 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR07 = book(m_name, "trk_nccls_maxDR07", "number of clusters within dR=0.7 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR08 = book(m_name, "trk_nccls_maxDR08", "number of clusters within dR=0.8 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR09 = book(m_name, "trk_nccls_maxDR09", "number of clusters within dR=0.9 from a selected track", nBinsCclN, minCclN, maxCclN); 
  m_trk_nccls_maxDR10 = book(m_name, "trk_nccls_maxDR10", "number of clusters within dR=1.0 from a selected track", nBinsCclN, minCclN, maxCclN); 

  // eoverp plots
  m_eop_neccl = book(m_name, "eop_neccls", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone = book(m_name, "eop_sum_ccls_dRcone", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);

  //// 2D histograms

  // plots of selected track and nearest cluster
  m_trk_neccl_dR_vs_neccl_e = book(m_name, "trk_neccl_dR_vs_neccl_e", "cluster e", 100, -5, 15, "dR between a track and the nearest cluster", nBinsDR, minDR, maxDR);
  m_trk_neccl_trk_eta_vs_neccl_eta = book(m_name, "trk_neccl_trk_eta_vs_neccl_eta", "cluster #eta", nBinsEta, minEta, maxEta, "track #eta", nBinsEta, minEta, maxEta);
  m_trk_neccl_trk_phi_vs_neccl_phi = book(m_name, "trk_neccl_trk_phi_vs_neccl_phi", "cluster #phi", nBinsPhi, minPhi, maxPhi, "track #phi", nBinsPhi, minPhi, maxPhi);

  // eoverp plots
  m_eop_neccl_vs_trk_pt = book(m_name, "eop_neccl_vs_trk_pt", "track p_{T}", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_vs_trk_eta = book(m_name, "eop_neccl_vs_trk_eta", "track p_{T}", nBinsEta, minEta, maxEta, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_vs_trk_phi = book(m_name, "eop_neccl_vs_trk_phi", "track p_{T}", nBinsPhi, minPhi, maxPhi, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_pt = book(m_name, "eop_sum_ccls_dRcone_vs_trk_pt", "track p_{T}", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_eta = book(m_name, "eop_sum_ccls_dRcone_vs_trk_eta", "track p_{T}", nBinsEta, minEta, maxEta, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_phi = book(m_name, "eop_sum_ccls_dRcone_vs_trk_phi", "track p_{T}", nBinsPhi, minPhi, maxPhi, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode TCMatchHists::execute( const xAOD::TrackParticleContainer* trks, const xAOD::CaloClusterContainer* ccls, float eventWeight )
{
  float trk_pt = 0;
  float trk_eta = 0; 
  float trk_phi = 0; 
  float trk_trk2_dR = 0; 
  float trk_ccl_dR = 0; 
  float trk_neccl_minDR = 1e8;
  float neccl_e = 0; 
  float neccl_eta = 0; 
  float neccl_phi = 0; 
  float ccl_e_sum_dRcone = 0;
  bool isolatedTrack = true;

  // loop over tracks and clusters around a selected track, plot how many tracks (clusters) fall within a certain dR
  float maxDR_ranges[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();
  xAOD::CaloClusterContainer::const_iterator ccl_itr = ccls->begin();
  xAOD::CaloClusterContainer::const_iterator ccl_end = ccls->end();

  // track-track and track-cluster (outer loop over all tracks)
  for( ; trk_itr != trk_end; ++trk_itr ) {
    const xAOD::TrackParticle* trk = (*trk_itr);
    trk_pt  = trk->pt()/1e3;
    trk_eta = trk->eta();
    trk_phi = trk->phi();

    // keep track of the number of tracks within a certain dR from the selected track
    int trk_ntrks_maxDR[10] = {0};
    isolatedTrack = true;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      const xAOD::TrackParticle* trk2 = (*trk2_itr);
      trk_trk2_dR = trk2->p4().DeltaR( trk->p4() );
      for (int i = 0; i < 10; ++i) {
        // std::cout << "maxDR_ranges[i]: " << maxDR_ranges[i] << std::endl;  
        // std::cout << "trk_trk2_dR: " << trk_trk2_dR << std::endl;  
        if (trk_trk2_dR > 1e-3 && trk_trk2_dR < maxDR_ranges[i]) { // set a threshold not to count the selected track
          trk_ntrks_maxDR[i]++;
        }
        // track isolation
        if (trk_trk2_dR > 1e-3 && trk_trk2_dR < m_cut_trk_trk2_maxDR) { // no other track allowed within dR < 0.4
          isolatedTrack = false;
        }
      }
    }

    // clusters 
    float trk_nccls_maxDR[10] = {0};
    ccl_e_sum_dRcone = 0;
    for( ; ccl_itr != ccl_end; ++ccl_itr ) {
      const xAOD::CaloCluster* ccl = (*ccl_itr);
      trk_ccl_dR = trk->p4().DeltaR( ccl->p4() );
      // figure out the dR of the nearest cluster
      if (trk_ccl_dR < trk_neccl_minDR) {
        neccl_e = ccl->e()/1e3;
        neccl_eta = ccl->eta();
        neccl_phi = ccl->phi();
        trk_neccl_minDR = trk_ccl_dR;   
      } 
      // keep track of the number of cluster within a certain dR from the selected track
      for (int i = 0; i < 10; ++i) {
        if (trk_ccl_dR < maxDR_ranges[i]) { 
          trk_nccls_maxDR[i]++;
        }
      }
      // sum the energy of all clusters in a cone around the selected track
      if (trk_ccl_dR < m_cut_trk_ccl_maxDR) {
        ccl_e_sum_dRcone += ccl->e()/1e3; 
      }
    }

    // track-cluster matching histograms
    m_trk_neccl_dR -> Fill(trk_neccl_minDR, eventWeight); 
    m_trk_neccl_dR_vs_neccl_e -> Fill(neccl_e, trk_neccl_minDR, eventWeight); 
    m_trk_neccl_trk_eta_vs_neccl_eta -> Fill(neccl_eta, trk_eta, eventWeight);
    m_trk_neccl_trk_phi_vs_neccl_phi -> Fill(neccl_phi, trk_phi, eventWeight);

    // histograms with information about closeby tracks
    m_trk_ntrks_maxDR01 -> Fill(trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02 -> Fill(trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03 -> Fill(trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04 -> Fill(trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05 -> Fill(trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06 -> Fill(trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07 -> Fill(trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08 -> Fill(trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09 -> Fill(trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10 -> Fill(trk_ntrks_maxDR[9], eventWeight);

    // histograms with information about closeby clusters 
    m_trk_nccls_maxDR01 -> Fill(trk_nccls_maxDR[0], eventWeight);
    m_trk_nccls_maxDR02 -> Fill(trk_nccls_maxDR[1], eventWeight);
    m_trk_nccls_maxDR03 -> Fill(trk_nccls_maxDR[2], eventWeight);
    m_trk_nccls_maxDR04 -> Fill(trk_nccls_maxDR[3], eventWeight);
    m_trk_nccls_maxDR05 -> Fill(trk_nccls_maxDR[4], eventWeight);
    m_trk_nccls_maxDR06 -> Fill(trk_nccls_maxDR[5], eventWeight);
    m_trk_nccls_maxDR07 -> Fill(trk_nccls_maxDR[6], eventWeight);
    m_trk_nccls_maxDR08 -> Fill(trk_nccls_maxDR[7], eventWeight);
    m_trk_nccls_maxDR09 -> Fill(trk_nccls_maxDR[8], eventWeight);
    m_trk_nccls_maxDR10 -> Fill(trk_nccls_maxDR[9], eventWeight);

    // eoverp histograms
    if (isolatedTrack) {
      float eop_neccl = neccl_e/trk_pt;
      float eop_sum_ccls_dRcone = ccl_e_sum_dRcone/trk_pt;
      // std::cout << "eop_neccl: " << eop_neccl << std::endl;
      // std::cout << "eop_sum_ccls_dRcone: " << eop_sum_ccls_dRcone << std::endl;

      m_eop_neccl                      -> Fill(eop_neccl, eventWeight); 
      m_eop_sum_ccls_dRcone            -> Fill(eop_sum_ccls_dRcone , eventWeight); 

      m_eop_neccl_vs_trk_pt            -> Fill(trk_pt,  eop_neccl, eventWeight); 
      m_eop_neccl_vs_trk_eta           -> Fill(trk_eta, eop_neccl, eventWeight); 
      m_eop_neccl_vs_trk_phi           -> Fill(trk_phi, eop_neccl, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_pt  -> Fill(trk_pt,  eop_sum_ccls_dRcone, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_eta -> Fill(trk_eta, eop_sum_ccls_dRcone, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_phi -> Fill(trk_phi, eop_sum_ccls_dRcone, eventWeight); 

      // if (abs(trk_eta) < 0.6) {
      // }
      // if (abs(trk_eta) > 0.6 && abs(trk_eta) < 1.2) {
      // }
    }

  } // END loop tracks 

  return StatusCode::SUCCESS;
}
