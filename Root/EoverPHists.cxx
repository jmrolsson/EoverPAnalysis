#include "EoverP/EoverPHists.h"
#include <math.h>
#include "xAODAnaHelpers/tools/ReturnCheck.h"

EoverPHists :: EoverPHists (std::string name, std::string detailStr) :
  HistogramManager(name, detailStr)
{
}

EoverPHists :: ~EoverPHists () {}

StatusCode EoverPHists::initialize()
{
  // number of bins and ranges for histograms
  int nBinsMu = 50;        float minMu = -0.5;            float maxMu = 49.5;
  int nBinsNPV = 50;       float minNPV = -0.5;           float maxNPV = 49.5;
  int nBinsDR = 100;       float minDR = 0;               float maxDR = 1;
  int nBinsEta = 80;       float minEta = -4.0;           float maxEta = 4.0;
  int nBinsPhi = 120;      float minPhi = -TMath::Pi();   float maxPhi = TMath::Pi();

  int nBinsTrkPt = 300;    float minTrkPt = 0;            float maxTrkPt = 30;
  int nBinsTrkN = 10;      float minTrkN = -0.5;          float maxTrkN = 9.5;
  int nBinsCclN = 100;     float minCclN = -0.5;          float maxCclN = 99.5;

  int nBinsEoP = 250;      float minEoP = -5;             float maxEoP = 20;

  //// 1D histograms
 
  // event level plots
  m_mu = book(m_name, "mu", "#mu", nBinsMu, minMu, maxMu); 
  m_avg_mu = book(m_name, "avg_mu", "<#mu>", nBinsMu, minMu, maxMu); 
  m_npv = book(m_name, "npv", "npv", nBinsNPV, minNPV, maxNPV); 

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
  m_eop_maxE_ccl_dRcone = book(m_name, "eop_maxE_ccl_dRcone", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_pG2  = book(m_name, "eop_neccls_pG2 ", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_pG2  = book(m_name, "eop_sum_ccls_dRcone_pG2 ", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_pG2  = book(m_name, "eop_maxE_ccl_dRcone_pG2 ", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  
  m_eop_neccl_etaL08 = book(m_name, "eop_neccls_etaL08", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_etaL08 = book(m_name, "eop_sum_ccls_dRcone_etaL08", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_etaL08 = book(m_name, "eop_maxE_ccl_dRcone_etaL08", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_etaL08_pG2 = book(m_name, "eop_neccls_etaL08_pG2 ", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_etaL08_pG2  = book(m_name, "eop_sum_ccls_dRcone_etaL08_pG2 ", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_etaL08_pG2  = book(m_name, "eop_maxE_ccl_dRcone_etaL08_pG2 ", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);

  m_eop_neccl_etaG08L17 = book(m_name, "eop_neccls_etaG08L17", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_etaG08L17 = book(m_name, "eop_sum_ccls_dRcone_etaG08L17", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_etaG08L17 = book(m_name, "eop_maxE_ccl_dRcone_etaG08L17", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_etaG08L17_pG2   = book(m_name, "eop_neccls_etaG08L17_pG2  ", "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_etaG08L17_pG2   = book(m_name, "eop_sum_ccls_dRcone_etaG08L17_pG2  ", "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_etaG08L17_pG2   = book(m_name, "eop_maxE_ccl_dRcone_etaG08L17_pG2  ", "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);

  //// 2D histograms

  // plots of selected track and nearest cluster
  m_trk_neccl_dR_vs_neccl_e = book(m_name, "trk_neccl_dR_vs_neccl_e", "cluster e", 100, -5, 15, "dR between a track and the nearest cluster", nBinsDR, minDR, maxDR);
  m_trk_neccl_trk_eta_vs_neccl_eta = book(m_name, "trk_neccl_trk_eta_vs_neccl_eta", "cluster #eta", nBinsEta, minEta, maxEta, "track #eta", nBinsEta, minEta, maxEta);
  m_trk_neccl_trk_phi_vs_neccl_phi = book(m_name, "trk_neccl_trk_phi_vs_neccl_phi", "cluster #phi", nBinsPhi, minPhi, maxPhi, "track #phi", nBinsPhi, minPhi, maxPhi);

  // eoverp plots
  m_eop_neccl_vs_avg_mu = book(m_name, "eop_neccl_vs_avg_mu", "<#mu>", nBinsMu, minMu, maxMu, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_vs_trk_p = book(m_name, "eop_neccl_vs_trk_p", "track |p|", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_vs_trk_eta = book(m_name, "eop_neccl_vs_trk_eta", "track #eta", nBinsEta, minEta, maxEta, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_neccl_vs_trk_phi = book(m_name, "eop_neccl_vs_trk_phi", "track #phi", nBinsPhi, minPhi, maxPhi, "E/p (nearest cluster, maxDR=0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_avg_mu = book(m_name, "eop_sum_ccls_dRcone_vs_avg_mu", "<#mu>", nBinsMu, minMu, maxMu, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_p = book(m_name, "eop_sum_ccls_dRcone_vs_trk_p", "track |p|", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_eta = book(m_name, "eop_sum_ccls_dRcone_vs_trk_eta", "track #eta", nBinsEta, minEta, maxEta, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_sum_ccls_dRcone_vs_trk_phi = book(m_name, "eop_sum_ccls_dRcone_vs_trk_phi", "track #phi", nBinsPhi, minPhi, maxPhi, "E/p (sum of clusters with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_vs_avg_mu = book(m_name, "eop_maxE_ccl_dRcone_vs_avg_mu", "<#mu>", nBinsMu, minMu, maxMu, "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_vs_trk_p = book(m_name, "eop_maxE_ccl_dRcone_vs_trk_p", "track |p|", nBinsTrkPt, minTrkPt, maxTrkPt, "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_vs_trk_eta = book(m_name, "eop_maxE_ccl_dRcone_vs_trk_eta", "track #eta", nBinsEta, minEta, maxEta, "E/p (highest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);
  m_eop_maxE_ccl_dRcone_vs_trk_phi = book(m_name, "eop_maxE_ccl_dRcone_vs_trk_phi", "track #phi", nBinsPhi, minPhi, maxPhi, "E/p (higest E cluster with dR<0.2)", nBinsEoP, minEoP, maxEoP);

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHists::execute( const xAOD::TrackParticleContainer* trks, const xAOD::CaloClusterContainer* ccls, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, float eventWeight )
{

  // get pileup
  float mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  float avg_mu(-1.);
  if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) ) {
    avg_mu = eventInfo->averageInteractionsPerCrossing();
  }
  m_mu -> Fill(mu, eventWeight);
  m_avg_mu -> Fill(avg_mu, eventWeight);
  
  // get number of primary vtxs 
  float npv = vtxs->size();
  m_npv -> Fill(npv, eventWeight);

  // track-cluster plots
  float trk_p = -1;
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
  float ccl_e_max_dRcone = 0;
  bool isolatedTrack = true;
  // double ccl_ENG_FRAC_EM

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
    if (fabs(trk->qOverP())>0.) trk_p = (1./fabs(trk->qOverP()))/1e3;
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
    ccl_e_max_dRcone = 0; 
    for( ; ccl_itr != ccl_end; ++ccl_itr ) {
      const xAOD::CaloCluster* ccl = (*ccl_itr);
      trk_ccl_dR = trk->p4().DeltaR( ccl->p4() );
      // std::cout << "\n\n----> ENG_FRAC_EM:\n";
      // std::cout << retrieveMoment(xAOD::CaloCluster::ENG_FRAC_EM, &ccl_ENG_FRAC_EM) << std::endl;
      // std::cout << ccl_ENG_FRAC_EM << std::endl;
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
        // max energy within a cone around the selected track
        if (ccl->e()/1e3 > ccl_e_max_dRcone) ccl_e_max_dRcone = ccl->e()/1e3;
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
      float eop_neccl = neccl_e/trk_p;
      float eop_sum_ccls_dRcone = ccl_e_sum_dRcone/trk_p;
      float eop_maxE_ccl_dRcone = ccl_e_max_dRcone/trk_p;

      m_eop_neccl                      -> Fill(eop_neccl, eventWeight); 
      m_eop_sum_ccls_dRcone            -> Fill(eop_sum_ccls_dRcone, eventWeight); 
      m_eop_maxE_ccl_dRcone            -> Fill(eop_maxE_ccl_dRcone, eventWeight); 

      m_eop_neccl_vs_avg_mu            -> Fill(avg_mu,  eop_neccl, eventWeight); 
      m_eop_neccl_vs_trk_p             -> Fill(trk_p,   eop_neccl, eventWeight); 
      m_eop_neccl_vs_trk_eta           -> Fill(trk_eta, eop_neccl, eventWeight); 
      m_eop_neccl_vs_trk_phi           -> Fill(trk_phi, eop_neccl, eventWeight); 

      m_eop_sum_ccls_dRcone_vs_avg_mu  -> Fill(avg_mu,  eop_sum_ccls_dRcone, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_p   -> Fill(trk_p,   eop_sum_ccls_dRcone, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_eta -> Fill(trk_eta, eop_sum_ccls_dRcone, eventWeight); 
      m_eop_sum_ccls_dRcone_vs_trk_phi -> Fill(trk_phi, eop_sum_ccls_dRcone, eventWeight); 

      m_eop_maxE_ccl_dRcone_vs_avg_mu  -> Fill(avg_mu,  eop_maxE_ccl_dRcone, eventWeight); 
      m_eop_maxE_ccl_dRcone_vs_trk_p   -> Fill(trk_p,   eop_maxE_ccl_dRcone, eventWeight); 
      m_eop_maxE_ccl_dRcone_vs_trk_eta -> Fill(trk_eta, eop_maxE_ccl_dRcone, eventWeight); 
      m_eop_maxE_ccl_dRcone_vs_trk_phi -> Fill(trk_phi, eop_maxE_ccl_dRcone, eventWeight); 

      if (fabs(trk_eta) < 0.8) {
        m_eop_neccl_etaL08 -> Fill(eop_neccl, eventWeight); 
        m_eop_sum_ccls_dRcone_etaL08 -> Fill(eop_sum_ccls_dRcone, eventWeight); 
        m_eop_maxE_ccl_dRcone_etaL08 -> Fill(eop_maxE_ccl_dRcone, eventWeight); 
        if (trk_p > 2.0) {
          m_eop_neccl_etaL08_pG2 -> Fill(eop_neccl, eventWeight); 
          m_eop_sum_ccls_dRcone_etaL08_pG2 -> Fill(eop_sum_ccls_dRcone, eventWeight); 
          m_eop_maxE_ccl_dRcone_etaL08_pG2 -> Fill(eop_maxE_ccl_dRcone, eventWeight); 
        }
      }
      if (fabs(trk_eta) > 0.8 && fabs(trk_eta < 1.7)) {
        m_eop_neccl_etaG08L17 -> Fill(eop_neccl, eventWeight); 
        m_eop_sum_ccls_dRcone_etaG08L17 -> Fill(eop_sum_ccls_dRcone, eventWeight); 
        m_eop_maxE_ccl_dRcone_etaG08L17 -> Fill(eop_maxE_ccl_dRcone, eventWeight); 
        if (trk_p > 2.0) {
          m_eop_neccl_etaG08L17_pG2 -> Fill(eop_neccl, eventWeight); 
          m_eop_sum_ccls_dRcone_etaG08L17_pG2 -> Fill(eop_sum_ccls_dRcone, eventWeight); 
          m_eop_maxE_ccl_dRcone_etaG08L17_pG2 -> Fill(eop_maxE_ccl_dRcone, eventWeight); 
        }
      }
    } // END eoverp histograms

  } // END loop tracks 

  return StatusCode::SUCCESS;
}
