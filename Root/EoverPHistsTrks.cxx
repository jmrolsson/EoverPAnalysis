#include "EoverP/EoverPHistsTrks.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <regex>
#include <math.h>

EoverPHistsTrks :: EoverPHistsTrks (std::string name, std::string detailStr, std::string trkExtrapol, float trkIsoDRmax, float trkIsoPfrac, bool doTrkPcut, float trkPmin, float trkPmax, bool doTrkEtacut, float trkEtamin, float trkEtamax, bool doTrkIsocut) :
  HistogramManager(name, detailStr)
{
  m_trkExtrapol = trkExtrapol; 
  m_trkIsoDRmax = trkIsoDRmax; 
  m_trkIsoPfrac = trkIsoPfrac; 
  m_doTrkPcut = doTrkPcut;
  m_trkPmin = trkPmin;
  m_trkPmax = trkPmax; 
  m_doTrkEtacut = doTrkEtacut;
  m_trkEtamin = trkEtamin;
  m_trkEtamax = trkEtamax; 
  m_doTrkIsocut = doTrkIsocut;
}

EoverPHistsTrks :: ~EoverPHistsTrks () {}

StatusCode EoverPHistsTrks::initialize()
{

  // number of bins and ranges for histograms
  int nBinsMu = 50;           float minMu = -0.5;    float maxMu = 49.5;
  int nBinsNPV = 50;          float minNPV = -0.5;   float maxNPV = 49.5;
  int nBinsTrkN = 200;        float minTrkN = -0.5;  float maxTrkN = 199.5;
  int nBinsE = 300;           float minE = 0;        float maxE = 30;
  unsigned int nBinsDR = 60;  float minDR = 0;       float maxDR = 3;

  //// Book histograms

  // event level plots
  m_mu = book(m_name, "mu", "#mu", nBinsMu, minMu, maxMu); 
  m_mu_avg = book(m_name, "mu_avg", "<#mu>", nBinsMu, minMu, maxMu); 
  m_npv = book(m_name, "npv", "NPV", nBinsNPV, minNPV, maxNPV); 

  // track plots
  m_trk_n_nocut = book(m_name, "trk_n_nocut", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_n = book(m_name, "trk_n", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_trk2_dR = book(m_name, "trk_trk2_dR", "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_vs_trk_p = book(m_name, "trk_trk2_dR_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);

  m_trk_ntrks_maxDR01 = book(m_name, "trk_ntrks_maxDR01", "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02 = book(m_name, "trk_ntrks_maxDR02", "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03 = book(m_name, "trk_ntrks_maxDR03", "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04 = book(m_name, "trk_ntrks_maxDR04", "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05 = book(m_name, "trk_ntrks_maxDR05", "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06 = book(m_name, "trk_ntrks_maxDR06", "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07 = book(m_name, "trk_ntrks_maxDR07", "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08 = book(m_name, "trk_ntrks_maxDR08", "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09 = book(m_name, "trk_ntrks_maxDR09", "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10 = book(m_name, "trk_ntrks_maxDR10", "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_maxDR01_vs_trk_p = book(m_name, "trk_ntrks_maxDR01_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02_vs_trk_p = book(m_name, "trk_ntrks_maxDR02_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03_vs_trk_p = book(m_name, "trk_ntrks_maxDR03_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04_vs_trk_p = book(m_name, "trk_ntrks_maxDR04_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05_vs_trk_p = book(m_name, "trk_ntrks_maxDR05_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06_vs_trk_p = book(m_name, "trk_ntrks_maxDR06_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07_vs_trk_p = book(m_name, "trk_ntrks_maxDR07_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08_vs_trk_p = book(m_name, "trk_ntrks_maxDR08_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09_vs_trk_p = book(m_name, "trk_ntrks_maxDR09_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10_vs_trk_p = book(m_name, "trk_ntrks_maxDR10_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p", "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p", "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsE, minE, maxE, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsE, minE, maxE); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsE, minE, maxE); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsE, minE, maxE, "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsE, minE, maxE); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{avg} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{leading} / p_{trk}", 100, 0, 5.0); 

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHistsTrks::execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer* vtxs, const xAOD::EventInfo* eventInfo, float eventWeight )
{

  // get pileup
  float mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  float mu_avg(-1.);
  if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) ) {
    mu_avg = eventInfo->averageInteractionsPerCrossing();
  }
  m_mu -> Fill(mu, eventWeight);
  m_mu_avg -> Fill(mu_avg, eventWeight);

  // get number of primary vtxs 
  float npv = HelperFunctions::countPrimaryVertices(vtxs, 2);
  m_npv -> Fill(npv, eventWeight);

  // number of tracks
  m_trk_n_nocut-> Fill(trks->size());

  // track-cluster plots
  float trk_n = 0.;
  float trk_p = 0.;
  float trk_etaID = 0; 
  float trk_phiID = 0; 
  float trk_etaCALO = 0; 
  float trk_phiCALO = 0; 
  float dR_CALO_ID = 0; 
  float dEta_CALO_ID = 0; 
  float dPhi_CALO_ID = 0; 
  float trk2_p = 0.;
  float trk2_etaCALO = 0; 
  float trk2_phiCALO = 0; 
  float trk_trk2_dR = 0; 
  float trk_trk2_dEta = 0; 
  float trk_trk2_dPhi = 0; 
  float surr_trk_sum_p = 0.;
  float surr_trk_avg_p = 0.;
  float surr_trk_leading_p = 0.;

  // loop over tracks and clusters around a selected track, plot how many tracks (clusters) fall within a certain dR
  float maxDR_ranges[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();

  // track-track and track-cluster (outer loop over all tracks)
  for( ; trk_itr != trk_end; ++trk_itr ) {

    // count the number of tracks passing all selections
    trk_n += 1;
    const xAOD::TrackParticle* trk = (*trk_itr);
    if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
    trk_etaID = trk->eta();
    trk_phiID = trk->phi();
    trk_etaCALO = trk->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
    trk_phiCALO = trk->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));

    dEta_CALO_ID = TMath::Abs(trk_etaCALO - trk_etaID);
    if (TMath::Abs(trk_phiCALO - trk_phiID) < TMath::Pi())
      dPhi_CALO_ID = TMath::Abs(trk_phiCALO - trk_phiID);
    else
      dPhi_CALO_ID = 2*TMath::Pi() - TMath::Abs(trk_phiCALO - trk_phiID);

    dR_CALO_ID = sqrt( pow(dEta_CALO_ID, 2) + pow(dPhi_CALO_ID, 2) );

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

    // keep track of the number of tracks within a certain dR from the selected track
    int trk_ntrks_maxDR[10] = {0};
    trk2_p = 0.;
    surr_trk_sum_p = 0.;
    surr_trk_leading_p = 0.;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      if (trk_itr != trk2_itr) { // do not double count the selected track 
        const xAOD::TrackParticle* trk2 = (*trk2_itr);
        trk2_etaCALO = trk2->auxdata<float>(std::string("CALO_trkEta_"+m_trkExtrapol));
        trk2_phiCALO = trk2->auxdata<float>(std::string("CALO_trkPhi_"+m_trkExtrapol));
        trk_trk2_dEta= TMath::Abs(trk2_etaCALO - trk_etaCALO);
        if (TMath::Abs(trk2_phiCALO - trk_phiCALO) < TMath::Pi())
          trk_trk2_dPhi = TMath::Abs(trk2_phiCALO - trk_phiCALO);
        else
          trk_trk2_dPhi = 2*TMath::Pi() - TMath::Abs(trk2_phiCALO - trk_phiCALO);
        trk_trk2_dR = sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );
        m_trk_trk2_dR -> Fill(trk_trk2_dR, eventWeight);
        m_trk_trk2_dR_vs_trk_p -> Fill(trk_p, trk_trk2_dR, eventWeight);

        // count the number of track wihin a given radius
        for (int i = 0; i < 10; ++i) {
          if (trk_trk2_dR < maxDR_ranges[i]) {
            trk_ntrks_maxDR[i]++;
          }
        }

        // track isolation
        if (trk_trk2_dR < m_trkIsoDRmax) { // check if trk2 falls within DRmax of trk 

          // calculate the leading and avg p of the surrounding tracks
          if (TMath::Abs(trk2->qOverP())>0.) trk2_p = (1./TMath::Abs(trk2->qOverP()))/1e3; 
          surr_trk_sum_p += trk2_p;
          if (trk2_p > surr_trk_leading_p) {
            surr_trk_leading_p = trk2_p;
          }
        }
      }
    } // END looping trk2

    // plots of the number of close-by tracks within a certain radius
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

    // plots of the number of close-by tracks vs p of the selected track
    m_trk_ntrks_maxDR01_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[9], eventWeight);

    // plots of leading p and avg p of the tracks surrounding the selected track 
    surr_trk_avg_p = surr_trk_sum_p/(trks->size()-1); // calculate avg p of the surrounding tracks 
    m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p -> Fill(surr_trk_sum_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p -> Fill(surr_trk_avg_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p -> Fill(surr_trk_leading_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_sum_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_avg_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_leading_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p -> Fill(surr_trk_sum_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p -> Fill(surr_trk_avg_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p -> Fill(surr_trk_leading_p/trk_p, eventWeight);

  } // END loop tracks 

  m_trk_n-> Fill(trk_n);

  return StatusCode::SUCCESS;
}
