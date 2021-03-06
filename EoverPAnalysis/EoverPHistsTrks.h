#ifndef EoverPAnalysis_EoverPHistsTrks_H
#define EoverPAnalysis_EoverPHistsTrks_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

class EoverPHistsTrks : public HistogramManager
{
  public:
    EoverPHistsTrks(std::string name, std::string detailStr, float trkIsoDRmax = .4, float trkIsoPfrac = 0., bool doTrkPcut = false, float trkPmin = 0., float trkPmax = 1e8, bool doTrkEtacut = false, float trkEtamin = 0., float trkEtamax = 1e8, bool doTrkIsocut = false);
    ~EoverPHistsTrks();

    StatusCode initialize();

    StatusCode execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer* vtxs, const xAOD::EventInfo* eventInfo, float eventWeight );
  
    float deltaR (float trk_eta, float trk_phi, float trk2_eta, float trk2_phi);

    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload

  protected:

    float m_trkIsoDRmax; //! track isolation max DR
    float m_trkIsoPfrac; //! track isolation max p fraction
    bool m_doTrkPcut; //!
    float m_trkPmin; //!
    float m_trkPmax;  //!
    bool m_doTrkEtacut; //!
    float m_trkEtamin; //!
    float m_trkEtamax;  //!
    bool m_doTrkIsocut; //!

  private:

    // event level plots
    TH1F* m_mu; //!
    TH1F* m_mu_avg; //!
    TH1F* m_mu_avg_many; //!
    TH1F* m_mu_avg_shift; //!
    TH1F* m_npv; //!
    TH1F* m_npv_trks; //!

    TH2F* m_mu_avg_vs_npv; //!
    TH2F* m_mu_avg_vs_trk_n_nocut; //!

    // track plots
    TH1F* m_trk_n_nocut; //!
    TH1F* m_trk_p_noiso; //!
    TH1F* m_trk_etaID_noiso; //!
    TH1F* m_trk_etaEMB2_noiso; //!
    TH1F* m_trk_etaEME2_noiso; //!
    TH1F* m_trk_phiID_noiso; //!
    TH1F* m_trk_phiEMB2_noiso; //!
    TH1F* m_trk_phiEME2_noiso; //!

    TH2F* m_trk_etaEMB2_vs_etaID_noiso; //!
    TH2F* m_trk_etaEME2_vs_etaID_noiso; //!
    TH2F* m_trk_etaEME2_vs_etaEMB2_noiso; //!

    TH1F* m_trk_trk2_dR_ID; //!
    TH2F* m_trk_trk2_dR_ID_vs_trk_p; //!
    TH1F* m_trk_trk2_dR_EMB2; //!
    TH2F* m_trk_trk2_dR_EMB2_vs_trk_p; //!
    TH1F* m_trk_trk2_dR_EME2; //!
    TH2F* m_trk_trk2_dR_EME2_vs_trk_p; //!
    TH1F* m_trk_trk2_dR_EMB2_EME2; //!
    TH2F* m_trk_trk2_dR_EMB2_EME2_vs_trk_p; //!
    TH1F* m_trk_trk2_dR_EME2_EMB2; //!
    TH2F* m_trk_trk2_dR_EME2_EMB2_vs_trk_p; //!

    TH1F* m_trk_ntrks_maxDR01; //!
    TH1F* m_trk_ntrks_maxDR02; //!
    TH1F* m_trk_ntrks_maxDR03; //!
    TH1F* m_trk_ntrks_maxDR04; //!
    TH1F* m_trk_ntrks_maxDR05; //!
    TH1F* m_trk_ntrks_maxDR06; //!
    TH1F* m_trk_ntrks_maxDR07; //!
    TH1F* m_trk_ntrks_maxDR08; //!
    TH1F* m_trk_ntrks_maxDR09; //!
    TH1F* m_trk_ntrks_maxDR10; //!

    TH2F* m_trk_ntrks_maxDR01_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR02_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR03_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR04_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR05_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR06_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR07_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR08_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR09_vs_trk_p; //!
    TH2F* m_trk_ntrks_maxDR10_vs_trk_p; //!

    TH2F* m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p; //!
    TH2F* m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p; //!
    TH1F* m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p; //!

};

#endif
