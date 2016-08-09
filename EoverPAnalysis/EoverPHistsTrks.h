#ifndef EoverPAnalysis_EoverPHistsTrks_H
#define EoverPAnalysis_EoverPHistsTrks_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

class EoverPHistsTrks : public HistogramManager
{
  public:
    EoverPHistsTrks(std::string name, std::string detailStr, std::string trkExtrapol = "EMB2", float trkIsoDRmax = .4, float trkIsoPfrac = 0., bool doTrkPcut = false, float trkPmin = 0., float trkPmax = 1e8, bool doTrkEtacut = false, float trkEtamin = 0., float trkEtamax = 1e8, bool doTrkIsocut = false);
    ~EoverPHistsTrks();

    StatusCode initialize();

    StatusCode execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer* vtxs, const xAOD::EventInfo* eventInfo, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload

  protected:
    
    std::string m_trkExtrapol; //! layer where tracks are extrapolated
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
    TH1F* m_npv; //!

    // track plots
    TH1F* m_trk_n_nocut; //!
    TH1F* m_trk_n; //!

    TH1F* m_trk_trk2_dR; //!
    TH2F* m_trk_trk2_dR_vs_trk_p; //!

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
