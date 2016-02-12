#ifndef EoverP_TCMatchHists_H
#define EoverP_TCMatchHists_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"

class TCMatchHists : public HistogramManager
{
  public:
    TCMatchHists(std::string name, std::string detailStr );
    ~TCMatchHists();

    StatusCode initialize();
    StatusCode execute( const xAOD::TrackParticleContainer* trks, const xAOD::CaloClusterContainer* ccls, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload

  protected:
    // bools to control which histograms are filled
    bool m_fillDebugging; //!

    // EoverP selections 
    const float m_cut_trk_ccl_maxDR = 0.2; //!
    const float m_cut_trk_trk2_maxDR = 0.4; //!

  private:

    //// 1D histograms

    //track-track plots
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
    
    // track-cluster plots
    TH1F* m_trk_neccl_dR; //!

    TH1F* m_trk_nccls_maxDR01; //!
    TH1F* m_trk_nccls_maxDR02; //!
    TH1F* m_trk_nccls_maxDR03; //!
    TH1F* m_trk_nccls_maxDR04; //!
    TH1F* m_trk_nccls_maxDR05; //!
    TH1F* m_trk_nccls_maxDR06; //!
    TH1F* m_trk_nccls_maxDR07; //!
    TH1F* m_trk_nccls_maxDR08; //!
    TH1F* m_trk_nccls_maxDR09; //!
    TH1F* m_trk_nccls_maxDR10; //!

    // eoverp plots
    TH1F* m_eop_neccl; //!
    TH1F* m_eop_sum_ccls_dRcone; //!

    // 2D histograms

    // plots of selected track and nearest cluster
    TH2F* m_trk_neccl_dR_vs_neccl_e; //! 
    TH2F* m_trk_neccl_trk_eta_vs_neccl_eta ; //! 
    TH2F* m_trk_neccl_trk_phi_vs_neccl_phi ; //! 

    // eoverp plots
    TH2F* m_eop_neccl_vs_trk_pt; //!
    TH2F* m_eop_neccl_vs_trk_eta; //!
    TH2F* m_eop_neccl_vs_trk_phi; //!
    TH2F* m_eop_sum_ccls_dRcone_vs_trk_pt; //!
    TH2F* m_eop_sum_ccls_dRcone_vs_trk_eta; //!
    TH2F* m_eop_sum_ccls_dRcone_vs_trk_phi; //!
};


#endif
