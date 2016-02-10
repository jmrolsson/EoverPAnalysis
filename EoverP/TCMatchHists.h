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
    const float m_cut_trk_ccl_dRmax = 0.2; //!

  private:
    // 1D histograms
    TH1F* m_trk_neccl_dR; //!
    TH1F* m_trk_allccl_dRmax; //!
    TH1F* m_eop_all; //!
    // 2D histograms
    TH2F* m_trk_neccl_dR_vs_neccl_e; //! 
    TH2F* m_trk_neccl_trk_eta_vs_neccl_eta ; //! 
    TH2F* m_trk_neccl_trk_phi_vs_neccl_phi ; //! 
    TH2F* m_eop_vs_trk_pt; //!
    TH2F* m_eop_vs_trk_eta; //!
    TH2F* m_eop_vs_trk_phi; //!
};


#endif
