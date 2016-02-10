# E/p analysis for run 2
# Joakim Olsson (joakim.olsson@cern.ch)

from xAH_config import xAH_config
c = xAH_config()

# input containers
# (abbreviations in this analysis: trk - track, ccl - calorimeter cluster)
trks = "InDetTrackParticles"
ccls = "CaloCalTopoClusters"

# selected version
trks_select = trks+"Select"

# selections
# https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPMoriond2016
trk_cutLevel = "LoosePrimary"
trk_maxD0 = 2.0
trk_maxZ0SinTheta = 3.0

''' Set up all the algorithms '''
c.setalg("BasicEventSelection", {"m_name": "BasicEventSelection",
                                 "m_applyGRLCut": False,
                                 "m_cleanPowheg": False,
                                 "m_doPUreweighting": False,
                                 "m_applyPrimaryVertexCut": False,
                                 "m_applyEventCleaningCut": False,
                                 "m_applyCoreFlagsCut": False,
                                 "m_applyTriggerCut": False,
                                 "m_useMetaData": False,
                                 "m_checkDuplicatesData": False,
                                 "m_checkDuplicatesMC": False})

''' Fill histograms with tracking details, passing only basic event selection'''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_JETM5",
                            "m_inContainerName": trks,
                            "m_detailStr": "IPDetails HitCounts",
                            "m_debug": False})

''' Fill histograms with cluster details, passing only basic event selection'''
c.setalg("ClusterHistsAlgo", {"m_name": "Clusters_JETM5",
                              "m_inContainerName": ccls,
                              "m_detailStr": "all",
                              "m_debug": False})

''' Select the Tracks passing the Tracking CP Recommendations (Moriond 2016)'''
c.setalg("TrackSelectorToolWrapper", {"m_name": "TrackSelector_"+trk_cutLevel,
                                      "m_inContainerName": trks,
                                      "m_decorateSelectedObjects": False,
                                      "m_createSelectedContainer": True,
                                      "m_cutLevel": trk_cutLevel,
                                      "m_maxD0": trk_maxD0,
                                      "m_maxZ0SinTheta": trk_maxZ0SinTheta,
                                      "m_outContainerName": trks_select,
                                      "m_useCutFlow": True,
                                      "m_debug": False})

''' Fill histograms with tracking details, after selection '''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_loosePrimary_mD02_mZ0ST3",
                             "m_inContainerName": trks_select,
                             "m_detailStr": "IPDetails HitCounts",
                             "m_debug": False})

''' Track-cluster matching and E/p histograms '''
c.setalg("EoverPAnalysis", {"m_name": "EoverP",
                            "m_inTrackContainerName": trks_select,
                            "m_inClusterContainerName": ccls,
                            "m_detailStr": "all",
                            "m_debug": False})
