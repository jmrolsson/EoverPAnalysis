# E/p analysis for run 2
# Joakim Olsson (joakim.olsson@cern.ch)

from xAH_config import xAH_config
c = xAH_config()

# input containers
trks = "InDetTrackParticles"

# selected version
trks_select = trks+"Select"

''' Set up all the algorithms '''
c.setalg("BasicEventSelection", {"m_name": "BasicEventSelection",
                                 "m_applyGRLCut": True,
                                 "m_GRLxml": "$ROOTCOREBIN/data/EoverP/data15_13TeV.periodB1_DetStatus-v62-pro18_DQDefects-00-01-02_PHYS_StandardModel_MinimuBias2010.xml",
                                 "m_doPUreweighting": False,
                                 #"m_PRWFileNames": "dev/PileupReweighting/mc15a_defaults.NotRecommended.prw.root",
                                 "m_applyPrimaryVertexCut": True,
                                 "m_PVNTrack": 4,
                                 "m_applyEventCleaningCut": True,
                                 "m_applyCoreFlagsCut": False,
                                 "m_applyTriggerCut": True,
                                 "m_triggerSelection": "HLT_noalg_mb_L1MBTS_1",
                                 "m_useMetaData": False,
                                 "m_checkDuplicatesData": False,
                                 "m_checkDuplicatesMC": False})

''' Fill histograms with tracking details, passing only basic event selection'''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_BasicEvtSel",
                            "m_inContainerName": trks,
                            "m_detailStr": "2D IPDetails HitCounts",
                            "m_debug": False})

# https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPMoriond2016
''' Select the Tracks passing the Tracking CP Recommendations (Moriond 2016)'''
c.setalg("TrackVertexSelection", {"m_name": "TrackSelector_LoosePrimary",
                                   "m_inContainerName": trks,
                                   "m_decorateSelectedObjects": False,
                                   "m_createSelectedContainer": True,
                                   "m_pass_min": 1.0,
                                   "m_cutLevel": "LoosePrimary",
                                   "m_maxD0": 2.0,
                                   "m_maxZ0SinTheta": 3.0,
                                   "m_minNTrtHits": 20,
                                   "m_outContainerName": trks_select,
                                   "m_useCutFlow": False,
                                   "m_debug": False})

''' Fill histograms with tracking details, after selection '''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_PassTrkSel",
                            "m_inContainerName": trks_select,
                            "m_detailStr": "2D IPDetails HitCounts",
                            "m_debug": False})

''' E/p histograms '''
c.setalg("EoverPAnalysis_eopxAOD", {"m_name": "EoverP_PassTrkSel",
                                    "m_inTrackContainerName": trks_select,
                                    "m_trkExtrapol": "EMB2",
                                    "m_doBgSubtr" : True,
                                    "m_doEMcalib": True,
                                    "m_doLCWcalib": True,
                                    "m_doCells": True,
                                    "m_doCaloEM": True,
                                    "m_doCaloHAD": True,
                                    "m_doCaloTotal": True,
                                    "m_trkIsoDRmax": .4,
                                    "m_doTrkPcut": True,
                                    "m_trkPmin": 0.,
                                    "m_trkPmax": 1e8,
                                    "m_doTrkEtacut": True,
                                    "m_trkEtamin": 0.,
                                    "m_trkEtamax": 1e8,
                                    "m_doTileCuts": True,
                                    "m_LarEmax": 1e8,
                                    "m_TileEfrac": -1,
                                    "m_doEtaPranges": True,
                                    "m_detailStr": "all",
                                    "m_debug": False})
