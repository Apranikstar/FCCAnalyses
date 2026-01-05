import json
import os
import copy
import urllib.request


processList = {
    "mgp8_pp_tt_cQd8_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQd8_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctu8_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctu8_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctq8_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctq8_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctd8_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctd8_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctG_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_ctG_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQu8_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQu8_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQq83_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQq83_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQq81_lin_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tt_cQq81_quad_PT5000_5f_84TeV": {"fraction": 1,'chunks': 10},

    "mgp8_pp_bb_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_cc_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_gg_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_thadthad_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_tleptlep_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_uuddss_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_whadwhad_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
    "mgp8_pp_zhadzhad_PTmin_5000_5f_84TeV": {"fraction": 1,'chunks': 10},
}

outputDir = "/eos/user/h/hfatehi/bounds/"
inputDir = "/eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II_trackCov"
nCPUS = -1


model_name = "topTagging-IR9"
url_model_dir = ""
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)
model_dir = "./"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

def get_file_path(url, filename):
    """Return local file path if exists else download from url and return basename."""
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        basename = os.path.basename(url)
        print(f"Downloading {basename} from {url}...")
        urllib.request.urlretrieve(url, basename)
        return os.path.abspath(basename)

weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)




from jetFlavourHelper import JetFlavourHelper

# helpers to be used inside analysers
jetFlavourHelper = None


class RDFanalysis:
    def analysers(df):



        global jetFlavourHelper

        coll = {
        "GenParticles": "Particle",
        "PFParticles": "ReconstructedParticles",
        "PFTracks": "EFlowTrack",
        "PFPhotons": "EFlowPhoton",
        "PFNeutralHadrons": "EFlowNeutralHadron",
        "TrackState": "_EFlowTrack_trackStates",
        "TrackerHits": "TrackerHits",
        "CalorimeterHits": "CalorimeterHits",
        "PathLength": "EFlowTrack_L",
        "Bz": "magFieldBz",
        }


        # Define RP kinematics
        ####################################################################################################
        df = df.Define("RP_px", "ReconstructedParticle::get_px(ReconstructedParticles)")
        df = df.Define("RP_py", "ReconstructedParticle::get_py(ReconstructedParticles)")
        df = df.Define("RP_pz", "ReconstructedParticle::get_pz(ReconstructedParticles)")
        df = df.Define("RP_e", "ReconstructedParticle::get_e(ReconstructedParticles)")
        df = df.Define("RP_m", "ReconstructedParticle::get_mass(ReconstructedParticles)")

        # Define pseudo-jets
        ####################################################################################################
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
        # Anti-kt clustering and jet constituents
        ####################################################################################################
        df = df.Define("_jet", "JetClustering::clustering_antikt(0.5, 0, 100, 0 , 1)(pjetc)")
        df = df.Define("jets","JetClusteringUtils::get_pseudoJets(_jet)" )
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)") 
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(ReconstructedParticles, _jetc)")
        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)" )
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
        df=df.Define("event_njet",   "JetConstituentsUtils::count_jets(jetc)")
        df.Filter("event_njet > 1")

        jetFlavourHelper = JetFlavourHelper(
            coll,
            "jets",
            "jetc",
        )
        df = jetFlavourHelper.define(df)
        df = df.Define("jet_nconst", "JetConstituentsUtils::count_consts(jetc)") 
        

        df = df.Define("jet_p", "JetClusteringUtils::get_p(jets)")
        df = df.Define("jet_px", "JetClusteringUtils::get_px(jets)")
        df = df.Define("jet_py", "JetClusteringUtils::get_py(jets)")
        df = df.Define("jet_pz", "JetClusteringUtils::get_pz(jets)")
        df = df.Define("jet_pT", "JetClusteringUtils::get_pt(jets)")

        df = df.Define("jet_e", "JetClusteringUtils::get_e(jets)")
        df = df.Define("jet_mass", "JetClusteringUtils::get_m(jets)")
        df = df.Define("jet_phi", "JetClusteringUtils::get_phi(jets)")
        df = df.Define("jet_theta", "JetClusteringUtils::get_theta(jets)")
        df = df.Define("jet_eta", "JetClusteringUtils::get_eta(jets)")

        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)


        return df

    def output():

        branchList = [
            "jet_e", "jet_mass", "jet_phi", "jet_pT","jet_theta","jet_eta", "jet_p",  "event_invariant_mass", "event_njet"
                ]
        branchList += jetFlavourHelper.outputBranches()
        return branchList
