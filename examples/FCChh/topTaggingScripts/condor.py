'''
Analysis for top tagging using jet flavour classification.
This analysis stage runs on HTCondor.
'''
from argparse import ArgumentParser
import copy
import urllib.request
import os


# Mandatory: Analysis class where the user defines the operations on the dataframe
class Analysis():
    '''
    Top tagging analysis using jet flavour classification.
    '''
    def __init__(self, cmdline_args):
        # Parse additional arguments
        parser = ArgumentParser(
            description='Additional analysis arguments',
            usage='Provided after "--"')
        # parser.add_argument('--jet-radius', default=0.5, type=float,
        #                     help='Anti-kt jet clustering radius.')
        # parser.add_argument('--pt-min', default=0, type=float,
        #                     help='Minimum pT for jets.')
        # parser.add_argument('--pt-max', default=100, type=float,
        #                     help='Maximum pT for jets.')
        
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        # Mandatory: List of processes used in the analysis
        self.process_list = {
            "mgp8_pp_tt_cQd8_lin_PT5000_5f_84TeV": {"fraction": 0.001, 'chunks': 1},
        #     "mgp8_pp_tt_cQd8_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctu8_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctu8_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctq8_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctq8_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctd8_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctd8_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctG_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_ctG_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQu8_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQu8_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQq83_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQq83_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQq81_lin_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tt_cQq81_quad_PT5000_5f_84TeV": {"fraction": 1, 'chunks': 10},

        #     "mgp8_pp_bb_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_cc_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_gg_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_thadthad_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_tleptlep_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_uuddss_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_whadwhad_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
        #     "mgp8_pp_zhadzhad_PTmin_5000_5f_84TeV": {"fraction": 1, 'chunks': 10},
         }

        # Mandatory: Production tag (input directory)
        self.input_dir = '/eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II_trackCov'

        # Optional: output directory
        self.output_dir = './SMEFT'

        # Optional: analysis name
        self.analysis_name = 'Top Tagging Analysis'

        # Optional: number of threads
        self.n_threads = 64
        self.run_batch = True
        self.eos_type = 'eospublic'

        # Optional: batch queue name
        self.batch_queue = 'testmatch'

        # Optional: computing account
        self.comp_group = 'group_u_FCC.local_gen'

        # Optional: output directory on eos
        #self.output_dir_eos = '/eos/user/h/hfatehi/SMEFT/output/'

        # Optional: include paths for custom functions
        self.include_paths = []

        # Jet flavor tagging model setup
# Jet flavor tagging model setup
        self.model_name = "topTagging-IR9"

        # Directory that actually contains the files
        #self.model_dir = "/eos/user/h/hfatehi/IR9/topTagging-IR9"
        
        self.model_dir = os.path.join(os.getcwd(), self.model_name)

        local_preproc = os.path.join(self.model_dir, f"{self.model_name}.json")
        local_model   = os.path.join(self.model_dir, f"{self.model_name}.onnx")

        self.weaver_preproc = self.get_file_path(local_preproc)
        self.weaver_model   = self.get_file_path(local_model)


        # Helper objects (will be initialized in analyzers)
        self.jetFlavourHelper = None

    @staticmethod
    def get_file_path(filename):
        """Return absolute path to a required local file."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Required file not found: {filename}")
        return os.path.abspath(filename)


    # Mandatory: analyzers function to define the analysis graph
    def analyzers(self, df):
        '''
        Analysis graph for top tagging.
        '''
        from jetFlavourHelper import JetFlavourHelper

        # # Get jet clustering parameters
        # jet_radius = self.ana_args.jet_radius
        # pt_min = self.ana_args.pt_min
        # pt_max = self.ana_args.pt_max

        # Define collections
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
        df = df.Define("jets", "JetClusteringUtils::get_pseudoJets(_jet)")
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)")
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(ReconstructedParticles, _jetc)")
        
        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)")
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
        df = df.Define("event_njet", "JetConstituentsUtils::count_jets(jetc)")
        df = df.Filter("event_njet > 1")

        self.jetFlavourHelper = JetFlavourHelper(
            coll,
            "jets",
            "jetc",
        )
        df = self.jetFlavourHelper.define(df)
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

        df = self.jetFlavourHelper.inference(self.weaver_preproc, self.weaver_model, df)

        return df

    # Mandatory: output function
    def output(self):
        '''
        Output variables which will be saved to output root file.
        '''
        branchList = [
            "jet_e", "jet_mass", "jet_phi", "jet_pT", "jet_theta", "jet_eta", "jet_p",
            "jet_px", "jet_py", "jet_pz", "jet_nconst",
            "event_invariant_mass", "event_njet"
        ]
        branchList += self.jetFlavourHelper.outputBranches()
        return branchList
