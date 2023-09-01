#!/bin/python3

# ----------- instructions -------------

# 1. compile gambit as normal

# 2. copy the gambit executable into files/
#    (overwrite the existing one)

# 3. setup the options below to you liking

# 4. set any other options in the yaml files
#    found in files/yaml_files/ but don't touch
#    any line containing '~'

# 5. setup your desired plots in files/plots.pip
#    but don't touch any line containing '~'

# 6. modify files/job.sh so that it works on
#    your HPC but don't touch any line containing '~'

# 7. run this python script

# 8. run the file: 'gen_path/runScans.sh'

# 9. run the file: 'gen_path/runPippi.sh'


# ----------- options -------------

print("\n-----------------\nrunning 2HDM generator\n")

# a GAMBIT is generated for all possible combinations below
# and a script is generated for running all simultaneously
# WARNING: be careful not to generate too many, otherwise you will run out of storage

yaml_dir = 'yaml_files_full'
# yaml_dir = 'yaml_files_med'
# yaml_dir = 'yaml_files_small'
use_speed_hacks = True # todo: maybe delete this ??
gen_path = 'gens'

conv_threshold = 1e-8
NP = 5000
required_printed_points = 5000000
required_points = -1
required_scan_duration = 1*35*60 # in seconds

NODE_COUNT = 1  # set to desired number of nodes per gambit
CORE_COUNT = 76 # set to number of cores per node

# either "DIRAC" or "BASH"
MODE = "DIRAC"

# allowed options: "THDM", "THDMI", "THDMII", "THDMLS", or "THDMflipped"
models = ["THDMI"]

# allowed options: "tree" or "loop"
runnings = ["loop"]

# allowed options: "generic", "hybrid_lambda_1", "hybrid_lambda_2", "hybrid_Higgs", "higgs", or "physical"

# instead of setting the basis, we set the file. This will allow us to run targetted scans in the same basis.
# not that we cant simply change the basis as the params will be wrong

# 22
bases = [
    
    # ("hybrid_Higgs", "test1"),
    # ("hybrid_Higgs", "test2"),
    # ("hybrid_Higgs", "test3"),
    # ("hybrid_Higgs", "test4"),

    # ("physical", "physical"),

    # ("generic", "generic"),
    # ("generic", "genericA"),
    # ("generic", "genericB"),
    # ("generic", "genericC"),
    # ("generic", "genericD"),

    # ("hybrid_Higgs", "hybrid1"),
    # ("hybrid_Higgs", "hybrid1A"),
    # ("hybrid_Higgs", "hybrid1B"),
    # ("hybrid_Higgs", "hybrid1C"),
    # ("hybrid_Higgs", "hybrid1D"),
    # ("hybrid_Higgs", "hybrid1E"),
    # ("hybrid_Higgs", "hybrid1F"),
    # ("hybrid_Higgs", "hybrid1G"),
    # ("hybrid_Higgs", "hybrid1H"),
    # ("hybrid_Higgs", "hybrid1I"),
    # ("hybrid_Higgs", "hybrid1J"),
    # ("hybrid_Higgs", "hybrid1K"),
    # ("hybrid_Higgs", "hybrid1L"),
    # ("hybrid_Higgs", "hybrid1M"),
    # ("hybrid_Higgs", "hybrid1N"),
    # ("hybrid_Higgs", "hybrid1O"),
    # ("hybrid_Higgs", "hybrid1P"),
    # ("hybrid_Higgs", "hybrid1Q"),
    # ("hybrid_Higgs", "hybrid1R"),

    # ("hybrid_Higgs2", "hybrid2"),
    # ("hybrid_Higgs2", "hybrid2A"),
    # ("hybrid_Higgs2", "hybrid2B"),
    # ("hybrid_Higgs2", "hybrid2C"),
    # ("hybrid_Higgs2", "hybrid2D"),
    # ("hybrid_Higgs2", "hybrid2E"),

]

# allowed options: all, theory, collider, electroweak, flavour
#                 (the name of any individual constraint)
#                 (in the case of perturbativity, or unitarity just name the function)


# 39-2

# ---- COMBINED SCANS ----

# 1
# constraints = [

#     (["all"], "all"),

# ]

# ---- THEORY SCANS ----

# 5
constraints = [

    # (["perturbativity_LogLikelihood_THDM"], "theoryG"),
    # (["theory"], "flavortest2_flat_grid36"),
    # (["scalar_mass_corrections_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM"], "NLO"),
    # (["scalar_mass_corrections_LogLikelihood_THDM", "stability_LogLikelihood_THDM"], "stability"),
    # (["runToScaleTest_LogLikelihood_THDM", "scalar_mass_corrections_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM"], "perturbativity"),

]
    # (["perturbativity_yukawas_LogLikelihood_THDM"], "pert_yukawas")

# ---- ELECTROWEAK SCANS ----

# 1
# constraints = [

#     (["theory", "electroweak"], "electroweak"),

# ]

# ---- COLLIDER SCANS ----

# # 6
# constraints = [

# #    #  (["theory", "electroweak"], "electroweak"),

# #    #  (["theory", "collider"], "collider"),
#     (["theory", "LEP_Higgs_LogLike"], "HB"),
# #    #  (["theory", "LHC_Higgs_LogLike", "HS_RUN1_SS"], "HSRUN1SS"),
# #    #  (["theory", "LHC_Higgs_LogLike", "HS_LATEST_SS"], "HSLATESTSS"),
# #    #  (["theory", "LHC_Higgs_LogLike", "HS_LATEST_STXS"], "HSLATESTSTXS"),

# ]
#     # (["theory", "LHC_Higgs_LogLike", "HS_ALL"], "HS"),

# ---- FLAVOR SCANS ----

# # 26-2
# constraints = [

#     # (["theory", "flavour"], "flavour"),
#     # (["theory", "deltaMB_LogLikelihood"], "deltaMB"),
#     (["theory", "deltaMBd_LogLikelihood"], "deltaMBd"),
#     # (["theory", "Bs2ll_LogLikelihood"], "Bs2ll"),
#     # (["theory", "B2Kll_LogLikelihood"], "B2Kll"),
#     # (["theory", "B2mumu_LogLikelihood_Atlas"], "B2mumuAtlas"),
#     # (["theory", "B2mumu_LogLikelihood_LHCb"], "B2mumuLHCb"),
#     # (["theory", "B2mumu_LogLikelihood_CMS"], "B2mumuCMS"),
#     # (["theory", "Bc_lifetime_LogLikelihood"], "Bclifetime"),
#     # (["theory", "B2Xsnunu_LogLikelihood"], "B2Xsnunu"),
#     # (["theory", "SL_LogLikelihood"], "SL"),
#     (["theory", "b2sgamma_LogLikelihood"], "b2sgamma"),
#     # (["theory", "B2Kstargamma_LogLikelihood"], "B2Kstargamma"),
#     # (["theory", "RK_LogLikelihood_LHCb"], "RK"),
#     # (["theory", "RKstar_LogLikelihood_LHCb"], "RKstar"),
#     # (["theory", "B2KstarmumuAng_LogLikelihood_Atlas"], "B2KstarmumuAngAtlas"),
#     # (["theory", "B2KstarmumuAng_LogLikelihood_CMS"], "B2KstarmumuAngCMS"),
#     # (["theory", "B2KstarmumuAng_LogLikelihood_LHCb_2020"], "B2KstarmumuAngLHCb2020"),
#     # (["theory", "B2KstarmumuAng_LogLikelihood_Belle"], "B2KstarmumuAngBelle"),
#     # (["theory", "B2KstarellellAng_LogLikelihood_Belle"], "B2KstarellellAng"),
#     # (["theory", "Bu2KstarmumuAng_LogLikelihood_LHCb_2020"], "Bu2KstarmumuAng"),
#     # (["theory", "B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020"], "B2KstareeAngLowq2"),
#     # (["theory", "B2KstarmumuBr_LogLikelihood_LHCb"], "B2KstarmumuBr"),
#     # (["theory", "B2KmumuBr_LogLikelihood_LHCb"], "B2KmumuBr"),
#     # (["theory", "Bs2phimumuBr_LogLikelihood"], "Bs2phimumuBr"),
#     # (["theory", "B2mumu_LogLikelihood_Atlas", "B2mumu_LogLikelihood_LHCb", "B2mumu_LogLikelihood_CMS"], "B2mumu"),
#     # (["theory", "B2KstarmumuAng_LogLikelihood_Atlas", "B2KstarmumuAng_LogLikelihood_CMS", "B2KstarmumuAng_LogLikelihood_LHCb_2020", "B2KstarmumuAng_LogLikelihood_Belle"], "B2KstarmumuAng"),

# ]

# note that all data for bases is combined for plotting
# whereas each [models,runnings,constraints] generate different sets of plots


# ----------- main script -------------

import copy
import re
import os
import shutil
import math
from distutils.dir_util import copy_tree
from distutils.dir_util import remove_tree
from sys import platform
from pathlib import Path
import yaml
import numpy as np
from math import log
from math import exp

len_bases = len(bases)

for (basis,file) in bases:
    name = "files/" + yaml_dir + "/" + file + ".yaml"
    file = open(name,'r').read()
    file = re.sub(r"^!import.*",r"",file,flags=re.MULTILINE)
    file = re.sub(r"!import",r"",file)
    yfile = yaml.safe_load(file)
    if "subscans" in yfile:
      len_bases += yfile["subscans"]["num_scans"] - 1

print("\nthese scans will eat ", len(models)*len(runnings)*len_bases*len(constraints)*required_scan_duration*CORE_COUNT/(60*60), " CPU hours\n\n")


def generate_gambit_name():
    generate_gambit_name.counter += 1
    name = "gambit_" + str(generate_gambit_name.counter)
    gambit_dirs.append(name)
    return name

generate_gambit_name.counter = 0
gambit_dirs = []

class Options:

    # the names of all constraints
    constraints_theory = ["runToScaleTest_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM", "LO_unitarity_LogLikelihood_THDM", "stability_LogLikelihood_THDM", "higgs_exp_mass_LogLikelihood_THDM", 
                         "scalar_mass_corrections_LogLikelihood_THDM", "higgs_scenario_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM", 
                         "perturbativity_yukawas_LogLikelihood_THDM"]
    constraints_collider = ["LEP_Higgs_LogLike", "LHC_Higgs_LogLike","HS_ALL","HS_RUN1_SS","HS_LATEST_SS","HS_LATEST_STXS"] #higgs_mass_LogLikelihood
    constraints_electroweak =  ["oblique_parameters_LogLikelihood"] # lnL_gm2
    constraints_flavour =  ["Bs2ll_LogLikelihood","B2Kll_LogLikelihood","B2mumu_LogLikelihood_Atlas","B2mumu_LogLikelihood_LHCb","B2mumu_LogLikelihood_CMS",
                            "dBRBDstartaunu_LogLikelihood","dBRBDtaunu_LogLikelihood","gmu_ge_LogLikelihood","FLDstar_LogLikelihood","Bc_lifetime_LogLikelihood",
                            "B2Xsnunu_LogLikelihood","h2taumu_LogLikelihood","t2ch_LogLikelihood","deltaMB_LogLikelihood","deltaMBd_LogLikelihood",
                            "SL_LogLikelihood","l2lgamma_LogLikelihood","l2lll_LogLikelihood","b2sgamma_LogLikelihood","B2Kstargamma_LogLikelihood",
                            "RK_LogLikelihood_LHCb","RKstar_LogLikelihood_LHCb","B2KstarmumuAng_LogLikelihood_Atlas","B2KstarmumuAng_LogLikelihood_CMS",
                            "B2KstarmumuAng_LogLikelihood_LHCb_2020","B2KstarmumuAng_LogLikelihood_Belle","B2KstarellellAng_LogLikelihood_Belle",
                            "Bu2KstarmumuAng_LogLikelihood_LHCb_2020","B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020","B2KstarmumuBr_LogLikelihood_LHCb",
                            "B2KmumuBr_LogLikelihood_LHCb","Bs2phimumuBr_LogLikelihood"]
    constraints_all = constraints_theory + constraints_collider + constraints_electroweak + constraints_flavour

    def __init__(self):

        # default scan options
        self.full_model_name = None
        self.model = None
        self.running = None
        self.basis = None

        # default convergence criteria
        self.conv_threshold = 1e-5
        self.NP = 40000
        self.required_printed_points = -1
        self.required_points = -1
        self.required_scan_duration = 7*24*60*60

        # set default constraints
        for c in self.constraints_all:
            setattr(self, c, False)

        # default paths
        self.results_folder = "../runs"
        self.plots_folder = "../plots"
        self.hdf5_name = "scan"

    def setModel(self, model, basis, running):

        self.model = model
        self.running = running
        self.basis = basis
        basis = "_" + basis

        if self.basis == "coupling" or self.basis == "generic" or self.basis == "general":
            basis = ""
        if running == "loop":
            running = "atQ"
        else:
            running = ""

        self.full_model_name = model + basis + running

    def setConstraint(self, name):

        setit = []

        if name == "all":
            setit = self.constraints_all
        elif name == "theory":
            setit = self.constraints_theory
        elif name == "collider":
            setit = self.constraints_collider
        elif name == "electroweak":
            setit = self.constraints_electroweak
        elif name == "flavour":
            setit = self.constraints_flavour
        else:
            setit = [name]
            if not hasattr(self, name):
                raise Exception("Error: no constraint called: " + name)

        for c in setit:
            setattr(self, c, True)

    def validate(self):

        # OFF FOR NOW
        self.B2KstarellellAng_LogLikelihood_Belle = False
        self.B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020 = False

        # DOES NOTHING likelihoods
        self.dBRBDstartaunu_LogLikelihood = False
        self.dBRBDtaunu_LogLikelihood = False
        self.gmu_ge_LogLikelihood = False
        self.FLDstar_LogLikelihood = False
        self.h2taumu_LogLikelihood = False
        self.t2ch_LogLikelihood = False
        self.l2lgamma_LogLikelihood = False
        self.l2lll_LogLikelihood = False

        # !!!!!!!!!!!!!
        # self.perturbativity_yukawas_LogLikelihood_THDM = False

        self.higgs_scenario_LogLikelihood_THDM = True
        self.higgs_exp_mass_LogLikelihood_THDM = True

        # # NOT WORKING likelihoods
        # self.deltaMB_LogLikelihood = False
        # self.deltaMBd_LogLikelihood = False

        # unitarity 
        if self.running == "tree" and self.NLO_unitarity_LogLikelihood_THDM:
            self.NLO_unitarity_LogLikelihood_THDM = False
            self.LO_unitarity_LogLikelihood_THDM = True

        if self.running == "loop" and self.LO_unitarity_LogLikelihood_THDM and self.NLO_unitarity_LogLikelihood_THDM:
            self.LO_unitarity_LogLikelihood_THDM = False
            self.NLO_unitarity_LogLikelihood_THDM = True

        # correction checks
        if self.running == "tree":
            self.scalar_mass_corrections_LogLikelihood_THDM = False
            self.runToScaleTest_LogLikelihood_THDM = False

        if self.HS_ALL:
            self.HS_RUN1_SS = False
            self.HS_LATEST_SS = False
            self.HS_LATEST_STXS = False

        if self.HS_RUN1_SS:
            self.HS_LATEST_SS = False
            self.HS_LATEST_STXS = False

        if self.HS_LATEST_SS:
            self.HS_RUN1_SS = False
            self.HS_LATEST_STXS = False

        if self.HS_LATEST_STXS:
            self.HS_RUN1_SS = False
            self.HS_LATEST_SS = False

def makeGambit(options, dir):

    print('making: ' + os.path.abspath(gen_path) + "/" + dir + " ... ")
    dir2 = gen_path+"/" + dir

    # delete contents if it already exists
    if os.path.exists(dir2):
        shutil.rmtree(dir2)

    # copy the files to a new gambit dir
    copy_tree("files/", dir2)
    # os.rename(gen_path+"/files", dir2)

    # figure out the yaml file name
    yaml_name = options.file + ".yaml"

    # patch the yaml file
    patchYaml(options, dir, yaml_name)
    patchYaml(options, dir, "THDM_constraints.yaml")
    remove_tree(dir2 + "/" + yaml_dir)

    # patch the run script
    if MODE == "BASH":
        patchRunScript(options, dir, yaml_name)
    elif MODE == "DIRAC":
        patchRunScriptJC(options, dir, yaml_name)
    else:
        raise Exception("unknown mode")
    
    os.remove(dir2 + "/job_jc.sh")
    os.remove(dir2 + "/job_dirac.sh")

    # create the output folders (otherwise hdf5_v1 will crash)
    Path(dir2 + '/yaml_files/' + options.results_folder + '/samples').mkdir(parents=True, exist_ok=True)

    print("done")
   
def patchYaml(options, dir, yaml_name):

    #  hack
    if yaml_name == "THDM_constraints.yaml":


        # read yaml file into string
        file = open(gen_path+"/" + dir + "/" + yaml_dir + "/" + yaml_name, 'r')
        s = file.read()
        file.close()

    else:
        s = options.subscan

    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_med_hhs")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_med_final_hhs")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_full")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_med")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_med_final")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_small")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_small_final")
    # shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_idm")

    # remove useless yaml dirs

    print("DEBUG: patching " + gen_path+"/" + dir + "/" + yaml_dir + "/" + yaml_name)

    # # set basis
    # s = s.replace("prior_type: tanb", "prior_type: " + "flat")

    # set the model name
    s = s.replace("TheModelName", options.full_model_name)

    # set the scan duration and point limit
    s = s.replace("12121212", str(int(options.required_points)))
    s = s.replace("23232323", str(int(options.required_printed_points)))
    s = s.replace("34343434", str(int(options.required_scan_duration)))
    s = s.replace("convthresh:", "convthresh: " + str(options.conv_threshold) + " #")
    s = s.replace("NP:", "NP: " + str(options.NP) + " #")

    # set the output hdf5 name
    s = s.replace("scan.hdf5", options.hdf5_name + ".hdf5")

    # set the output folder
    s = s.replace("the_output_folder", options.results_folder)

    # add parameters for loop-level models
    if options.running == "loop":
        s = s.replace("#~~Qin:","Qin:")

    # comment other scale for tree level
    if options.running == "tree":
        s = s.replace("check_other_scale:", "check_other_scale: -1 # ")
    
    # set the speed hacks
    if use_speed_hacks:
        s = s.replace("only_perturbativity: false","only_perturbativity: false")
    else:
        s = s.replace("only_perturbativity: false","only_perturbativity: true")


    # uncomment required constraints
    for c in  options.constraints_all:
        if getattr(options,c):

            # do the HS Analyses separately
            if c == "HS_ALL":
                s = s.replace("HS_analysis:", "HS_analysis: -1 #")
                continue
            if c == "HS_RUN1_SS":
                s = s.replace("HS_analysis:", "HS_analysis: 0 #")
                continue
            if c == "HS_LATEST_SS":
                s = s.replace("HS_analysis:", "HS_analysis: 1 #")
                continue
            if c == "HS_LATEST_STXS":
                s = s.replace("HS_analysis:", "HS_analysis: 2 #")
                continue

            # get the marker
            i = s.find("capability: " + c)
            if i == -1:
                i = s.find("function: " + c)
            if i != -1:
                marker = s[i-1]
                t = s[i-2]
                if t == "~":
                    s = s.replace("#~"+marker, "")

    # write patched file to disk
    file = open(gen_path+"/" + dir + "/yaml_files/" + yaml_name, 'w')
    file.write(s)
    file.close()

def patchRunScript(options, dir, yaml_name):

    # load run script into string
    file = open(gen_path+"/" + dir + "/job.sh", 'r')
    s = file.read()
    file.close()
    
    # set the yaml name
    s = s.replace("THDM_physical.yaml", yaml_name)

    # write patched file to disk
    file = open(gen_path+"/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

def patchRunScriptJC(options, dir, yaml_name):

    # load run script into string
    file = open(gen_path+"/" + dir + "/job_dirac.sh", 'r')
    s = file.read()
    file.close()

    # set the job name
    s = s.replace("SBATCH -J gambit_thdm", "SBATCH -J thdm_" + options.hdf5_name)

    # set the stdout and stderr paths (todo)
    dir_abs = os.path.abspath(gen_path+"/" + dir)
    # s = s.replace("stdo_", dir + "/work/stdo_")
    # s = s.replace("stde_", dir + "/work/stde_")

    # set gambit dir
    s = s.replace("GAMBIT_BASE_DIR", dir_abs)

    # set yaml name
    s = s.replace("THDM_physical.yaml", yaml_name)

    # set scan name
    # s = s.replace("scan.tar.gz", options.outputName + ".tar.gz")
    # s = s.replace("scan.hdf5", options.outputName + ".hdf5")

    # set scan duration

    secs = 60 + 1.03*options.required_scan_duration
    mins = secs / 60.
    hours = mins / 60.
    secs = math.floor(secs) % 60
    mins = math.floor(mins) % 60
    hours = math.floor(hours)
    
    print("HH:MM:SS " + str(hours) + ":" + str(mins) + ":" + str(secs))

    s = s.replace("878787", str(math.floor(options.required_scan_duration + 40)))
    s = s.replace("MMMMM", str(mins))
    s = s.replace("SSSSS", str(secs))
    s = s.replace("HHHHH", str(hours))

    # write patched file to disk
    file = open(gen_path+"/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

# load a yaml file with (optionally) a subscan node and convert to string
def load_subscan_yaml(name):

    subscans = []

    file = open(name,'r').read()

    # quick and dirty hack to deal with imports
    file = re.sub(r"^!import.*",r"",file,flags=re.MULTILINE)
    file = re.sub(r"!import",r"",file)

    yfile = yaml.safe_load(file)

    # check if we have a subscan node

    num_scans = 1
    overlap = 0.0

    if "subscans" in yfile:
        num_scans = yfile["subscans"]["num_scans"]
        overlap = yfile["subscans"]["overlap"]

    # get param names
    params = [s for s in yfile["Parameters"]["TheModelName"]]

    tmp = yfile["Parameters"]["TheModelName"]

    # get parameter subscan weights
    param_weights = { }

    # float, [f, f], [f, f, f, f]
    param_ranges = { }

    # fixed_value (x2), flat/log, double_log_flat_join
    param_priors = { }

    # get the weights, ranges and priors
    for param in params:
        
        # get a simple fixed value prior
        if type(tmp[param]) is not dict:
            param_weights[param] = 0
            param_ranges[param] = tmp[param]
            param_priors[param] = "fixed_value"
        
        # get a normal map-type prior
        else:
            # a slightly less simple fixed value
            if "fixed_value" in tmp[param]:
                param_weights[param] = 0
                param_ranges[param] = tmp[param]["fixed_value"]
                param_priors[param] = "fixed_value"

            # otherwise, it should have a prior_type and range/ranges
            else:
                param_priors[param] = tmp[param]["prior_type"]

                if "range" in tmp[param]:
                    param_ranges[param] = tmp[param]["range"]
                else:
                    param_ranges[param] = tmp[param]["ranges"]

            # get the subscan_weight if any
            if "subscan_weight" in tmp[param]:
                param_weights[param] = tmp[param]["subscan_weight"]
            else:
                param_weights[param] = 0

    # delete the unwanted nodes (which are not compatible with gambit yet)

    if "subscans" in yfile:
        del yfile["subscans"]

    for param in params:
        if type(yfile["Parameters"]["TheModelName"][param]) is dict:
            if "subscan_weight" in yfile["Parameters"]["TheModelName"][param]:
                del yfile["Parameters"]["TheModelName"][param]["subscan_weight"]
    
    # figure out number of divisions

    pieces = { param:1 for param in params }

    sum_scans_so_far = 1

    while sum_scans_so_far < num_scans:

        weighted_count = [ param_weights[param]/pieces[param] for param in params ]
        sort_indices = np.argsort(weighted_count)
        next_param = params[sort_indices[-1]]
        pieces[next_param] += 1
        sum_scans_so_far = np.prod(list(pieces.values()))

    num_scans = sum_scans_so_far

    # loop over the number of subscans
    for i in range(0,num_scans):

        prod_pieces = 1

        # update the parameter ranges
        for param in params:

            # don't divide up fixed values!
            if (param_priors[param]) == "fixed_value": continue
            
            # get number of parameter range pieces & index
            nPieces = pieces[param]
            if (nPieces) == 1: continue
            slice_index = (i // prod_pieces) % nPieces
            prod_pieces *= nPieces

            # divide up current param range
            param_range = copy.deepcopy(param_ranges[param])
            param_range[0] = float(param_range[0])
            param_range[1] = float(param_range[1])

            # print(param_range)

            # for now, double_log_flat_join is not supported

            # deal with log scale
            if (param_priors[param]) == "log":
              param_range[0] = log(param_range[0])
              param_range[1] = log(param_range[1])

            # calculate the new parameter range
            rangee = (param_range[1] - param_range[0])/nPieces
            param_range_new = [param_range[0]+rangee*slice_index,param_range[0]+rangee*(1+slice_index)]

            # add the overlap
            overlap_range = overlap*rangee
            param_range_new[0] -= overlap_range
            param_range_new[1] += overlap_range

            # deal with log scale
            if (param_priors[param]) == "log":
                param_range_new[0] = exp(param_range_new[0])
                param_range_new[1] = exp(param_range_new[1])
                param_range = copy.deepcopy(param_ranges[param])
                param_range[0] = float(param_range[0])
                param_range[1] = float(param_range[1])

            # don't go outside the full range
            param_range_new[0] = max(param_range_new[0], param_range[0])
            param_range_new[1] = min(param_range_new[1], param_range[1])

            yfile["Parameters"]["TheModelName"][param]["range"] = param_range_new
            # yfile["Parameters"]["TheModelName"][param]["shift"] = -param_range_new[0]
            # yfile["Parameters"]["TheModelName"][param]["shift"] = -max(0.0,param_range[0]-0.1)

        # convert yaml file to string
        contents = yaml.dump(yfile)

        # add Qin back
        contents += "    #~~Qin: 91.1876 # = mZ\n"

        # add the imports back
        contents += "!import ../yaml_files/THDM_constraints.yaml\n"
        contents = contents.replace("StandardModel_SLHA2:", "StandardModel_SLHA2: !import")

        # add to list of yamls
        subscans.append(contents)

    # print results

    # for i in range(0,len(subscans)):
    #     with open("dump"+str(i)+".yaml","w") as f:
    #         f.write(subscans[i])

    return subscans


def main():

    # ad postfix to gen_path so that it is unique
    global gen_path
    postfix = 0
    while os.path.exists(gen_path+"_"+str(postfix)):
        postfix += 1
    postfix = str(postfix)
    gen_path = gen_path + "_" + postfix


    folders_to_merge = {}

    # --- loop over all options ---

    # loop over all models (different output folder)
    for model in models:

        # loop over all runnings (different output folder)
        for running in runnings:

            # loop over all constraints (different output folder)
            for (constraint, constraint_name) in constraints:

                # uniquely identifies results folders (to be merged later)
                resultsSuffix = 0

                # loop over all bases (merged output folder)
                for (basis,file) in bases:

                    # get the list of subscan files
                    subscans = load_subscan_yaml("files/" + yaml_dir + "/" + file + ".yaml")

                    # loop over all subscans (merged output folder)
                    for subscan in subscans:

                        resultsSuffix += 1
                        resultsSuffixStr = "_" + str(resultsSuffix) if MODE != "BASH" else ""

                        # setup the scan-specific options
                        options = Options()
                        options.setModel(model, basis, running)
                        options.subscan = subscan
                        options.file = file
                        for c in constraint:
                            options.setConstraint(c)

                        # setup convergence criteria
                        options.conv_threshold = conv_threshold
                        options.NP = NP
                        options.required_printed_points = required_printed_points
                        options.required_points = required_points
                        if required_printed_points != -1:
                            options.required_printed_points /= (len_bases*CORE_COUNT)
                        if required_points != -1:
                            options.required_points /= (len_bases*CORE_COUNT)
                        options.required_scan_duration = required_scan_duration

                        # make sure everything is valid
                        options.validate()

                        # setup paths (also stored in options)
                        fullName = model + running + "_" + constraint_name
                        options.results_folder = "../../../runs/" + fullName + resultsSuffixStr + "/"
                        options.plots_folder = "../plots/" + fullName + "/"
                        options.hdf5_name = file + "_" + postfix + "_" + resultsSuffixStr
                        
                        # dict that tells us which folders to merge (essentially we will just get rid of resultsSuffixStr)
                        tmp = "../runs/" + fullName + resultsSuffixStr + "/"
                        folders_to_merge[tmp] = "../runs/" + fullName + "/"

                        # for i in range(0,CORE_COUNT):
                        #     folders_to_merge[options.results_folder + 'samples/' + options.hdf5_name + ".hdf5_temp_" + str(i)] = "../../../runs/" + model + running + constraint_name_full + "/"

                        # make a new gambit with the options specified above
                        makeGambit(options, generate_gambit_name())

    parent_abs = os.path.abspath(gen_path)
    if MODE != "BASH":
        with open(parent_abs + "/merge.py", "w") as f:
            f.write("#!/bin/python3\n")
            f.write("import os, shutil, pathlib, fnmatch\n")
            f.write("def move_dir(src: str, dst: str, pattern: str = '*'):\n")
            f.write("    if not os.path.isdir(src):\n")
            f.write("        return\n")
            f.write("    if not os.path.isdir(dst):\n")
            f.write("        pathlib.Path(dst).mkdir(parents=True, exist_ok=True)\n")
            f.write("    for f in fnmatch.filter(os.listdir(src), pattern):\n")
            f.write("        shutil.move(os.path.join(src, f), os.path.join(dst, f))\n\n")
            

            # f.write("cd gambit_1/yaml_files\n")
            # f.write( "import os\n")
            # f.write( "import shutil\n")
            # f.write( "from pathlib import Path\n")
            for src,dst in folders_to_merge.items():
                f.write('move_dir(\"{0}samples\",\"{1}samples\")\n'.format(src,dst))
                f.write("if os.path.isdir(\"{0}\"):\n".format(src))
                f.write("    os.system('rm -rf \"{0}\"')\n\n".format(src))
                # f.write( 'Path("{0}samples").mkdir(parents=True, exist_ok=True)\n'.format(dst))
                # f.write( 'if os.path.exists(\"{0}samples\"):\n'.format(src))
                # f.write( '    files = os.listdir(\"{0}samples\")\n'.format(src))
                # f.write( '    for f in files:\n')
                # f.write( '        shutil.move(f.name,\"{0}samples\")\n'.format(dst))

                # f.write("mv \"{0}\" \"{1}samples\"\n".format(k,v))
                # f.write("rsync -a \"{0}\" \"{1}\"\n".format(k,v))
                # f.write("rm -rf \"{0}\"\n".format(k))

        if platform == "linux" or platform == "linux2":
            os.system("chmod +x " + parent_abs + "/merge.py")

    # --- create the runner script --- 

    # write it
    print("generating runner script...")
    file = open(parent_abs + "/runScans.sh", "w")
    for i,dir in enumerate(gambit_dirs):
        file.write("echo \"------------------------------\"\n")
        file.write("echo \"---------- "+str(i)+" of " + str(len(gambit_dirs)) + " ----------\"\n")
        file.write("echo \"------------------------------\"\n")
        file.write('cd "' + parent_abs + "/" + dir + '"\n')
        if MODE == "BASH":
            file.write("./job.sh\n")
        else:
            file.write("sbatch job.sh\n")
        # file.write('ccc_msub "' + parent_abs + "/" + dir + "/job.sh" + '"\n') # Jolit-Curie
    file.close()
    print("done")

    # make it executable
    if platform == "linux" or platform == "linux2":
        os.system("chmod +x " + parent_abs + "/runScans.sh")

# run the main function
main()

