#!/bin/python3

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
from system import *
from input import *

# check that system is setup as expected
if not os.path.isdir(ROOT_DIR):
  print(f"File tree not setup correctly")

# get actual PP file location; make sure it exists
if not os.path.isfile(POSTPROCESSOR_FILE):
  POSTPROCESSOR_FILE = f"{ROOT_DIR}/data/scansPP/{POSTPROCESSOR_FILE}"
for b in BASES:
  if "pp" == b[1]:
    for model in MODELS:
      for running in SPECTRUM_ORDER:
        file = POSTPROCESSOR_FILE.replace("TheModelPP", model+running)
        if not os.path.isfile(file):
          raise Exception(f"PP file: {file} does not exist")

# add the additional constraints and output prefix
CONSTRAINTS = [(c[0]+ADDITIONAL_CONSTRAINTS, OUTPUT_PREFIX + c[1]) for c in CONSTRAINTS] 

if THE_SYSTEM == "dirac_icelake":
  NODE_COUNT = 1  # set to desired number of nodes per gambit
  CORE_COUNT = 76 # set to number of cores per node
elif THE_SYSTEM == "UbuntuDesktop":
  NODE_COUNT = 1  # set to desired number of nodes per gambit
  CORE_COUNT = 16 # set to number of cores per node
elif THE_SYSTEM == "UbuntuLegion":
  NODE_COUNT = 1  # set to desired number of nodes per gambit
  CORE_COUNT = 8 # set to number of cores per node
  
# path for generated scans
GEN_PATH = 'gens'

print("\n-----------------\nrunning 2HDM generator\n")

len_bases = len(BASES)

for (basis,file) in BASES:
    name = "files/" + BASIS_YAML_DIRECTORY + "/" + file + ".yaml"
    file = open(name,'r').read()
    file = re.sub(r"^!import.*",r"",file,flags=re.MULTILINE)
    file = re.sub(r"!import",r"",file)
    yfile = yaml.safe_load(file)
    if "subscans" in yfile:
      len_bases += yfile["subscans"]["num_scans"] - 1

print("\nthese scans will eat ", len(MODELS)*len(SPECTRUM_ORDER)*len_bases*len(CONSTRAINTS)*MAX_SCAN_DURATION*CORE_COUNT/(60*60), " CPU hours\n\n")


def generate_gambit_name():
    generate_gambit_name.counter += 1
    name = "gambit_" + str(generate_gambit_name.counter)
    gambit_dirs.append(name)
    return name

generate_gambit_name.counter = 0
gambit_dirs = []

class Options:

    # the names of all constraints
    constraints_theory = ["VS_likelihood", "runToScaleTest_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM", "LO_unitarity_LogLikelihood_THDM", "stability_LogLikelihood_THDM", "higgs_exp_mass_LogLikelihood_THDM", 
                         "scalar_mass_corrections_LogLikelihood_THDM", "higgs_scenario_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM", 
                         "perturbativity_yukawas_LogLikelihood_THDM"]
    constraints_collider = ["LEP_Higgs_LogLike", "LHC_Higgs_LogLike","HS_ALL","HS_RUN1_SS","HS_LATEST_SS","HS_LATEST_STXS"] #higgs_mass_LogLikelihood
    constraints_electroweak =  ["oblique_parameters_LogLikelihood"] # lnL_gm2
    constraints_flavour =  ["b2sgamma_LogLikelihood","B2Kstargamma_LogLikelihood","B2mumu_LogLikelihood_Atlas","B2mumu_LogLikelihood_LHCb","B2mumu_LogLikelihood_CMS","B2mumu_LogLikelihood_CMS_ATLAS_LHCb","Bd2KmumuBr_LogLikelihood_LHCb","Bd2KmumuBr_LogLikelihood_Belle","B2KmumuBr_LogLikelihood_LHCb","B2KmumuBr_LogLikelihood_CMS","B2KmumuBr_LogLikelihood_Belle","B2KeeBr_LogLikelihood_Belle","Bd2KeeBr_LogLikelihood_Belle","B2KstarmumuAng_LogLikelihood_Atlas","B2KstarmumuAng_LogLikelihood_CMS","B2KstarmumuAng_LogLikelihood_Belle","B2KstarmumuAng_LogLikelihood_LHCb_2020","Bu2KstarmumuAng_LogLikelihood_LHCb_2020","B2KstarmumuBr_LogLikelihood_LHCb","Bs2phimumuBr_LogLikelihood","B2KstarmumuAng_CPAssym_LogLikelihood_LHCb","B2KstarellellAng_LogLikelihood_Belle","B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020","RKRKstar_LogLikelihood_LHCb","RK_LogLikelihood_CMS","RK_LogLikelihood_Belle","BKnunu_LogLikelihood_Belle_sl","BKnunu_LogLikelihood_Belle_had","BuKnunu_LogLikelihood_Belle_sl","BuKnunu_LogLikelihood_Belle_had","BuKnunu_LogLikelihood_BelleII","BKnunu_LogLikelihood_BaBar","BuKnunu_LogLikelihood_BaBar","BKstarnunu_LogLikelihood_Belle_sl","BKstarnunu_LogLikelihood_Belle_had","BuKstarnunu_LogLikelihood_Belle_sl","BuKstarnunu_LogLikelihood_Belle_had","BKstarnunu_LogLikelihood_BaBar","BuKstarnunu_LogLikelihood_BaBar","SL_LogLikelihood","FLDstar_LogLikelihood","dBRBDstartaunu_LogLikelihood","dBRBDtaunu_LogLikelihood","Bc_lifetime_LogLikelihood","Delta_MBs_LogLikelihood","Delta_MBd_LogLikelihood","l2lgamma_LogLikelihood","l2lll_LogLikelihood","h2ltau_LogLikelihood","t2ch_LogLikelihood","t2bbc_LogLikelihood","t2mutauc_LogLikelihood","Bc2taunu_LogLikelihood","Bs2ll_LogLikelihood","B2Kll_LogLikelihood","B2Xsnunu_LogLikelihood","gmu_ge_LogLikelihood"]
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
        self.results_folder = None
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

        self.higgs_scenario_LogLikelihood_THDM = True
        self.higgs_exp_mass_LogLikelihood_THDM = True

        # unitarity 
        if self.running == "tree" and self.NLO_unitarity_LogLikelihood_THDM:
            self.NLO_unitarity_LogLikelihood_THDM = False
            self.LO_unitarity_LogLikelihood_THDM = True

        if self.running == "loop" and self.LO_unitarity_LogLikelihood_THDM and self.NLO_unitarity_LogLikelihood_THDM:
            self.LO_unitarity_LogLikelihood_THDM = False
            self.NLO_unitarity_LogLikelihood_THDM = True

        if self.B2mumu_LogLikelihood_CMS_ATLAS_LHCb:
            self.B2mumu_LogLikelihood_CMS = False
            self.B2mumu_LogLikelihood_ATLAS = False
            self.B2mumu_LogLikelihood_LHCb = False

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

    print('making: ' + os.path.abspath(GEN_PATH) + "/" + dir + " ... ")
    dir2 = GEN_PATH+"/" + dir

    # delete contents if it already exists
    if os.path.exists(dir2):
        shutil.rmtree(dir2)

    # copy the files to a new gambit dir
    copy_tree("files/copyme/", dir2)
    # os.rename(GEN_PATH+"/files", dir2)

    # figure out the yaml file name
    yaml_name = options.file + ".yaml"

    # patch the yaml file
    patchYaml(options, dir, yaml_name)
    patchYaml(options, dir, "THDM_constraints.yaml")
    # remove_tree(dir2 + "/" + BASIS_YAML_DIRECTORY)

    # patch the run script
    if THE_SYSTEM in ["UbuntuDesktop", "UbuntuLegion"]:
        patchRunScript(options, dir, yaml_name)
    elif THE_SYSTEM == "dirac_icelake":
        patchRunScriptJC(options, dir, yaml_name)
    else:
        raise Exception("unknown mode")
    
    # os.remove(dir2 + "/job_jc.sh")
    os.remove(dir2 + "/job_dirac.sh")

    # create the output folders (otherwise hdf5_v1 will crash)
    Path(options.results_folder + '/samples').mkdir(parents=True, exist_ok=True)

    print("done")
   
def patchYaml(options, dir, yaml_name):

    #  hack
    if yaml_name == "THDM_constraints.yaml":


        # read yaml file into string
        file = open("files/" + BASIS_YAML_DIRECTORY + "/" + yaml_name, 'r')
        s = file.read()
        file.close()

    else:
        s = options.subscan

    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_med_hhs")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_med_final_hhs")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_full")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_med")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_med_final")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_small")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_small_final")
    # shutil.rmtree(GEN_PATH+"/" + dir + "/" + "yaml_files_idm")

    # remove useless yaml dirs

    print("DEBUG: patching " + GEN_PATH+"/" + dir + "/" + BASIS_YAML_DIRECTORY + "/" + yaml_name)

    # # set basis
    # s = s.replace("prior_type: tanb", "prior_type: " + "flat")

    # set the model name
    s = s.replace("TheModelName", options.full_model_name)

    # enable postprocessor scanner for PP basis
    if options.file == "pp":
        s = s.replace("use_scanner:", "use_scanner: postprocessor # ")

    # replace model name for PP
    s = s.replace("ppfile:", "file: " + '"' + POSTPROCESSOR_FILE + '"' + " # ")
    s = s.replace("TheModelPP", options.model+options.running)

    # set the scan duration and point limit
    if options.required_points: s = s.replace("12121212", str(int(options.required_points)))
    if options.required_printed_points: s = s.replace("23232323", str(int(options.required_printed_points)))
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
    if USE_SPEED_HACKS:
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
    file = open(GEN_PATH+"/" + dir + "/yaml_files/" + yaml_name, 'w')
    file.write(s)
    file.close()

def patchRunScript(options, dir, yaml_name):

    # load run script into string
    file = open(GEN_PATH+"/" + dir + "/job.sh", 'r')
    s = file.read()
    file.close()
    
    # set the yaml name
    s = s.replace("THDM_physical.yaml", yaml_name)

    # write patched file to disk
    file = open(GEN_PATH+"/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

def patchRunScriptJC(options, dir, yaml_name):

    # load run script into string
    file = open(GEN_PATH+"/" + dir + "/job_dirac.sh", 'r')
    s = file.read()
    file.close()

    # set the job name
    s = s.replace("SBATCH -J gambit_thdm", "SBATCH -J thdm_" + options.hdf5_name)

    # set the stdout and stderr paths (todo)
    dir_abs = os.path.abspath(GEN_PATH+"/" + dir)
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
    file = open(GEN_PATH+"/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

# load a yaml file with (optionally) a subscan node and convert to string
def load_subscan_yaml(running, name):

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

                if param_priors[param] == "none":
                    param_weights[param] = 0
                    param_ranges[param] = tmp[param]["prior_type"]
                elif "range" in tmp[param]:
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
            if (param_priors[param]) == "none": continue
            
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
            if (param_priors[param]) == "pow":
              yfile["Parameters"]["TheModelName"][param]["shift"] = -param_range_new[0]
              # yfile["Parameters"]["TheModelName"][param]["shift"] = -max(0.0,param_range[0]-0.1)

        # convert yaml file to string
        contents = yaml.dump(yfile)

        # # add Qin back
        if running == "loop":
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

    # ad postfix to GEN_PATH so that it is unique
    global GEN_PATH
    postfix = 0
    while os.path.exists(GEN_PATH+"_"+str(postfix)):
        postfix += 1
    postfix = str(postfix)
    GEN_PATH = GEN_PATH + "_" + postfix


    folders_to_merge = {}

    # --- loop over all options ---

    # loop over all models (different output folder)
    for model in MODELS:

        # loop over all runnings (different output folder)
        for running in SPECTRUM_ORDER:

            # loop over all constraints (different output folder)
            for (constraint, constraint_name) in CONSTRAINTS:

                # uniquely identifies results folders (to be merged later)
                resultsSuffix = 0

                # loop over all bases (merged output folder)
                for (basis,file) in BASES:

                    # get the list of subscan files
                    subscans = load_subscan_yaml(running, "files/" + BASIS_YAML_DIRECTORY + "/" + file + ".yaml")

                    # loop over all subscans (merged output folder)
                    for x,subscan in enumerate(subscans):

                        resultsSuffix += 1
                        resultsSuffixStr = "_" + str(resultsSuffix) if THE_SYSTEM == "dirac_icelake" else ""

                        # setup the scan-specific options
                        options = Options()
                        options.setModel(model, basis, running)
                        options.subscan = subscan
                        options.file = file
                        if len(subscans) > 1:
                            options.file += "x"+str(x)
                        for c in constraint:
                            options.setConstraint(c)

                        # setup convergence criteria
                        options.conv_threshold = DIVER_CONV_THRESHOLD
                        options.NP = DIVER_NP
                        options.required_printed_points = MAX_PRINTED_POINTS
                        options.required_points = MAX_POINTS
                        if MAX_PRINTED_POINTS != None:
                            options.required_printed_points /= (len_bases*CORE_COUNT)
                        if MAX_POINTS != None:
                            options.required_points /= (len_bases*CORE_COUNT)
                        options.required_scan_duration = MAX_SCAN_DURATION

                        # make sure everything is valid
                        options.validate()

                        # setup paths (also stored in options)
                        fullName = model + running + "_" + constraint_name
                        options.results_folder = f"{ROOT_DIR}/data/scans/{fullName}{resultsSuffixStr}/"
                        options.hdf5_name = options.file + "_" + postfix + "_" + resultsSuffixStr
                        
                        # dict that tells us which folders to merge (essentially we will just get rid of resultsSuffixStr)
                        tmp = "{ROOT_DIR}/data/scans/{fullName}{resultsSuffixStr}/"
                        folders_to_merge[tmp] = "{ROOT_DIR}{fullName}/"

                        # make a new gambit with the options specified above
                        makeGambit(options, generate_gambit_name())

    parent_abs = os.path.abspath(GEN_PATH)
    if THE_SYSTEM == "dirac_icelake":
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
        if THE_SYSTEM in ["UbuntuDesktop", 'UbuntuLegion']:
            file.write("./job.sh\n")
        else:
            file.write("sbatch job.sh\n")
        # file.write('ccc_msub "' + parent_abs + "/" + dir + "/job.sh" + '"\n') # Jolit-Curie
    file.write("rm ../runScans.sh")
    file.close()
    print("done")

    # make it executable
    if platform == "linux" or platform == "linux2":
        os.system("chmod +x " + parent_abs + "/runScans.sh")

# run the main function
main()

