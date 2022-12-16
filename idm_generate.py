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

# a GAMBIT is generated for all possible combinations below
# and a script is generated for running all simultaneously
# WARNING: be careful not to generate too many, otherwise you will run out of storage

yaml_dir = 'yaml_files_idm'
use_speed_hacks = False
postfix = '_1'

gen_path = 'gens' + postfix

conv_threshold = 1e-6
required_printed_points = 400000
required_points = -1
required_scan_duration = 10*60/60 # in seconds

NODE_COUNT = 1  # set to desired number of nodes per gambit
CORE_COUNT = 72 # set to number of cores per node

# either "DIRAC" or "BASH"
MODE = "DIRAC"

# allowed options: "THDM", "THDMI", "THDMII", "THDMLS", or "THDMflipped"
models = ["Inert2"]

# allowed options: "tree" or "loop"
runnings = ["tree"]

# allowed options: "generic", "hybrid_lambda_1", "hybrid_lambda_2", "hybrid_Higgs", "higgs", or "physical"

# instead of setting the basis, we set the file. This will allow us to run targetted scans in the same basis.
# not that we cant simply change the basis as the params will be wrong

# 20
bases = [

    ("physical", "physical_high_lam3456"),
    ("physical", "physical_high_lam3458"),
    ("physical", "physical_high_lam34510"),
    ("physical", "physical_logmH"),
    ("physical", "physical_low_lam345"),
    ("physical", "physical_low_mH140"),
    ("physical", "physical_low_mH400"),
    ("physical", "physical"),

]


# allowed options: "flat", or "log"
tanb_types = ["flat"]

# allowed options: all, theory, collider, electroweak, flavour
#                 (the name of any individual constraint)
#                 (in the case of perturbativity, or unitarity just name the function)


# 39-2

# ---- COMBINED SCANS ----

# 1
constraints = [

    (["all"], "all"),

]

# ---- THEORY SCANS ----

# 5
# constraints = [

#     (["theory"], "theory"),
#     (["check_unitarity"], "LO"),
#     (["check_vacuum_stability"], "stability"),
#     (["check_perturbativity"], "perturbativity"),
#     (["check_charged_vacuum"], "charged_vacuum"),
#     (["check_vacuum_metastability"], "vacuum_metastability"),

# ]

# ---- ELECTROWEAK SCANS ----

# 1
# constraints = [

#     (["theory", "electroweak"], "electroweak"),

# ]

# ---- COLLIDER SCANS ----

# 6
# constraints = [

#     (["theory", "collider"], "collider"),
#     (["theory", "LEP_Higgs_LogLike"], "HB"),
#     (["theory", "LHC_Higgs_LogLike", "HS_ALL"], "HS"),
#     (["theory", "LHC_Higgs_LogLike", "HS_RUN1_SS"], "HSRUN1SS"),
#     (["theory", "LHC_Higgs_LogLike", "HS_LATEST_SS"], "HSLATESTSS"),
#     (["theory", "LHC_Higgs_LogLike", "HS_LATEST_STXS"], "HSLATESTSTXS"),

# ]


# ---- DARK MATTER SCANS ----

# constraints = [

#     # (["theory", "dark_matter"], "dark_matter"),
#     # (["theory", "XENON1T_2018_LogLikelihood", "XENON1T_2017_LogLikelihood", "XENON100_2012_LogLikelihood", "DARWIN_LogLikelihood", "LUX_2013_LogLikelihood", "LUX_2015_LogLikelihood", "LUX_2016_LogLikelihood", "LZ_LogLikelihood", "PandaX_2016_LogLikelihood", "PandaX_2017_LogLikelihood", "DarkSide_50_LogLikelihood", "DarkSide_50_S2_LogLikelihood", "CRESST_II_LogLikelihood", "CRESST_III_LogLikelihood", "SuperCDMS_2014_LogLikelihood", "CDMSlite_LogLikelihood", "SIMPLE_2014_LogLikelihood", "PICO_2L_LogLikelihood", "PICO_60_LogLikelihood", "PICO_60_2017_LogLikelihood", "PICO_60_2019_LogLikelihood", "PICO_500_LogLikelihood"], "direct_detection"),
#     # (["theory", "lnL_FermiLATdwarfs", "lnL_FermiGC", "lnL_CTAGC", "lnL_HESSGC"], "indirect_detection"), # excluding neutrinos
#     # (["theory", "lnL_oh2"], "lnL_oh2"),
#     # (["theory", "lnL_FermiLATdwarfs"], "lnL_FermiLATdwarfs"),
#     # (["theory", "lnL_FermiGC", "lnL_CTAGC", "lnL_HESSGC"], "galactic_center"),
#     # (["theory", "XENON1T_2018_LogLikelihood", "XENON1T_2017_LogLikelihood", "XENON100_2012_LogLikelihood"], "XENON"),
#     # (["theory", "DARWIN_LogLikelihood"], "DARWIN"),
#     # (["theory", "LUX_2013_LogLikelihood", "LUX_2015_LogLikelihood", "LUX_2016_LogLikelihood"], "LUX"),
#     # (["theory", "LZ_LogLikelihood"], "LZ"),
#     # (["theory", "PandaX_2016_LogLikelihood", "PandaX_2017_LogLikelihood"], "PandaX"),
#     # (["theory", "DarkSide_50_LogLikelihood", "DarkSide_50_S2_LogLikelihood"], "DarkSide_50"),
#     # (["theory", "CRESST_II_LogLikelihood", "CRESST_III_LogLikelihood"], "CRESST"),
#     # (["theory", "SuperCDMS_2014_LogLikelihood", "CDMSlite_LogLikelihood"], "CDMS"),
#     # (["theory", "SIMPLE_2014_LogLikelihood"], "SIMPLE"),
#     # (["theory", "PICO_2L_LogLikelihood", "PICO_60_LogLikelihood", "PICO_60_2017_LogLikelihood", "PICO_60_2019_LogLikelihood", "PICO_500_LogLikelihood"], "PICO"),

#     # NOT WORKING
#     # (["theory", "IceCube_likelihood"], "IceCube_likelihood"),

# ]

# ---- FLAVOR SCANS ----

# NONE

# note that all data for [bases,tanb_types] is combined for plotting
# whereas each [models,runnings,constraints] generate different sets of plots


# ----------- main script -------------

import re
import os
import shutil
import math
from distutils.dir_util import copy_tree
from sys import platform
from pathlib import Path


def generate_gambit_name():
    generate_gambit_name.counter += 1
    name = "gambit_" + str(generate_gambit_name.counter)
    gambit_dirs.append(name)
    return name

generate_gambit_name.counter = 0
gambit_dirs = []

class Options:

    # the names of all constraints
    constraints_theory = ["check_unitarity", "check_vacuum_stability", "check_perturbativity","check_charged_vacuum","check_vacuum_metastability"]
    constraints_collider = ["LEP_Higgs_LogLike", "LHC_Higgs_LogLike","HS_ALL","HS_RUN1_SS","HS_LATEST_SS","HS_LATEST_STXS"] #higgs_mass_LogLikelihood
    constraints_electroweak =  ["oblique_parameters_LogLikelihood"] # lnL_gm2
    constraints_dark_matter = ["XENON1T_2018_LogLikelihood", "XENON1T_2017_LogLikelihood", "XENON100_2012_LogLikelihood", "DARWIN_LogLikelihood", 
                               "LUX_2013_LogLikelihood", "LUX_2015_LogLikelihood", "LUX_2016_LogLikelihood", "LZ_LogLikelihood", "PandaX_2016_LogLikelihood",
                               "PandaX_2017_LogLikelihood", "DarkSide_50_LogLikelihood", "DarkSide_50_S2_LogLikelihood", "CRESST_II_LogLikelihood", 
                               "CRESST_III_LogLikelihood", "SuperCDMS_2014_LogLikelihood", "CDMSlite_LogLikelihood", "SIMPLE_2014_LogLikelihood", 
                               "PICO_2L_LogLikelihood", "PICO_60_LogLikelihood", "PICO_60_2017_LogLikelihood", "PICO_60_2019_LogLikelihood", 
                               "PICO_500_LogLikelihood","lnL_FermiLATdwarfs", "lnL_FermiGC", "lnL_CTAGC", "lnL_HESSGC", "lnL_oh2", "IceCube_likelihood"]
    constraints_flavour =  ["Bs2ll_LogLikelihood","B2Kll_LogLikelihood","B2mumu_LogLikelihood_Atlas","B2mumu_LogLikelihood_LHCb","B2mumu_LogLikelihood_CMS",
                            "dBRBDstartaunu_LogLikelihood","dBRBDtaunu_LogLikelihood","gmu_ge_LogLikelihood","FLDstar_LogLikelihood","Bc_lifetime_LogLikelihood",
                            "B2Xsnunu","h2taumu_LogLikelihood","t2ch_LogLikelihood","deltaMB_LogLikelihood","deltaMBd_LogLikelihood",
                            "SL_LogLikelihood","l2lgamma_LogLikelihood","l2lll_LogLikelihood","b2sgamma_LogLikelihood","B2Kstargamma_LogLikelihood",
                            "RK_LogLikelihood_LHCb","RKstar_LogLikelihood_LHCb","B2KstarmumuAng_LogLikelihood_Atlas","B2KstarmumuAng_LogLikelihood_CMS",
                            "B2KstarmumuAng_LogLikelihood_LHCb_2020","B2KstarmumuAng_LogLikelihood_Belle","B2KstarellellAng_LogLikelihood_Belle",
                            "Bu2KstarmumuAng_LogLikelihood_LHCb_2020","B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020","B2KstarmumuBr_LogLikelihood_LHCb",
                            "B2KmumuBr_LogLikelihood_LHCb","Bs2phimumuBr_LogLikelihood"]
    constraints_trivial = ["lnL_v0","lnL_SI_nuclear_parameters",]

    constraints_all = constraints_theory + constraints_collider + constraints_electroweak + constraints_dark_matter + constraints_trivial

    def __init__(self):

        # default scan options
        self.tanb_type = None
        self.full_model_name = None
        self.model = None
        self.running = None
        self.basis = None

        # default convergence criteria
        self.conv_threshold = 1e-5
        self.required_printed_points = -1
        self.required_points = -1
        self.required_scan_duration = 7*24*60*60

        # set default constraints
        for c in self.constraints_all:
            setattr(self, c, False)

        # default paths
        self.results_folder = "../runs"
        self.plots_folder = "../plots"
        self.scan_name = "scan"

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
        elif name == "dark_matter":
            setit = self.constraints_dark_matter
        else:
            setit = [name]
            if not hasattr(self, name):
                raise Exception("Error: no constraint called: " + name)

        for c in setit:
            setattr(self, c, True)

    def validate(self):

        # trivial likelihoods

        DD = ["XENON1T_2018_LogLikelihood", "XENON1T_2017_LogLikelihood", "XENON100_2012_LogLikelihood", "DARWIN_LogLikelihood", "LUX_2013_LogLikelihood", "LUX_2015_LogLikelihood", "LUX_2016_LogLikelihood", "LZ_LogLikelihood", "PandaX_2016_LogLikelihood", "PandaX_2017_LogLikelihood", "DarkSide_50_LogLikelihood", "DarkSide_50_S2_LogLikelihood", "CRESST_II_LogLikelihood", "CRESST_III_LogLikelihood", "SuperCDMS_2014_LogLikelihood", "CDMSlite_LogLikelihood", "SIMPLE_2014_LogLikelihood", "PICO_2L_LogLikelihood", "PICO_60_LogLikelihood", "PICO_60_2017_LogLikelihood", "PICO_60_2019_LogLikelihood", "PICO_500_LogLikelihood"]
        ID = ["lnL_FermiGC", "lnL_CTAGC", "lnL_HESSGC"]

        for c in DD:
            if getattr(self,c):
                self.lnL_SI_nuclear_parameters = True
                self.lnL_v0 = True
                break

        for c in ID:
            if getattr(self,c):
                self.lnL_v0 = True
                break

        # not ready likelihoods
        self.lnL_CTAGC = False
        self.LZ_LogLikelihood = False

        # DOES NOTHING likelihoods
        self.dBRBDstartaunu_LogLikelihood = False
        self.dBRBDtaunu_LogLikelihood = False
        self.gmu_ge_LogLikelihood = False
        self.FLDstar_LogLikelihood = False
        self.h2taumu_LogLikelihood = False
        self.t2ch_LogLikelihood = False
        self.l2lgamma_LogLikelihood = False
        self.l2lll_LogLikelihood = False

        # NOT WORKING likelihoods
        self.deltaMB_LogLikelihood = False
        self.deltaMBd_LogLikelihood = False
        self.IceCube_likelihood = False

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

    # patch the run script
    if MODE == "BASH":
        patchRunScript(options, dir, yaml_name)
    elif MODE == "DIRAC":
        patchRunScriptJC(options, dir, yaml_name)
    else:
        raise Exception("unknown mode")

    # create the output folders (otherwise hdf5_v1 will crash)
    Path(dir2 + '/yaml_files/' + options.results_folder + '/samples').mkdir(parents=True, exist_ok=True)

    print("done")
   
def patchYaml(options, dir, yaml_name):

    # read yaml file into string
    file = open(gen_path+"/" + dir + "/" + yaml_dir + "/" + yaml_name, 'r')
    s = file.read()
    file.close()

    shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_full")
    shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_med")
    shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_small")
    shutil.rmtree(gen_path+"/" + dir + "/" + "yaml_files_idm")

    # remove useless yaml dirs

    print("DEBUG: patching " + gen_path+"/" + dir + "/" + yaml_dir + "/" + yaml_name)

    # set basis
    s = s.replace("prior_Type: tanb", "prior_type: " + options.tanb_type)

    # set the model name
    s = s.replace("TheModelName", options.full_model_name)

    # set the scan duration and point limit
    s = s.replace("12121212", str(int(options.required_points)))
    s = s.replace("23232323", str(int(options.required_printed_points)))
    s = s.replace("34343434", str(int(options.required_scan_duration)))
    s = s.replace("convthresh:", "convthresh: " + str(options.conv_threshold) + " #")

    # set the output hdf5 name
    s = s.replace("scan.hdf5", options.scan_name + ".hdf5")

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
        s = s.replace("#~o","")
    else:
        s = s.replace("#~o","use_speedhacks: false # ")


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
            if i == -1:
                i = s.find(c)
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
    s = s.replace("SBATCH -J gambit_thdm", "SBATCH -J thdm_" + options.scan_name)

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

def main():

    folders_to_merge = {}

    # --- loop over all options ---

    for model in models:
        for running in runnings:
            for (constraint, constraint_name) in constraints:
                counter = 0
                for (basis,file) in bases:
                    counter += 1
                    for tanb in tanb_types:

                        options = Options()

                        # setup the scan options
                        options.setModel(model, basis, running)
                        options.file = file
                        options.tanb_type = tanb
                        for c in constraint:
                            options.setConstraint(c)

                        # setup convergence criteria
                        options.conv_threshold = conv_threshold
                        options.required_printed_points = required_printed_points
                        options.required_points = required_points
                        if required_printed_points != -1:
                            options.required_printed_points /= CORE_COUNT
                        if required_points != -1:
                            options.required_points /= CORE_COUNT
                        options.required_scan_duration = required_scan_duration

                        # make sure everything is valid
                        options.validate()

                        scanner_suffix = "_" + str(counter) if MODE != "BASH" else ""

                        # setup paths
                        constraint_name_full = "_" + constraint_name
                        options.results_folder = "../../../runs/" + model + running + constraint_name_full + scanner_suffix + "/"
                        options.plots_folder = "../plots/" + model + running + constraint_name_full + "/"
                        options.scan_name = file + "_tb"  + tanb + postfix

                        tmp = "../runs/" + model + running + constraint_name_full + scanner_suffix + "/"
                        folders_to_merge[tmp] = "../runs/" + model + running + constraint_name_full + "/"

                        # for i in range(0,CORE_COUNT):
                        #     folders_to_merge[options.results_folder + 'samples/' + options.scan_name + ".hdf5_temp_" + str(i)] = "../../../runs/" + model + running + constraint_name_full + "/"

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


