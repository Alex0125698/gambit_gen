# ----------- instructions -------------

# 1. compile gambit as normal

# 2. copy the gambit executable into files/

# 3. setup the options below to you liking

# 4. set any other options in the yaml files
#    found in files/yaml_files/ but don't touch
#    any line containing '~'

# 5. setup your desired plots in files/plots.pip
#    but don't touch any line containing '~'

# 6. modify files/job.sh so that it works on
#    your HPC but don't touch any line containing '~'

# 7. run this python script

# 8. run the file: 'gens/runScans.sh'

# 9. run the file: 'gens/runPippi.sh'


# ----------- options -------------

# a GAMBIT is generated for all possible combinations below
# and a script is generated for running all simultaneously
# WARNING: be careful not to generate too many, otherwise you will run out of memory

conv_threshold = 1e-6
required_printed_points = 3500000
required_points = -1
required_scan_duration = 90*60

NODE_COUNT = 1  # set to desired number of nodes per gambit
CORE_COUNT = 16 # set to number of cores per node

# allowed options: "THDM", "THDMI", "THDMII", "THDMLS", or "THDMflipped"
models = ["THDMI"]

# allowed options: "tree" or "loop"
runnings = ["tree"]

# allowed options: "generic", "hybrid_lambda_1", "hybrid_lambda_2", "hybrid_Higgs", "higgs", or "physical"

# instead of setting the basis, we set the file. This will allow us to run targetted scans in the same basis.
# not that we cant simply change the basis as the params will be wrong
bases = [#("hybrid_Higgs", "THDMI_hybrid_Higgs"),
         ("hybrid_Higgs", "THDMI_high_cosba"), 
         ("hybrid_Higgs", "THDMI_low_cosba")]

# allowed options: "flat", or "log"
tanb_types = ["flat"]

# allowed options: all, theory, collider, electroweak, flavour
#                 (the name of any individual constraint)
#                 (in the case of perturbativity, or unitarity just name the function)
constraints = [["theory"], 
               ["theory", "collider"], 
               ["theory", "electroweak"], 
               ["theory", "flavour"],
               ["all"]]

# note that all data for [bases,tanb_types] is combined for plotting
# whereas each [models,runnings,constraints] generate different sets of plots


# ----------- main script -------------

import re
import os
import shutil
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
    constraints_theory = ["NLO_unitarity_LogLikelihood_THDM", "LO_unitarity_LogLikelihood_THDM", "stability_LogLikelihood_THDM", "light_scalar_mass_corrections_LogLikelihood_THDM", "heavy_scalar_mass_corrections_LogLikelihood_THDM", "scalar_mass_range_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM", "perturbativity_lambdas_LogLikelihood_THDM", "perturbativity_yukawas_LogLikelihood_THDM"]
    constraints_collider = ["LEP_Higgs_LogLike", "LHC_Higgs_LogLike"]
    constraints_electroweak =  ["oblique_parameters_LogLikelihood_THDM", "lnL_gm2"]
    constraints_flavour =  ["BDstartaunu_LogLikelihood", "BDtaunu_LogLikelihood", "gmu_ge_LogLikelihood", "FLDstar_LogLikelihood", "Bc_lifetime_LogLikelihood", "Bs2llp_LogLikelihood", "B2Kllp_LogLikelihood", "RK_RKstarnunu_LogLikelihood", "h2taumu_LogLikelihood", "t2ch_LogLikelihood", "deltaMB_LogLikelihood", "deltaMBd_LogLikelihood", "SL_LogLikelihood", "l2lgamma_LogLikelihood", "l2lll_LogLikelihood", "RDRDstar_LogLikelihood", "b2sgamma_LogLikelihood", "B2Kstargamma_LogLikelihood", "B2mumu_LogLikelihood_LHCb", "B2mumu_LogLikelihood_CMS", "B2mumu_LogLikelihood_Atlas", "B2KstarmumuAng_LogLikelihood_Atlas", "B2KstarmumuAng_LogLikelihood_CMS", "B2KstarmumuAng_LogLikelihood_LHCb_2020", "B2KstarmumuAng_LogLikelihood_Belle", "B2KstarmumuAng_LogLikelihood_LHCb", "B2KstarellellAng_LogLikelihood_Belle", "Bu2KstarmumuAng_LogLikelihood_LHCb_2020", "B2KstarmumuAng_CPAssym_LogLikelihood_LHCb", "B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020", "B2KstarmumuBr_LogLikelihood_LHCb", "B2KmumuBr_LogLikelihood_LHCb", "Bs2phimumuBr_LogLikelihood"]
    constraints_all = constraints_theory + constraints_collider + constraints_electroweak + constraints_flavour

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

        if basis == "coupling":
            basis = ""
        if running == "loop":
            running = "atQ"
        else:
            running = ""

        self.full_model_name = model + "_" + basis + running

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

        # unitarity 
        if self.running == "tree" and self.NLO_unitarity_LogLikelihood_THDM:
            self.NLO_unitarity_LogLikelihood_THDM = False
            self.LO_unitarity_LogLikelihood_THDM = True

        if self.running == "loop" and self.LO_unitarity_LogLikelihood_THDM:
            self.LO_unitarity_LogLikelihood_THDM = False
            self.NLO_unitarity_LogLikelihood_THDM = True

        # correction checks
        if self.running == "tree":
            self.light_scalar_mass_corrections_LogLikelihood_THDM = False
            self.heavy_scalar_mass_corrections_LogLikelihood_THDM = False

        if  self.perturbativity_LogLikelihood_THDM:
            self.perturbativity_lambdas_LogLikelihood_THDM = False

def makeGambit(options, dir):

    print('making: ' + os.path.abspath("gens") + "/" + dir + " ... ")
    dir2 = "gens/" + dir

    # delete contents if it already exists
    if os.path.exists(dir2):
        shutil.rmtree(dir2)

    # copy the files to a new gambit dir
    copy_tree("files/", dir2)
    # os.rename("gens/files", dir2)

    # figure out the yaml file name
    yaml_name = options.file + ".yaml"

    # patch the yaml file
    patchYaml(options, dir, yaml_name)

    # patch the run script
    patchRunScript(options, dir, yaml_name)

    # create the output folders
    # Path(dir2 + '/../runs/samples/' + options.results_folder).mkdir(parents=True, exist_ok=True)

    print("done")
   
def patchYaml(options, dir, yaml_name):

    # read yaml file into string
    file = open("gens/" + dir + "/yaml_files/" + yaml_name, 'r')
    s = file.read()
    file.close()

    # set basis
    s = s.replace("prior_type: tanb", "prior_type: " + options.tanb_type)

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
    
    # uncomment required constraints
    for c in  options.constraints_all:
        if getattr(options,c):
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
    file = open("gens/" + dir + "/yaml_files/" + yaml_name, 'w')
    file.write(s)
    file.close()

def patchRunScript(options, dir, yaml_name):

    # load run script into string
    file = open("gens/" + dir + "/job.sh", 'r')
    s = file.read()
    file.close()
    
    # set the yaml name
    s = s.replace("THDM_physical.yaml", yaml_name)

    # write patched file to disk
    file = open("gens/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

def patchRunScriptJC(options, dir, yaml_name):

    # load run script into string
    file = open("gens/" + dir + "/job.sh", 'r')
    s = file.read()
    file.close()

    # set the job name
    s = s.replace("#MSUB -r gambit_thdm", "#MSUB -r thdm_" + options.outputName)

    # set the stdout and stderr paths
    dir_abs = os.path.abspath("gens/" + dir)
    # s = s.replace("stdo_", dir + "/work/stdo_")
    # s = s.replace("stde_", dir + "/work/stde_")

    # set gambit dir
    s = s.replace("${CCCSCRATCHDIR}/gambit", dir_abs)

    # set yaml name
    s = s.replace("THDM_physical.yaml", yaml_name)

    # set scan name
    s = s.replace("scan.tar.gz", options.outputName + ".tar.gz")
    s = s.replace("scan.hdf5", options.outputName + ".hdf5")

    # set scan duration
    dur = 0.8 # hours
    scan_percent = 0.93 # use 94% of time for scan
    s = s.replace("9999", str(int(60*60*dur)))
    s = s.replace("8888", str(int(60*60*dur*scan_percent-120)))

    # write patched file to disk
    file = open("gens/" + dir + "/job.sh", 'w')
    file.write(s)
    file.close()

def main():

    counter = 0

    # --- loop over all options ---

    for model in models:
        for running in runnings:
            for constraint in constraints:

                constraint_name = "_set" + str(counter) + "_" + constraint[0]
                counter += 1

                for (basis,file) in bases:
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
                        if required_printed_points is not -1:
                            options.required_printed_points /= CORE_COUNT
                        if required_points is not -1:
                            options.required_points /= CORE_COUNT
                        options.required_scan_duration = required_scan_duration

                        # setup paths
                        options.results_folder = "../../runs/" + model + running + constraint_name + "/"
                        options.plots_folder = "../plots/" + model + running + constraint_name + "/"
                        options.scan_name = file + "_tb"  + tanb

                        # make sure everything is valid
                        options.validate()

                        # make a new gambit with the options specified above
                        makeGambit(options, generate_gambit_name())

    # --- create the runner script --- 

    # write it
    print("generating runner script...")
    parent_abs = os.path.abspath("gens")
    file = open(parent_abs + "/runScans.sh", "w")
    for dir in gambit_dirs:
        file.write('cd "' + parent_abs + "/" + dir + '"\n')
        file.write("./job.sh\n")
        file.write("echo \"------------------------------\"\n")
        # file.write('ccc_msub "' + parent_abs + "/" + dir + "/job.sh" + '"\n') # Jolit-Curie
    file.close()
    print("done")

    # make it executable
    if platform == "linux" or platform == "linux2":
        os.system("chmod +x " + parent_abs + "/runScans.sh")

# run the main function
main()
