import re
import os
import shutil
from distutils.dir_util import copy_tree
from sys import platform

CORE_COUNT = 16

def generate_gambit_name():
    generate_gambit_name.counter += 1
    name = "gambit_" + str(generate_gambit_name.counter)
    gambit_dirs.append(name)
    return name

generate_gambit_name.counter = 0
gambit_dirs = []

class Options:

    def __init__(self):

        self.required_points = -1
        self.required_valid_points = -1
        self.required_scan_duration = -1

        # name of the hdf5 file
        self.outputName = "scan"

        # --- basis ---

        # either: "coupling", "physical", "physical-zoom"
        self.basis = "physical-zoom" 
        self.tanb_type = "flat"

        # --- models ---

        # either: "THDMI", "THDMIatQ", "THDMI_physical", "THDMI_physicalatQ"
        self.model = "THDMI_physical"

        # --- constraints ---

        # theory
        self.negativeMass = True # TODO
        self.stability = False
        self.metaStability = False
        self.unitarity = False
        self.perturbativity_couplings = False
        self.perturbativity_h0_mass = False
        self.perturbativity_scalar_mass = False
        # Collider
        self.higgsBounds = False
        self.higgsSignals = False
        # (todo) XX -> h
        # (todo) h -> YY
        # Precision
        self.obliqueParameters = False
        # Flavor
        self.b2sgamma = False
        self.deltaMB = False
        self.deltaMBd = False
        self.b2ll = False
        self.SL = False
        self.LUV = False

    def allTheory(self):
        self.negativeMass = True
        self.stability = True
        self.metaStability = True
        self.unitarity = True
        self.perturbativity_couplings = True
        if self.isLoopLevel():
            self.perturbativity_h0_mass = True
            self.perturbativity_scalar_mass = True

    def allCollider(self):
        self.higgsBounds = True
        self.higgsSignals = True

    def allPrecision(self):
        self.obliqueParameters = True

    def allFlavour(self):
        self.b2sgamma = True
        self.deltaMB = True
        self.deltaMBd = True
        self.b2ll = True
        self.SL = True
        self.LUV = True

    def isLoopLevel(self):
        return self.model == "THDMIatQ" or self.model == "THDMI_physicalatQ"

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
    if options.basis == "coupling":
        yaml_name = "THDM_coupling.yaml"
    if options.basis == "physical":
        yaml_name = "THDM_physical.yaml"
    if options.basis == "physical-zoom":
        yaml_name = "THDM_physical_zoom.yaml"

    # patch the yaml file
    patchYaml(options, dir, yaml_name)

    # patch the run script
    patchRunScript(options, dir, yaml_name)

    print("done")
   
def patchYaml(options, dir, yaml_name):

    # read yaml file into string
    file = open("gens/" + dir + "/yaml_files/" + yaml_name, 'r')
    s = file.read()
    file.close()

    # set basis
    s = s.replace("prior_type: tanb", "prior_type: " + options.tanb_type)

    # set the model name
    s = s.replace("TheModelName", options.model)

    # set the scan duration and point limit
    s = s.replace("12121212", str(int(options.required_points / CORE_COUNT)))
    s = s.replace("23232323", str(int(options.required_valid_points / CORE_COUNT)))
    s = s.replace("34343434", str(int(options.required_scan_duration)))

    # set the output hdf5 name
    s = s.replace("scan.hdf5", options.outputName + ".hdf5")

    # add parameters for loop-level models
    if options.isLoopLevel():
        s = s.replace("#~0","")

    # set the constraints
    if options.negativeMass:
        pass
    if options.stability:
        s = s.replace("#~1","")
    if options.metaStability:
        pass
    if options.unitarity:
        if options.isLoopLevel():
            s = s.replace("#~2","")
        else:
            s = s.replace("#~F", "")
    if options.perturbativity_couplings:
        s = s.replace("#~3","")
    if options.perturbativity_h0_mass:
        s = s.replace("#~4","")
    if options.perturbativity_scalar_mass:
        s = s.replace("#~5","")
    if options.higgsBounds:
        s = s.replace("#~6","")
    if options.higgsSignals:
        s = s.replace("#~7","")
    if options.obliqueParameters:
        s = s.replace("#~8","")
    if options.b2sgamma:
        s = s.replace("#~9","")
    if options.deltaMB:
        s = s.replace("#~A","")
    if options.deltaMBd:
        s = s.replace("#~B","")
    if options.b2ll:
        s = s.replace("#~C","")
    if options.SL:
        s = s.replace("#~D","")
    if options.LUV:
        s = s.replace("#~E","")

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

    runnings = ["tree", "loop"]
    bases = ["coupling", "physical", "physical-zoom"]
    tanb_types = ["log", "flat"]

    # loop over all bases + runnings
    for running in runnings:
        for basis in bases:
            for tanb_type in tanb_types:

                # figure out the model name
                if running == "loop":
                    if basis == "coupling":
                        model = "THDMIatQ"
                    else:
                        model = "THDMI_physicalatQ"
                else:
                    if basis == "coupling":
                        model = "THDMI"
                    else:
                        model = "THDMI_physical"

                # --- generate scans for each set of constraints
                constraint_types = [Options.allTheory, Options.allPrecision, Options.allFlavour, Options.allCollider]
                constraint_names = ["theory", "precision", "flavour", "collider"]

                for i in range(0, len(constraint_types)):

                    fcn = constraint_types[i]
                    fcn_name = constraint_names[i]

                    # setup the scan options
                    options = Options()

                    options.required_points = 10000000
                    options.required_valid_points = 200000
                    options.required_scan_duration = 3*60*60
                    if tanb_types == "log" and ( bases == "coupling" or bases == "physical-zoom" ):
                        options.required_scan_duration = 6*60*60

                    # for our preliminary scan...
                    options.required_scan_duration /= 3

                    outputFolder = running + "_set_" + constraint_names[i] + "/"
                    options.outputName = outputFolder + "scan_" + basis + "_" + tanb_type
                    options.model = model
                    options.basis = basis
                    options.tanb_type = tanb_type
                    fcn(options)

                    # make a new gambit with the options specified above
                    makeGambit(options, generate_gambit_name())


                # --- generate scans for individual theory constraints

                theory_names = ["negativeMass", "stability", "metaStability", "unitarity", "perturbativity_couplings", "perturbativity_h0_mass", "perturbativity_scalar_mass"]

                for i in range(0, len(theory_names)):

                    theory = theory_names[i]

                    # setup the scan options
                    options = Options()

                    options.required_points = 10000000
                    options.required_valid_points = 200000
                    options.required_scan_duration = 3*60*60
                    if tanb_types == "log" and ( bases == "coupling" or bases == "physical-zoom" ):
                        options.required_scan_duration = 6*60*60

                    # for our preliminary scan...
                    options.required_scan_duration /= 3

                    outputFolder = running + "_th_" + theory_names[i] + "/"
                    options.outputName = outputFolder + "scan_" + basis + "_" + tanb_type
                    options.model = model
                    options.basis = basis
                    options.tanb_type = tanb_type
                    setattr(options,theory,True)

                    # make a new gambit with the options specified above
                    makeGambit(options, generate_gambit_name())


                # generate scans for individual collider channels (XX -> h)

                # generate scans for individual collider channels (h -> YY)

                # generate scans for individual flavour constraints



    # --- create the runner script --- 

    # write it
    print("generating runner script...")
    parent_abs = os.path.abspath("gens")
    file = open(parent_abs + "/runScans.sh", "w")
    for dir in gambit_dirs:
        file.write('cd "' + parent_abs + "/" + dir + '"\n')
        file.write("./job.sh\n")
        # file.write('ccc_msub "' + parent_abs + "/" + dir + "/job.sh" + '"\n') # Jolit-Curie
    file.close()
    print("done")

    # make it executable
    if platform == "linux" or platform == "linux2":
        os.system("chmod +x " + parent_abs + "/runScans.sh")

# run the main function
main()