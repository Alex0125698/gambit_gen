# Check some theory LLs prior to spectrum generation
USE_SPEED_HACKS = False

# Convergence options (used for ALL scans)
DIVER_CONV_THRESHOLD = 1e-8
DIVER_NP = 10000
MAX_PRINTED_POINTS = None
MAX_POINTS = None
MAX_SCAN_DURATION = 1*60*60 # in seconds; REQUIRED

# The models to generate scans for ("THDM", "THDMI", "THDMII", "THDMLS", or "THDMflipped")
MODELS = ["THDMII"]

# The spectrum orders to generate scans for ("tree" or "loop")
SPECTRUM_ORDER = ["loop"]

# Prefix for output directory; needed if constraint name is not enough to differentiate different scans
OUTPUT_PREFIX = ""

# Use this postprocessor file for the pp basis (see below); ignored for other bases
POSTPROCESSOR_FILE = "PP_2025_02_26_TheModelPP_theory_40000p.hdf5"

# The directory of the YAML files (including the bases and constraints files)
BASIS_YAML_DIRECTORY = 'yaml_files_full'

# The basis/prior combinations to use
BASES = [
  
    # ---- post-processor

    ("generic", "pp"),

    # ---- full basis

    # ("physical", "physical"),       # a standard physical basis
    # ("generic", "generic"),         # b standard generic basis
    # ("hybrid_Higgs", "hybrid1"),    # c standard hybrid1 basis
    # ("hybrid_Higgs2", "hybrid2"),   # d standard hybrid2 basis
    
    # ---- grids
    
    # ("physical", "physicalA"),       # s tanb grid (log)

    # ("generic", "genericA"),        # e tanb grid (log)
    # ("generic", "genericB"),        # f lam2, lam4 grid
    # ("generic", "genericC"),        # g lam1, lam2 grid
    # ("generic", "genericD"),        # r lam4, tanb grid (flat-log)
    # ("generic", "genericE"),        # t tanb grid (log) (high m122)
    # ("generic", "genericF"),        # u tanb grid (log) (vlow m122)
    # ("generic", "genericG"),        # v tanb grid (log) (low m122)

    # ("hybrid_Higgs", "hybrid1A"),   # h Lam4, Lam5 grid
    # ("hybrid_Higgs", "hybrid1B"),   # i Lam4, Lam5, Lam7 grid
    # ("hybrid_Higgs", "hybrid1H"),   # q tanb grid (log)
    # ("hybrid_Higgs", "hybrid1C"),   # j cba, tanb grid (flat-log)
    # ("hybrid_Higgs", "hybrid1D"),   # k mH, tanb grid (log-log)
    # ("hybrid_Higgs", "hybrid1E"),   # l mH, tanb grid (log-log) (higher mH)
    # ("hybrid_Higgs", "hybrid1F"),   # m cba, tanb grid (flat-log) Lower cba
    # ("hybrid_Higgs", "hybrid1G"),   # p mH, tanb (lower mH) -> better off doing full tb grid

    # ("hybrid_Higgs2", "hybrid2A"),  # n mHp,tanb grid (log-log) -> better off doing full tb grid
    # ("hybrid_Higgs2", "hybrid2B"),  # o mHp,tanb grid (log-log) (bottom right corner)
]


# ---- Likelihood conbinations to generate scans for -----

# Combined Scans
CONSTRAINTS_A = [

    (["theory"], "theory"),
    (["theory", "electroweak"], "electroweak"),
    (["theory", "collider"], "collider"),
    (["theory", "collider", "electroweak"], "most"),
]

# Theory Scans
CONSTRAINTS_B = [
  
    (["scalar_mass_corrections_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM"], "NLO"),
    (["theory"], "theoryZ"),
    (["theory", "electroweak"], "electroweakZ"),
    (["theory", "collider"], "colliderZ"),
    (["theory", "electroweak", "collider"], "mostZ"),
    (["scalar_mass_corrections_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM"], "NLO"),
    (["scalar_mass_corrections_LogLikelihood_THDM", "stability_LogLikelihood_THDM"], "stability"),
    (["runToScaleTest_LogLikelihood_THDM", "scalar_mass_corrections_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM"], "perturbativity_hybrid2X"),
    (["runToScaleTest_LogLikelihood_THDM", "LO_unitarity_LogLikelihood_THDM"  , "stability_LogLikelihood_THDM", 
      "higgs_exp_mass_LogLikelihood_THDM", "higgs_scenario_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM", "perturbativity_yukawas_LogLikelihood_THDM"], "theoryTA"),
    (["runToScaleTest_LogLikelihood_THDM", "NLO_unitarity_LogLikelihood_THDM" , "stability_LogLikelihood_THDM", 
      "higgs_exp_mass_LogLikelihood_THDM", "higgs_scenario_LogLikelihood_THDM", "perturbativity_LogLikelihood_THDM", "perturbativity_yukawas_LogLikelihood_THDM"], "bayesianA"),
    (["perturbativity_yukawas_LogLikelihood_THDM"], "pert_yukawas")
]

# Electroweak Scans
CONSTRAINTS_C = [
  
    (["theory", "electroweak"], "electroweak"),
]

# Collider Scans
CONSTRAINTS_D = [
  
    (["theory", "electroweak"], "electroweak"),
    (["theory", "collider"], "collider"),
    (["theory", "LEP_Higgs_LogLike"], "HB"),
    (["theory", "LHC_Higgs_LogLike", "HS_RUN1_SS"], "HSRUN1SS"),
    (["theory", "LHC_Higgs_LogLike", "HS_LATEST_SS"], "HSLATESTSS"),
    (["theory", "LHC_Higgs_LogLike", "HS_LATEST_STXS"], "HSLATESTSTXS"),
    (["theory", "LHC_Higgs_LogLike", "HS_ALL"], "HS"),
]

# Flavour Scans (Individual)
CONSTRAINTS_E = [

    # (["b2sgamma_LogLikelihood"], "flavor/B2Xsgamma/b2sgamma"),
    # (["B2Kstargamma_LogLikelihood"], "flavor/B2Xsgamma/B2Kstargamma"),
    # (["B2Kstargamma_LogLikelihood", "b2sgamma_LogLikelihood"], "flavor/B2Xsgamma/comb"),

    # (["B2mumu_LogLikelihood_Atlas"], "flavor/B2mumu/B2mumu_Atlas"),
    # (["B2mumu_LogLikelihood_LHCb"], "flavor/B2mumu/B2mumu_LHCb"), # missing -> fixed
    # (["B2mumu_LogLikelihood_CMS"], "flavor/B2mumu/B2mumu_CMS"),
    # (["B2mumu_LogLikelihood_Atlas", "B2mumu_LogLikelihood_LHCb", "B2mumu_LogLikelihood_CMS"], "flavor/B2mumu/combA"), # missing -> fixed
    # (["B2mumu_LogLikelihood_CMS_ATLAS_LHCb"], "flavor/B2mumu/combB"), # wrong -> fixed

    # (["Bd2KmumuBr_LogLikelihood_LHCb"], "flavor/B2KmumuBr/Bd2KmumuBr_LHCb"), # wrong?
    # (["Bd2KmumuBr_LogLikelihood_Belle"], "flavor/B2KmumuBr/Bd2KmumuBr_Belle"), # wrong?
    # (["B2KmumuBr_LogLikelihood_LHCb"], "flavor/B2KmumuBr/B2KmumuBr_LHCb"),
    # (["B2KmumuBr_LogLikelihood_CMS"], "flavor/B2KmumuBr/B2KmumuBr_CMS"), # wrong?
    # (["B2KmumuBr_LogLikelihood_Belle"], "flavor/B2KmumuBr/B2KmumuBr_Belle"), # wrong
    # (["Bd2KmumuBr_LogLikelihood_LHCb", "Bd2KmumuBr_LogLikelihood_Belle", "B2KmumuBr_LogLikelihood_LHCb", "B2KmumuBr_LogLikelihood_CMS", "B2KmumuBr_LogLikelihood_Belle"], "flavor/B2KmumuBr/comb"),

    # (["B2KeeBr_LogLikelihood_Belle"], "flavor/B2Kee/B2KeeBr_Belle"),
    # (["Bd2KeeBr_LogLikelihood_Belle"], "flavor/B2Kee/Bd2KeeBr_Belle"),
    # (["Bd2KeeBr_LogLikelihood_Belle", "B2KeeBr_LogLikelihood_Belle"], "flavor/B2Kee/comb"),

    # (["B2KstarmumuAng_LogLikelihood_Atlas"], "flavor/B2Kstarmumu/B2KstarmumuAng_Atlas"), # wrong?
    # (["B2KstarmumuAng_LogLikelihood_Belle"], "flavor/B2Kstarmumu/B2KstarmumuAng_Belle"), # wrong?
    (["B2KstarmumuAng_LogLikelihood_LHCb_2020"], "flavor/B2Kstarmumu/B2KstarmumuAng_LHCb_2020"), # wrong?
    # (["Bu2KstarmumuAng_LogLikelihood_LHCb_2020"], "flavor/B2Kstarmumu/Bu2KstarmumuAng_LHCb_2020"), # wrong?
    # (["B2KstarmumuBr_LogLikelihood_LHCb"], "flavor/B2Kstarmumu/B2KstarmumuBr_LHCb"),
    # (["Bs2phimumuBr_LogLikelihood"], "flavor/B2Kstarmumu/Bs2phimumuBr"),
    # (["B2KstarmumuAng_CPAssym_LogLikelihood_LHCb"], "flavor/B2Kstarmumu/B2KstarmumuAng_CPAssym_LHCb"),
    # (["B2KstarmumuAng_LogLikelihood_CMS"], "flavor/B2Kstarmumu/B2KstarmumuAng_CMS"),
    (["B2KstarellellAng_LogLikelihood_Belle"], "flavor/B2Kstarmumu/B2KstarellellAng_Belle"), # wrong
    # (["B2KstarmumuAng_LogLikelihood_Atlas", "B2KstarmumuAng_LogLikelihood_Belle", "B2KstarmumuAng_LogLikelihood_LHCb_2020", "Bu2KstarmumuAng_LogLikelihood_LHCb_2020", "B2KstarmumuBr_LogLikelihood_LHCb", "Bs2phimumuBr_LogLikelihood", "B2KstarmumuAng_CPAssym_LogLikelihood_LHCb","B2KstarmumuAng_LogLikelihood_CMS", "B2KstarellellAng_LogLikelihood_Belle"], "flavor/B2Kstarmumu/comb"),

    # (["B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020"], "flavor/B2KstareeAng_Lowq2_LHCb_2020/B2KstareeAng_Lowq2_LHCb_2020"), # wrong?

    # (["RKRKstar_LogLikelihood_LHCb"], "flavor/RKRKstar/RKRKstar_LHCb"), # wrong?
    # (["RK_LogLikelihood_CMS"], "flavor/RKRKstar/RK_CMS"), # wrong?
    # (["RK_LogLikelihood_Belle"], "flavor/RKRKstar/RK_Belle"), # wrong?
    # (["RKRKstar_LogLikelihood_LHCb", "RK_LogLikelihood_Belle", "RK_LogLikelihood_CMS"], "flavor/RKRKstar/comb"),

    # # (["BKnunu_LogLikelihood_Belle_sl"], "flavor/B2Knunu/BKnunu_Belle_sl"),
    # # (["BKnunu_LogLikelihood_Belle_had"], "flavor/B2Knunu/BKnunu_Belle_had"),
    # # (["BuKnunu_LogLikelihood_Belle_sl"], "flavor/B2Knunu/BuKnunu_Belle_sl"),
    # # (["BuKnunu_LogLikelihood_Belle_had"], "flavor/B2Knunu/BuKnunu_Belle_had"),
    # # (["BuKnunu_LogLikelihood_BelleII"], "flavor/B2Knunu/BuKnunu_BelleII"),
    # # (["BKnunu_LogLikelihood_BaBar"], "flavor/B2Knunu/BKnunu_BaBar"),
    # # (["BuKnunu_LogLikelihood_BaBar"], "flavor/B2Knunu/BuKnunu_BaBar"),
    # # (["BKnunu_LogLikelihood_Belle_sl", "BKnunu_LogLikelihood_Belle_had", "BuKnunu_LogLikelihood_Belle_sl", "BuKnunu_LogLikelihood_Belle_had", "BuKnunu_LogLikelihood_BelleII", "BKnunu_LogLikelihood_BaBar", "BuKnunu_LogLikelihood_BaBar"], "flavor/B2Knunu/comb"),

    # # (["BKstarnunu_LogLikelihood_Belle_sl"], "flavor/BKstarnunu_Belle_sl"), # not ready
    # # (["BKstarnunu_LogLikelihood_Belle_had"], "flavor/BKstarnunu_Belle_had"), # not ready
    # # (["BuKstarnunu_LogLikelihood_Belle_sl"], "flavor/BuKstarnunu_Belle_sl"), # not ready
    # # (["BuKstarnunu_LogLikelihood_Belle_had"], "flavor/BuKstarnunu_Belle_had"), # not ready
    # # (["BKstarnunu_LogLikelihood_BaBar"], "flavor/BKstarnunu_BaBar"), # not ready
    # # (["BuKstarnunu_LogLikelihood_BaBar"], "flavor/BuKstarnunu_BaBar"), # not ready

    # # ([" SL_LogLikelihood"], "flavor/SL_FCCC/RD_RDstar"),
    # (["SL_LogLikelihood"], "flavor/SL_FCCC/SL"),
    # (["FLDstar_LogLikelihood"], "flavor/SL_FCCC/FLDstar"),
    # (["dBRBDstartaunu_LogLikelihood"], "flavor/SL_FCCC/dBRBDstartaunu"),
    # (["dBRBDtaunu_LogLikelihood"], "flavor/SL_FCCC/dBRBDtaunu"),
    # (["SL_LogLikelihood", "FLDstar_LogLikelihood", "dBRBDstartaunu_LogLikelihood", "dBRBDtaunu_LogLikelihood"], "flavor/SL_FCCC/comb"),

    # (["Bc_lifetime_LogLikelihood"], "flavor/DeltaMB/Bc_lifetime"), # wrong?
    # (["Delta_MBs_LogLikelihood"], "flavor/DeltaMB/Delta_MBs"),
    # (["Delta_MBd_LogLikelihood"], "flavor/DeltaMB/Delta_MBd"),
    # (["Bc_lifetime_LogLikelihood", "Delta_MBs_LogLikelihood", "Delta_MBd_LogLikelihood"], "flavor/DeltaMB/comb"),

    # # (["l2lgamma_LogLikelihood"], "flavor/l2lgamma"), # only for g2hdm
    # # (["l2lll_LogLikelihood"], "flavor/l2lll"), # only for g2hdm
    # # (["h2ltau_LogLikelihood"], "flavor/h2ltau"), # only for g2hdm
    # # (["l2lgamma_LogLikelihood", "l2lll_LogLikelihood", "h2ltau_LogLikelihood"], "flavor/comb_LFV"), # only for g2hdm

    # # (["t2ch_LogLikelihood"], "flavor/t2ch"), # only for g2hdm
    # # (["t2bbc_LogLikelihood"], "flavor/t2bbc"), # only for g2hdm
    # # (["t2mutauc_LogLikelihood"], "flavor/t2mutauc"), # only for g2hdm
    # # (["Bc2taunu_LogLikelihood"], "flavor/Bc2taunu"), # only for g2hdm
    # # (["Bs2ll_LogLikelihood"], "flavor/Bs2ll"), # only for g2hdm
    # # (["B2Kll_LogLikelihood"], "flavor/B2Kll"), # only for g2hdm
    # # (["t2ch_LogLikelihood", "t2bbc_LogLikelihood", "t2mutauc_LogLikelihood", "Bc2taunu_LogLikelihood", "Bs2ll_LogLikelihood", "B2Kll_LogLikelihood"], "flavor/comb_FV_top"), # only for g2hdm

    # (["B2Xsnunu_LogLikelihood"], "flavor/B2Xsnunu"), # missing WCs
    # # (["gmu_ge_LogLikelihood"], "flavor/gmu_ge"), # only for g2hdm

]

# Flavour Scans (Combined)
CONSTRAINTS_F = [

    ([
        "B2Kstargamma_LogLikelihood",
        "b2sgamma_LogLikelihood",
        # -----
        "B2mumu_LogLikelihood_Atlas", 
        "B2mumu_LogLikelihood_LHCb",
        # "B2mumu_LogLikelihood_CMS", # OFF - broken; two lines
        # "B2mumu_LogLikelihood_CMS_ATLAS_LHCb", # OFF - broken
        # -----
        "Bd2KmumuBr_LogLikelihood_LHCb", 
        "Bd2KmumuBr_LogLikelihood_Belle", 
        "B2KmumuBr_LogLikelihood_LHCb", 
        # "B2KmumuBr_LogLikelihood_CMS", # OFF - broken; spikes; strong for Type-II
        "B2KmumuBr_LogLikelihood_Belle",
        # -----
        "B2KeeBr_LogLikelihood_Belle",
        "Bd2KeeBr_LogLikelihood_Belle", 
        # -----
        "B2KstarmumuAng_LogLikelihood_Atlas",
        # "B2KstarmumuAng_LogLikelihood_CMS", # OFF - broken; data looks weird; missing point issue??
        "B2KstarmumuAng_LogLikelihood_Belle",
        "B2KstarmumuAng_LogLikelihood_LHCb_2020",
        "Bu2KstarmumuAng_LogLikelihood_LHCb_2020", 
        "B2KstarmumuBr_LogLikelihood_LHCb", 
        "Bs2phimumuBr_LogLikelihood", 
        "B2KstarmumuAng_CPAssym_LogLikelihood_LHCb",
        "B2KstarellellAng_LogLikelihood_Belle", # OFF - strong for Type-II
        # -----
        "B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020",
        # -----
        "RKRKstar_LogLikelihood_LHCb",
        "RK_LogLikelihood_Belle", 
        "RK_LogLikelihood_CMS",
        # -----
        "BKnunu_LogLikelihood_Belle_sl", 
        "BKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_Belle_sl", 
        "BuKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_BelleII", 
        "BKnunu_LogLikelihood_BaBar", 
        "BuKnunu_LogLikelihood_BaBar",
        # -----
        "SL_LogLikelihood", 
        # "FLDstar_LogLikelihood", 
        # "dBRBDstartaunu_LogLikelihood", 
        # "dBRBDtaunu_LogLikelihood",
        # -----
        # "Bc_lifetime_LogLikelihood", 
        "Delta_MBs_LogLikelihood"], "flavor/combA"),

    ([
        "B2Kstargamma_LogLikelihood",
        "b2sgamma_LogLikelihood",
        # -----
        "B2mumu_LogLikelihood_Atlas", 
        "B2mumu_LogLikelihood_LHCb",
        # "B2mumu_LogLikelihood_CMS", # OFF - broken; two lines
        # "B2mumu_LogLikelihood_CMS_ATLAS_LHCb", # OFF - broken
        # -----
        "Bd2KmumuBr_LogLikelihood_LHCb", 
        "Bd2KmumuBr_LogLikelihood_Belle", 
        "B2KmumuBr_LogLikelihood_LHCb", 
        # "B2KmumuBr_LogLikelihood_CMS", # OFF - broken; spikes; strong for Type-II
        "B2KmumuBr_LogLikelihood_Belle",
        # -----
        "B2KeeBr_LogLikelihood_Belle",
        "Bd2KeeBr_LogLikelihood_Belle", 
        # -----
        "B2KstarmumuAng_LogLikelihood_Atlas",
        # "B2KstarmumuAng_LogLikelihood_CMS", # OFF - broken; data looks weird; missing point issue??
        "B2KstarmumuAng_LogLikelihood_Belle",
        "B2KstarmumuAng_LogLikelihood_LHCb_2020",
        "Bu2KstarmumuAng_LogLikelihood_LHCb_2020", 
        "B2KstarmumuBr_LogLikelihood_LHCb", 
        "Bs2phimumuBr_LogLikelihood", 
        "B2KstarmumuAng_CPAssym_LogLikelihood_LHCb",
        # "B2KstarellellAng_LogLikelihood_Belle", # OFF - strong for Type-II
        # -----
        "B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020",
        # -----
        "RKRKstar_LogLikelihood_LHCb",
        "RK_LogLikelihood_Belle", 
        "RK_LogLikelihood_CMS",
        # -----
        "BKnunu_LogLikelihood_Belle_sl", 
        "BKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_Belle_sl", 
        "BuKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_BelleII", 
        "BKnunu_LogLikelihood_BaBar", 
        "BuKnunu_LogLikelihood_BaBar",
        # -----
        "SL_LogLikelihood", 
        "FLDstar_LogLikelihood", 
        "dBRBDstartaunu_LogLikelihood", 
        "dBRBDtaunu_LogLikelihood",
        # -----
        "Bc_lifetime_LogLikelihood", 
        "Delta_MBs_LogLikelihood"], "flavor/combB"),

    ([
        "B2Kstargamma_LogLikelihood",
        "b2sgamma_LogLikelihood",
        # -----
        "B2mumu_LogLikelihood_Atlas", 
        "B2mumu_LogLikelihood_LHCb",
        # "B2mumu_LogLikelihood_CMS", # OFF - broken; two lines
        # "B2mumu_LogLikelihood_CMS_ATLAS_LHCb", # OFF - broken
        # -----
        # "Bd2KmumuBr_LogLikelihood_LHCb", 
        # "Bd2KmumuBr_LogLikelihood_Belle", 
        # "B2KmumuBr_LogLikelihood_LHCb", 
        # # "B2KmumuBr_LogLikelihood_CMS", # OFF - broken; spikes; strong for Type-II
        # "B2KmumuBr_LogLikelihood_Belle",
        # -----
        # "B2KeeBr_LogLikelihood_Belle",
        # "Bd2KeeBr_LogLikelihood_Belle", 
        # -----
        # "B2KstarmumuAng_LogLikelihood_Atlas",
        # # "B2KstarmumuAng_LogLikelihood_CMS", # OFF - broken; data looks weird; missing point issue??
        # "B2KstarmumuAng_LogLikelihood_Belle",
        # "B2KstarmumuAng_LogLikelihood_LHCb_2020",
        # "Bu2KstarmumuAng_LogLikelihood_LHCb_2020", 
        # "B2KstarmumuBr_LogLikelihood_LHCb", 
        # "Bs2phimumuBr_LogLikelihood", 
        # "B2KstarmumuAng_CPAssym_LogLikelihood_LHCb",
        # # "B2KstarellellAng_LogLikelihood_Belle", # OFF - strong for Type-II
        # -----
        # "B2KstareeAng_Lowq2_LogLikelihood_LHCb_2020",
        # -----
        # "RKRKstar_LogLikelihood_LHCb",
        # "RK_LogLikelihood_Belle", 
        # "RK_LogLikelihood_CMS",
        # -----
        "BKnunu_LogLikelihood_Belle_sl", 
        "BKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_Belle_sl", 
        "BuKnunu_LogLikelihood_Belle_had", 
        "BuKnunu_LogLikelihood_BelleII", 
        "BKnunu_LogLikelihood_BaBar", 
        "BuKnunu_LogLikelihood_BaBar",
        # -----
        "SL_LogLikelihood", 
        # "FLDstar_LogLikelihood", 
        # "dBRBDstartaunu_LogLikelihood", 
        # "dBRBDtaunu_LogLikelihood",
        # -----
        # "Bc_lifetime_LogLikelihood", 
        "Delta_MBs_LogLikelihood"], "flavor/combC"),
]

# -----

# Set the final list of likelihood combinations
CONSTRAINTS = CONSTRAINTS_E

# These are appended to all above combinations
ADDITIONAL_CONSTRAINTS = [] # ("theory", "electroweak", or "collider")
