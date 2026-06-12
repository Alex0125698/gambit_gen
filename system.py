import os

_known_systems = [
  
  "dirac_icelake",
  "UbuntuDesktop",
  "UbuntuLegion",
  
]

# try to load system string

THE_SYSTEM = ""
ROOT_DIR = ""

try:
  with open("../system.txt","r") as f:
    THE_SYSTEM = f.read()
    ROOT_DIR = os.path.abspath("..")
except:
  pass
try:
  with open("../../system.txt","r") as f:
    THE_SYSTEM = f.read()
    ROOT_DIR = os.path.abspath("../..")
except:
  pass

if not THE_SYSTEM in _known_systems:
  raise Exception(f"Unknown/missing system: {THE_SYSTEM}")
