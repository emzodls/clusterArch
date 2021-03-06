import sys
from glob import glob
from cx_Freeze import setup, Executable



# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"
    includefiles = glob('winDeps/*') + ['testHMM.fasta','htmlvis']
if sys.platform == "darwin":
    includefiles = glob('macDeps/*') + ['testHMM.fasta','htmlvis']
if sys.platform == "linux":
    includefiles = glob('linuxDeps/*') + ['testHMM.fasta','htmlvis']
# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": [],
                     "excludes": ["tkinter","BioSQL","Bio.PDB","Bio.Nexus","Bio.SwissProt","BioSQL.BioSeq","Bio.AlignIO"],
                     'include_files':includefiles,"optimize":2}
if sys.platform == "win32":
    build_exe_options['include_msvcr'] = True
setup(  name = "clusterTools",
        version = "0.2",
        description = "BGC cluster analysis",
        options = {"build_exe": build_exe_options},
        executables = [Executable("main.py", base=base,icon='icon.ico')])