#!/usr/bin/env python3
import os

from wisdem import run_wisdem

mydir = os.path.dirname(os.path.realpath(__file__))  # get path to this file
fname_wt_input = mydir + os.sep + "simpleTurbine_tower.yaml"
fname_modeling_options = mydir + os.sep + "modeling.yaml"
fname_analysis_options = mydir + os.sep + "analysis.yaml"

wt_opt, analysis_options, opt_options = run_wisdem(fname_wt_input, fname_modeling_options, fname_analysis_options)

# print results from the analysis or optimization
print("mass (kg) =", wt_opt["towerse.tower_mass"])
print("cg (m) =", wt_opt["towerse.tower_center_of_mass"])
print("freq (Hz) =", wt_opt["towerse.tower.structural_frequencies"])
print("Fore-aft mode shapes =", wt_opt["towerse.tower.fore_aft_modes"])
print("Side-side mode shapes =", wt_opt["towerse.tower.side_side_modes"])
