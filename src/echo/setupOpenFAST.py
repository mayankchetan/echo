from echo.openfast_io_overload import InputReader_OpenFAST, InputWriter_OpenFAST
from openfast_io.FAST_output_reader import FASTOutputFile
from scipy import fft
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy


# Sweep variables
RPM = np.arange(0, 25, 1)

# Simulation location
sim_location = "/Users/mchetan/Desktop/nrel/projects/3-STABLE/workFolder/timeSeriesModeShapes/blade_forcing_x"


# Read in the base case
base_case = InputReader_OpenFAST()
base_case.FAST_InputFile = "simpleTurbine.fst"
base_case.FAST_directory = "/Users/mchetan/Desktop/nrel/projects/3-STABLE/tools/echo/tests/simplifiedTurbine/rot-blade-x"
base_case.execute()




for i in range(len(RPM)):

    print(f"Running case {i+1} of {len(RPM)}, setting RPM to {RPM[i]}")

    # Copy the base case
    case = copy.deepcopy(base_case)

    # Set the RPM
    case.fst_vt["ElastoDyn"]["RotSpeed"] = RPM[i]

    # create the write object
    writer = InputWriter_OpenFAST()
    writer.FAST_namingOut = f'bladeForcing_{RPM[i]:02d}_RPM'
    writer.FAST_runDirectory = sim_location #Output directory

    writer.fst_vt = case.fst_vt
    writer.execute()
