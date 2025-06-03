import os
import numpy as np
from openfast_io.FAST_reader import InputReader_OpenFAST as OrigInputReader_OpenFAST
from openfast_io.FAST_writer import InputWriter_OpenFAST as OrigInputWriter_OpenFAST
from openfast_io.FAST_reader import bool_read, quoted_read, float_read, read_array

class InputReader_OpenFAST(OrigInputReader_OpenFAST):
    """Extended version of InputReader_OpenFAST with custom read_MainInput method"""
    
    def read_MainInput(self):
        """Override of the original read_MainInput method to include extra lines"""

        # Main FAST Input File
        fst_file = os.path.join(self.FAST_directory, self.FAST_InputFile)
        f = open(fst_file)

        # Header of .fst file
        f.readline()
        self.fst_vt['description'] = f.readline().rstrip()

        # Simulation Control (fst_sim_ctrl)
        f.readline()
        self.fst_vt['Fst']['Echo'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['AbortLevel'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['TMax'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['DT']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['InterpOrder']  = int(f.readline().split()[0])
        self.fst_vt['Fst']['NumCrctn']  = int(f.readline().split()[0])
        
        ########## Custom lines added here for tight coupling

        self.fst_vt['Fst']['RhoInf']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['ConvTol']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['MaxConvIter']  = int(f.readline().split()[0])

        ###########
        self.fst_vt['Fst']['DT_UJac']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['UJacSclFact']  = float_read(f.readline().split()[0])

        # Feature Switches and Flags (ftr_swtchs_flgs)
        f.readline()
        self.fst_vt['Fst']['CompElast'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompInflow'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompAero'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompServo'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompSeaSt'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompHydro'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompSub'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompMooring'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['CompIce'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['MHK'] = int(f.readline().split()[0])

        # Environmental conditions
        f.readline()
        self.fst_vt['Fst']['Gravity']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['AirDens']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['WtrDens']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['KinVisc']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['SpdSound']  = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['Patm']      = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['Pvap']      = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['WtrDpth']   = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['MSL2SWL']   = float_read(f.readline().split()[0])

        # Input Files (input_files)
        f.readline()
        self.fst_vt['Fst']['EDFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['BDBldFile(1)'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['BDBldFile(2)'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['BDBldFile(3)'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['InflowFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['AeroFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['ServoFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['SeaStFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['HydroFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['SubFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['MooringFile'] = quoted_read(f.readline().split()[0])
        self.fst_vt['Fst']['IceFile'] = quoted_read(f.readline().split()[0])

        # FAST Output Parameters (fst_output_params)
        f.readline()
        self.fst_vt['Fst']['SumPrint'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['SttsTime'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['ChkptTime'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['DT_Out'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['TStart'] = float_read(f.readline().split()[0])
        self.fst_vt['Fst']['OutFileFmt'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['TabDelim'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['OutFmt'] = quoted_read(f.readline().split()[0])

        # Fst
        f.readline()
        self.fst_vt['Fst']['Linearize']  = f.readline().split()[0]
        self.fst_vt['Fst']['CalcSteady'] = f.readline().split()[0]
        self.fst_vt['Fst']['TrimCase']   = f.readline().split()[0]
        self.fst_vt['Fst']['TrimTol']    = f.readline().split()[0]
        self.fst_vt['Fst']['TrimGain']   = f.readline().split()[0]
        self.fst_vt['Fst']['Twr_Kdmp']   = f.readline().split()[0]
        self.fst_vt['Fst']['Bld_Kdmp']   = f.readline().split()[0]
        self.fst_vt['Fst']['NLinTimes']  = int(f.readline().split()[0])
        self.fst_vt['Fst']['LinTimes']   = read_array(f, self.fst_vt['Fst']['NLinTimes'], array_type=float) 
        self.fst_vt['Fst']['LinInputs']  = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutputs'] = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutJac']  = f.readline().split()[0]
        self.fst_vt['Fst']['LinOutMod']  = f.readline().split()[0]

        # Visualization ()
        f.readline()
        self.fst_vt['Fst']['WrVTK'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_type'] = int(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_fields'] = bool_read(f.readline().split()[0])
        self.fst_vt['Fst']['VTK_fps'] = float_read(f.readline().split()[0])
        
        f.close()
        
        # File paths
        self.fst_vt['Fst']['EDFile_path']       = os.path.split(self.fst_vt['Fst']['EDFile'])[0]
        self.fst_vt['Fst']['BDBldFile(1_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(1)'])[0]
        self.fst_vt['Fst']['BDBldFile(2_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(2)'])[0]
        self.fst_vt['Fst']['BDBldFile(3_path)'] = os.path.split(self.fst_vt['Fst']['BDBldFile(3)'])[0]
        self.fst_vt['Fst']['InflowFile_path']   = os.path.split(self.fst_vt['Fst']['InflowFile'])[0]
        self.fst_vt['Fst']['AeroFile_path']     = os.path.split(self.fst_vt['Fst']['AeroFile'])[0]
        self.fst_vt['Fst']['ServoFile_path']    = os.path.split(self.fst_vt['Fst']['ServoFile'])[0]
        self.fst_vt['Fst']['HydroFile_path']    = os.path.split(self.fst_vt['Fst']['HydroFile'])[0]
        self.fst_vt['Fst']['SubFile_path']      = os.path.split(self.fst_vt['Fst']['SubFile'])[0]
        self.fst_vt['Fst']['MooringFile_path']  = os.path.split(self.fst_vt['Fst']['MooringFile'])[0]
        self.fst_vt['Fst']['IceFile_path']      = os.path.split(self.fst_vt['Fst']['IceFile'])[0]

class InputWriter_OpenFAST(OrigInputWriter_OpenFAST):
    """Extended version of InputWriter_OpenFAST with custom write_MainInput method"""
    
    def write_MainInput(self):
        """Override of the original write_MainInput method to include extra lines"""
        
        # Main FAST Input File
        self.FAST_InputFileOut = os.path.join(self.FAST_runDirectory, self.FAST_namingOut+'.fst')

        # Keep simple for now:
        f = open(self.FAST_InputFileOut, 'w')

        # ===== .fst Input File =====

        f.write('------- OpenFAST INPUT FILE -------------------------------------------\n')
        f.write('Generated with OpenFAST_IO - Extended version\n')
        f.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Echo'], 'Echo', '- Echo input data to <RootName>.ech (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['AbortLevel']+'"', 'AbortLevel', '- Error level when simulation should abort (string) {"WARNING", "SEVERE", "FATAL"}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TMax'], 'TMax', '- Total run time (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT'], 'DT', '- Recommended module time step (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['InterpOrder'], 'InterpOrder', '- Interpolation order for input/output time history (-) {1=linear, 2=quadratic}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['NumCrctn'], 'NumCrctn', '- Number of correction iterations (-) {0=explicit calculation, i.e., no corrections}\n'))

        ### Custom lines added here for tight coupling

        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['RhoInf'], 'RhoInf', '- Convergence iteration error tolerance for tight coupling generalized alpha integrator (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['ConvTol'], 'ConvTol', '- Maximum number of convergence iterations for tight coupling generalized alpha integrator (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['MaxConvIter'], 'MaxConvIter', '- Number of correction iterations (-) {0=explicit calculation, i.e., no corrections}\n'))

        ####

        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT_UJac'], 'DT_UJac', '- Time between calls to get Jacobians (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['UJacSclFact'], 'UJacSclFact', '- Scaling factor used in Jacobians (-)\n'))
        f.write('---------------------- FEATURE SWITCHES AND FLAGS ------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompElast'], 'CompElast', '- Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades; 3=Simplified ElastoDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompInflow'], 'CompInflow', '- Compute inflow wind velocities (switch) {0=still air; 1=InflowWind; 2=external from ExtInflow}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompAero'], 'CompAero', '- Compute aerodynamic loads (switch) {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompServo'], 'CompServo', '- Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompSeaSt'], 'CompSeaSt', '- Compute sea state information (switch) {0=None; 1=SeaState}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompHydro'], 'CompHydro', '- Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompSub'], 'CompSub', '- Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=External Platform MCKF}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompMooring'], 'CompMooring', '- Compute mooring system (switch) {0=None; 1=MAP++; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CompIce'], 'CompIce', '- Compute ice loads (switch) {0=None; 1=IceFloe; 2=IceDyn}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['MHK'], 'MHK', '- MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}\n'))
        f.write('---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Gravity'], 'Gravity', '- Gravitational acceleration (m/s^2)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['AirDens'], 'AirDens', '- Air density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WtrDens'], 'WtrDens', '- Water density (kg/m^3)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['KinVisc'], 'KinVisc', '- Kinematic viscosity of working fluid (m^2/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SpdSound'], 'SpdSound', '- Speed of sound in working fluid (m/s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Patm'], 'Patm', '- Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Pvap'], 'Pvap', '- Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WtrDpth'], 'WtrDpth', '- Water depth (m)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['MSL2SWL'], 'MSL2SWL', '- Offset between still-water level and mean sea level (m) [positive upward]\n'))        
        f.write('---------------------- INPUT FILES ---------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['EDFile']+'"', 'EDFile', '- Name of file containing ElastoDyn input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(1)']+'"', 'BDBldFile(1)', '- Name of file containing BeamDyn input parameters for blade 1 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(2)']+'"', 'BDBldFile(2)', '- Name of file containing BeamDyn input parameters for blade 2 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['BDBldFile(3)']+'"', 'BDBldFile(3)', '- Name of file containing BeamDyn input parameters for blade 3 (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['InflowFile']+'"', 'InflowFile', '- Name of file containing inflow wind input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['AeroFile']+'"', 'AeroFile', '- Name of file containing aerodynamic input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['ServoFile']+'"', 'ServoFile', '- Name of file containing control and electrical-drive input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['SeaStFile']+'"', 'SeaStFile', '- Name of file containing sea state input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['HydroFile']+'"', 'HydroFile', '- Name of file containing hydrodynamic input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['SubFile']+'"', 'SubFile', '- Name of file containing sub-structural input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['MooringFile']+'"', 'MooringFile', '- Name of file containing mooring system input parameters (quoted string)\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['IceFile']+'"', 'IceFile', '- Name of file containing ice input parameters (quoted string)\n'))
        f.write('---------------------- OUTPUT --------------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SumPrint'], 'SumPrint', '- Print summary data to "<RootName>.sum" (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['SttsTime'], 'SttsTime', '- Amount of time between screen status messages (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['ChkptTime'], 'ChkptTime', '- Amount of time between creating checkpoint files for potential restart (s)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['DT_Out'], 'DT_Out', '- Time step for tabular output (s) (or "default")\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TStart'], 'TStart', '- Time to begin tabular output (s)\n'))
        f.write('{:<22d} {:<11} {:}'.format(self.fst_vt['Fst']['OutFileFmt'], 'OutFileFmt', '- Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both 1 and 2, 4: uncompressed binary [<RootName>.outb], 5: both 1 and 4}\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TabDelim'], 'TabDelim', '- Use tab delimiters in text tabular output file? (flag) {uses spaces if false}\n'))
        f.write('{:<22} {:<11} {:}'.format('"'+self.fst_vt['Fst']['OutFmt']+'"', 'OutFmt', '- Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)\n'))
        f.write('---------------------- LINEARIZATION -------------------------------------------\n')
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Linearize'],   'Linearize',    '- Linearization analysis (flag)\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['CalcSteady'],  'CalcSteady',   '- Calculate a steady-state periodic operating point before linearization? [unused if Linearize=False] (flag)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimCase'],      'TrimCase',     '- Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only if CalcSteady=True] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimTol'],       'TrimTol',      '- Tolerance for the rotational speed convergence [used only if CalcSteady=True] (-)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['TrimGain'],      'TrimGain',     '- Proportional gain for the rotational speed error (>0) [used only if CalcSteady=True] (rad/(rad/s) for yaw or pitch; Nm/(rad/s) for torque)\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Twr_Kdmp'],      'Twr_Kdmp',     '- Damping factor for the tower [used only if CalcSteady=True] (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['Bld_Kdmp'],      'Bld_Kdmp',     '- Damping factor for the blades [used only if CalcSteady=True] (N/(m/s))\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['NLinTimes'],     'NLinTimes',    '- Number of times to linearize (-) [>=1] [unused if Linearize=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(', '.join(['%f'%i for i in np.array(self.fst_vt['Fst']['LinTimes'], dtype=float)]), 'LinTimes', '- List of times at which to linearize (s) [1 to NLinTimes] [used only when Linearize=True and CalcSteady=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinInputs'],     'LinInputs',    '- Inputs included in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)} [unused if Linearize=False]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutputs'],    'LinOutputs',   '- Outputs included in linearization (switch) {0=none; 1=from OutList(s); 2=all module outputs (debug)} [unused if Linearize=False]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutJac'],   'LinOutJac',    '- Include full Jacobians in linearization output (for debug) (flag) [unused if Linearize=False; used only if LinInputs=LinOutputs=2]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['LinOutMod'],   'LinOutMod',    '- Write module-level linearization output files in addition to output for full system? (flag) [unused if Linearize=False]\n'))
        f.write('---------------------- VISUALIZATION ------------------------------------------\n')
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['WrVTK'], 'WrVTK', '- VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_type'], 'VTK_type', '- Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]\n'))
        f.write('{!s:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_fields'], 'VTK_fields', '- Write mesh fields to VTK data files? (flag) {true/false} [unused if WrVTK=0]\n'))
        f.write('{:<22} {:<11} {:}'.format(self.fst_vt['Fst']['VTK_fps'], 'VTK_fps', '-Frame rate for VTK output (frames per second){will use closest integer multiple of DT} [used only if WrVTK=2 or WrVTK=3]\n'))

        f.flush()
        os.fsync(f)
        f.close()