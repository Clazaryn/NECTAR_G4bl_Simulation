# NECTAR G4bl Simulation

The NECTAR G4 Beamline (G4bl) simulation is a 3-step simulation pipeline that simulates the NECTAR surrogate reaction experiments conducted at the Experimental Storage Ring at GSI, Darmstadt. This repository hosts the most generalised version of the simulation, but specialised analysis formats can and should be forked.

Package in current form assembled by Guy Leckenby. `event_generator.py`, `nuc_reaction.py` originally written by Ana Henriques, Michele Sguazzin 
from the MatLab event generator by Manfred Grieser, with relativistic updates from Camille Berthelot. `ESR_ring.in` simulation originally written by Michele Sguazzin and Manfred Grieser. `analysis.cxx` originally written by Michele Sguazzin and updated by Camille Berthelot.

## Structure of Simulation

The simulation uses 3 bash scripts to run three important but separate components of the simulation successively:

1. **`bscript-gen_events.sh`** - Uses a Python script implementing relativistic two-body kinematics with realistic CN decay to generate a coupled pair of ejectile and heavy recoil event lists.
2. **`bscript-G4bl_run.sh`** - Runs G4 beamline using the events produced in step 1.
3. **`bash_Analysis_nChamb.sh`** - Not yet organised but will analyse the G4bl output and construct a ROOT tree of all the simulated events.

The outputs of these bash scripts will be stored in a reaction directory in the *parent directory of the code files* (the bash scripts automatically construct `reaction_sim` directories). Thus the code should be hosted in a directory called `./Code` to create the following structure:

```
NECTAR_G4bl_Sim/
└───Code/
│   │   (all code files)
│
└───238U_dp_sim/  
│   └───Detector_output/
│   └───Event_output/
│   └───Hist_output/
```

### Structure via Excitation Energy

To parallelise the calculation, the simulation is broken down into 1 MeV chunks. The number of cores used by each bash script can be edited at the top; the default is 14 to leave two spare on borlin305.

Whilst E* is randomly sampled across these chunks with the exact value tracked in the final ROOT tree, the output files are separated to allow easy tracking. Some features of the sampling are also E*-dependent (e.g. theta sampling to match data rates), so these are tracked with 1 MeV bins.

### `reac_info.txt` Input File

All aspects of the simulation refer to `reac_info.txt` where the user input is kept, so the user should only need to modify one file to simulate a new reaction. Whilst the fields are self-explanatory, some require specific attention:

- The magnetic optic values can be calculated using the `Qi_dipole` Excel spreadsheets stored in `./ReacInfoFiles`
- The masses used in the input file are the ones sourced from the NuclideDataMaster in `event_generator.py`. These are printed if `event_generator.py` is run in Verbose mode.

If you change the beam type or energy, you MUST update the ring parameters Brho, Bfield, and Qi.

### Limiting angular range sampled
In `nuc_reaction.py`, you can limit the angular range and kinematic solution that is sampled:
 - on line 37: you can define the telescope theta range. The range is defined from the CDF to ensure isotropic sampling.
 - on line 112: you can define the telescope phi range in radians.
 - on line 91: you can define the kinematic solution. Default is random.

## Detailed Description of Files

### Main Simulation Scripts

- **`bscript-gen_events.sh`** - Bash script for generating events. The complicated structure is to enable the script to run processes in parallel automatically. Calls:
  - **`event_generator.py`** - Runs a while loop that generates events for a given channel. Calls a two-body collision and then decays the compound nucleus. Currently uses custom coded 4-vectors to keep things correctly relativistic.
  - **`nuc_reaction.py`** - Where the two-body collision and CN decay are implemented. Lorentz boosts are correctly implemented.

- **`bscript-G4bl_run.sh`** - Runs the G4bl based on the input files. The heavy recoils take quite long due to simulating the magnetic fields, so it can be efficient to just rerun the light ejectiles if only the telescope geometry has changed. This must be chosen manually at the top of the script. Note that all output is suppressed and sent to the output file, so be sure to check the output file for potential bugs and errors. Takes:
  - **`ESR_ring-mode.in`** - Two input setups are available, the Proof-of-Principle detectors in full detail, and the new telescopes in a mostly schematic way.

- **`bash_Analysis_nChamb.sh`** - Currently not fixed... Calls:
  - **`analysis_nChamb.C`** - The C++ macro to analyse the G4bl output and convert to a ROOT tree

### Supporting Directories

#### NuclideDataMaster

- **`nuclide_data.py`** - Implements a complicated dictionary to access NIST data for all isotopes
- **`Electron_Binding_Energies-2007.dat`** - Yuri's special binding energy table for calculating the nuclear masses of highly-charged ions.

#### ReacInfoFiles

- **`Qi_dipole-mode.xlsx`** - These Excel sheets replicate Manfred's calculations to adjust the quadrupole strengths given a certain beam + energy. The kq coefficients change based on the ESR operating mode, so two are available.
- Some other precalculated `reac_info.txt` files are available.

#### UtilityScripts

- **`calc_hr_ranges.py`** - A helper script to automatically calculate the iteration ranges in the bash script for each channel based off the separation energies.
- **`funcs_plotting.h`** / **`funcs_plotting.C`** - Scripts including many functions used to create plots for the new telescope simulations.
- **`ini_parser.h`** - A header containing functions to allow C++ to read the INI format in general.
- **`phys_functions.py`** - Some helper functions used by `event_generator.py`.
- **`reaction_info.h`** - A header containing a function to allow C++ to specifically read our `reac_info.txt` input file.
