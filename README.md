# CRADLE++

Customisable RAdioactive Decay for Low Energy Particle Physics: A C++ event generator

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

* [BOOST](https://www.boost.org/) - Linear algebra rotuines, and the `program_options` part of the library
* [GSL](https://www.gnu.org/software/gsl/) - Used for the calculation of the calculation of the complex Gamma function in the Fermi function
* [Geant4 Data](http://geant4.web.cern.ch/support/download) - G4PhotonEvaporation and G4RadioactiveDecay files are used internally
* [ROOT](https://root.cern/install/) - Used to generate ROOT file 

### Installing

Installing is done using [CMake](https://cmake.org/) using the supplied `CMakeLists.txt` file

Perform an out of source build, for example:

```
mkdir build
cd build && cmake ..
make
```

Setting some environnement variables : 
```
setenv Radiationdata ../RadiationData
setenv Gammadata ../GammaData
setenv AMEdata ../Nuclear_Databases/AMEdata.txt
```
## Usage 
### Config file
This file is an input of CRADLE to setting the generator. You can configure verbosity, coupling constant *C<sub>i*, cuts, $\beta$ decay type.

* Verbosity: 0 or 1 or 2 to see more or less information in the terminal.
* Verbosity_file: 0 to save only particles in the output file, 1 to save all the generated particles.
* *C<sub>i*: set all the coupling constant of the Lee Yang Lagragian.
* *a*, *b*, *A*, *B*, *D*, *c*: if not NaN overwritting the calculated parameters with*C<sub>i*
* Cuts: distance, time and energy limit of calculation
* Beta Decay: *Fermi* or *Gamow-Teller* to impose $\beta$ decay type. *Auto* take into account real $\beta$ decay type extracted from Geant4 *GammaData*, only available for pure transition (set on Gamow-Teller if mixed transition)
* FermiFunction: Simple or Advanced
* Alignement: setting its value
* Polarisation: setting its value and direction
* InFlightDecay: *true*/*false* (if *false* the kinetic energy of the recoil will be set to 0.)
* NuclearLevelWidth: *true*/*false* (if *true* Breit-Wigner distribution will be used for *proton and $\alpha$ decay if a lifetime is given in the RadiationData file)
* GammaGammaCorrelation: *true*/*false* (if *true* the angular correlation between 2 following $\gamma$ in a cascade is calculated)

### Run 
Running command example for 32Ar : 
```bash
./CRADLE++ nucleus --name 32Ar -Z 18 -A 32 general -l [NumberOfEvents] -o [OuputFileName] -c [ConfigFileName] -t [NumberOfThread]
```
Note that you can generate .txt file or .root file (TTree)
In the case of a ROOT file, input files as the Radioactive/Evaporation data and the config file will be saved using a TObjString.

## Details
See *Generator_Documentation.pdf* for more details about the generator.

## Reading output
In the folder Reader, one python and one C++ (ROOT macro based) script can read the ROOT output file.

## Version 
### v1.0 
- All features working
### v2.0
- New file writter, possibility to save .root or .txt file
- Auto mode for $\beta$ decay type
### v2.1
- Quick implementation of nuclear level width (Breit-Wigner) for proton and alpha decay (taking into account only the width of the current decaying level.)
### v2.2
- Implementation of Internal Conversion 
- Implementation of Electron Capture
### v2.3
- Improved Multi-Threading Mode
### v2.4 
- Cleaning code 
- Using PDG instead particle names
- Fixing errors
- Printing paramaters at the begining
- Writting input decay data in the final ROOT file
- Using Von Neumann rejection when possible
- Oriented nuclei linked $\beta$ decay correlations 
- Include $\gamma$-$\gamma$ angular correlation (HAL Id: hal-04964841)

### TODO 
- Four body decay for radiative correction
- Include matrix nuclear data for mixed decay
- Using NUDAT data
- Include R-Matrix 
- Include $\beta$-decay mixing ratio
- Inlcude $\gamma$ angular distribution (relevant for polarised nuclei)
- Decay scheme generator

## Authors

* **Leendert Hayen** - *KU Leuven*
* **Lecanuet Samuel** - *LP2iB*
* **Canovas Pablo** - **
* **Dumenil Victor** - *LPC Caen*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
