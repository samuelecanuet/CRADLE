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

* Verbosity : 0 or 1 to see more or less information in the terminal.
* Verbosity_file : 0 to save only particles in the ouput file, 1 to save also the nucleus kinematics.
* *C<sub>i* : set all the coupling constant of the Lee Yang Lagragian.
* *a*, *b* : if not NaN overwritting the *C<sub>i*
* Cuts : distance, time and energy limit of calculation
* Beta Decay : *Fermi* or *Gamow-Teller* to impose $\beta$ decay type. *Auto* take into account real $\beta$ decay type extracted from Geant4 *GammaData*, only available for pure transition (set on Gamow-Teller if none pure transition)
* FermiFunction : Simple or Advanced

### Run 
Running command example for 32Ar : 
```bash
./CRADLE++ nucleus --name 32Ar -Z 18 -A 32 -c [ConfigFileName] general -l [EventNumber] -o [OuputFileName]
```
Note that you can generate .txt file or .root file (TTree)

## Details
See *Generator_Documentation.pdf* for more details about the generator.

## Version 
### v1.0 
- All feature working
### v2.0
- New file writter, possibility to save .root or .txt file
- Auto mode for $\beta$ decay type
### v2.1
- Quick implementation of nuclear level width (Breit-Wigner) for proton decay (taking into account only the width of the initial level.)
### v2.2
- Implementation of Internal Conversion 
- Implementation of Electron Capture

### TODO 
- Four body decay for radiative correction
- Include matrix nuclear data for mixed decay
- Using NUDAT data
- Polarized nuclei
- Include R-Matrix 

## Authors

* **Leendert Hayen** - *KU Leuven*
* **Lecanuet Samuel** - *LP2iB*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
