# CRADLE++

Customisable RAdioactive Decay for Low Energy Particle Physics: A C++ event generator

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

* [BOOST](https://www.boost.org/) - Linear algebra routines, and the `program_options` part of the library
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
* Verbosity_file: 0 to save only Geant4 releveant particles in the output file, 1 to save all the generated particles.
* *C<sub>i*: set all the coupling constant of the Lee Yang Lagragian.
* *a*, *b*, *A*, *B*, *D*, *c*: if not NaN overwritting the calculated parameters with *C<sub>i*
* Cuts: distance, time and energy limit of calculation
* Beta Decay: *Fermi* or *Gamow-Teller* to impose $\beta$ decay type. *Auto* take into account real $\beta$ decay type deduced from $J^{\pi}$ state included in Geant4 *GammaData*, only available for pure transition (set on Gamow-Teller if Mixed transition)
* FermiFunction: Simple or Advanced
* BetaSpectrumCorrections: *true*/*false* (if *false* only Fermi function and phase space factor are used for the $\beta$ spectrum shape else all the correction of [Rev. Mod. Phys. 90, 015008 (2018)](https://doi.org/10.1103/RevModPhys.90.015008) are included)
* Alignement: setting its value
* Polarisation: setting its value and direction
* InFlightDecay: *true*/*false* (if *false* the kinetic energy of the recoil will be set to 0.)
* NuclearLevelWidth: *true*/*false* (if *true* Breit-Wigner distribution will be used for proton and $\alpha$ decay if a lifetime is given in the RadiationData file)
* GammaGammaCorrelation: *true*/*false* (if *true* the angular correlation between 2 following $\gamma$ in a cascade is calculated)

### Run 
Running command example for $^{32}\mathrm{Ar}$ : 
```bash
./CRADLE++ nucleus --n 32Ar -Z 18 -A 32 general -l [NumberOfEvents] -o [OuputFileName] -c [ConfigFileName] -t [NumberOfThread]
```
Note that you can generate .txt file or .root file (TTree) by specify the format of the ouput filename.

## OUTPUT
### ROOT
In the case of a ROOT file, input files as the Radioactive/Evaporation data and the config file will be saved using a TObjString.
Additionnaly, a TTree is created with the following branches:
- time (Creation time of the particle, *vector\<double>*)
- code (PDG code of the particle, *vector\<int>*)
- energy (Kinetic energy of the particle, *vector\<double>*)
- excitation_energy (Excitation energy of the particle, *vector\<double>*)
- p (Momentum of the particle, *vector\<double>*)
- px (x component of the direction, *vector\<double>*)
- py (y component of the direction, *vector\<double>*)
- pz (z component of the direction, *vector\<double>*)

Each element of the vector correspond to a particle constituting the event. The number of particle in the event is given by the size of the vector.

### TXT
Exemple (Verbosity_file = 2)
```
0		6
0       0		1
0		0.0000	32Ar	0	0	2.98056e+07	0	0	0
0       1		3
0		0.0138	enu	0	1637.93	1637.93	849.525	-966.08	1013.81
0		0.0138	e+	0	4449.64	4960.63	938.992	3849.79	2940.11
0		0.0138	32Cl	5046.3	0.455511	2.97995e+07	-1788.52	-2883.71	-3953.91
0       2		2
0		0.0138	p	0	3358.21	941630	-14713.3	-68411.4	37636.6
0		0.0138	31S	0	107.264	2.88574e+07	12924.8	65527.7	-41590.5
1		6
1       0		1
1		0.0000	32Ar	0	0	2.98056e+07	0	0	0
1       1		3
1		0.0138	enu	0	1668.53	1668.53	251.982	987.288	-1321.27
1		0.0138	e+	0	4419.23	4930.23	531.898	533.121	4845.51
1		0.0138	32Cl	5046.3	0.257493	2.97995e+07	-783.881	-1520.41	-3524.24
1       2		2
1		0.0138	p	0	3357.21	941629	-42398	60299.5	-29624.5
1		0.0138	31S	0	108.026	2.88574e+07	41614.1	-61819.9	26100.3
```
The first row correspond to event number. Second is either total number of particle in the event, sub event number  or time of creation. Third is the total number of particle in the sub event or the particle name. Following row correspond to excitation energy, kinetic energy, momentum, momentum vector componant in this order.

Exemple (Verbosity_file = 1) 
```
0		0.0138	e+	0	4449.64	4960.63	938.992	3849.79	2940.11
0		0.0138	p	0	3358.21	941630	-14713.3	-68411.4	37636.6
1		0.0138	e+	0	4419.23	4930.23	531.898	533.121	4845.51
1		0.0138	p	0	3357.21	941629	-42398	60299.5	-29624.5
```

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
- Oriented nuclei linked to $\beta$ decay correlations 
- Include $\gamma-\gamma$ angular correlation [T. Lauritsen et al.](https://hal.science/hal-04964841v1)
### v2.5
- Fixing errors
- Adding 4-body decay for radiative correction ([F.Glück 1998](https://doi.org/10.1016/S0010-4655(96)00168-3))
- Adding a unit test folder to crosscheck with Litt.
- Include $\beta$-decay mixing ratio, data coming from [Phys. Rev. C 107, 015502](https://doi.org/10.1103/PhysRevC.107.015502)
- TTree Threaded filling (one TTree per thread) to improve the speed of the writing process in multi-threading mode.

### TODO 
- Using NUDAT data
- Include R-Matrix input
- Include $\gamma$ angular distribution (relevant for polarised nuclei)
- Decay scheme generator
- Include some Matrix element calculation for Gamow-Teller shape factor correction.

## Authors

* **Leendert Hayen** - *KU Leuven*
* **Lecanuet Samuel** - *LP2iB*
* **Canovas Pablo** - **
* **Dumenil Victor** - *LPC Caen*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
