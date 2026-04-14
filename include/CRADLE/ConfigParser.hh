#ifndef CRADLE_CONFIG_CONTAINER_H
#define CRADLE_CONFIG_CONTAINER_H

#include "CLI11.hpp"
#include "CRADLE/Messenger.hh"
#include <string>
#include <complex>

namespace CRADLE {

struct NuclearOptions {
  std::string Name = "";
  int Charge = 1;
  int Nucleons = 0;
  double Energy = 0.0;
  double WeakMagnetism = 0.0;
  double Alignment = 0.0;
  double PolarisationMag = 0.0;
  double PolarisationX = 0.0;
  double PolarisationY = 0.0;
  double PolarisationZ = 0.0;
};

struct General {
  int Verbosity = 1;
  int Verbosity_file = 1;
  int Loop = 0;
  int Threads = 5;
  std::string Output = "output.txt";
};

struct CouplingConstants {
  std::complex<double> CS = std::complex<double>(0.0,0.0);
  std::complex<double> CSP = std::complex<double>(0.0,0.0);
  std::complex<double> CV = std::complex<double>(1.0,0.0);
  std::complex<double> CVP = std::complex<double>(1.0,0.0);
  std::complex<double> CT = std::complex<double>(0.0,0.0);
  std::complex<double> CTP = std::complex<double>(0.0,0.0);
  std::complex<double> CA = std::complex<double>(1.2754,0.0);
  std::complex<double> CAP = std::complex<double>(1.2754,0.0);

  double iCS = 0.0;
  double iCSP = 0.0;
  double iCV = 0.0;
  double iCVP = 0.0;
  double iCT = 0.0;
  double iCTP = 0.0;
  double iCA = 0.0;
  double iCAP = 0.0;
  
  double b = std::nan("");
  double a = std::nan("");
  double A = std::nan("");
  double B = std::nan("");
  double D = std::nan("");
  double c = std::nan("");
};

struct Cuts {
  double Distance = 1.E10;
  double Lifetime = 1.E10;
  double Energy = 1.E10;
};

struct BetaDecay {
  std::string Default = "Auto";
  std::string FermiFunction = "Advanced";   
  bool BetaSpectrumCorrections = true;
  bool RadiativeCorrections = true;
  double Cs = 1e-3;
};

struct Decay{
  bool InFlightDecay = true;
  bool NuclearLevelWidth = true;
  bool GammaGammaCorrelation = true;
};

struct EnvOptions {
  std::string AMEdata;
  std::string Gammadata;
  std::string Radiationdata;
  std::string BetaMixingRatios;
};

struct ConfigOptions{
  NuclearOptions nuclear;
  General general;
  CouplingConstants couplingConstants;
  Cuts cuts;
  BetaDecay betaDecay;
  Decay decay;
  EnvOptions envOptions;
};

ConfigOptions ParseOptions(std::string, int argc = 0, const char** argv = nullptr);

void SetGeneralOptions(CLI::App&, General&);
void SetNuclearOptions(CLI::App&, NuclearOptions&);
void SetCouplingConstants(CLI::App&, CouplingConstants&);
void SetCuts(CLI::App&, Cuts&);
void SetBetaDecayOptions(CLI::App&, BetaDecay&);
void SetDecayOptions(CLI::App&, Decay&);
void SetEnvironmentOptions(CLI::App&, EnvOptions&);
void PrintingAllOptions(const ConfigOptions&);

}

#endif
