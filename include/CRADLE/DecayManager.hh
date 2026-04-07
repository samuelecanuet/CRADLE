#ifndef DECAYMANAGER
#define DECAYMANAGER

#include <vector>
#include <map>
#include <string>
#include <random>

#include "CRADLE/ConfigParser.hh"
#include "CRADLE/Messenger.hh"
#include "CRADLE/PDGcode.hh"

#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"

namespace CRADLE {

class Particle;
class DecayMode;
class SpectrumGenerator;

struct ParticleData {
    int code = 0;
    double time = 0.;
    double excitation_energy = 0.;
    double kinetic_energy = 0.;
    double p = 0.;
    double px = 0.; 
    double py = 0.;
    double pz = 0.;
};

struct ChannelProperties {
    std::vector<std::vector<double>>* distribution;
    double MAX_distribution = 0.;
    int betaType; // 0: Fermi, 1: Gamow-Teller, 2: Mixed 
    double j_i = 0.;
    double j_f = 0.;
};

class DecayManager {
  public:
    static DecayManager& GetInstance() {
      static DecayManager instance;
      return instance;
    }

    ~DecayManager();

    bool Initialise(std::string, int argc = 0, const char** argv = nullptr);
    bool Initialise(ConfigOptions);
    bool Initialise(std::string, int, int, double, std::string, int);
    bool MainLoop();
    bool GenerateNucleus(std::string, int, int);
    void RegisterBasicParticles();
    void RegisterBasicDecayModes();
    void RegisterDecayMode(const std::string, DecayMode&);
    void RegisterParticle(Particle*);
    
    void RegisterChannelPropreties(const std::string, std::vector<std::vector<double> >*, double, int, double, double);
    ChannelProperties GetChannelPropreties(const std::string);
    std::vector<std::vector<double> >* GetChannelDistribution(const std::string);
    double GetChannelDistributionMax(const std::string);
    int GetChannelBetaType(const std::string);
    double GetChannelJi(const std::string);
    double GetChannelJf(const std::string);

    
    void RegisterSpectrumGenerator(const std::string, SpectrumGenerator&);
    void RegisterBasicSpectrumGenerators();
    void ListRegisteredParticles();
    std::vector<ParticleData> GenerateEvent_ROOT(int, int);
    std::string GenerateEvent_TXT(int, int);
    Particle* GetNewParticle(const int, int Z=0, int A=0, bool temp = false);
    DecayMode& GetDecayMode(const std::string);
    ConfigOptions configOptions;
    std::mt19937 generator;

    void WriteDecayData(std::string, std::string);
    void WriteConfigData(std::string);

  private:
    DecayManager() {};
    DecayManager(DecayManager const&);
    void operator=(DecayManager const&);

    std::map<const std::string, ChannelProperties> registeredChannelProperties;
    std::map<const std::string, DecayMode&> registeredDecayModes;
    std::vector<Particle*> particleStack;
    std::map<const int, Particle*> registeredParticles;
    std::string outputName;
    std::string ConfigFilename;
    std::string initStateName;
    int initStatePDG;
    double initExcitationEn;
    int NRTHREADS;

    TFile *outputFile;
};

}//End of CRADLE namespace
#endif
