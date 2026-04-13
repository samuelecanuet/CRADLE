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
    std::vector<std::vector<double>>* distribution; // 2D vector for the distribution of the channel
    double MAX_distribution = 0.; // Maximum of the distribution for rejection sampling
    int betaType; // 0: Fermi, 1: Gamow-Teller, 2: Mixed 
    double j_i = 0.; // Initial State Spin
    double j_f = 0.; // Final State Spin
    double j_m = 0.; // Intermediate State Spin (for gamma-gamma correlation)

    // In the case of 4body decay 
    double W_max_H = 0.; // Maximum of the distribution for the Hard Bremsstrahlung
    double W_max_S = 0.; // Maximum of the distribution for the Virtual/Soft Bremsstrahlung
    double PH = 0.; // Probability of Hard Bremsstrahlung
    //

    // Matrix element for Beta decay
    double mf = -1;
    double mgt = -1;
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
    
    void RegisterChannelPropreties(const std::string, std::vector<std::vector<double> >*, double, int, double, double, double = 0, double = 0., double = 0., double = 0.);
    ChannelProperties GetChannelPropreties(const std::string);
    std::vector<std::vector<double> >* GetChannelDistribution(const std::string);
    double GetChannelDistributionMax(const std::string);
    int GetChannelBetaType(const std::string);
    double GetChannelJi(const std::string);
    double GetChannelJf(const std::string);
    double GetChannelJm(const std::string);
    std::pair<double, double> GetChannelDistributionMaxs(const std::string);
    double GetChannelPH(const std::string);
    void SetChannelMf(const std::string, double);
    double GetChannelMf(const std::string);
    void SetChannelMgt(const std::string, double);
    double GetChannelMgt(const std::string);
    
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

    std::map<const std::string, ChannelProperties> registeredChannelProperties;

  private:
    DecayManager() {};
    DecayManager(DecayManager const&);
    void operator=(DecayManager const&);

    
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
