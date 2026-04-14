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
#include "TKey.h"

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

    // bool MergeParticleTreesInPlace(const std::string& filename,
    //                            const std::string& inputPrefix = "ParticleTree_",
    //                            const std::string& mergedTreeName = "ParticleTree");

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


inline bool MergeParticleTreesInPlace(const std::string& filename,
                               const std::string& inputPrefix = "ParticleTree_",
                               const std::string& mergedTreeName = "ParticleTree")
{
  DecayManager& dm = DecayManager::GetInstance();
    Start("Merging Thread TTree");

    TFile* file = TFile::Open(filename.c_str(), "UPDATE");
    if (!file || file->IsZombie()) {
        Error("Error: output ROOT file is corrupted or not finalized: " + filename);
        if (file) {
            file->Close();
            delete file;
        }
        return false;
    }

    std::vector<std::string> treeNames;

    TIter next(file->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = dynamic_cast<TKey*>(next()))) {
        // Do not instantiate objects just to check the class
        if (std::string(key->GetClassName()) != "TTree")
            continue;

        const std::string name = key->GetName();
        if (name.rfind(inputPrefix, 0) == 0)
            treeNames.push_back(name);
    }

    if (treeNames.empty()) {
        Error("No trees found with prefix " + inputPrefix + " in file " + filename);
        file->Close();
        delete file;
        return false;
    }

    std::sort(treeNames.begin(), treeNames.end());

    // Remove previous merged tree if present
    if (file->GetListOfKeys()->FindObject(mergedTreeName.c_str())) {
        file->Delete((mergedTreeName + ";*").c_str());
        file->SaveSelf(kTRUE);
    }

    TTree* firstTree = nullptr;
    file->GetObject(treeNames.front().c_str(), firstTree);
    if (!firstTree) {
        Warning("Cannot read first tree: " + treeNames.front() + ". Aborting merge.");
        file->Close();
        delete file;
        return false;
    }

    file->cd();

    // IMPORTANT: raw pointer, ROOT/file owns it after SetDirectory
    TTree* mergedTree = firstTree->CloneTree(0);
    if (!mergedTree) {
        Warning("Failed to clone tree structure for merging. Aborting merge.");
        file->Close();
        delete file;
        return false;
    }

    mergedTree->SetName(mergedTreeName.c_str());
    mergedTree->SetTitle(mergedTreeName.c_str());
    mergedTree->SetDirectory(file);

    Long64_t totalEntries = 0;

    for (const auto& name : treeNames) {
        TTree* inTree = nullptr;
        file->GetObject(name.c_str(), inTree);
        if (!inTree) {
            Warning("Cannot read tree: " + name + ". Skipping.");
            continue;
        }

        const Long64_t copied = mergedTree->CopyEntries(inTree, -1, "fast");
        totalEntries += copied;

        if (dm.configOptions.general.Verbosity >= 2) {
            Info("Merged tree " + name + " with " + std::to_string(inTree->GetEntries()) + " entries.", 1);
        }
    }

    file->cd();
    mergedTree->Write("", TObject::kOverwrite);
    file->SaveSelf(kTRUE);
    file->Flush();

    if (dm.configOptions.general.Verbosity >= 2) {
        Info("Merged tree " + mergedTreeName + " written with " +
             std::to_string(mergedTree->GetEntries()) + " entries.", 1);
    }

    // Delete original thread trees only after merged tree is safely written
    for (const auto& name : treeNames) {
        file->Delete((name + ";*").c_str());
        if (dm.configOptions.general.Verbosity >= 2) {
            Info("Deleted original tree " + name + " from file after merging.", 1);
        }
    }

    file->SaveSelf(kTRUE);
    file->Write("", TObject::kOverwrite);
    file->Flush();
    file->Close();
    delete file;

    Success("TTree has been merged successfully");
    return true;
}

}//End of CRADLE namespace
#endif
