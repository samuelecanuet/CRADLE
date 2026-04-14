#include "CRADLE/DecayManager.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/DecayChannel.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/SpectrumGenerator.hh"
#include "CRADLE/ThreadPool.hh"
#include "CRADLE/ECShell.hh"
#include "CRADLE/RadiativeCorrections.hh"

#include <ROOT/TBufferMerger.hxx>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/RDataFrame.hxx>

// #include <boost/progress.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <future>

#include <iomanip>

template <typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B> &p)
{
  return std::pair<B, A>(p.second, p.first);
}

template <typename A, typename B>
std::multimap<B, A> flip_map(const std::map<A, B> &src)
{
  std::multimap<B, A> dst;
  std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                 flip_pair<A, B>);
  return dst;
}

namespace CRADLE
{
  using std::cout;
  using std::endl;
  using std::map;
  using std::pair;
  using std::string;
  using std::vector;

  DecayManager::~DecayManager()
  {
    for (map<const int, Particle *>::iterator it = registeredParticles.begin();
         it != registeredParticles.end(); ++it)
    {
      for (vector<DecayChannel *>::iterator it2 =
               (it->second)->GetDecayChannels().begin();
           it2 != (it->second)->GetDecayChannels().end(); ++it2)
      {
        delete *it2;
      }
      delete it->second;
    }
    registeredParticles.clear();

    for (map<const string, ChannelProperties>::iterator it = registeredChannelProperties.begin();
         it != registeredChannelProperties.end(); ++it)
    {
      delete it->second.distribution;
    }
  }

  void DecayManager::RegisterDecayMode(const string name, DecayMode &dm)
  {
    registeredDecayModes.insert(pair<string, DecayMode &>(name, dm));
    if (configOptions.general.Verbosity >= 2)
      Info("Registered decay mode " + name);
  }

  DecayMode &DecayManager::GetDecayMode(const string name)
  {
    if (registeredDecayModes.count(name) == 0)
    {
      throw std::invalid_argument("DecayMode " + name + " not registered. Aborting.");
    }
    return registeredDecayModes.at(name);
  }

  void DecayManager::RegisterParticle(Particle *p)
  {
    registeredParticles.insert(pair<int, Particle *>(p->GetPDG(), p));
    if (configOptions.general.Verbosity >= 2)
      Info("Registered particle " + p->GetName() + " with PDG code " + std::to_string(p->GetPDG()));
  }

  Particle *DecayManager::GetNewParticle(const int pdg, int Z, int A, bool temp)
  {
    if (registeredParticles.count(pdg) == 0)
    {
      GenerateNucleus(PDGtoName(pdg), Z, A);
    }
    Particle *p = new Particle(*(registeredParticles.at(pdg)));
    if (configOptions.general.Verbosity >= 2 && temp == false)
      Info("Generated new particle " + p->GetName() + " with PDG code " + std::to_string(pdg));
    return p;
  }

  void DecayManager::RegisterChannelPropreties(const string name, vector<vector<double>> *dist, double Max, int betaType, double j_i, double j_f, double j_m, double W_max_H, double W_max_VS, double PH)
  {
    ChannelProperties cp;
    cp.distribution = dist;
    cp.MAX_distribution = Max;
    cp.betaType = betaType;
    cp.j_i = j_i;
    cp.j_f = j_f;

    // only for gamma cascade
    cp.j_m = j_m;

    // only for radiative beta decay
    cp.W_max_H = W_max_H;
    cp.W_max_S = W_max_VS;
    cp.PH = PH;

    registeredChannelProperties.insert(pair<string, ChannelProperties>(name, cp));
    if (configOptions.general.Verbosity >= 2)
      Info(Form("Registered channel properties for %s with beta type %d, j_i = %.1f and j_f = %.1f", name.c_str(), betaType, j_i, j_f));
  }

  ChannelProperties DecayManager::GetChannelPropreties(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name);
  }

  int DecayManager::GetChannelBetaType(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).betaType;
  }

  double DecayManager::GetChannelJi(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).j_i;
  }

  double DecayManager::GetChannelJf(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).j_f;
  }

  double DecayManager::GetChannelJm(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).j_m;
  }

  vector<vector<double>> *DecayManager::GetChannelDistribution(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).distribution;
  }

  double DecayManager::GetChannelDistributionMax(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).MAX_distribution;
  }

  std::pair<double, double> DecayManager::GetChannelDistributionMaxs(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return std::make_pair(registeredChannelProperties.at(name).W_max_H, registeredChannelProperties.at(name).W_max_S);
  }

  double DecayManager::GetChannelPH(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).PH;
  }

  void DecayManager::SetChannelMf(const string name, double mf)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    registeredChannelProperties.at(name).mf = mf;
  }

  void DecayManager::SetChannelMgt(const string name, double mgt)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    registeredChannelProperties.at(name).mgt = mgt;
  }

  double DecayManager::GetChannelMf(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).mf;
  }

  double DecayManager::GetChannelMgt(const string name)
  {
    if (registeredChannelProperties.count(name) == 0)
    {
      throw std::invalid_argument("Channel properties not registered.");
    }
    return registeredChannelProperties.at(name).mgt;
  }

  void DecayManager::RegisterBasicParticles()
  {
    RegisterParticle(new Particle(NametoPDG("e-"), utilities::EMASSC2, -1, 0, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("e+"), utilities::EMASSC2, 1, 0, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("p"), utilities::PMASSC2, 1, 0, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("n"), utilities::NMASSC2, 0, 1, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("alpha"), utilities::ALPHAMASSC2, 2, 2, 0., 0.));
    RegisterParticle(new Particle(NametoPDG("enu"), 0., 0, 0, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("enubar"), 0., 0, 0, 0.5, 0.));
    RegisterParticle(new Particle(NametoPDG("gamma"), 0., 0, 0, 0., 0.));
  }

  void DecayManager::RegisterBasicDecayModes()
  {
    RegisterDecayMode("Beta", Beta::GetInstance());
    RegisterDecayMode("Beta_RC", BetaRadiative::GetInstance());
    RegisterDecayMode("ConversionElectron", ConversionElectron::GetInstance());
    RegisterDecayMode("Proton", Proton::GetInstance());
    RegisterDecayMode("Alpha", Alpha::GetInstance());
    RegisterDecayMode("Gamma", Gamma::GetInstance());
    RegisterDecayMode("IT", Gamma::GetInstance());
    RegisterDecayMode("EC", ElectronCapture::GetInstance());
  }

  void DecayManager::RegisterSpectrumGenerator(const string decayMode, SpectrumGenerator &sg)
  {
    try
    {
      DecayMode &dm = GetDecayMode(decayMode);
      dm.SetSpectrumGenerator(&sg);
      if (configOptions.general.Verbosity >= 2)
        Info("Registered " + decayMode + " Spectrum Generator " + (std::string)(typeid(sg).name()));
    }
    catch (const std::invalid_argument &e)
    {
      Warning("Cannot register" + (std::string)(typeid(sg).name()) + "spectrum generator. Decay mode " + decayMode + " not registered.");
    }
  }

  void DecayManager::RegisterBasicSpectrumGenerators()
  {
    RegisterSpectrumGenerator("Proton", DeltaSpectrumGenerator::GetInstance());
    RegisterSpectrumGenerator("Alpha", DeltaSpectrumGenerator::GetInstance());
    RegisterSpectrumGenerator("Gamma", DeltaSpectrumGenerator::GetInstance());
    RegisterSpectrumGenerator("IT", DeltaSpectrumGenerator::GetInstance());
    RegisterSpectrumGenerator("EC", DeltaSpectrumGenerator::GetInstance());
    RegisterSpectrumGenerator("ConversionElectron", DeltaSpectrumGenerator::GetInstance());
    if (configOptions.betaDecay.FermiFunction == "Advanced")
    {
      RegisterSpectrumGenerator("Beta", AdvancedBetaDecay::GetInstance());
      RegisterSpectrumGenerator("Beta_RC", AdvancedBetaDecay::GetInstance());
    }
    else
    {
      RegisterSpectrumGenerator("Beta", SimpleBetaDecay::GetInstance());
      RegisterSpectrumGenerator("Beta_RC", AdvancedBetaDecay::GetInstance());
    }
  }

  void DecayManager::ListRegisteredParticles()
  {
    cout << "--------------------------------------------------------\n";
    cout << " List of registered particles\n";
    cout << "--------------------------------------------------------\n\n";
    for (map<const int, Particle *>::iterator it = registeredParticles.begin();
         it != registeredParticles.end(); ++it)
    {
      it->second->ListInformation();
      cout << "\n";
    }
    cout << "--------------------------------------------------------\n\n"
         << endl;
  }

  void DecayManager::WriteDecayData(string filename, string type)
  {
    std::ifstream DataFile(filename.c_str());
    if (!DataFile.is_open())
    {
      Warning("Could not open file " + filename);
      return;
    }

    std::stringstream buffer;
    buffer << DataFile.rdbuf();
    std::string content = buffer.str();

    TObjString *stringObject_data = new TObjString(content.c_str());
    if (outputFile == nullptr)
      return;

    if (outputFile->GetDirectory("DecayData") == nullptr)
      outputFile->mkdir("DecayData");
    if (outputFile->GetDirectory(Form("DecayData/%s", type.c_str())) == nullptr)
      outputFile->mkdir(Form("DecayData/%s", type.c_str()));

    filename = filename.substr(filename.find_last_of("/\\") + 1);
    outputFile->cd(Form("DecayData/%s", type.c_str()));
    stringObject_data->Write(filename.c_str(), TObject::kOverwrite);
    outputFile->cd();
  }

  void DecayManager::WriteConfigData(string filename)
  {
    std::ifstream DataFile(filename.c_str());
    if (!DataFile.is_open())
    {
      Warning("Could not open file " + filename);
      return;
    }

    std::stringstream buffer;
    buffer << DataFile.rdbuf();
    std::string content = buffer.str();

    TObjString *stringObject_data = new TObjString(content.c_str());
    if (outputFile == nullptr)
      return;

    outputFile->mkdir("Configuration");
    outputFile->cd("Configuration");
    filename = filename.substr(filename.find_last_of("/\\") + 1);
    stringObject_data->Write(filename.c_str(), TObject::kOverwrite);
    outputFile->cd();
  }

  bool DecayManager::GenerateNucleus(string name, int Z, int A)
  {
    if (configOptions.general.Verbosity >= 2)
      Info("Generating nucleus " + name + " with Z = " + std::to_string(Z) + " and A = " + std::to_string(A));

    std::ostringstream filename;
    filename << configOptions.envOptions.Radiationdata;
    filename << "/z" << Z << ".a" << A;
    std::ifstream radDataFile((filename.str()).c_str());

    string line;
    double excitationEnergy = 0.;
    double lifetime;
    double atomicMass = utilities::GetAMEMass(configOptions.envOptions.AMEdata, Z, A);

    if (atomicMass == 0)
    {
      Warning("No AME data found for " + std::to_string(Z) + " " + std::to_string(A));
      atomicMass = utilities::GetApproximateMass(Z, A);
    }

    // cout << "Generating nucleus " << name << " with Z = " << Z << " and A = " << A << endl;

    Particle *p = new Particle(GetPDG(Z, A), atomicMass, Z, (A - Z), 0., 0);

    // cout << filename.str() << endl;

    while (getline(radDataFile, line))
    {
      // cout<<line<<endl;

      if (!line.compare(0, 1, "#"))
      {
        // Comment line
        continue;
      }
      else if (!line.compare(0, 1, "P"))
      {
        // Parent line
        std::istringstream iss(line);
        string p;
        string flag;

        iss >> p >> excitationEnergy >> flag >> lifetime;
        continue;
      }
      // cout << "Lifetime: " << lifetime << endl;
      string mode;
      double daughterExcitationEnergy = 0;
      double intensity = 0;
      double Q = 0;
      string modifier;
      string flag;

      std::istringstream iss(line);
      iss >> mode >> daughterExcitationEnergy >> flag >> intensity >> Q >> modifier;
      // cout << "Mode : " << mode << endl;
      // cout << "Daughter Energy :" << daughterExcitationEnergy << endl;
      // cout << "Flag :" << flag <<endl;
      // cout << "Intensity : " << intensity << endl;
      // cout << "Q : " << Q << endl;
      // cout << "Modifier : " << modifier << endl;
      // cout << "\n" <<endl;

      if (Q > 0.)
      {
        /*cout << "Adding DecayChannel " << mode << " Excitation Energy " <<
        excitationEnergy << " to " << daughterExcitationEnergy << endl;*/
        DecayChannel *dc;
        if (mode.find("shellEC") != string::npos)
        {
          int NbShell = ecshell::fNumberOfShells[Z];
          if (mode.find("K") != string::npos)
          {
            dc = new DecayChannel(mode, &GetDecayMode("EC"), Q - ecshell::GetBindingEnergy(Z, ecshell::K), intensity, lifetime, excitationEnergy,
                                  daughterExcitationEnergy);
          }
          else if (mode.find("L") != string::npos)
          {
            // L1
            dc = new DecayChannel(mode, &GetDecayMode("EC"), Q - ecshell::GetBindingEnergy(Z, ecshell::L1), intensity * ecshell::ProbabilityL1(Z), lifetime, excitationEnergy,
                                  daughterExcitationEnergy);

            // L2
            if (NbShell > 2)
            {
              dc = new DecayChannel(mode, &GetDecayMode("EC"), Q - ecshell::GetBindingEnergy(Z, ecshell::L2), intensity * (1 - ecshell::ProbabilityL1(Z)), lifetime, excitationEnergy,
                                    daughterExcitationEnergy);
            }
          }
          else
          {
            // M1
            dc = new DecayChannel(mode, &GetDecayMode("EC"), Q - ecshell::GetBindingEnergy(Z, ecshell::M1), intensity * ecshell::ProbabilityM1(Z), lifetime, excitationEnergy,
                                  daughterExcitationEnergy);

            // M2
            if (NbShell > 4)
            {
              dc = new DecayChannel(mode, &GetDecayMode("EC"), Q - ecshell::GetBindingEnergy(Z, ecshell::M2), intensity * (1 - ecshell::ProbabilityM1(Z)), lifetime, excitationEnergy,
                                    daughterExcitationEnergy);
            }
          }
        }
        else if (mode.find("Beta") != string::npos)
        {
          if (mode.find("Plus") != string::npos)
            Q = -Q;
          string mode = "Beta";
          if (configOptions.betaDecay.RadiativeCorrections)
            mode += "_RC";
          dc = new DecayChannel(mode, &GetDecayMode(mode), Q, intensity, lifetime, excitationEnergy,
                                daughterExcitationEnergy);
        }
        else
        {
          dc = new DecayChannel(mode, &GetDecayMode(mode), Q, intensity, lifetime, excitationEnergy,
                                daughterExcitationEnergy);
        }
        p->AddDecayChannel(dc);
      }
    }
    radDataFile.close();

    std::ostringstream gammaFileSS;
    gammaFileSS << configOptions.envOptions.Gammadata;
    gammaFileSS << "/z" << Z << ".a" << A;
    std::ifstream gammaDataFile(gammaFileSS.str().c_str());
    if (gammaDataFile.is_open())
    {
      while (getline(gammaDataFile, line))
      {
        int levelNr;
        double initEnergy, E;
        double intensity;
        double convIntensity;
        double kCoeff, lCoeff1, lCoeff2, lCoeff3, mCoeff1, mCoeff2, mCoeff3, mCoeff4, mCoeff5;
        double lifetime;
        string angMom;
        string polarity;
        string flag;

        int nGammas;

        std::istringstream iss(line);
        iss >> levelNr >> flag >> initEnergy >> lifetime >> angMom >> nGammas;

        double other_process_intensity = p->GetTotalIntensity(initEnergy);
        double feeding_intensity = 0.;
        // looking for decay feeding the level in registeredParticles
        for (map<const int, Particle *>::iterator it = registeredParticles.begin();
             it != registeredParticles.end(); ++it)
        {
          for (vector<DecayChannel *>::iterator it2 = (it->second)->GetDecayChannels().begin();
               it2 != (it->second)->GetDecayChannels().end(); ++it2)
          {
            if (std::abs(((*it2)->GetDaughterExcitationEnergy()) - initEnergy) < 1e-3)
            {
              feeding_intensity += (*it2)->GetIntensity();
            }
          }
        }

        double factor = feeding_intensity - other_process_intensity;

        for (int i = 0; i < nGammas; ++i)
        {
          getline(gammaDataFile, line);
          int daughterLevelNr;
          int multipolarity;
          double multipolarityMixing;
          std::istringstream issLevel(line);
          issLevel >> daughterLevelNr >> E >> intensity >> multipolarity >> multipolarityMixing >> convIntensity >> kCoeff >> lCoeff1 >> lCoeff2 >> lCoeff3 >> mCoeff1 >>
              mCoeff2 >> mCoeff3 >> mCoeff4 >> mCoeff5;

          intensity *= factor / 100.; // correcting to get the right gamma decay branching ratio

          // cout << "Adding gamma decay level " << initEnergy << " " << E << endl;
          if ((initEnergy - E) >= 0)
          {

            // Multipolarity
            std::pair<int, int> possibleMultipolarities;
            if (multipolarity > 100)
            {
              possibleMultipolarities.first = int(multipolarity / 100) / 2.;
              possibleMultipolarities.second = int(multipolarity - int(multipolarity / 100) * 100) / 2.;
            }
            else
            {
              possibleMultipolarities.first = int(multipolarity / 2.);
              possibleMultipolarities.second = 0;
            }
            //

            DecayChannel *dcGamma =
                new DecayChannel("Gamma", &GetDecayMode("Gamma"), E, intensity / (1. + convIntensity),
                                 lifetime, initEnergy, initEnergy - E, possibleMultipolarities, multipolarityMixing);
            p->AddDecayChannel(dcGamma);

            if (convIntensity == 0)
              continue;

            DecayChannel *dcConvK =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::K), intensity * convIntensity * (1. + convIntensity) * kCoeff,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvK);

            DecayChannel *dcConvL1 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::L1), intensity * convIntensity * (1. + convIntensity) * lCoeff1,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvL1);

            DecayChannel *dcConvL2 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::L2), intensity * convIntensity * (1. + convIntensity) * lCoeff2,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvL2);

            DecayChannel *dcConvL3 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::L3), intensity * convIntensity * (1. + convIntensity) * lCoeff3,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvL3);

            DecayChannel *dcConvM1 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::M1), intensity * convIntensity * (1. + convIntensity) * mCoeff1,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvM1);

            DecayChannel *dcConvM2 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::M2), intensity * convIntensity * (1. + convIntensity) * mCoeff2,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvM2);

            DecayChannel *dcConvM3 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::M3), intensity * convIntensity * (1. + convIntensity) * mCoeff3,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvM3);

            DecayChannel *dcConvM4 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::M4), intensity * convIntensity * (1. + convIntensity) * mCoeff4,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvM4);

            DecayChannel *dcConvM5 =
                new DecayChannel("ConversionElectron", &GetDecayMode("ConversionElectron"), E - ecshell::GetBindingEnergy(p->GetCharge(), ecshell::M5), intensity * convIntensity * (1. + convIntensity) * mCoeff5,
                                 lifetime, initEnergy, initEnergy - E);
            p->AddDecayChannel(dcConvM5);
          }

          else
          {
            Warning("Attempted to add gamma branch to a final state with negative excitation energy. Please check you are using the correct version of PhotonEvaporation.\nCurrent filename: " + gammaFileSS.str());
          }
        }
      }
      gammaDataFile.close();
    }
    RegisterParticle(p);

    if (configOptions.general.Verbosity >= 2)
      Info("Nucleus " + name + " generated with " + std::to_string(p->GetDecayChannels().size()) + " decay channels.");
    return true;
  }

  bool DecayManager::Initialise(std::string configFilename, int argc, const char **argv)
  {
    ConfigOptions configOptions = ParseOptions(configFilename, argc, argv);
    ConfigFilename = configFilename;
    return Initialise(configOptions);
  }

  bool DecayManager::Initialise(ConfigOptions _configOptions)
  {
    // cout << "Initialising..." << endl;
    configOptions = _configOptions;
    initStateName = configOptions.nuclear.Name;
    initStatePDG = GetPDG(configOptions.nuclear.Charge, configOptions.nuclear.Nucleons);
    initExcitationEn = configOptions.nuclear.Energy;
    outputName = configOptions.general.Output;
    NRTHREADS = configOptions.general.Threads;

    if (initStateName != "" && configOptions.nuclear.Nucleons > 0)
    {
      struct stat infoRD;
      struct stat infoG;
      int i = stat(
          configOptions.envOptions.Radiationdata.c_str(),
          &infoRD);
      int j = stat(
          configOptions.envOptions.Gammadata.c_str(),
          &infoG);

      if (i == 0 && j == 0 && S_ISDIR(infoRD.st_mode) && S_ISDIR(infoG.st_mode))
      {
        RegisterBasicParticles();
        RegisterBasicDecayModes();
        RegisterBasicSpectrumGenerators();
        return GenerateNucleus(initStateName, configOptions.nuclear.Charge, configOptions.nuclear.Nucleons);
      }
      else
      {
        Error("Data files not found. Set Radiationdata and Gammadata to their correct folders.");
        return false;
      }
    }
    else
    {
      Error("Initial nucleus is not defined.");
      return false;
    }
  }

  std::vector<ParticleData> DecayManager::GenerateEvent_ROOT(int eventNr, int verbosity)
  {
    double time = 0.;
    double checkTime = 0.;

    std::vector<ParticleData> vec;

    std::vector<Particle *> particleStack;
    Particle *ini = GetNewParticle(initStatePDG);
    ini->SetExcitationEnergy(initExcitationEn);

    // SET INITIAL KIONETIC ENERGY TO 10keV and the momentum only on the z axis
    // ini->SetKinEnergy(10);
    // ublas::vector<double> momentum = ini->GetMomentum();
    // momentum[3] = sqrt(ini->GetKinEnergy() * 2.0 * ini->GetMass());
    // ini->SetMomentum(momentum);
    /////////////////////////////////////////////////////////////////////////////

    if (this->configOptions.general.Verbosity >= 2)
      Start("## Event n°" + std::to_string(eventNr));

    particleStack.push_back(ini);
    while (!particleStack.empty())
    {
      double mom;
      ParticleData ParticleData_ini;

      Particle *p = particleStack.back();
      vector<Particle *> finalStates;
      double decayTime = p->GetDecayTime();
      bool filling = true;
      // cout << "\n Decaying particle " << p->GetName() << endl;
      // std::cout << eventNr << "\t" << subEventNr << std::endl;
      // std::cout << "     Time =\t" << time    if (particleDefinition) {
      //           << "CheckTime =\t" << checkTime << "\n "
      //           << "decayTime =\t" << decayTime << std::endl;
      if (verbosity == 0)
      {
        if (find(RestrictedParticleList.begin(), RestrictedParticleList.end(), p->GetPDG()) == RestrictedParticleList.end())
        {
          filling = false;
        }
      }

      if (filling)
      {
        ParticleData_ini.code = p->GetPDG();
        ParticleData_ini.time = time;
        ParticleData_ini.excitation_energy = p->GetExcitationEnergy();
        ParticleData_ini.kinetic_energy = p->GetKinEnergy();
        mom = sqrt(pow(p->GetMomentum()[1], 2) + pow(p->GetMomentum()[2], 2) + pow(p->GetMomentum()[3], 2));
        ParticleData_ini.p = p->GetMomentum()[0];
        ParticleData_ini.px = p->GetMomentum()[1] / mom;
        ParticleData_ini.py = p->GetMomentum()[2] / mom;
        ParticleData_ini.pz = p->GetMomentum()[3] / mom;
        vec.push_back(ParticleData_ini);
      }

      if ((time + decayTime) <= configOptions.cuts.Lifetime)
      {
        try
        {
          time += decayTime;
          finalStates = p->Decay();
        }
        catch (const std::invalid_argument &e)
        {
          std::cerr << "Decay Mode for particle " << p->GetName() << " not found. Aborting." << endl;
        }
      }
      delete particleStack.back();
      particleStack.pop_back();
      if (!finalStates.empty())
      {
        particleStack.insert(particleStack.end(), finalStates.begin(),
                             finalStates.end());
      }
    }
    return vec;
  }

  std::string DecayManager::GenerateEvent_TXT(int eventNr, int verbosity)
  {
    double time = 0.;
    double checkTime = 0.;
    int subEventNr = 0;
    int totSubEvents = 0;
    int totEvents = 0;
    std::ostringstream eventData;
    std::ostringstream subHeader;
    std::ostringstream subEventData;
    std::vector<Particle *> particleStack;
    Particle *ini = GetNewParticle(initStatePDG);
    ini->SetExcitationEnergy(initExcitationEn);

    // SET INITIAL KIONETIC ENERGY TO 10keV and the momentuml only on the z axis
    ini->SetKinEnergy(0.0);
    ublas::vector<double> momentum = ini->GetMomentum();
    momentum[3] = sqrt(ini->GetKinEnergy() * 2.0 * ini->GetMass());
    ini->SetMomentum(momentum);
    /////////////////////////////////////////////////////////////////////////////

    particleStack.push_back(ini);

    if (verbosity == 0)
    {
      while (!particleStack.empty())
      {
        Particle *p = particleStack.back();
        vector<Particle *> finalStates;
        double decayTime = p->GetDecayTime();
        // cout << "\n Decaying particle " << p->GetRawName() << endl;
        //  std::cout << eventNr << "\t" << subEventNr << std::endl;
        //  std::cout << "     Time =\t" << time      << "\n "
        //   << "CheckTime =\t" << checkTime << "\n "
        //   << "decayTime =\t" << decayTime << std::endl;

        if (decayTime > configOptions.cuts.Lifetime && (find(RestrictedParticleList.begin(), RestrictedParticleList.end(), p->GetPDG()) != RestrictedParticleList.end()))
        {
          ++totSubEvents;
          subEventData << eventNr << "\t\t" << std::fixed << std::setprecision(4) << roundf(time * 10000) / 10000. << "\t" << p->GetInfoForFile() << "\n";
        }
        if ((time + decayTime) <= configOptions.cuts.Lifetime)
        {
          try
          {
            time += decayTime;
            finalStates = p->Decay();
          }
          catch (const std::invalid_argument &e)
          {
            std::cout << "Decay Mode for particle " << p->GetName() << " not found. Aborting." << endl;
            return "";
          }
        }

        delete particleStack.back();
        particleStack.pop_back();
        if (!finalStates.empty())
        {
          particleStack.insert(particleStack.end(), finalStates.begin(),
                               finalStates.end());
        }
      }
    }
    else
    {
      while (!particleStack.empty())
      {
        Particle *p = particleStack.back();
        vector<Particle *> finalStates;
        double decayTime = p->GetDecayTime();
        // cout << "\n Decaying particle " << p->GetRawName() << endl;
        //  std::cout << eventNr << "\t" << subEventNr << std::endl;
        //  std::cout << "     Time =\t" << time      << "\n "
        //   << "CheckTime =\t" << checkTime << "\n "
        //   << "decayTime =\t" << decayTime << std::endl;

        ++totSubEvents;
        subEventData << eventNr << "\t\t" << std::fixed << std::setprecision(4) << roundf(time * 10000) / 10000. << "\t" << p->GetInfoForFile() << "\n";

        if ((time + decayTime) <= configOptions.cuts.Lifetime)
        {
          try
          {
            time += decayTime;
            finalStates = p->Decay();
            subHeader << eventNr << std::setw(8) << subEventNr << "\t\t" << totSubEvents << "\n"
                      << subEventData.str();
            subEventData.str(std::string());
            totEvents += totSubEvents;
            totSubEvents = 0;
            ++subEventNr;
          }
          catch (const std::invalid_argument &e)
          {
            std::cout << "Decay Mode for particle " << p->GetName() << " not found. Aborting." << endl;
            return "";
          }
        }

        delete particleStack.back();
        particleStack.pop_back();
        if (!finalStates.empty())
        {
          particleStack.insert(particleStack.end(), finalStates.begin(),
                               finalStates.end());
        }
      }
    }
    totEvents += totSubEvents;
    subHeader << eventNr << std::setw(8) << subEventNr << "\t\t" << totSubEvents << "\n"
              << subEventData.str();
    eventData << eventNr << "\t\t" << totEvents << "\n"
              << subHeader.str();

    return eventData.str();
  }

  bool DecayManager::MainLoop()
  {
    int nrParticles = configOptions.general.Loop;
    int verbosity = configOptions.general.Verbosity_file;
    if (nrParticles < 1)
    {
      Error("ERROR: Incorrect number of events (" + std::to_string(nrParticles) + ")");
      return false;
    }
    Start("Generating " + std::to_string(nrParticles) + " events...");

    std::ios::sync_with_stdio(false);
    // int show_progress = 0;
    int steppingProgress = std::max(1, nrParticles / 10000);
    clock_t start = clock();

    if (outputName.find("root") != std::string::npos)
    {
      // outputFile = new TFile(outputName.c_str(), "RECREATE");
      // if (!outputFile->IsOpen())
      // {
      //   Error("Could not open output file " + outputName);
      //   return false;
      // }
      // TTree *tree = new TTree("ParticleTree", "Tree for Particle Data");
      // if (!tree)
      // {
      //   Error("Could not create TTree in output file " + outputName);
      //   return false;
      // }
      // vector<double> Time;
      // vector<int> Code;
      // vector<double> Kinetic_energy;
      // vector<double> Excitation_energy;
      // vector<double> p;
      // vector<double> Px;
      // vector<double> Py;
      // vector<double> Pz;

      // tree->Branch("time", &Time);
      // tree->Branch("code", &Code);
      // tree->Branch("energy", &Kinetic_energy);
      // tree->Branch("excitation_energy", &Excitation_energy);
      // tree->Branch("p", &p);
      // tree->Branch("px", &Px);
      // tree->Branch("py", &Py);
      // tree->Branch("pz", &Pz);

      //////// OPTIMIZATION WITH THREADPOOL //////////////
      // ThreadPool pool(NRTHREADS);
      // std::mutex result_mutex;

      // for (int i = 0; i < nrParticles; ++i)
      // {
      //   pool.enqueue([&, i]
      //                {
      //     auto particles = this->GenerateEvent_ROOT(i, verbosity);

      //     std::lock_guard<std::mutex> lock(result_mutex);

      //     for (const auto &particle : particles) {
      //         Time.push_back(particle.time);
      //         Code.push_back(particle.code);
      //         Kinetic_energy.push_back(particle.kinetic_energy);
      //         Excitation_energy.push_back(particle.excitation_energy);
      //         p.push_back(particle.p);
      //         Px.push_back(particle.px);
      //         Py.push_back(particle.py);
      //         Pz.push_back(particle.pz);
      //     }
      //     tree->Fill();
      //     Time.clear(); Code.clear(); Kinetic_energy.clear();
      //     Excitation_energy.clear(); p.clear();
      //     Px.clear(); Py.clear(); Pz.clear();
      //     show_progress++;
      //     ProgressBar(show_progress, nrParticles, start, "", steppingProgress, NRTHREADS);
      //     });
      // }
      // pool.wait_all();

      // ...

      ROOT::EnableThreadSafety();
      // ROOT::EnableImplicitMT(NRTHREADS); // not needed with TThreadExecutor

      std::atomic<int> show_progress{0};

      {
        // Keep merger inside a scope so it is destroyed before reopening the file
        ROOT::TBufferMerger merger(outputName.c_str(), "RECREATE");
        ROOT::TThreadExecutor executor(NRTHREADS);

        const int chunkSize = 1000000;
        std::vector<std::pair<int, int>> ranges;
        ranges.reserve((nrParticles + chunkSize - 1) / chunkSize);

        for (int first = 0; first < nrParticles; first += chunkSize)
        {
          ranges.emplace_back(first, std::min(first + chunkSize, nrParticles));
        }

        std::atomic<int> treeCounter{0};

        executor.Foreach([&](const std::pair<int, int> &range)
                         {
        auto file = merger.GetFile();
        file->cd();

        // Unique tree name per task/chunk
        const int treeId = treeCounter++;
        const std::string treeName = "ParticleTree_" + std::to_string(treeId);

        TTree tree(treeName.c_str(), "ParticleTree");
        tree.SetDirectory(file.get());

        std::vector<double> Time, Kinetic_energy, Excitation_energy, p, Px, Py, Pz;
        std::vector<int> Code;

        tree.Branch("time",              &Time);
        tree.Branch("code",              &Code);
        tree.Branch("energy",    &Kinetic_energy);
        tree.Branch("excitation_energy", &Excitation_energy);
        tree.Branch("p",                 &p);
        tree.Branch("px",                &Px);
        tree.Branch("py",                &Py);
        tree.Branch("pz",                &Pz);

        tree.SetAutoFlush(-20 * 1024 * 1024);

        std::vector<ParticleData> particles;
        particles.reserve(16);

        Time.reserve(16);
        Code.reserve(16);
        Kinetic_energy.reserve(16);
        Excitation_energy.reserve(16);
        p.reserve(16);
        Px.reserve(16);
        Py.reserve(16);
        Pz.reserve(16);

        for (int i = range.first; i < range.second; ++i) {
            particles = this->GenerateEvent_ROOT(i, verbosity);
            const std::size_t n = particles.size();

            Time.clear();
            Code.clear();
            Kinetic_energy.clear();
            Excitation_energy.clear();
            p.clear();
            Px.clear();
            Py.clear();
            Pz.clear();

            if (Time.capacity() < n)              Time.reserve(n);
            if (Code.capacity() < n)              Code.reserve(n);
            if (Kinetic_energy.capacity() < n)    Kinetic_energy.reserve(n);
            if (Excitation_energy.capacity() < n) Excitation_energy.reserve(n);
            if (p.capacity() < n)                 p.reserve(n);
            if (Px.capacity() < n)                Px.reserve(n);
            if (Py.capacity() < n)                Py.reserve(n);
            if (Pz.capacity() < n)                Pz.reserve(n);

            for (const auto& particle : particles) {
                Time.push_back(particle.time);
                Code.push_back(particle.code);
                Kinetic_energy.push_back(particle.kinetic_energy);
                Excitation_energy.push_back(particle.excitation_energy);
                p.push_back(particle.p);
                Px.push_back(particle.px);
                Py.push_back(particle.py);
                Pz.push_back(particle.pz);
            }

            tree.Fill();

            const int done = ++show_progress;
            if (done % steppingProgress == 0 || done == nrParticles) {
                ProgressBar(done, nrParticles, start, "", steppingProgress, NRTHREADS);
            }
        }

        file->cd();
        tree.Write("", TObject::kOverwrite);
        file->Write(); }, ranges);

      }

      MergeParticleTreesInPlace(outputName);
      
      outputFile = new TFile(outputName.c_str(), "UPDATE");
      if (!outputFile->IsOpen())      {
        Error("Could not open output file " + outputName);
        return false;
      }
      // Writting Data File
      for (const auto &p : registeredParticles)
      {
        if (p.second->GetDecayChannels().size() > 0)
        {
          WriteDecayData(configOptions.envOptions.Radiationdata + "/z" + std::to_string(p.second->GetCharge()) + ".a" + std::to_string(p.second->GetNeutrons()+p.second->GetCharge()), "Radiation");
          WriteDecayData(configOptions.envOptions.Gammadata + "/z" + std::to_string(p.second->GetCharge()) + ".a" + std::to_string(p.second->GetNeutrons()+p.second->GetCharge()), "Gamma");
        }
      }
      // Writting config file
      WriteConfigData(ConfigFilename);
      outputFile->Close();

      ///////////////////////////////////////////////////
    }

    else if (outputName.find("txt") != std::string::npos)
    {
      int show_progress = 0;
      std::ofstream fileStream;
      fileStream.open(outputName.c_str());
      for (int i = 0; i < nrParticles; i += NRTHREADS)
      {
        // cout << "LOOP NR " << i+1 << endl;
        const int threads = std::min(NRTHREADS, nrParticles - i);
        std::vector<std::future<std::string>> f(threads);
        for (int t = 0; t < threads; t++)
        {
          f[t] = std::async(std::launch::async, &DecayManager::GenerateEvent_TXT, this, i + t, verbosity);
        }

        for (int t = 0; t < threads; t++)
        {
          fileStream << f[t].get();
          show_progress++;
          ProgressBar(show_progress, nrParticles, start, "", steppingProgress, NRTHREADS);
        }
      }
      fileStream.flush();
      fileStream.close();
    }
    else
    {
      Error("Choose .txt or .root for your output file");
      return false;
    }

    Success(Form("Done! Generated in %.1f seconds.", (double)(clock() - start) / CLOCKS_PER_SEC / NRTHREADS));
    return true;
  }

} // End of CRADLE namespace