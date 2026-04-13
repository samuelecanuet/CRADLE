#ifndef DECAYMODE
#define DECAYMODE

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <string>
#include "CRADLE/Messenger.hh"

namespace CRADLE {

class Particle;
class SpectrumGenerator;

namespace ublas = boost::numeric::ublas;

class DecayMode {
  public:
    virtual std::vector<Particle*> Decay(Particle*, double, double) = 0;
    DecayMode();
    virtual ~DecayMode();

    void SetSpectrumGenerator(SpectrumGenerator*);

  protected:
    static void FourBodyDecay(ublas::vector<double>&, Particle*, Particle*, Particle*, Particle*);
    static void ThreeBodyDecay(ublas::vector<double>&, Particle*, Particle*, Particle*, ublas::vector<double>&, double);
    static void TwoBodyDecay(ublas::vector<double>&, Particle*, Particle*, double);
    static void TwoBodyDecay(ublas::vector<double>&, Particle*, Particle*, double, ublas::vector<double> dir);
    SpectrumGenerator* spectrumGen;
};

class Beta: public DecayMode {
  public:
    static Beta& GetInstance() {
      static Beta instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Beta();
    Beta(Beta const& copy);
    Beta& operator=(Beta const& copy);
};

class BetaRadiative: public DecayMode {
  public:
    static BetaRadiative& GetInstance() {
      static BetaRadiative instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaRadiative();
    BetaRadiative(BetaRadiative const& copy);
    BetaRadiative& operator=(BetaRadiative const& copy);
};

class ConversionElectron: public DecayMode {
    public:
    static ConversionElectron& GetInstance() {
      static ConversionElectron instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    ConversionElectron();
    ConversionElectron(ConversionElectron const& copy);
    ConversionElectron& operator=(ConversionElectron const& copy);
};

class Proton: public DecayMode {
    public:
    static Proton& GetInstance() {
      static Proton instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Proton();
    Proton(Proton const& copy);
    Proton& operator=(Proton const& copy);
};

class Alpha: public DecayMode {
    public:
    static Alpha& GetInstance() {
      static Alpha instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Alpha();
    Alpha(Alpha const& copy);
    Alpha& operator=(Alpha const& copy);
};

class Gamma: public DecayMode {
    public:
    static Gamma& GetInstance() {
      static Gamma instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Gamma();
    Gamma(Gamma const& copy);
    Gamma& operator=(Gamma const& copy);
};

class ElectronCapture: public DecayMode {
    public:
    static ElectronCapture& GetInstance() {
      static ElectronCapture instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    ElectronCapture();
    ElectronCapture(ElectronCapture const& copy);
    ElectronCapture& operator=(ElectronCapture const& copy);
};

}//End of CRADLE namespace
#endif
