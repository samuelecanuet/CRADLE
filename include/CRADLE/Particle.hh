#ifndef PARTICLE
#define PARTICLE

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <string>
#include <sstream>
#include <random>
#include "CRADLE/Messenger.hh"
#include "CRADLE/PDGcode.hh"

namespace CRADLE
{

  class DecayChannel;

  namespace ublas = boost::numeric::ublas;

  class Particle
  {
  private:
    double mass;
    int charge;
    int neutrons;
    double spin;
    double currentExcitationEnergy;
    double currentDeltaExcitationEnergy;
    int PDG;

    std::string name;

    ublas::vector<double> fourMomentum;
    std::vector<DecayChannel *> decayChannels;

    static std::default_random_engine randomGen;

    double LevelEnergyUncertainty = 2; // keV, threshold for considering two levels as degenerate

  public:
    Particle();
    Particle(const std::string &, double, int, int, double, double);
    Particle(const int &, double, int, int, double, double);
    Particle(const Particle &);
    Particle &operator=(const Particle &rhs);

    inline void ListInformation() const
    {
      Info(GetName());
      Info("Mass: " + std::to_string(GetMass()) + " keV", 1);
      Info("Charge: " + std::to_string(GetCharge()), 1);
      Info("Neutrons: " + std::to_string(GetNeutrons()), 1);
      Info("Spin: " + std::to_string(spin), 1);
      Info("Excitation energy: " + std::to_string(GetExcitationEnergy()) + " keV", 1);
      Info("G.S. Lifetime: " + std::to_string(GetLifetime()) + " s", 1);
    }

    std::vector<Particle *> Decay();
    double GetLifetime() const;
    double GetDecayTime();
    double GetMixingRatio(double, double) const;
    std::pair<int, int> GetMultipolarities(double, double) const;

    inline std::string GetName() const
    {
      std::string n;
      std::ostringstream oss(n);
      oss << name << " (" << currentExcitationEnergy << " keV) Kin. En.: " << GetKinEnergy();
      return oss.str();
    };
    inline std::string GetRawName() const { return name; };
    inline int GetPDG() const { return PDG; };
    inline std::string GetInfoForFile() const
    {
      std::string n;
      std::ostringstream oss(n);
      oss << name << "\t" << currentExcitationEnergy << "\t" << GetKinEnergy() << "\t" << fourMomentum(0) << "\t" << fourMomentum(1) << "\t" << fourMomentum(2) << "\t" << fourMomentum(3);
      return oss.str();
    };

    inline void SetMomentum(ublas::vector<double> v) { fourMomentum = v; };
    inline void SetExcitationEnergy(double e)
    {
      currentExcitationEnergy = e;
      fourMomentum(0) += e;
    };
    inline void SetDeltaExcitationEnergy(double e)
    {
      currentDeltaExcitationEnergy = e;
      fourMomentum(0) += e;
    };
    
    inline int GetCharge() const { return charge; };
    inline int GetNeutrons() const { return neutrons; };
    inline int GetNucleons() const { return charge + neutrons; };
    inline double GetMass() const { return mass + currentExcitationEnergy; };
    inline ublas::vector<double> GetMomentum() const { return fourMomentum; };
    inline double GetExcitationEnergy() const { return currentExcitationEnergy; };
    inline void AddDecayChannel(DecayChannel *dc) { decayChannels.push_back(dc); };
    ublas::vector<double> GetVelocity() const;
    inline double GetKinEnergy() const { return fourMomentum(0) - GetMass(); };
    inline void SetKinEnergy(double e) { fourMomentum(0) = GetMass() + e; };
    inline ublas::vector<double> Get3Momentum() const
    {
      ublas::vector<double> v(3);
      v[0] = fourMomentum[1];
      v[1] = fourMomentum[2];
      v[2] = fourMomentum[3];
      return v;
    };
    inline std::vector<DecayChannel *> &GetDecayChannels() { return decayChannels; };
    double GetTotalIntensity(double ) const; 

    Particle* lastGamma = nullptr;
    inline void SetLastGamma(Particle* g) { delete lastGamma; lastGamma = g; };
    inline Particle* GetLastGamma() const { return lastGamma; };
  };

} // End of CRADLE namespace

#endif
