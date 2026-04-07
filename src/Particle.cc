#include "CRADLE/Particle.hh"
#include "CRADLE/DecayChannel.hh"

#include <stdlib.h>
#include <stdexcept>
#include <cmath>

namespace CRADLE {

std::default_random_engine Particle::randomGen;

Particle::Particle(const std::string& _name, double _mass, int _charge, int _neutrons, double _spin, double _excitationEnergy): name(_name), mass(_mass), charge(_charge), neutrons(_neutrons), spin(_spin), currentExcitationEnergy(_excitationEnergy) {
  //std::cout << "Creating new particle " << name << std::endl;
  Warning("Creating new particle with string name " + name);
  PDG = NametoPDG(name);
  fourMomentum.resize(4);
  fourMomentum(0) = mass + currentExcitationEnergy;
}

Particle::Particle(const int &_PDG, double _mass, int _charge, int _neutrons, double _spin, double _excitationEnergy): PDG(_PDG), mass(_mass), charge(_charge), neutrons(_neutrons), spin(_spin), currentExcitationEnergy(_excitationEnergy) {
  //std::cout << "Creating new particle with PDG code " << PDG << std::endl;
  name = PDGtoName(PDG);
  fourMomentum.resize(4);
  fourMomentum(0) = mass + currentExcitationEnergy;
}

Particle::Particle(const Particle& orig) {
  //std::cout << "Copy constructor for particle " << orig.name << std::endl;
  fourMomentum.resize(4);
  name = orig.name;
  mass = orig.mass;
  charge = orig.charge;
  neutrons = orig.neutrons;
  spin = orig.spin;
  currentExcitationEnergy = orig.currentExcitationEnergy;
  currentDeltaExcitationEnergy = 0;
  fourMomentum(0) = mass + currentExcitationEnergy;

  PDG = orig.PDG;

  //std::cout << charge << " " << neutrons << std::endl;

  decayChannels.insert(decayChannels.end(), orig.decayChannels.begin(), orig.decayChannels.end());
}

double Particle::GetLifetime() const {
  double t = 1.e46;
  for(int i = 0; i < decayChannels.size(); ++i) {
    // Look for decay channels from current excitation state
    if (std::abs(decayChannels[i]->GetParentExcitationEnergy()-currentExcitationEnergy) < LevelEnergyUncertainty) {
      t = decayChannels[i]->GetLifetime();
       break;
    }
  }
  return t;
}

double Particle::GetDecayTime() {
  double lifetime = GetLifetime();
  //std::cout << name << " Lifetime " << lifetime << std::endl;
  if (lifetime!= 1.e46) {
    std::exponential_distribution<double> distribution(1./lifetime);
    double decayTime = distribution(randomGen);
    return decayTime;
  }
  return lifetime;
}

double Particle::GetMixingRatio(double InitExcistationEnergy, double FinalExcitationEnergy) const {
  for(int i = 0; i < decayChannels.size(); ++i) {
    if (std::abs(decayChannels[i]->GetParentExcitationEnergy()-InitExcistationEnergy) < LevelEnergyUncertainty && std::abs(decayChannels[i]->GetDaughterExcitationEnergy()-FinalExcitationEnergy) < LevelEnergyUncertainty) {
      return decayChannels[i]->GetMixingRatio();
    }
  }
  Warning(Form("Mixing ratio not found for transition from %.1f keV to %.1f keV for particle %s. Returning 0.", InitExcistationEnergy, FinalExcitationEnergy, name.c_str()));
  return 0.;
}

std::pair<int, int> Particle::GetMultipolarities(double InitExcistationEnergy, double FinalExcitationEnergy) const {
  for(int i = 0; i < decayChannels.size(); ++i) {
    if (std::abs(decayChannels[i]->GetParentExcitationEnergy()-InitExcistationEnergy) < LevelEnergyUncertainty && std::abs(decayChannels[i]->GetDaughterExcitationEnergy()-FinalExcitationEnergy) < LevelEnergyUncertainty) {
      return decayChannels[i]->GetMultipolarities();
    }
  }
  Warning(Form("Multipolarities not found for transition from %.1f keV to %.1f keV for particle %s. Returning (0, 0).", InitExcistationEnergy, FinalExcitationEnergy, name.c_str()));
  return std::make_pair<int, int>(0, 0);
}

ublas::vector<double> Particle::GetVelocity() const {
  ublas::vector<double> velocity = Get3Momentum()/fourMomentum(0);
  return velocity;
}

std::vector<Particle*> Particle::Decay() {
  double totalIntensity = 0.;
  for(std::vector<DecayChannel*>::size_type i = 0; i != decayChannels.size(); i++) {
    // Look for decay channels at current excitation level
    if (std::abs(currentExcitationEnergy - decayChannels[i]->GetParentExcitationEnergy()) < LevelEnergyUncertainty) {
      totalIntensity+=decayChannels[i]->GetIntensity();
    }
  }
  double r = rand()/(double)RAND_MAX*totalIntensity;
  double intensity = 0.;
  double index = 0.;
  // Sample randomly from the different decay channels
  for(std::vector<DecayChannel*>::size_type i = 0; i != decayChannels.size(); i++) {
    if (std::abs(currentExcitationEnergy - decayChannels[i]->GetParentExcitationEnergy()) < LevelEnergyUncertainty) {
      intensity+=decayChannels[i]->GetIntensity();
      if (r <= intensity) {
        break;
      }
    }
    index++;
  }
  try {
    return decayChannels[index]->Decay(this);
  }
  catch (const std::invalid_argument& e) {
    Warning(Form("Error in decay channel %d for particle %s. Check the decay channel properties and the decay mode implementation.", (int)index, name.c_str()));
    throw e;
  }
}

double Particle::GetTotalIntensity(double excitationEnergy) const {
  double intensity = 0.;
  for(std::vector<DecayChannel*>::size_type i = 0; i != decayChannels.size(); i++) {
    if (std::abs(excitationEnergy - decayChannels[i]->GetParentExcitationEnergy()) < LevelEnergyUncertainty) {
      intensity+=decayChannels[i]->GetIntensity();
    }
  }
  return intensity;
}

}//End of CRADLE namespace
