#include "CRADLE/DecayChannel.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/Particle.hh"

#include <iostream>
#include <stdexcept>

namespace CRADLE {

  DecayChannel::DecayChannel(std::string md, DecayMode* dm, double q, double i, double t, double pExEn, double dExEn):
    modeName(md), decayMode(dm), Q(q), intensity(i), lifetime(t), parentExcitationEnergy(pExEn), daughterExcitationEnergy(dExEn){
  }

  // DecayChannel::DecayChannel(std::string md, DecayMode* dm, double q, double i, double t, double pExEn, double dExEn, double ph):
  //   modeName(md), decayMode(dm), Q(q), intensity(i), lifetime(t), parentExcitationEnergy(pExEn), daughterExcitationEnergy(dExEn), PH(ph){
  // }

  DecayChannel::DecayChannel(std::string md, DecayMode* dm, double q, double i, double t, double pExEn, double dExEn, std::pair<int, int> multipolarities, double mixingratio):
    modeName(md), decayMode(dm), Q(q), intensity(i), lifetime(t), parentExcitationEnergy(pExEn), daughterExcitationEnergy(dExEn), Multipolarities(multipolarities), mixingRatio(mixingratio){
  }

  std::vector<Particle*> DecayChannel::Decay (Particle* initState) {
    try {
      DecayManager &dm = DecayManager::GetInstance();
      if (dm.configOptions.general.Verbosity >= 2)
        Info(Form("Performing decay %s with Q = %.1f, intensity = %.1f, lifetime = %.1f, parent excitation energy = %.1f and daughter excitation energy = %.1f", modeName.c_str(), Q, intensity, lifetime, parentExcitationEnergy, daughterExcitationEnergy));
      
      if (!dm.configOptions.decay.InFlightDecay) 
      {
        initState->SetMomentum(ublas::vector<double>(4, 0.));
        initState->SetKinEnergy(0.);
      }
      
      return decayMode->Decay(initState, Q, daughterExcitationEnergy);
    }
    catch (const std::invalid_argument &e) {
      std::cout << "Forwarding exception up the stack" << std::endl;
      throw e;
    }
  }
}//End of CRADLE namespace
