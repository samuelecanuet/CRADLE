#ifndef DECAYCHANNEL
#define DECAYCHANNEL

#include <vector>
#include <string>
#include "CRADLE/Messenger.hh"

namespace CRADLE {

class Particle;
class DecayMode;

class DecayChannel {
  public:
    DecayChannel(std::string, DecayMode*, double, double, double, double, double);
    DecayChannel(std::string, DecayMode*, double, double, double, double, double, std::pair<int, int>, double);

    inline double GetIntensity() { return intensity; };
    inline double GetQValue() { return Q; };
    inline double GetLifetime() { return lifetime; };
    inline double GetDaughterExcitationEnergy() { return daughterExcitationEnergy; };
    inline double GetParentExcitationEnergy() { return parentExcitationEnergy; };
    inline std::string GetModeName() { return modeName; };
    inline std::pair<int, int> GetMultipolarities() { return Multipolarities; };
    inline double GetMixingRatio() { return mixingRatio; };

    std::vector<Particle*> Decay(Particle*);
  private:
    double daughterExcitationEnergy;
    double parentExcitationEnergy;
    double intensity;
    double Q;
    double lifetime;
    std::pair<int, int> Multipolarities;
    double mixingRatio;
    std::string modeName;
    DecayMode* decayMode;
};

}//End of CRADLE namespace
#endif
