#include "CRADLE/DecayMode.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/SpectrumGenerator.hh"

#include <string>
#include <sstream>


namespace CRADLE {

void DecayMode::ThreeBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, ublas::vector<double>& dir2, double Q) {
  //Perform decay in CoM frame
  ublas::vector<double> momentum1 = finalState1->GetMomentum();
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);

  ublas::vector<double> p2 (3);
  double p2Norm = 0.;

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  ublas::vector<double> p1 = finalState1->Get3Momentum();

  double a = mass2*mass2;
  double b = mass3*mass3;
  double c = utilities::GetNorm(p1);
  double d = Q + mass1 + mass2 + mass3 - momentum1(0);
  double e = inner_prod(p1, dir2)/c;

  double first = 1./2./(c*c*e*e-d*d);
  double second = a*a*d*d-2*a*b*d*d+4.*a*c*c*d*d*e*e-2.*a*c*c*d*d-2.*a*d*d*d*d+b*b*d*d+2*b*c*c*d*d-2.*b*d*d*d*d+c*c*c*c*d*d-2.*c*c*d*d*d*d+d*d*d*d*d*d;
  double third = a*c*e-b*c*e-c*c*c*e+c*d*d*e;

  p2Norm = first*(-std::sqrt(second)+third);
  p2 = p2Norm*dir2;

  ublas::vector<double> p3 = -(p1+p2);
  double p3Norm = utilities::GetNorm(p3);

  momentum2(0) = std::sqrt(a+p2Norm*p2Norm);
  momentum2(1) = p2(0);
  momentum2(2) = p2(1);
  momentum2(3) = p2(2);

  momentum3(0) = std::sqrt(b+p3Norm*p3Norm);
  momentum3(1) = p3(0);
  momentum3(2) = p3(1);
  momentum3(3) = p3(2);

  //std::cout << "\t" << inner_prod(p1, p2)/p2Norm/c << std::endl;

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));

}


void DecayMode::TwoBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, double Q, ublas::vector<double>& dir) {
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();

  double M = Q + mass1 + mass2;

  double p = 1./(2.*M)*std::sqrt((M*M-std::pow(mass1-mass2, 2.))*(M*M-std::pow(mass1+mass2, 2.)));

  double energy1 = std::sqrt(mass1*mass1+p*p);
  double energy2 = std::sqrt(mass2*mass2+p*p);

  momentum1(0) = energy1;
  momentum1(1) = -p*dir[0];
  momentum1(2) = -p*dir[1];
  momentum1(3) = -p*dir[2];

  momentum2(0) = energy2;
  momentum2(1) = p*dir[0];
  momentum2(2) = p*dir[1];
  momentum2(3) = p*dir[2];

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
}

void DecayMode::TwoBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, double Q) {
  ublas::vector<double> dir = utilities::RandomDirection();
  TwoBodyDecay(velocity, finalState1, finalState2, Q, dir);
}

std::vector<Particle*> Beta::Decay(Particle* initState, double Q, double daughterExEn) {
  
  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();

  if (dm.configOptions.general.Verbosity >= 2)
    Info(Form("Beta Decay with Q = %.1f, parent excitation energy = %.1f and daughter excitation energy = %.1f", Q, initState->GetExcitationEnergy(), daughterExEn), 1);

  // Getting Beta +/- from sign of Q
  int BetaSign = 0;
  if (Q < 0)
    BetaSign = 1;
  else
    BetaSign = -1;
  Q = abs(Q);

  // Modifying the endpoint energy for Beta+ decay
  double E0 = Q;
  if (BetaSign == 1) { // Beta+ case
    E0 -= 2*utilities::EMASSC2;
  }

  // Init Final States
  int Recoil_PDG = GetPDG(initState->GetCharge()-BetaSign, initState->GetCharge()+initState->GetNeutrons());
  Particle* Recoil = dm.GetNewParticle(Recoil_PDG, initState->GetCharge()-BetaSign, initState->GetCharge()+initState->GetNeutrons());
  Recoil->SetExcitationEnergy(daughterExEn);
  Particle* ChargedLepton = dm.GetNewParticle(- BetaSign * 11);
  Particle* NeutralLepton = dm.GetNewParticle(BetaSign * 12);
  int Recoil_Z = Recoil->GetCharge();
  double Recoil_ExEn = Recoil->GetExcitationEnergy();

  // name as key for channel properties
  std::ostringstream oss;
  oss << "Beta:" << "Sign" << BetaSign << "Z" << Recoil_Z << "A" << Recoil_Z + Recoil->GetNeutrons() << "Q" << Q;

  // Build or Get info of the channel using oss
  double mf = 1.;
  double mgt = 0.;

  std::vector<std::vector<double> >* dist;
  double dist_max = 0.;
  double j_i;
  double j_f;
  int Type;
  try
  {
    dist = dm.GetChannelDistribution(oss.str());
    dist_max = dm.GetChannelDistributionMax(oss.str());
    j_i = dm.GetChannelJi(oss.str());
    j_f = dm.GetChannelJf(oss.str());
    Type = dm.GetChannelBetaType(oss.str());
  }
  catch (const std::invalid_argument &e)
  {
    dist = spectrumGen->GenerateSpectrum(initState, Recoil, E0, Type);
    double b = correlation::CalculateFierz(mf, mgt, initState->GetCharge(), BetaSign);
    for (int i = 0; i < dist->size(); i++)
    {
      double E = ((*dist)[i])[0] + utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH * (1 + b * utilities::EMASSC2 / E + (-BetaSign) * 4. / 3. * E / (initState->GetMass()) * dm.configOptions.nuclear.WeakMagnetism);
    }
    dist_max = utilities::CalculateMax(*dist);
    j_i = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy());
    j_f = utilities::GetJpi(Recoil_Z + Recoil->GetNeutrons(), Recoil_Z, Recoil->GetExcitationEnergy());
    Type = utilities::FindBetaType(initState, Recoil);
    dm.RegisterChannelPropreties(oss.str(), dist, dist_max, Type, j_i, j_f);

    if (Type == FERMI)
    {
      mf = 1.;
      mgt = 0.;
    }
    else if (Type == GAMOW_TELLER)
    {
      mf = 0.;
      mgt = 1;
    }
    else 
    {
      Error("MIXED IS NOT HANDLE FOR NOW");
    }

  }

  // Angle correlation
  double ChargedLepton_Energy = utilities::RandomFromDistribution(*dist, dist_max) + utilities::EMASSC2;
  ublas::vector<double> ChargedLepton_FourMomentum(4);
  double ChargedLepton_Momentum = std::sqrt(ChargedLepton_Energy*ChargedLepton_Energy-std::pow(utilities::EMASSC2, 2.));
  ublas::vector<double> NeutralLepton_Dir;
  ublas::vector<double> ChargedLepton_Dir;
  
  if (dm.configOptions.nuclear.Alignment == 0 && dm.configOptions.nuclear.PolarisationMag == 0)
  {
    ///// IF THE NUCLEUS IS NOT POLARISED ////
    NeutralLepton_Dir = utilities::RandomDirection();
    double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, ChargedLepton_Energy, Recoil_Z, BetaSign);
    double b = correlation::CalculateFierz(mf, mgt, Recoil_Z, BetaSign);
    ChargedLepton_Dir = utilities::GetParticleDirection(NeutralLepton_Dir, 1 + b * utilities::EMASSC2/ChargedLepton_Energy, a * ChargedLepton_Momentum / ChargedLepton_Energy);
  }
  else
  {
    ///// IF THE NUCLEUS IS POLARISED ////
    double b = correlation::CalculateFierz(mf, mgt, Recoil_Z, BetaSign);
    double align = dm.configOptions.nuclear.Alignment;
    ublas::vector<double> polDir(3);
    polDir(0) = dm.configOptions.nuclear.PolarisationX;
    polDir(1) = dm.configOptions.nuclear.PolarisationY;
    polDir(2) = dm.configOptions.nuclear.PolarisationZ;
    polDir = utilities::NormaliseVector(polDir);
    double polMag = dm.configOptions.nuclear.PolarisationMag;

    double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, ChargedLepton_Energy, Recoil_Z, BetaSign);
    double c = 0;
    double A = 0;
    double B = 0;
    double D = 0;

    if (j_i > 0)
    {
      A = correlation::CalculateBetaAssymetry(mf, mgt, j_i, j_f, BetaSign, Recoil_Z, ChargedLepton_Energy);
      B = correlation::CalculateNeutrinoAssymetry(mf, mgt, j_i, j_f, BetaSign, Recoil_Z, ChargedLepton_Energy);
      D = correlation::CalculateDTripleCorrelation(mf, mgt, j_i, j_f, BetaSign, Recoil_Z, ChargedLepton_Energy);
      if (j_i > 0.5)
      {
        c = correlation::CalculateAlignmentCorrelation(mf, mgt, j_i, j_f, BetaSign, Recoil_Z, ChargedLepton_Energy);
      }
    }
    c *= align * (-1);
    A *= polMag;
    B *= polMag;
    D *= polMag;

    double F_max = correlation::AnalyticalMaximumAngCorrFactor(a, b, c, A, B, D, ChargedLepton_Energy);
    double F_point = 0;
    std::uniform_real_distribution<double> distribution_F(0.0, F_max);
    double F = F_max;
    while (F > F_point)
    {
      F = distribution_F(dm.generator);
      ChargedLepton_Dir = utilities::RandomDirection();
      NeutralLepton_Dir = utilities::RandomDirection();
      F_point = correlation::CalculateAngularCorrelationFactor(a, b, c, A, B, D, ChargedLepton_Energy, ChargedLepton_Dir, NeutralLepton_Dir, polDir);
    }
  }

  // Setting Charged Lepton 4-momentum
  ChargedLepton_FourMomentum(0) = ChargedLepton_Energy;
  ChargedLepton_FourMomentum(1) = ChargedLepton_Momentum*ChargedLepton_Dir[0];
  ChargedLepton_FourMomentum(2) = ChargedLepton_Momentum*ChargedLepton_Dir[1];
  ChargedLepton_FourMomentum(3) = ChargedLepton_Momentum*ChargedLepton_Dir[2];
  ChargedLepton->SetMomentum(ChargedLepton_FourMomentum);

  // 3-body decay kinematics
  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, ChargedLepton, NeutralLepton, Recoil, NeutralLepton_Dir, E0);

  // Adding final states to vector
  finalStates.push_back(Recoil);
  finalStates.push_back(ChargedLepton);
  finalStates.push_back(NeutralLepton);
  
 
  return finalStates;
}
  

std::vector<Particle*> ConversionElectron::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  ublas::vector<double> velocity = -initState->GetVelocity();
  Particle* Recoil =  DecayManager::GetInstance().GetNewParticle(initState->GetPDG());
  Particle* e =  DecayManager::GetInstance().GetNewParticle(11);
  Recoil->SetExcitationEnergy(daughterExEn);

  TwoBodyDecay(velocity, Recoil, e, Q);

  finalStates.push_back(Recoil);
  finalStates.push_back(e);

  return finalStates;
}

std::vector<Particle*> Proton::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();
  //// nuclear level width
  if (dm.configOptions.decay.Nuclear_Level_Width) 
    Q = utilities::BreitWigner(Q, initState->GetLifetime());
  
  int daughter_PDG = GetPDG(initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons()-1);
  Particle* Recoil = dm.GetNewParticle(daughter_PDG, initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons()-1);
  Recoil->SetExcitationEnergy(daughterExEn);
  Particle* p = dm.GetNewParticle(2212);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, Recoil, p, Q);

  finalStates.push_back(Recoil);
  finalStates.push_back(p);

  return finalStates;
}

std::vector<Particle*> Alpha::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();
  //// nuclear level width
  if (dm.configOptions.decay.Nuclear_Level_Width) 
    Q = utilities::BreitWigner(Q, initState->GetLifetime());

  int daughter_PDG = GetPDG(initState->GetCharge()-2, initState->GetCharge()+initState->GetNeutrons()-4);
  Particle* Recoil = dm.GetNewParticle(daughter_PDG, initState->GetCharge()-2, initState->GetCharge()+initState->GetNeutrons()-4);
  Particle* alpha = dm.GetNewParticle(1000020040);
  Recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, Recoil, alpha, Q);

  finalStates.push_back(Recoil);
  finalStates.push_back(alpha);

  return finalStates;
}

std::vector<Particle*> Gamma::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  DecayManager& dm = DecayManager::GetInstance();
  ublas::vector<double> velocity = -initState->GetVelocity();
  Particle* Recoil = DecayManager::GetInstance().GetNewParticle(initState->GetPDG());
  Particle* gamma = DecayManager::GetInstance().GetNewParticle(22);
  Recoil->SetExcitationEnergy(daughterExEn);

  if (initState->GetLastGamma() != nullptr && dm.configOptions.decay.GammaGammaCorrelation)
  {
    // Gamma - Gamma correlation (E_i --> E --> E_f)
    Particle *gamma_1 = initState->GetLastGamma();
    ublas::vector<double> gamma_1_dir = gamma_1->Get3Momentum();

    double j_i = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy() + gamma_1->GetKinEnergy());
    double j = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy());
    double j_f = utilities::GetJpi(Recoil->GetCharge() + Recoil->GetNeutrons(), Recoil->GetCharge(), Recoil->GetExcitationEnergy());
    double delta_1 = initState->GetMixingRatio(initState->GetExcitationEnergy() + gamma_1->GetKinEnergy(), initState->GetExcitationEnergy());
    std::pair<int, int> l1 = initState->GetMultipolarities(initState->GetExcitationEnergy() + gamma_1->GetKinEnergy(), initState->GetExcitationEnergy());
    double delta_2 = Recoil->GetMixingRatio(initState->GetExcitationEnergy(), daughterExEn);
    std::pair<int, int> l2 = Recoil->GetMultipolarities(initState->GetExcitationEnergy(), daughterExEn);

    // Info("State i:");
    // Info(Form("E = %.1f keV", initState->GetExcitationEnergy() + gamma_1->GetKinEnergy()), 1);
    // Info(Form("J = %.1f", j_i), 1);
    // Info("State m:");
    // Info(Form("E = %.1f keV", initState->GetExcitationEnergy()), 1);
    // Info(Form("J = %.1f", j), 1);
    // Info("State f:");
    // Info(Form("E = %.1f keV", Recoil->GetExcitationEnergy()), 1);
    // Info(Form("J = %.1f", j_f), 1);
    // Info("Gamma 1:");
    // Info(Form("E = %.1f keV", gamma_1->GetKinEnergy()), 1);
    // Info(Form("Multipolarity = %d, %d", l1.first, l1.second), 1);
    // Info(Form("Mixing Ratio = %.4f", delta_1), 1);
    // Info("Gamma 2:");
    // Info(Form("E = %.1f keV", gamma->GetKinEnergy()), 1);
    // Info(Form("Multipolarity = %d, %d", l2.first, l2.second), 1);
    // Info(Form("Mixing Ratio = %.4f", delta_2), 1);

    // test 0->2->0 (ok)
    // j_i = 0.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(2, 0);
    // delta_1 = 0.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 1->2->0 (close)
    // j_i = 1.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(1, 2);
    // delta_1 = 5.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 3->2->0 (ok)
    // j_i = 3.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(1, 0);
    // delta_1 = 0.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 2->2->0 (ok)
    // j_i = 2.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(1, 2);
    // delta_1 = 0.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 2->2->0 (ok)
    // j_i = 2.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(1, 2);
    // delta_1 = -1.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 60Co 4->2->0
    // j_i = 4.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(2, 0);
    // delta_1 = 0.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;
    //
    // 88Y 3->2->0
    // j_i = 3.;
    // j = 2.;
    // j_f = 0.;
    // l1 = std::make_pair<int, int>(1, 0);
    // delta_1 = 0.;
    // l2 = std::make_pair<int, int>(2, 0);
    // delta_2 = 0.;

    // coef 
    std::vector<double> ak = correlation::CaluclateGammaCoefficient_a(std::abs(j_i), std::abs(j), std::abs(j_f), l1, delta_1, l2, delta_2);

    // Info(Form("Gamma-Gamma correlation coefficients: a0 = %.4f, a2 = %.4f, a4 = %.4f", ak[0], ak[1], ak[2]), 1);
    double W_max = correlation::MaxAnalyticalGammaCorrelation(ak);
    double W = W_max;
    double W_point = 0.;
    std::uniform_real_distribution<double> distribution(0.0, W_max);
    ublas::vector<double> gamma_2_dir;
    while (W > W_point)
    {
      W = distribution(dm.generator);
      gamma_2_dir = utilities::RandomDirection();
      double cos_theta = inner_prod(gamma_1_dir, gamma_2_dir) / utilities::GetNorm(gamma_1_dir) / utilities::GetNorm(gamma_2_dir);
      W_point = correlation::GammaGammaW(cos_theta, ak);
    }

    TwoBodyDecay(velocity, Recoil, gamma, Q, gamma_2_dir);
  }
  else
    TwoBodyDecay(velocity, Recoil, gamma, Q);

  // Saving the last gamma with a temporary particle pointer to let the particlestack deletable
  Particle *temp_gamma = DecayManager::GetInstance().GetNewParticle(22, 0, 0, true);
  temp_gamma->SetMomentum(gamma->GetMomentum());
  Recoil->SetLastGamma(temp_gamma);
  //
  finalStates.push_back(Recoil);
  finalStates.push_back(gamma);

  return finalStates;
}

std::vector<Particle*> ElectronCapture::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  int daughter_PDG = GetPDG(initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  Particle* Recoil = DecayManager::GetInstance().GetNewParticle(daughter_PDG, initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  Particle* enu = DecayManager::GetInstance().GetNewParticle(12);
  Recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, Recoil, enu, Q);

  finalStates.push_back(Recoil);
  finalStates.push_back(enu);

  return finalStates;
}

DecayMode::DecayMode() { }

DecayMode::~DecayMode() { }

void DecayMode::SetSpectrumGenerator(SpectrumGenerator* sg) {
  spectrumGen = sg;
}

Beta::Beta() { }

ConversionElectron::ConversionElectron() { }

Proton::Proton () { }

Alpha::Alpha () { }

Gamma::Gamma () { }

ElectronCapture::ElectronCapture () { }

}//End of CRADLE namespace
