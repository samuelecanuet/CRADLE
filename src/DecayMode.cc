#include "CRADLE/DecayMode.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/SpectrumGenerator.hh"
#include "CRADLE/RadiativeCorrections.hh"

#include <string>
#include <sstream>


namespace CRADLE {

void DecayMode::FourBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, Particle* finalState4) //, ublas::vector<double>& dir1, ublas::vector<double>& dir2, ublas::vector<double>& dirg, double EnergyElectron, double EnergyNeutrino, double EnergyBrPhoton) 
{
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);
  ublas::vector<double> momentumg (4);

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  double massg = finalState4->GetMass();

  ublas::vector<double> p1 = finalState1->Get3Momentum();
  ublas::vector<double> p2 = finalState2->Get3Momentum();
  ublas::vector<double> pg = finalState4->Get3Momentum();
  ublas::vector<double> p3 = -(p1+p2+pg) ;

  double p1Norm = utilities::GetNorm(p1);
  double p2Norm = utilities::GetNorm(p2);
  double p3Norm = utilities::GetNorm(p3);
  double pgNorm = utilities::GetNorm(pg);

  momentum1(0) = std::sqrt(mass1*mass1 + std::pow(p1Norm, 2)) ;
  momentum1(1) = p1[0] ;
  momentum1(2) = p1[1] ;
  momentum1(3) = p1[2] ;

  momentum2(0) = std::sqrt(mass2*mass2 + std::pow(p2Norm, 2)) ;
  momentum2(1) = p2[0] ;
  momentum2(2) = p2[1] ;
  momentum2(3) = p2[2] ;

  momentumg(0) = std::sqrt(massg*massg + std::pow(pgNorm, 2)) ;
  momentumg(1) = pg[0] ;
  momentumg(2) = pg[1] ;
  momentumg(3) = pg[2] ;

  momentum3(0) = std::sqrt(mass3*mass3 + std::pow(p3Norm, 2)) ;
  momentum3(1) = p3[0] ;
  momentum3(2) = p3[1] ;
  momentum3(3) = p3[2] ;

  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));
  finalState4->SetMomentum(utilities::LorentzBoost(velocity, momentumg));
}

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


void DecayMode::TwoBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, double Q, ublas::vector<double> dir) {
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

std::vector<Particle*> BetaRadiative::Decay(Particle* initState, double Q, double daughterExEn) {
  // Beta decay using 4-body decay radiative correction
  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();

  if (dm.configOptions.general.Verbosity >= 2)
    Info(Form("Beta Radiative Decay with Q = %.1f, parent excitation energy = %.1f and daughter excitation energy = %.1f", Q, initState->GetExcitationEnergy(), daughterExEn), 1);

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
  double RecoilMass = Recoil->GetMass();
  double RecoilRadius = utilities::ApproximateRadius(Recoil->GetCharge()+Recoil->GetNeutrons());
  Particle* ChargedLepton = dm.GetNewParticle(- BetaSign * 11);
  Particle* NeutralLepton = dm.GetNewParticle(BetaSign * 12);
  double InitialMass = initState->GetMass();
  
  //Name as key for channel properties
  std::ostringstream oss;
  oss << "BetaRadiative:" << "Sign" << BetaSign << "Z" << Recoil->GetCharge() << "A" << Recoil->GetCharge() + Recoil->GetNeutrons() << "Q" << Q;
  std::string ChannelName = oss.str();

  // Build or Get info of the channel using oss
  double mf;
  double mgt;
  double PH;
  double W_max_H;
  double W_max_VS;
  int Type;

  try
  {
    W_max_H = dm.GetChannelDistributionMaxs(ChannelName).first;
    W_max_VS = dm.GetChannelDistributionMaxs(ChannelName).second;
    PH = dm.GetChannelPH(ChannelName);
    Type = dm.GetChannelBetaType(ChannelName);
    mf = dm.GetChannelMf(ChannelName);
    mgt = dm.GetChannelMgt(ChannelName);
  }
  catch(const std::exception& e)
  {
    double mixing_ratio;
    Type = utilities::FindMatrixElement(initState, Recoil, mf, mgt, mixing_ratio);
    double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, E0/3.+utilities::EMASSC2, Recoil->GetCharge(), -BetaSign);
    PH = radiativecorrections::PH(dm.configOptions.betaDecay.Cs, mf, mgt, a, InitialMass, RecoilMass, Recoil->GetCharge(), RecoilRadius, BetaSign);
    W_max_H = radiativecorrections::WH_max(1e6, dm.configOptions.betaDecay.Cs, mf, mgt, a, InitialMass, RecoilMass, Recoil->GetCharge(), RecoilRadius, BetaSign);
    W_max_VS = radiativecorrections::W0VS_max(1e6, dm.configOptions.betaDecay.Cs, mf, mgt, a, InitialMass, RecoilMass, Recoil->GetCharge(), RecoilRadius, BetaSign);
    dm.RegisterChannelPropreties(ChannelName, nullptr, 0., Type, 0., 0., 0., W_max_H, W_max_VS, PH);
    dm.SetChannelMf(ChannelName, mf);
    dm.SetChannelMgt(ChannelName, mgt);
  }

  //
  ublas::vector<double> NeutralLepton_FourMomentum(4);
  ublas::vector<double> NeutralLepton_Dir;
  ublas::vector<double> ChargedLepton_FourMomentum(4);
  ublas::vector<double> ChargedLepton_Dir;
  //

  // Random for Hard vs Soft/Virtual Bremsstrahlung
  std::uniform_real_distribution<double> distribution(0., 1.);
  double ph = distribution(dm.generator);
  if (ph < PH)
  {
    if (dm.configOptions.general.Verbosity >= 2)
      Info(Form("Hard Bremsstrahlung"), 2);
    // Hard Bremsstrahlung
    Particle *Gamma = DecayManager::GetInstance().GetNewParticle(22);
    ublas::vector<double> Gamma_FourMomentum(4);
    std::uniform_real_distribution<double> distribution_wHmax(0.0, W_max_H);
    double W_point_H = 0;
    double W_H = W_max_H;

    double E2;
    double E1;
    double K;

    ublas::vector<double> n_ELECTRON(3);
    ublas::vector<double> n_GAMMA(3);
    ublas::vector<double> n_NEUTRINO(3);

    while (W_H > W_point_H)
    {
      W_H = distribution_wHmax(dm.generator);

      double U[8] = {distribution(dm.generator), distribution(dm.generator), distribution(dm.generator), distribution(dm.generator), distribution(dm.generator), distribution(dm.generator), distribution(dm.generator), distribution(dm.generator)};
      E2 = 1. + (radiativecorrections::delta(InitialMass, RecoilMass, BetaSign) - 1.) * U[0];
      double E10 = radiativecorrections::delta(InitialMass, RecoilMass, BetaSign) - E2;
      double omega = dm.configOptions.betaDecay.Cs * E10;
      K = omega * exp(-U[1] * log(dm.configOptions.betaDecay.Cs));
      E1 = E10 - K;

      double BETA = std::sqrt(1. - 1. / std::pow(E2, 2));
      double N = 0.5 * log((1. + BETA) / (1. - BETA));

      double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2])) / BETA;
      double COS_NEUTRINO = 2. * U[3] - 1.;
      double COS_ELECTRON = 2. * U[4] - 1.;

      double PHI_GAMMA = 2. * utilities::PI * U[5];
      double PHI_NEUTRINO = 2. * utilities::PI * U[6];
      double PHI_ELECTRON = 2. * utilities::PI * U[7];

      double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2)));
      double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
      double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2)));


      n_ELECTRON[0] = SIN_ELECTRON * cos(PHI_ELECTRON);
      n_ELECTRON[1] = SIN_ELECTRON * sin(PHI_ELECTRON);
      n_ELECTRON[2] = COS_ELECTRON;
      double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0};
      double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

      double n_PERPENDICULAIRE_GAMMA[3];
      n_NEUTRINO[0] = SIN_NEUTRINO * cos(PHI_NEUTRINO);
      n_NEUTRINO[1] = SIN_NEUTRINO * sin(PHI_NEUTRINO);
      n_NEUTRINO[2] = COS_NEUTRINO;
      for (int j = 0; j < 3; j++)
      {
        n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA);
        n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA;
      }

      double N1_N2 = n_NEUTRINO[0] * n_ELECTRON[0] + n_NEUTRINO[1] * n_ELECTRON[1] + n_NEUTRINO[2] * n_ELECTRON[2];
      double N1_K = n_NEUTRINO[0] * n_GAMMA[0] + n_NEUTRINO[1] * n_GAMMA[1] + n_NEUTRINO[2] * n_GAMMA[2];

      double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, E2, Recoil->GetCharge(), -BetaSign);
      W_point_H = radiativecorrections::WH(E2, K, COS_GAMMA, N1_K, N1_N2, mf, mgt, a, InitialMass, RecoilMass, Recoil->GetCharge(), RecoilRadius, BetaSign);
    }

    ublas::vector<double> velocity = -initState->GetVelocity();
    double eMomentum = std::sqrt(std::pow(E2 * utilities::EMASSC2, 2) - std::pow(utilities::EMASSC2, 2.));
    double enubarMomentum = E1 * utilities::EMASSC2;
    double gammaMomentum = K * utilities::EMASSC2;

    ChargedLepton_FourMomentum(0) = E2 * utilities::EMASSC2;
    ChargedLepton_FourMomentum(1) = eMomentum * n_ELECTRON[0];
    ChargedLepton_FourMomentum(2) = eMomentum * n_ELECTRON[1];
    ChargedLepton_FourMomentum(3) = eMomentum * n_ELECTRON[2];

    NeutralLepton_FourMomentum(0) = enubarMomentum;
    NeutralLepton_FourMomentum(1) = enubarMomentum * n_NEUTRINO[0];
    NeutralLepton_FourMomentum(2) = enubarMomentum * n_NEUTRINO[1];
    NeutralLepton_FourMomentum(3) = enubarMomentum * n_NEUTRINO[2];

    Gamma_FourMomentum(0) = gammaMomentum;
    Gamma_FourMomentum(1) = gammaMomentum * n_GAMMA[0];
    Gamma_FourMomentum(2) = gammaMomentum * n_GAMMA[1];
    Gamma_FourMomentum(3) = gammaMomentum * n_GAMMA[2];

    ChargedLepton->SetMomentum(ChargedLepton_FourMomentum);
    NeutralLepton->SetMomentum(NeutralLepton_FourMomentum);
    Gamma->SetMomentum(Gamma_FourMomentum);
    FourBodyDecay(velocity, ChargedLepton, NeutralLepton, Recoil, Gamma);

    finalStates.push_back(Recoil);
    finalStates.push_back(ChargedLepton);
    finalStates.push_back(NeutralLepton);
    finalStates.push_back(Gamma);

    return finalStates;
  }
  else
  {
    if (dm.configOptions.general.Verbosity >= 2)
      Info(Form("Soft/Virtual Bremsstrahlung"), 2);
    // Soft/Virtual Bremsstrahlung
    std::uniform_real_distribution<double> distribution_wVSmax(0.0, W_max_VS);
    double W_point_VS = 0;
    double W_VS = W_max_VS;
    double E2;
    double COS_NEUTRINO;
    double U[5];

    while (W_VS > W_point_VS)
    {
      W_VS = distribution_wVSmax(dm.generator);
      
      U[0] = distribution(dm.generator);
      U[1] = distribution(dm.generator);
      U[2] = distribution(dm.generator);
      U[3] = distribution(dm.generator);
      U[4] = distribution(dm.generator);

      E2 = 1. + (radiativecorrections::delta(InitialMass, RecoilMass, BetaSign) - 1.) * U[0];
      COS_NEUTRINO = 2. * U[1] - 1.;

      double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, E2, Recoil->GetCharge(), -BetaSign);
      W_point_VS = radiativecorrections::W0VS(E2, COS_NEUTRINO, dm.configOptions.betaDecay.Cs, mf, mgt, a, InitialMass, RecoilMass, Recoil->GetCharge(), RecoilRadius, BetaSign);
    }

    double E10 = radiativecorrections::delta(InitialMass, RecoilMass, BetaSign) - E2;
    double BETA = std::sqrt(1. - 1. / std::pow(E2, 2));

    double COS_ELECTRON = 2. * U[2] - 1.;
    double PHI_NEUTRINO = 2. * utilities::PI * U[3];
    double PHI_ELECTRON = 2. * utilities::PI * U[4];
    double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
    double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2)));

    ublas::vector<double> n_ELECTRON(3);
    n_ELECTRON[0] = SIN_ELECTRON * cos(PHI_ELECTRON);
    n_ELECTRON[1] = SIN_ELECTRON * sin(PHI_ELECTRON);
    n_ELECTRON[2] = COS_ELECTRON;
    double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0};
    double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

    double n_PERPENDICULAIRE_NEUTRINO[3];
    ublas::vector<double> n_NEUTRINO(3);
    for (int j = 0; j < 3; j++)
    {
      n_PERPENDICULAIRE_NEUTRINO[j] = n_ELECTRON_PRIME[j] * cos(PHI_NEUTRINO) + n_ELECTRON_SECOND[j] * sin(PHI_NEUTRINO);
      n_NEUTRINO[j] = n_ELECTRON[j] * COS_NEUTRINO + n_PERPENDICULAIRE_NEUTRINO[j] * SIN_NEUTRINO;
    }

    ublas::vector<double> velocity = -initState->GetVelocity();
    double ChargedLeptonMomentum = std::sqrt(std::pow(E2 * utilities::EMASSC2, 2) - std::pow(utilities::EMASSC2, 2.));
    double NeutralLeptonMomentum = E10 * utilities::EMASSC2;

    ChargedLepton_FourMomentum(0) = E2 * utilities::EMASSC2;
    ChargedLepton_FourMomentum(1) = ChargedLeptonMomentum * n_ELECTRON[0];
    ChargedLepton_FourMomentum(2) = ChargedLeptonMomentum * n_ELECTRON[1];
    ChargedLepton_FourMomentum(3) = ChargedLeptonMomentum * n_ELECTRON[2];
    ChargedLepton->SetMomentum(ChargedLepton_FourMomentum);

    ThreeBodyDecay(velocity, ChargedLepton, NeutralLepton, Recoil, n_NEUTRINO, Q);

    finalStates.push_back(Recoil);
    finalStates.push_back(ChargedLepton);
    finalStates.push_back(NeutralLepton);

    return finalStates;
  }
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

  // Name as key for channel properties
  std::ostringstream oss;
  oss << "Beta:" << "Sign" << BetaSign << "Z" << Recoil_Z << "A" << Recoil_Z + Recoil->GetNeutrons() << "Q" << Q;
  std::string ChannelName = oss.str();
  // Build or Get info of the channel using oss
  double mf;
  double mgt;

  std::vector<std::vector<double> >* dist;
  double dist_max = 0.;
  double j_i;
  double j_f;
  int Type;
  try
  {
    dist = dm.GetChannelDistribution(ChannelName);
    dist_max = dm.GetChannelDistributionMax(ChannelName);
    j_i = dm.GetChannelJi(ChannelName);
    j_f = dm.GetChannelJf(ChannelName);
    Type = dm.GetChannelBetaType(ChannelName);
    mf = dm.GetChannelMf(ChannelName);
    mgt = dm.GetChannelMgt(ChannelName);
  }
  catch (const std::invalid_argument &e)
  {
    double mixing_ratio;
    Type = utilities::FindMatrixElement(initState, Recoil, mf, mgt, mixing_ratio);
    dist = spectrumGen->GenerateSpectrum(initState, Recoil, E0, Type, mf, mgt, mixing_ratio);
    double b = correlation::CalculateFierz(mf, mgt, initState->GetCharge(), -BetaSign);
    for (int i = 0; i < dist->size(); i++)
    {
      double E = ((*dist)[i])[0] + utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH * (1 + b * utilities::EMASSC2 / E + (-BetaSign) * 4. / 3. * E / (initState->GetMass()) * dm.configOptions.nuclear.WeakMagnetism);
    }
    dist_max = utilities::CalculateMax(*dist);
    j_i = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy());
    j_f = utilities::GetJpi(Recoil_Z + Recoil->GetNeutrons(), Recoil_Z, Recoil->GetExcitationEnergy());    
    dm.RegisterChannelPropreties(ChannelName, dist, dist_max, Type, j_i, j_f);
    dm.SetChannelMf(ChannelName, mf);
    dm.SetChannelMgt(ChannelName, mgt);
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
    double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, ChargedLepton_Energy, Recoil_Z, -BetaSign);
    double b = correlation::CalculateFierz(mf, mgt, Recoil_Z, -BetaSign);
    ChargedLepton_Dir = utilities::GetParticleDirection(NeutralLepton_Dir, 1 + b * utilities::EMASSC2/ChargedLepton_Energy, a * ChargedLepton_Momentum / ChargedLepton_Energy);
  }
  else
  {
    ///// IF THE NUCLEUS IS POLARISED ////
    double b = correlation::CalculateFierz(mf, mgt, Recoil_Z, -BetaSign);
    double align = dm.configOptions.nuclear.Alignment;
    ublas::vector<double> polDir(3);
    polDir(0) = dm.configOptions.nuclear.PolarisationX;
    polDir(1) = dm.configOptions.nuclear.PolarisationY;
    polDir(2) = dm.configOptions.nuclear.PolarisationZ;
    polDir = utilities::NormaliseVector(polDir);
    double polMag = dm.configOptions.nuclear.PolarisationMag;

    double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, ChargedLepton_Energy, Recoil_Z, -BetaSign);
    double c = 0;
    double A = 0;
    double B = 0;
    double D = 0;

    if (j_i > 0)
    {
      A = correlation::CalculateBetaAssymetry(mf, mgt, j_i, j_f, -BetaSign, Recoil_Z, ChargedLepton_Energy);
      B = correlation::CalculateNeutrinoAssymetry(mf, mgt, j_i, j_f, -BetaSign, Recoil_Z, ChargedLepton_Energy);
      D = correlation::CalculateDTripleCorrelation(mf, mgt, j_i, j_f, -BetaSign, Recoil_Z, ChargedLepton_Energy);
      if (j_i > 0.5)
      {
        c = correlation::CalculateAlignmentCorrelation(mf, mgt, j_i, j_f, -BetaSign, Recoil_Z, ChargedLepton_Energy);
      }
    }
    c *= align * (-1);
    A *= polMag;
    B *= polMag;
    D *= polMag;

    if (dm.configOptions.general.Verbosity >= 2)
    {
      Info(Form("Correlation Coefficients: a = %.3f, b = %.3f, c = %.3f, A = %.3f, B = %.3f, D = %.3f", a, b, c, A, B, D), 3);
    }

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
  if (dm.configOptions.decay.NuclearLevelWidth) 
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
  if (dm.configOptions.decay.NuclearLevelWidth) 
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
  int Recoil_Z = Recoil->GetCharge();
  int Recoil_A = Recoil_Z + Recoil->GetNeutrons();

  if (initState->GetLastGamma() != nullptr && dm.configOptions.decay.GammaGammaCorrelation)
  {

    std::ostringstream oss;
    oss << "GammaGamma:" << "Z" << Recoil_Z << "A" << Recoil_A << "Ei" << initState->GetExcitationEnergy() + initState->GetLastGamma()->GetKinEnergy() << "Em" << initState->GetExcitationEnergy() << "Ef" << Recoil->GetExcitationEnergy();
    // std::cout << oss.str() << std::endl;
    // Gamma - Gamma correlation (E_i --> E --> E_f)
    Particle *gamma_1 = initState->GetLastGamma();
    ublas::vector<double> gamma_1_dir = gamma_1->Get3Momentum();

    //
    std::vector<double> ak;
    double W_max;
    try
    {
      ak = dm.GetChannelDistribution(oss.str())->at(0);
      W_max = dm.GetChannelDistributionMax(oss.str());  
    }
    catch (const std::invalid_argument &e)
    {
      // Calculation
      double j_i = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy() + gamma_1->GetKinEnergy());
      double j = utilities::GetJpi(initState->GetCharge() + initState->GetNeutrons(), initState->GetCharge(), initState->GetExcitationEnergy());
      double j_f = utilities::GetJpi(Recoil->GetCharge() + Recoil->GetNeutrons(), Recoil->GetCharge(), Recoil->GetExcitationEnergy());
    
      double delta_1 = initState->GetMixingRatio(initState->GetExcitationEnergy() + gamma_1->GetKinEnergy(), initState->GetExcitationEnergy());
      std::pair<int, int> l1 = initState->GetMultipolarities(initState->GetExcitationEnergy() + gamma_1->GetKinEnergy(), initState->GetExcitationEnergy());
      double delta_2 = Recoil->GetMixingRatio(initState->GetExcitationEnergy(), daughterExEn);
      std::pair<int, int> l2 = Recoil->GetMultipolarities(initState->GetExcitationEnergy(), daughterExEn);

      ak = correlation::CaluclateGammaCoefficient_a(std::abs(j_i), std::abs(j), std::abs(j_f), l1, delta_1, l2, delta_2);
      W_max = correlation::MaxAnalyticalGammaCorrelation(ak);

      std::vector<std::vector<double>>* dist = new std::vector<std::vector<double>>();
      dist->push_back(ak);
      dm.RegisterChannelPropreties(oss.str(), dist, W_max, j_i, j_f, j); 
    }

    if (dm.configOptions.general.Verbosity >= 2)
      Info(Form("Gamma-Gamma correlation coefficients: a0 = %.4f, a2 = %.4f, a4 = %.4f", ak[0], ak[1], ak[2]), 1);
    
    ublas::vector<double> gamma_2_dir(3);
    double W = W_max;
    double W_point = 0;
    std::uniform_real_distribution<double> distribution(0.0, W_max);
    std::uniform_real_distribution<double> theta_dist(0, M_PI);
    double theta = 0;
    while (W_point < W)
    {
      W = distribution(dm.generator);    
      theta = theta_dist(dm.generator);
      W_point = correlation::GammaGammaW(cos(theta), ak);
    }

    gamma_1_dir = utilities::NormaliseVector(gamma_1_dir);
    ublas::vector<double> perp = utilities::CrossProduct(gamma_1_dir, utilities::RandomDirection());
    perp = utilities::NormaliseVector(perp);
    gamma_2_dir = utilities::RotateAroundVector(gamma_1_dir, perp, theta);

    TwoBodyDecay(velocity, Recoil, gamma, Q, gamma_2_dir);
  }
  else
  {
    // Info(Form("Single gamma decay"), 1);
    TwoBodyDecay(velocity, Recoil, gamma, Q);
  }
 

  // Saving the last gamma with a temporary particle pointer to let the particlestack deletable
  Particle *temp_gamma = DecayManager::GetInstance().GetNewParticle(22, 0, 0, true);
  temp_gamma->SetMomentum(gamma->GetMomentum());
  Recoil->SetLastGamma(temp_gamma);
  //
  finalStates.push_back(Recoil);
  finalStates.push_back(gamma);

  // std::cout << gamma->Get3Momentum() << std::endl;

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

BetaRadiative::BetaRadiative() { }

ConversionElectron::ConversionElectron() { }

Proton::Proton () { }

Alpha::Alpha () { }

Gamma::Gamma () { }

ElectronCapture::ElectronCapture () { }

}//End of CRADLE namespace
