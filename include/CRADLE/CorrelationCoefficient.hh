#ifndef CRADLE_CORRELATIONCOEFFICIENT_H
#define CRADLE_CORRELATIONCOEFFICIENT_H

#include "CLI11.hpp"
#include "CRADLE/Messenger.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/ConfigParser.hh"
#include <complex>
#include <string>

namespace CRADLE
{
  namespace correlation
  {
    using namespace boost::numeric::ublas;

    const double PI = 3.14159265359;
    const double C = 299792458;                     // m/s
    const double EMASSC2 = 510.9989461;             // keV
    const double PMASSC2 = 938272.046;              // keV
    const double NMASSC2 = 939565.4133;             // keV
    const double UMASSC2 = 931494.10242;            // keV
    const double ALPHAMASSC2 = 3727379.4066;        // keV
    const double FINESTRUCTUREMASSC2 = 3727379.508; // keV
    const double FINESTRUCTURE = 0.0072973525664;
    const double E = 2.718281828459045;
    const double HBAR = 6.58211889e-16;                         // ev*s
    const double NATURALLENGTH = HBAR * C / EMASSC2 / 1000.;    // m
    const double EULER_MASCHERONI_CONSTANT = 0.577215664901532; /**< the Euler-Mascheroni constant */

    // doubly coded no choice ?
    inline vector<double> CrossProduct(vector<double> first, vector<double> second)
    {
      vector<double> v(3);

      v(0) = first(1) * second(2) - first(2) * second(1);
      v(1) = first(2) * second(0) - first(0) * second(2);
      v(2) = first(0) * second(1) - first(1) * second(0);

      return v;
    }
    //

    // Polarisation Variables (from Nuclear Physics 4 (1957) 206-212; by J. D. Jackson)
    inline double SmallLambdaJiJfFactor(double j_i, double j_f)
    {
      if (j_f == (j_i - 1))
      {
        return 1;
      }
      else if (j_f == j_i)
      {
        return 1 / (j_i + 1);
      }
      else if (j_f == (j_i + 1))
      {
        return -j_i / j_f;
      }
      else
        return 0;
    }
    inline double BigLambdaJiJfFactor(double j_i, double j_f)
    {
      if (j_f == (j_i - 1))
      {
        return 1;
      }
      else if (j_f == j_i)
      {
        return -(2 * j_i - 1) / (j_i + 1);
      }
      else if (j_f == (j_i + 1))
      {
        return j_i * (2 * j_i - 1) / j_f / (2 * j_f + 1);
      }
      else
        return 0;
    }
    //

    //////////// *** BETA DECAY *** ////////////
    // Xi
    inline double CalculateXiBetaDecay(double mf, double mgt)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      return mf * mf * (norm(c.CS) + norm(c.CV) + norm(c.CSP) + norm(c.CVP)) + mgt * mgt * (norm(c.CT) + norm(c.CTP) + norm(c.CA) + norm(c.CAP));
    }

    // b (Fierz Interference Term)
    inline double CalculateFierz(double mf, double mgt, int Z, int betaType)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (c.b != std::nan(""))
        return c.b; // set by user

      double gamma = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.));
      return 2. * gamma * betaType * (mf * mf * (c.CS * conj(c.CV) + c.CSP * conj(c.CVP)).real() + mgt * mgt * (c.CT * conj(c.CA) + c.CTP * conj(c.CAP)).real()) / CalculateXiBetaDecay(mf, mgt);
    }

    // a (Beta-Neutrino Angular Correlation)
    inline double CalculateBetaNeutrinoAsymmetry(double mf, double mgt, double energy, int Z, int betaType)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (c.a != std::nan(""))
        return c.a; // set by user

      double coulombCorr = FINESTRUCTURE * Z / std::sqrt(energy * energy / EMASSC2 / EMASSC2 - 1);
      double a = mf * mf * (-norm(c.CS) - norm(c.CSP) + norm(c.CV) + norm(c.CVP) - betaType * coulombCorr * 2. * (c.CS * conj(c.CV) + c.CSP * conj(c.CVP)).imag()) +
                 mgt * mgt / 3. * (-norm(c.CA) - norm(c.CAP) + norm(c.CT) + norm(c.CTP) + betaType * coulombCorr * 2. * (c.CT * conj(c.CA) + c.CTP * conj(c.CAP)).imag());
      return a / CalculateXiBetaDecay(mf, mgt);
    }

    // A (Beta asymmetry)
    inline double CalculateBetaAssymetry(double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (c.A != std::nan("")) // set by user
        return c.A;

      double coulombCorr = FINESTRUCTURE * Z / std::sqrt(energy * energy / EMASSC2 / EMASSC2 - 1); // alpha*Z*m_e/p_e
      double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.)) * EMASSC2 / energy;     // gamma*m_e/E
      double A = mgt * mgt * SmallLambdaJiJfFactor(j_i, j_f) * (2. * betaType * (c.CT * conj(c.CTP) - c.CA * conj(c.CAP)).real() + 2. * gamma_ratio * (c.CT * conj(c.CAP) + c.CTP * conj(c.CA)).imag());
      if (j_i == j_f)
      { // Mixed decay factors
        A += mf * mgt * std::sqrt(j_i / (j_i + 1)) * (2. * (c.CS * conj(c.CTP) + c.CSP * conj(c.CT) - c.CV * conj(c.CAP) - c.CVP * conj(c.CA)).real() + 2 * betaType * coulombCorr * (c.CS * conj(c.CAP) + c.CSP * conj(c.CA) - c.CV * conj(c.CTP) - c.CVP * conj(c.CT)).imag());
      }
      A /= CalculateXiBetaDecay(mf, mgt);
      return A;
    }

    // B (Neutrino asymmetry)
    inline double CalculateNeutrinoAssymetry(double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (c.B != std::nan(""))
        return c.B; // set by user

      double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.)) * EMASSC2 / energy; // gamma*m_e/E
      double B = 2. * mgt * mgt * SmallLambdaJiJfFactor(j_i, j_f) * (gamma_ratio * (c.CT * conj(c.CAP) + c.CTP * conj(c.CA)).real() + betaType * (c.CT * conj(c.CTP) + c.CA * conj(c.CAP)).real());
      if (j_i == j_f)
      { // Mixed decay factors
        B -= mf * mgt * std::sqrt(j_i / (j_i + 1)) * (2. * (c.CS * conj(c.CTP) + c.CSP * conj(c.CT) + c.CV * conj(c.CAP) + c.CVP * conj(c.CA)).real() + 2. * betaType * gamma_ratio * (c.CS * conj(c.CAP) + c.CSP * conj(c.CA) + c.CV * conj(c.CTP) + c.CVP * conj(c.CT)).real());
      }
      B /= CalculateXiBetaDecay(mf, mgt);
      return B;
    }

    // D (Triple correlation)
    inline double CalculateDTripleCorrelation(double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy)
    {
      CouplingConstants c = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (c.D != std::nan(""))
        return c.D; // set by user
      // Only mixed decay factors
      if (j_i == j_f)
      {
        double coulombCorr = FINESTRUCTURE * Z / std::sqrt(energy * energy / EMASSC2 / EMASSC2 - 1);
        double D = mf * mgt * std::sqrt(j_i / (j_i + 1)) * (-2. * coulombCorr * betaType * (c.CS * conj(c.CA) + c.CSP * conj(c.CAP) - c.CV * conj(c.CT) - c.CVP * conj(c.CTP)).real() + 2. * (c.CS * conj(c.CT) + c.CSP * conj(c.CTP) - c.CV * conj(c.CA) - c.CVP * conj(c.CAP)).imag());
        return D / CalculateXiBetaDecay(mf, mgt);
      }
      else
        return 0;
    }

    // c (Alignment correlation)
    inline double CalculateAlignmentCorrelation(double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy)
    {
      CouplingConstants coupling = (DecayManager::GetInstance()).configOptions.couplingConstants;
      if (coupling.c != std::nan(""))
        return coupling.c; // set by user

      double coulombCorr = FINESTRUCTURE * Z / std::sqrt(energy * energy / EMASSC2 / EMASSC2 - 1);
      double c = mgt * mgt * BigLambdaJiJfFactor(j_i, j_f) * (norm(coupling.CT) + norm(coupling.CTP) - norm(coupling.CA) - norm(coupling.CAP) + 2 * betaType * coulombCorr * (coupling.CT * conj(coupling.CA) + coupling.CTP * conj(coupling.CAP)).imag());
      return c / CalculateXiBetaDecay(mf, mgt);
    }

    inline double CalculateAngularCorrelationFactor(double a, double A, double B, double D, double E, double cosTheta_e, double cosTheta_enu, double phi)
    {
      /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. No c term, kept for compatibility purposes*/
      ublas::vector<double> elDir(3);
      ublas::vector<double> enuDir(3);

      // electron in XZ plane, using axial symmetry in Z direction (direction of J)
      double sinTheta_e = std::sqrt(1 - cosTheta_e * cosTheta_e);
      elDir(0) = sinTheta_e;
      elDir(1) = 0;
      elDir(2) = cosTheta_e;

      double sinTheta_enu = std::sqrt(1 - cosTheta_enu * cosTheta_enu);
      enuDir(0) = sinTheta_enu * std::cos(phi);
      enuDir(1) = sinTheta_enu * std::sin(phi);
      enuDir(2) = cosTheta_enu;

      double beta_e = std::sqrt(1 - EMASSC2 * EMASSC2 / E / E); // p_e/E_e; p_nu/E_nu = 1

      double angCorrFactor = 1.;
      angCorrFactor += a * beta_e * (elDir(0) * enuDir(0) + elDir(1) * enuDir(1) + elDir(2) * enuDir(2));
      angCorrFactor += A * beta_e * elDir(2);
      angCorrFactor += B * enuDir(2);
      angCorrFactor += D * beta_e * (elDir(0) * enuDir(1) - elDir(1) * enuDir(0));

      return angCorrFactor;
    }

    inline double CalculateAngularCorrelationFactor(double a, double b, double c, double A, double B, double D, double E, double cosTheta_e, double cosTheta_enu, double phi)
    {
      /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. Includes b and c term, and assumes perfect alignment. Assumes J/|J| is a unit vector in the z component. Used in tests*/
      vector<double> elDir(3);
      vector<double> enuDir(3);

      // electron in XZ plane, using axial symmetry in Z direction (direction of J)
      double sinTheta_e = std::sqrt(1 - cosTheta_e * cosTheta_e);
      elDir(0) = sinTheta_e;
      elDir(1) = 0;
      elDir(2) = cosTheta_e;

      double sinTheta_enu = std::sqrt(1 - cosTheta_enu * cosTheta_enu);
      enuDir(0) = sinTheta_enu * std::cos(phi);
      enuDir(1) = sinTheta_enu * std::sin(phi);
      enuDir(2) = cosTheta_enu;

      double beta_e = std::sqrt(1 - EMASSC2 * EMASSC2 / E / E); // p_e/E_e; p_nu/E_nu = 1

      double angCorrFactor = 1.;
      angCorrFactor += b * EMASSC2 / E;
      angCorrFactor += a * beta_e * (elDir(0) * enuDir(0) + elDir(1) * enuDir(1) + elDir(2) * enuDir(2));
      angCorrFactor += A * beta_e * elDir(2);
      angCorrFactor += B * enuDir(2);
      angCorrFactor += D * beta_e * (elDir(0) * enuDir(1) - elDir(1) * enuDir(0));
      angCorrFactor -= c * beta_e * (elDir(0) * enuDir(0) + elDir(1) * enuDir(1) - 2 * elDir(2) * enuDir(2)) / 3;

      return angCorrFactor;
    }

    inline double CalculateAngularCorrelationFactor(double a, double b, double c, double A, double B, double D, double E, vector<double> elDir, vector<double> enuDir, vector<double> polDir)
    {
      /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. Includes b and c term, noting that c, A, B and D already account for the polarisation and alignment dependent factors. Here starting from vectors themselves, so j isn't bound to z axis. polDir is unit vector*/

      double beta_e = std::sqrt(1 - EMASSC2 * EMASSC2 / E / E); // p_e/E_e; p_nu/E_nu = 1

      double angCorrFactor = 1.;
      angCorrFactor += b * EMASSC2 / E;
      angCorrFactor += a * beta_e * inner_prod(elDir, enuDir);
      angCorrFactor += A * beta_e * inner_prod(elDir, polDir);
      angCorrFactor += B * inner_prod(enuDir, polDir);
      angCorrFactor += D * beta_e * inner_prod(polDir, CrossProduct(elDir, enuDir));
      angCorrFactor += c * beta_e * (inner_prod(elDir, enuDir) / 3 - inner_prod(elDir, polDir) * inner_prod(enuDir, polDir));

      return angCorrFactor;
    }

    // ########## Analytical Maximum of all the Angular Correlation Factor ########## //

    inline double MaximumF(double a, double A, double B, double K, double znu)
    {
      return std::sqrt(std::pow(a * znu + A, 2) + K * K * (1 - znu * znu)) + B * znu;
    }

    inline double AnalyticalMaximumAngCorrFactor(double a, double b, double c, double A, double B, double D, double E)
    {
      /*Search of the maximum value analitically*/
      double beta = std::sqrt(1 - EMASSC2 * EMASSC2 / E / E);
      // scale a, c, A and D by beta. Note that c is multiplied already by the alignment, and A, B and D are multiplied by J as well
      a *= beta;
      c *= beta;
      A *= beta;
      D *= beta;

      double K = std::sqrt(D * D + (a + c / 3) * (a + c / 3));
      double a_st = a - 2. * c / 3;

      double A_m = (a_st * a_st - K * K) * (a_st * a_st - K * K - B * B);
      double B_m = 2 * a_st * A * (a_st * a_st - K * K - B * B);
      double C_m = a_st * a_st * A * A - A * A * B * B - B * B * K * K;
      double znu_m, znu_m2; // candidates for maximum
      vector<double> F_cand(4, 0);

      // computing the values at the extrema of the interval
      F_cand(0) = MaximumF(a_st, A, B, K, 1);
      F_cand(1) = MaximumF(a_st, A, B, K, -1);

      if (A_m == 0)
      {
        znu_m = -C_m / B_m;
        if ((znu_m > -1) && (znu_m < 1))
        {
          F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
        }
      }
      else if (B_m == 0 && C_m == 0)
      {
        F_cand(2) = MaximumF(a_st, A, B, K, 0); // cover some edge cases like only D non-zero that should never happen in reality
      }
      else
      {
        double det = B_m * B_m - 4 * A_m * C_m;
        if ((det < 0) && (det > -1e-16))
          det = 0; // avoid spurious cases
        if (det >= 0)
        {
          znu_m = (-B_m + std::sqrt(det)) / 2 / A_m;
          if ((znu_m > -1) && (znu_m < 1))
          {
            F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
          }
          znu_m2 = (-B_m - std::sqrt(det)) / 2 / A_m;
          if ((znu_m2 > -1) && (znu_m2 < 1))
          {
            F_cand(3) = MaximumF(a_st, A, B, K, znu_m2);
          }
        }
      }
      // std::cout << F_cand << std::endl;
      double F_max = *std::max_element(F_cand.begin(), F_cand.end());
      F_max += 1 + b * EMASSC2 / E; // adding the constant terms
      return F_max;
    }

    inline vector<double> MaximumAngCorrFactorPos(double a, double b, double c, double A, double B, double D, double E)
    {
      /*Search of the position of the maximum analitically*/
      double beta = std::sqrt(1 - EMASSC2 * EMASSC2 / E / E);
      // scale a, c, A and D by beta. Note that c is multiplied already by the alignment, and A, B and D are multiplied by J as well
      a *= beta;
      c *= beta;
      A *= beta;
      D *= beta;

      vector<double> max_pos(3);

      double K = std::sqrt(D * D + (a + c / 3) * (a + c / 3));
      double a_st = a - 2. * c / 3;

      max_pos(2) = std::atan2(D, a + c / 3); // phi

      double A_m = (a_st * a_st - K * K) * (a_st * a_st - K * K - B * B);
      double B_m = 2 * a_st * A * (a_st * a_st - K * K - B * B);
      double C_m = a_st * a_st * A * A - A * A * B * B - B * B * K * K;
      double znu_m, znu_m2; // candidates for maximum
      vector<double> F_cand(4, 0);
      double F_max;
      // computing the values at the extrema of the interval
      F_cand(0) = MaximumF(a_st, A, B, K, 1);
      F_cand(1) = MaximumF(a_st, A, B, K, -1);

      if (A_m == 0 && B_m != 0)
      {
        znu_m = -C_m / B_m;
        if ((znu_m > -1) && (znu_m < 1))
        {
          F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
        }
      }
      else if (B_m == 0 && C_m == 0)
      {
        F_cand(2) = MaximumF(a_st, A, B, K, 0); // cover some edge cases like only D non-zero that should never happen in reality
        znu_m = 0;
      }
      else
      {
        double det = B_m * B_m - 4 * A_m * C_m;
        if ((det < 0) && (det > -1e-16))
          det = 0; // avoid spurious cases
        if (det >= 0)
        {
          znu_m = (-B_m + std::sqrt(det)) / 2 / A_m;
          if ((znu_m > -1) && (znu_m < 1))
          {
            F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
          }
          znu_m2 = (-B_m - std::sqrt(det)) / 2 / A_m;
          if ((znu_m2 > -1) && (znu_m2 < 1))
          {
            F_cand(3) = MaximumF(a_st, A, B, K, znu_m2);
          }
        }
      }

      int indexMax = std::distance(F_cand.begin(), std::max_element(F_cand.begin(), F_cand.end()));
      switch (indexMax)
      {
      case 0:
        max_pos(1) = 1;                           // znu
        max_pos(0) = std::copysign(1., A + a_st); // ze
        break;
      case 1:
        max_pos(1) = -1;
        max_pos(0) = std::copysign(1., A - a_st);
        break;
      case 2:
        max_pos(1) = znu_m;
        max_pos(0) = (a_st * znu_m + A) / (MaximumF(a_st, A, 0, K, znu_m));
        break;
      case 3:
        max_pos(1) = znu_m2;
        max_pos(0) = (a_st * znu_m + A) / (MaximumF(a_st, A, 0, K, znu_m2));
        break;
      }
      return max_pos;
    }

//////////// *** GAMMA DECAY *** ////////////
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>

    using std::abs;
    using std::max;
    using std::min;
    using std::pow;
    using std::sqrt;
    using std::vector;

    //------------------------------------------------------------
    // Utility: factorial for integer or half-integer arguments
    // x must be integer or half-integer >= 0
    //------------------------------------------------------------
    inline double FactHalf(double x)
    {
      if (x < 0.0)
        return 0.0;

      // Gamma(x+1) = x!
      return std::tgamma(x + 1.0);
    }

    //------------------------------------------------------------
    // Triangle condition
    //------------------------------------------------------------
    inline bool Triangle(double a, double b, double c)
    {
      return (c <= a + b + 1e-12) && (c + 1e-12 >= abs(a - b));
    }

    //------------------------------------------------------------
    // Check if j,m values are integer/half-integer compatible
    //------------------------------------------------------------
    inline bool IsHalfInteger(double x)
    {
      double twox = 2.0 * x;
      return abs(twox - std::round(twox)) < 1e-10;
    }

    //------------------------------------------------------------
    // Delta coefficient for angular momentum algebra
    //------------------------------------------------------------
    inline double DeltaCoeff(double a, double b, double c)
    {
      if (!Triangle(a, b, c))
        return 0.0;

      double num = FactHalf(a + b - c) *
                   FactHalf(a - b + c) *
                   FactHalf(-a + b + c);
      double den = FactHalf(a + b + c + 1.0);

      if (den == 0.0)
        return 0.0;
      return sqrt(num / den);
    }

    //------------------------------------------------------------
    // Wigner 3j symbol
    // ( j1 j2 j3 )
    // ( m1 m2 m3 )
    // Racah formula
    //------------------------------------------------------------
    inline double Wigner3j(double j1, double j2, double j3,
                           double m1, double m2, double m3)
    {
      if (!IsHalfInteger(j1) || !IsHalfInteger(j2) || !IsHalfInteger(j3) ||
          !IsHalfInteger(m1) || !IsHalfInteger(m2) || !IsHalfInteger(m3))
        throw std::runtime_error("Wigner3j: arguments must be integer or half-integer.");

      if (abs(m1 + m2 + m3) > 1e-12)
        return 0.0;
      if (!Triangle(j1, j2, j3))
        return 0.0;
      if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
        return 0.0;

      double pref_phase = std::pow(-1.0, std::llround(j1 - j2 - m3));
      double pref = pref_phase *
                    DeltaCoeff(j1, j2, j3) *
                    sqrt(FactHalf(j1 + m1) * FactHalf(j1 - m1) *
                         FactHalf(j2 + m2) * FactHalf(j2 - m2) *
                         FactHalf(j3 + m3) * FactHalf(j3 - m3));

      double zmin = max({0.0, j2 - j3 - m1, j1 - j3 + m2});
      double zmax = min({j1 + j2 - j3, j1 - m1, j2 + m2});

      double sum = 0.0;
      for (int iz = (int)std::ceil(zmin - 1e-12); iz <= (int)std::floor(zmax + 1e-12); ++iz)
      {
        double z = (double)iz;

        double den =
            FactHalf(z) *
            FactHalf(j1 + j2 - j3 - z) *
            FactHalf(j1 - m1 - z) *
            FactHalf(j2 + m2 - z) *
            FactHalf(j3 - j2 + m1 + z) *
            FactHalf(j3 - j1 - m2 + z);

        if (den == 0.0)
          continue;

        sum += std::pow(-1.0, iz) / den;
      }

      return pref * sum;
    }

    //------------------------------------------------------------
    // Wigner 6j symbol using Racah formula
    //------------------------------------------------------------
    inline double Wigner6j(double j1, double j2, double j3,
                           double l1, double l2, double l3)
    {
      if (!IsHalfInteger(j1) || !IsHalfInteger(j2) || !IsHalfInteger(j3) ||
          !IsHalfInteger(l1) || !IsHalfInteger(l2) || !IsHalfInteger(l3))
        throw std::runtime_error("Wigner6j: arguments must be integer or half-integer.");

      if (!Triangle(j1, j2, j3) ||
          !Triangle(j1, l2, l3) ||
          !Triangle(l1, j2, l3) ||
          !Triangle(l1, l2, j3))
        return 0.0;

      double pref =
          DeltaCoeff(j1, j2, j3) *
          DeltaCoeff(j1, l2, l3) *
          DeltaCoeff(l1, j2, l3) *
          DeltaCoeff(l1, l2, j3);

      double zmin = max({j1 + j2 + j3,
                         j1 + l2 + l3,
                         l1 + j2 + l3,
                         l1 + l2 + j3});

      double zmax = min({j1 + j2 + l1 + l2,
                         j2 + j3 + l2 + l3,
                         j3 + j1 + l3 + l1});

      double sum = 0.0;
      for (int iz = (int)std::ceil(zmin - 1e-12); iz <= (int)std::floor(zmax + 1e-12); ++iz)
      {
        double z = (double)iz;

        double den =
            FactHalf(z - j1 - j2 - j3) *
            FactHalf(z - j1 - l2 - l3) *
            FactHalf(z - l1 - j2 - l3) *
            FactHalf(z - l1 - l2 - j3) *
            FactHalf(j1 + j2 + l1 + l2 - z) *
            FactHalf(j2 + j3 + l2 + l3 - z) *
            FactHalf(j3 + j1 + l3 + l1 - z);

        if (den == 0.0)
          continue;

        sum += std::pow(-1.0, iz) * FactHalf(z + 1.0) / den;
      }

      return pref * sum;
    }

    //------------------------------------------------------------
    // F_k coefficient for one gamma transition
    //
    // Transition: J_initial -> J_final
    // Multipole pair: L, Lp
    //------------------------------------------------------------
    inline double GammaFk(double L, double Lp, double I1, double I2, double k)
    {
      if (k < 0 || ((int)k % 2 != 0))
        return 0.0;

      double first = std::pow(-1.0, std::llround(I1 + I2 - 1));

      double second = sqrt((2.0 * L + 1.0) *
                           (2.0 * Lp + 1.0) *
                           (2.0 * I2 + 1.0) *
                           (2.0 * k + 1.0));

      double third = Wigner3j(L, Lp, k, 1., -1., 0.);

      double fourth = Wigner6j(L, Lp, k, I2, I2, I1);

      return first * second * third * fourth;
    }

    //------------------------------------------------------------
    // A_k for one mixed transition
    //
    // If delta = 0, pure L
    // If transition is mixed, use L and L+1
    //------------------------------------------------------------
    inline double GammaAkSingleTransition(double L, double Lp, double I1, double I2, double delta, double k)
    {
      // (eq.3)
      double F_LLI1I2 = GammaFk(L, L, I1, I2, k);
      // Info(Form("GammaAkSingleTransition: F(%d, %d, %f, %f) = %f", L, L, I1, I2, F_LLI1I2), 2);
      double F_LLpI1I2 = GammaFk(L, Lp, I1, I2, k);
      // Info(Form("GammaAkSingleTransition: F(%d, %d, %f, %f) = %f", L, Lp, I1, I2, F_LLpI1I2), 2);
      double F_LpLpI1I2 = GammaFk(Lp, Lp, I1, I2, k);
      // Info(Form("GammaAkSingleTransition: F(%d, %d, %f, %f) = %f", Lp, Lp, I1, I2, F_LpLpI1I2), 2);

      return (F_LLI1I2 + 2.0 * delta * F_LLpI1I2 + delta * delta * F_LpLpI1I2) /
             (1.0 + delta * delta);
    }

    //------------------------------------------------------------
    // Maximum allowed k for a pair of multipoles
    //------------------------------------------------------------
    inline int MaxEvenKForTransition(std::pair<int, int> L)
    {
      return 2 * std::max(L.first, L.second);
    }

    //------------------------------------------------------------
    // Requested final function
    //
    // Ji  -> J  by gamma1 with lowest multipole L1 and mixing delta1
    // J   -> Jf by gamma2 with lowest multipole L2 and mixing delta2
    //
    // Returns: [a0, a2, a4, ...]
    // with a0 = 1
    //------------------------------------------------------------
    inline vector<double> CaluclateGammaCoefficient_a(double Ji, double J, double Jf,
                                                      std::pair<int, int> L1, double delta1,
                                                      std::pair<int, int> L2, double delta2)
    {

      vector<double> a;
      a.push_back(1.0); // a0

      // HAL Id: hal-04964841
      if (!IsHalfInteger(Ji) || !IsHalfInteger(J) || !IsHalfInteger(Jf))
        throw std::runtime_error("CaluclateGammaCoefficient_a: spins must be integer or half-integer.");

      if (L1.first < 1 || L2.first < 1)
        return a;

      double maxK = std::min(MaxEvenKForTransition(L1), MaxEvenKForTransition(L2));

      // Info(Form("CaluclateGammaCoefficient_a: maxK1 = %d", maxK));

      for (double k = 2.; k <= maxK; k += 2.)
      {
        // Info("k = " + std::to_string(k), 1);
        double A1k = GammaAkSingleTransition(L1.first, L1.second, Ji, J, delta1, k);
        double A2k = GammaAkSingleTransition(L2.first, L2.second, Jf, J, delta2, k);

        double ak = A1k * A2k;
        a.push_back(ak);
      }

      return a;
    }

    //------------------------------------------------------------
    // Legendre polynomials needed for W(theta)
    //------------------------------------------------------------
    inline double P2(double x)
    {
      return 0.5 * (3.0 * x * x - 1.0);
    }

    inline double P4(double x)
    {
      double x2 = x * x;
      return (35.0 * x2 * x2 - 30.0 * x2 + 3.0) / 8.0;
    }

    //------------------------------------------------------------
    // Evaluate W(theta) from the a coefficients
    // a = [a0, a2, a4, ...]
    //------------------------------------------------------------
    inline double GammaGammaW(double x, const vector<double> &a)
    {
      double w = 0.0;

      if (a.size() >= 1)
        w += a[0];
      if (a.size() >= 2)
        w += a[1] * P2(x);
      if (a.size() >= 3)
        w += a[2] * P4(x);

      // Extend here for P6, P8 if needed
      return w;
    }

    //------------------------------------------------------------
    // Sampling theta analytically from the a coefficients
    // a = [a0, a2, a4, ...]
    //------------------------------------------------------------
    inline double MaxAnalyticalGammaCorrelation(const vector<double> &a)
    {
      // For now only implemented for a0 and a2, extend here for a4 if needed
      if (a.size() < 2)
        return 1.0; // isotropic

      else if (a.size() < 3)
      {
        double a0 = a[0];
        double a2 = a[1];

        // W(theta) = a0 + a2*P2(cos(theta))
        // Max at cos(theta) = 1 or -1 or -a0/a2*3/5 if in [-1,1]
        double w_max = a0 + abs(a2); // max at cos(theta) = 1 or -1
        double cos_theta_extremum = -a0 / a2 * 3.0 / 5.0;
        if (cos_theta_extremum >= -1.0 && cos_theta_extremum <= 1.0)
        {
          double w_extremum = a0 + a2 * P2(cos_theta_extremum);
          w_max = std::max(w_max, w_extremum);
        }
        return w_max;
      }

      else if (a.size() >= 3)
      {
        // For a0, a2 and a4, the maximum can be found by solving the derivative of W(theta) = a0 + a2*P2(cos(theta)) + a4*P4(cos(theta)) equal to zero. This leads to a cubic equation in cos(theta) that can be solved analytically.
        double a0 = a[0];
        double a2 = a[1];
        double a4 = a[2];
        // Coefficients of the cubic equation a4*P4'(x) + a2*P2'(x) = 0
        double c3 = 35.0 * a4 / 8.0;
        double c2 = -30.0 * a4 / 8.0;
        double c1 = 3.0 * a4 / 8.0 + 3.0 * a2 / 2.0;
        double c0 = -3.0 * a2 / 2.0;

        // Solve the cubic equation c3*x^3 + c2*x^2 + c1*x + c0 = 0 for x in [-1,1]
        // Using Cardano's formula for cubic equations
        double p = (3.0 * c3 * c1 - c2 * c2) / (3.0 * c3 * c3);
        double q = (2.0 * c2 * c2 * c2 - 9.0 * c3 * c2 * c1 + 27.0 * c3 * c3 * c0) / (27.0 * c3 * c3 * c3);
        double discriminant = (q * q / 4.0) + (p * p * p / 27.0);
        vector<double> candidates;

        if (discriminant > 0)
        {
          // One real root
          double sqrt_disc = std::sqrt(discriminant);
          double u = std::cbrt(-q / 2.0 + sqrt_disc);
          double v = std::cbrt(-q / 2.0 - sqrt_disc);
          double x = u + v - c2 / (3.0 * c3);
          if (x >= -1.0 && x <= 1.0)
            candidates.push_back(x);
        }
        else if (abs(discriminant) < 1e-12)
        {
          // Multiple real roots
          double u = std::cbrt(-q / 2.0);
          double x1 = 2.0 * u - c2 / (3.0 * c3);
          double x2 = -u - c2 / (3.0 * c3);
          if (x1 >= -1.0 && x1 <= 1.0)
            candidates.push_back(x1);
          if (x2 >= -1.0 && x2 <= 1.0)
            candidates.push_back(x2);
        }
        else
        {
          // Three real roots
          double r = std::sqrt(-p * p * p / 27.0);
          double phi = std::acos(-q / (2.0 * r));
          double m = 2.0 * std::cbrt(r);
          for (int k = 0; k < 3; ++k)
          {
            double x = m * std::cos((phi + 2.0 * M_PI * k) / 3.0) - c2 / (3.0 * c3);
            if (x >= -1.0 && x <= 1.0)
              candidates.push_back(x);
          }
        }

        // Evaluate W(theta) at the candidates and at the endpoints
        double w_max = a0 + abs(a2) + abs(a4); // max at endpoints
        for (double x : candidates)
        {
          double w = a0 + a2 * P2(x) + a4 * P4(x);
          w_max = std::max(w_max, w);
        }
        return w_max;
      }
      else
      {
        // Not implemented
        Warning("MaxAnalyticalGammaCorrelation: not implemented for a size >= 4, returning 1.0");
        return 1.0;
      }
    }
  }
}

#endif