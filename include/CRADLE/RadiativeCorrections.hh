#ifndef RADIATIVECORRECTIONS
#define RADIATIVECORRECTIONS

#include "CRADLE/Utilities.hh"
#include "CRADLE/DecayManager.hh"

namespace CRADLE
{

    namespace radiativecorrections
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
        const double e = std::sqrt((4 * PI * FINESTRUCTURE));
        const double HBAR = 6.58211889e-16;                         // ev*s
        const double NATURALLENGTH = HBAR * C / EMASSC2 / 1000.;    // m
        const double EULER_MASCHERONI_CONSTANT = 0.577215664901532; /**< the Euler-Mascheroni constant */
        const double LAMBDA = -1.2754;
        const double Vud = 0.97435;
        const double GF = 1.166378 * std::pow(10, -5);
        const double FERMICONSTANT = GF; // GF * Vud ;

        const std::string atoms[] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"};

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        inline double Spence_function(double x)
        {
            return -gsl_sf_dilog(x);
        }

        inline double delta(double MIMASSC2, double MFMASSC2, int mode)
        {
            double E_transition = (MIMASSC2 - MFMASSC2);
            double Er_max = (std::pow((E_transition - EMASSC2) / EMASSC2, 2) + 2 * (E_transition - EMASSC2) / EMASSC2) / (2 * MFMASSC2 / EMASSC2);
            /*std::cout << "masse i " << MIMASSC2 << "\n";
            std::cout << "masse f " << MFMASSC2 << "\n";
            std::cout << "e transi : " << E_transition << "\n" ;
            std::cout << "delta : " << (E_transition - Er_max*EMASSC2)/EMASSC2 << "\n";
            std::cout << "er : " << Er_max*EMASSC2 << "\n";*/
            double d;
            // Er_max = 0 ;
            if (mode == -1)
            {
                d = (E_transition - Er_max * EMASSC2) / EMASSC2 + 1;
            }
            else if (mode == 1)
            {
                d = (E_transition - Er_max * EMASSC2) / EMASSC2 - 1;
            }
            // std::cout << "delta : " << d << "\n";
            return d;
        }

        // inline double xhi(double MF, double MGT) {
        //     double x = std::pow(MF, 2) + std::pow(LAMBDA, 2) * std::pow(MGT, 2) ;
        //     return x ;
        // }

        // inline double xhi_a(double MF, double MGT) {
        //     double x_a = (std::pow(MF, 2) - std::pow(LAMBDA, 2) * std::pow(MGT, 2)/3.)/xhi(MF, MGT) ;
        //     return x_a ;
        // }

        inline double H0(double E2, double K, double COS_GAMMA, double MIMASSC2, double MFMASSC2, int mode)
        {
            double E1 = delta(MIMASSC2, MFMASSC2, mode) - K - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double p2_k = E2 * K - BETA * E2 * K * COS_GAMMA;
            double P_2 = 1 / std::pow(K, 2) + 1 / std::pow(p2_k, 2) - (2 * E2) / (K * p2_k);
            return E1 * (-(E2 + K) * P_2 + K / (p2_k));
        }

        inline double H1(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double MIMASSC2, double MFMASSC2, int mode)
        {
            double E1 = delta(MIMASSC2, MFMASSC2, mode) - K - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double p2_k = E2 * K - BETA * E2 * K * COS_GAMMA;
            double p1_p2 = BETA * E1 * E2 * N1_N2;
            double p1_k = E1 * K * N1_K;
            double P_2 = 1 / std::pow(K, 2) + 1 / std::pow(p2_k, 2) - (2 * E2) / (K * p2_k);
            return p1_p2 * (-P_2 + 1 / p2_k) + p1_k * ((E2 + K) / (K * p2_k) - 1 / std::pow(p2_k, 2));
        }

        inline double MBR(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            /*double EPS = (MF_2) + (std::pow(LAMBDA, 2)) * (MGT_2) ;
            double EPSa = (MF_2) - (std::pow(LAMBDA, 2)) * (MGT_2)/3. ;
            double a = EPSa/EPS ;*/
            // double a = xhi_a(MF, MGT) ;
            return 16 * std::pow(FERMICONSTANT, 2) * std::pow((MIMASSC2 / EMASSC2), 2) * std::pow(e, 2) * (1 * H0(E2, K, COS_GAMMA, MIMASSC2, MFMASSC2, mode) + a * H1(E2, K, COS_GAMMA, N1_K, N1_N2, MIMASSC2, MFMASSC2, mode)) * utilities::FermiFunction(Z, E2, R, -mode);
        }

        inline double M0(double E2, double COS, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            /*double EPS = (MF_2) + (std::pow(LAMBDA, 2)) * (MGT_2) ;
            double EPSa = (MF_2) - (std::pow(LAMBDA, 2)) * (MGT_2)/3. ;
            double a = EPSa/EPS ;*/
            // double xi = xhi(MF, MGT) ;
            double xi = correlation::CalculateXiBetaDecay(MF, MGT);
            // double a = xhi_a(MF, MGT) ;
            //  double xi = 1;
            return 16. * std::pow(FERMICONSTANT, 2) * xi * std::pow(MIMASSC2 / EMASSC2, 2) * E10 * E2 * (1 + a * BETA * COS) * utilities::FermiFunction(Z, E2, R, -mode);
        }

        inline double Mtilde(double E2, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double N = 0.5 * log((1. + BETA) / (1. - BETA));
            // double EPS = (MF_2) + (std::pow(LAMBDA, 2)) * (MGT_2) ;
            //  double xi = xhi(MF, MGT) ; /// MODIF, BEFORE ==1
            double xi = 1;
            return -(FINESTRUCTURE / PI) * (16. * std::pow(FERMICONSTANT, 2) * ((1 - std::pow(BETA, 2)) / BETA) * (N * std::pow((MIMASSC2 / EMASSC2), 2) * E10 * E2 * xi)) * utilities::FermiFunction(Z, E2, R, -mode);
        }

        inline double MVS(double E2, double COS, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double N = 0.5 * log((1. + BETA) / (1. - BETA));
            double omega = E10 * Cs;
            double zVS = (FINESTRUCTURE / PI) * (1.5 * log(PMASSC2 / EMASSC2) + 2. * ((N / BETA) - 1.) * log((2. * omega)) + 2 * (N / BETA) * (1. - N) + (2. / BETA) * (Spence_function(2. * BETA / (1. + BETA))) - (3. / 8.));
            return zVS * M0(E2, COS, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) + Mtilde(E2, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
        }

        inline double g_weight(double E2, double K, double COS_GAMMA)
        {
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double N = 0.5 * log((1. + BETA) / (1. - BETA));
            double p2_k = E2 * K - BETA * E2 * K * COS_GAMMA;
            return (BETA * E2) / (2 * N * p2_k);
        }

        inline double WH(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            if (E2 == 1.0) Error("WH : E2 is 1, returning 0 to avoid division by zero in beta calculation (Kinetic energy of the charged lepton is 0)");
            double E1 = delta(MIMASSC2, MFMASSC2, mode) - E2 - K;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            return (K * BETA * E1 * E2 * MBR(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode));
        }

        inline double rho_H(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            DecayManager &dm = DecayManager::GetInstance();
            if (dm.configOptions.general.Verbosity >= 2)
                Info("Caluculation of rho_H", 2);
            
            double somme_rhoH = 0.;
            double Vg = -32. * std::pow(PI, 3) * (delta(MIMASSC2, MFMASSC2, mode) - 1) * log(Cs);
            std::random_device rd;
            std::mt19937 generator(rd());
            int nout = 0;

            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            for (int i = 0; i < n; i++)
            {
                double u1 = distribution(generator);
                double u2 = distribution(generator);
                double u3 = distribution(generator);
                double u4 = distribution(generator);
                double u5 = distribution(generator);
                double u6 = distribution(generator);
                double u7 = distribution(generator);
                double u8 = distribution(generator);

                double E2 = 1. + (delta(MIMASSC2, MFMASSC2, mode) - 1.) * u1;
                double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
                double omega = Cs * E10;
                double K = omega * exp(-u2 * log(Cs));

                // if (E10 < 0.) {std::cout << "E10 : " << E10 << "\n";}
                // std::cout << "omega : " << omega << "\n";
                // std::cout << "K : " << K << "\n";

                double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
                double N = 0.5 * log((1. + BETA) / (1. - BETA));

                double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * u3)) / BETA;
                double COS_NEUTRINO = 2. * u4 - 1.;
                double COS_ELECTRON = 2. * u5 - 1.;

                double PHI_GAMMA = 2. * PI * u6;
                double PHI_NEUTRINO = 2. * PI * u7;
                double PHI_ELECTRON = 2. * PI * u8;

                double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2)));
                double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
                double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2)));

                double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON};
                double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.};
                double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

                double n_PERPENDICULAIRE_GAMMA[3];
                double n_PERPENDICULAIRE_NEUTRINO[3];
                double n_GAMMA[3];
                double n_NEUTRINO[3] = {SIN_NEUTRINO * cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};
                for (int j = 0; j < 3; j++)
                {
                    n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA);
                    n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA;
                }

                double N1_N2 = n_NEUTRINO[0] * n_ELECTRON[0] + n_NEUTRINO[1] * n_ELECTRON[1] + n_NEUTRINO[2] * n_ELECTRON[2];
                double N1_K = n_NEUTRINO[0] * n_GAMMA[0] + n_NEUTRINO[1] * n_GAMMA[1] + n_NEUTRINO[2] * n_GAMMA[2];
                double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / (g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2 / EMASSC2), 2));

                if (std::isnan(wh) || std::isinf(wh))
                {
                    /*double p2_k = E2 * K - BETA * E2 * K * COS_GAMMA ;
                    double E1 = delta(MIMASSC2, MFMASSC2, mode) - K - E2 ;
                    double P_2 = 1/std::pow(K, 2) + 1/std::pow(p2_k, 2) - (2 * E2)/(K * p2_k) ;
                    double p1_p2 = BETA * E1 * E2 * N1_N2 ;
                    double p1_k = E1 * K *  N1_K ;

                    double h0 = E1 * ( -(E2 + K)*P_2 + K/(p2_k) ) ;
                    double h1 = p1_p2 * ( -P_2 + 1/p2_k ) + p1_k * ( (E2 + K)/(K * p2_k) - 1/std::pow(p2_k, 2) ) ;

                    double EPSa = MF_2 - std::pow(LAMBDA, 2)*MGT_2/3 ;
                    double EPS = MF_2 + std::pow(LAMBDA, 2)*MGT_2 ;
                    double mbr = 16 * std::pow(FERMICONSTANT, 2) * std::pow( (MIMASSC2/EMASSC2) , 2 ) * std::pow(e, 2) * (EPS * h0 + EPSa * h1) ;

                    double g = (BETA * E2)/(2 * N * p2_k) ;

                    double whhh = ( K * BETA * E1 * E2 * mbr )/( std::pow(2, 13) * std::pow(utilities::PI, 8) * g * std::pow( (MIMASSC2/EMASSC2) , 2 )) ;


                    std::cout << "WH = nan" << "\n";
                    std::cout << "E1 : " << E1 << "\n";
                    std::cout << "E2 : " << E2 << "\n";
                    std::cout << "K : " << K << "\n";
                    std::cout << "delta : " << delta(MIMASSC2, MFMASSC2, mode) << "\n" ;

                    std::cout << "p2_k : " << p2_k << "\n";
                    std::cout << "P_2 : " << P_2 << "\n";
                    std::cout << "p1_p2 : " << p1_p2 << "\n";
                    std::cout << "p1_k : " << p1_k << "\n";

                    std::cout << "h0 : " << h0 << "\n";
                    std::cout << "h1 : " << h1 << "\n";
                    std::cout << "mbr : " << mbr << "\n";
                    std::cout << "wh : " << wh << "\n";
                    std::cout << "whhh : " << whhh << "\n";*/
                    nout += 1;
                }
                else
                {

                    somme_rhoH += wh;
                }

                // std::cout << "Wh one value : " << WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF_2, MGT_2, MIMASSC2, MFMASSC2, Z, R, mode) << "\n";
            }
            // std::cout << "nout : " << nout << "\n";
            double RHOH = (Vg * somme_rhoH) / (n - nout);
            return RHOH;
        }

        inline double delta_rho_H(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double somme_rhoH2 = 0.;
            double somme_rhoH = 0.;
            double Vg = -32. * std::pow(PI, 3) * (delta(MIMASSC2, MFMASSC2, mode) - 1) * log(Cs);
            std::random_device rd;
            std::mt19937 generator(rd());
            double nout = 0;

            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            for (int i = 0; i < n; i++)
            {
                double u1 = distribution(generator);
                double u2 = distribution(generator);
                double u3 = distribution(generator);
                double u4 = distribution(generator);
                double u5 = distribution(generator);
                double u6 = distribution(generator);
                double u7 = distribution(generator);
                double u8 = distribution(generator);

                double E2 = 1. + (delta(MIMASSC2, MFMASSC2, mode) - 1.) * u1;
                double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
                double omega = Cs * E10;
                double K = omega * exp(-u2 * log(Cs));

                // if (E10 < 0.) {std::cout << "E10 : " << E10 << "\n";}
                // std::cout << "omega : " << omega << "\n";
                // std::cout << "K : " << K << "\n";

                double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
                double N = 0.5 * log((1. + BETA) / (1. - BETA));

                double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * u3)) / BETA;
                double COS_NEUTRINO = 2. * u4 - 1.;
                double COS_ELECTRON = 2. * u5 - 1.;

                double PHI_GAMMA = 2. * PI * u6;
                double PHI_NEUTRINO = 2. * PI * u7;
                double PHI_ELECTRON = 2. * PI * u8;

                double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2)));
                double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
                double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2)));

                double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON};
                double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.};
                double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

                double n_PERPENDICULAIRE_GAMMA[3];
                double n_PERPENDICULAIRE_NEUTRINO[3];
                double n_GAMMA[3];
                double n_NEUTRINO[3] = {SIN_NEUTRINO * cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};
                for (int j = 0; j < 3; j++)
                {
                    n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA);
                    n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA;
                }

                double N1_N2 = n_NEUTRINO[0] * n_ELECTRON[0] + n_NEUTRINO[1] * n_ELECTRON[1] + n_NEUTRINO[2] * n_ELECTRON[2];
                double N1_K = n_NEUTRINO[0] * n_GAMMA[0] + n_NEUTRINO[1] * n_GAMMA[1] + n_NEUTRINO[2] * n_GAMMA[2];
                double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / (g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2 / EMASSC2), 2));

                if (!std::isnan(wh) && !std::isinf(wh))
                {
                    nout += 1;
                    somme_rhoH2 += pow(wh, 2);
                    somme_rhoH += wh;
                }

                // std::cout << "Wh one value : " << WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF_2, MGT_2, MIMASSC2, MFMASSC2, Z, R, mode) << "\n";
            }
            // std::cout << "nout : " << nout << "\n";
            double DELTA_RHO_H = Vg / (nout)*sqrt(somme_rhoH2 - pow(somme_rhoH, 2) / (nout));
            return DELTA_RHO_H;
        }

        inline double WH_max(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {

            DecayManager &dm = DecayManager::GetInstance();
            if (dm.configOptions.general.Verbosity >= 2)
                Info("Calculating maximum of WH for radiative corrections...", 2);
            double max = 0.;

            std::random_device rd;
            std::mt19937 generator(rd());

            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            for (int i = 0; i < n; i++)
            {
                double u1 = distribution(generator);
                double u2 = distribution(generator);
                double u3 = distribution(generator);
                double u4 = distribution(generator);
                double u5 = distribution(generator);
                double u6 = distribution(generator);
                double u7 = distribution(generator);
                double u8 = distribution(generator);

                double E2 = 1. + (delta(MIMASSC2, MFMASSC2, mode) - 1.) * u1;
                double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
                double omega = Cs * E10;
                double K = omega * exp(-u2 * log(Cs));

                double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
                double N = 0.5 * log((1. + BETA) / (1. - BETA));

                double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * u3)) / BETA;
                double COS_NEUTRINO = 2. * u4 - 1.;
                double COS_ELECTRON = 2. * u5 - 1.;

                double PHI_GAMMA = 2. * PI * u6;
                double PHI_NEUTRINO = 2. * PI * u7;
                double PHI_ELECTRON = 2. * PI * u8;

                double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2)));
                double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
                double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2)));

                double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON};
                double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.};
                double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

                double n_PERPENDICULAIRE_GAMMA[3];
                double n_PERPENDICULAIRE_NEUTRINO[3];
                double n_GAMMA[3];
                double n_NEUTRINO[3] = {SIN_NEUTRINO * cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};

                for (int j = 0; j < 3; j++)
                {
                    n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA);
                    n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA;
                }

                double N1_N2 = n_NEUTRINO[0] * n_ELECTRON[0] + n_NEUTRINO[1] * n_ELECTRON[1] + n_NEUTRINO[2] * n_ELECTRON[2];
                double N1_K = n_NEUTRINO[0] * n_GAMMA[0] + n_NEUTRINO[1] * n_GAMMA[1] + n_NEUTRINO[2] * n_GAMMA[2];

                // double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) ;
                double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / (g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2 / EMASSC2), 2));

                if (!std::isnan(wh) && !std::isinf(wh))
                {
                    if (wh > max)
                    {
                        max = wh;
                    }
                }
            }

            if (std::isnan(max) || std::isinf(max))
            {
                Error("WH_max : max is nan or inf, returning 0");
                return 0.;
            }

            return max;
        }

        inline double WH_mean(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double somme_wh = 0.;
            double nout = 0.;
            std::random_device rd;
            std::mt19937 generator(rd());

            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            for (int i = 0; i < n; i++)
            {
                double u1 = distribution(generator);
                double u2 = distribution(generator);
                double u3 = distribution(generator);
                double u4 = distribution(generator);
                double u5 = distribution(generator);
                double u6 = distribution(generator);
                double u7 = distribution(generator);
                double u8 = distribution(generator);

                double E2 = 1. + (delta(MIMASSC2, MFMASSC2, mode) - 1.) * u1;
                double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
                double omega = Cs * E10;
                double K = omega * exp(-u2 * log(Cs));

                double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
                double N = 0.5 * log((1. + BETA) / (1. - BETA));

                double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * u3)) / BETA;
                double COS_NEUTRINO = 2. * u4 - 1.;
                double COS_ELECTRON = 2. * u5 - 1.;

                double PHI_GAMMA = 2. * PI * u6;
                double PHI_NEUTRINO = 2. * PI * u7;
                double PHI_ELECTRON = 2. * PI * u8;

                double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2)));
                double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2)));
                double SIN_ELECTRON = std::sqrt((1 - std::pow(COS_ELECTRON, 2)));

                double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON};
                double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.};
                double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

                double n_PERPENDICULAIRE_GAMMA[3];
                double n_PERPENDICULAIRE_NEUTRINO[3];
                double n_GAMMA[3];
                double n_NEUTRINO[3] = {SIN_NEUTRINO * cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};

                for (int j = 0; j < 3; j++)
                {
                    n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA);
                    n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA;
                }

                double N1_N2 = n_NEUTRINO[0] * n_ELECTRON[0] + n_NEUTRINO[1] * n_ELECTRON[1] + n_NEUTRINO[2] * n_ELECTRON[2];
                double N1_K = n_NEUTRINO[0] * n_GAMMA[0] + n_NEUTRINO[1] * n_GAMMA[1] + n_NEUTRINO[2] * n_GAMMA[2];

                // double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) ;
                double wh = WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / (g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2 / EMASSC2), 2));

                if (!std::isnan(wh) && !std::isinf(wh))
                {
                    somme_wh += wh;
                    nout += 1;
                }
            }
            return somme_wh / nout;
        }

        // Equation 5.18
        inline double w0(double E2, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            // double xi = xhi(MF, MGT) ;
            double xi = 1;
            // double EPS = (MF_2) + (std::pow(LAMBDA, 2)) * (MGT_2) ;
            double weight0 = (std::pow(FERMICONSTANT, 2) * xi * BETA * std::pow(E10, 2) * std::pow(E2, 2)) / (2 * std::pow(PI, 3));

            // Info("Print all parameters for w0 calculation : ", 2);
            // Info("E2 : " + std::to_string(E2), 2);
            // Info("MF : " + std::to_string(MF), 2);
            // Info("MGT : " + std::to_string(MGT), 2);
            // Info("MIMASSC2 : " + std::to_string(MIMASSC2), 2);
            // Info("MFMASSC2 : " + std::to_string(MFMASSC2), 2);
            // Info("Z : " + std::to_string(Z), 2);
            // Info(Form("R : %e", R*1e15), 2);
            // Info(Form("E10 : %e", E10), 2);
            // Info(Form("BETA : %e", BETA), 2);
            // Info(Form("weight0 : %e", weight0), 2);
            // Info(Form("Fermi function : %e", utilities::FermiFunction(Z, E2, R, -mode)), 2);

            double res = weight0 * utilities::FermiFunction(Z, E2, R, -mode);

            if (std::isnan(res) || std::isinf(res))
                return 0;
        
            return res;
        }

        // Calcul du max de w0 pour Neumann Rejection
        inline double w0_max(double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double w0max = 0;
            double n = 5000;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode) - 1) / (n);
            std::vector<double> tableau_E2(n);
            for (int i = 0; i < n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                double w = w0(tableau_E2[i], MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);

                if (!std::isnan(w) && !std::isinf(w))
                {
                    if (w0max < w)
                    {
                        w0max = w;
                    }
                }
            }

            if (std::isnan(w0max) || std::isinf(w0max))
            {
                Error("w0_max : w0max is nan or inf, returning 0");
                return 0;
            }
            return w0max;
        }

        // Calcul rho0 avec la méthode de Neumann Rejection
        /*
        inline double rho0 (int n, double MF_2, double MGT_2, double MIMASSC2, double MFMASSC2, int Z, double R) {
            double somme_w0 = 0 ;
            double EPS = (MF_2) + (std::pow(LAMBDA, 2)) * (MGT_2) ;

            double w0max = w0_max(MF_2, MGT_2, MIMASSC2, MFMASSC2, Z, R) ;
            std::random_device rd;
            std::mt19937 generator(rd());
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            std::uniform_real_distribution<double> distribution_w0(0.0, w0max);

            for (int i=0 ; i <  n ; i++) {
                double weight0 = distribution_w0(generator);

                double u1 = distribution(generator) ;

                double E2 = 1. + (delta(MIMASSC2, MFMASSC2) - 1.) * u1 ;
                double E10 = delta(MIMASSC2, MFMASSC2) - E2 ;
                double BETA = std::sqrt( (1. - 1./std::pow(E2, 2)) ) ;

                if (weight0 <= w0(E2, MF_2, MGT_2, MIMASSC2, MFMASSC2, Z, R) ) {
                    somme_w0 += 1 ;
                }
            }
            double result = (delta(MIMASSC2, MFMASSC2) - 1) * w0max * somme_w0 / n ;
            return result ;
        }
        */
        inline double rho0(int n, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            DecayManager &dm = DecayManager::GetInstance();
            if (dm.configOptions.general.Verbosity >= 2)
                Info("Caluculation of rho0", 2);

            double b = delta(MIMASSC2, MFMASSC2, mode) - 1;
            double h = b / (n + 1);
            double result = 0;
            for (int i = 0; i < n; i++)
            {
                if (!std::isnan(w0(1 + h * i, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode)))
                {
                    double dx1 = 1 + i * h;
                    double dx2 = dx1 + h;
                    result += (h)*w0(0.5 * (dx1 + dx2), MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
                }
            }
            // std::cout << "res : " << w0(1+0.01, a, MIMASSC2, MFMASSC2, Z, R, mode) << "\n";
            if (std::isnan(result) || std::isinf(result))
            {
                Warning("rho0 : result is nan or inf, returning 0");
                return 0;
            }
            return result;
        }

        // Equation 5.20
        inline double wVS(double E2, double Cs, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            if (E2 == 1.0) return 0; // to avoid division by zero in beta calculation (Kinetic energy of the charged lepton is 0)
            
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double N = 0.5 * log((1. + BETA) / (1. - BETA));
            double omega = E10 * Cs;
            double zVS = (FINESTRUCTURE / PI) * (1.5 * log(PMASSC2 / EMASSC2) + 2. * ((N / BETA) - 1.) * log(2. * omega) + 2 * (N / BETA) * (1. - N) + (2. / BETA) * (Spence_function(2. * BETA / (1. + BETA))) - (3. / 8.));
            double weightVS = w0(E2, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode) * (zVS - (FINESTRUCTURE * N / PI) * (1. - std::pow(BETA, 2)) / BETA);

            if (std::isnan(weightVS) || std::isinf(weightVS))
            {
                Warning("wVS : weightVS is nan or inf, returning 0");
                return 0;
            }
            return weightVS;
        }

        // Calcul du max de wVS pour Neumann Rejection
        inline double wVS_max(double Cs, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double wVSmax = 0;
            double n = 5000;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode) - 1) / (n);
            std::vector<double> tableau_E2(n);
            for (int i = 0; i < n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                double w = wVS(tableau_E2[i], Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);

                if (!std::isnan(w) && !std::isinf(w))
                {
                    if (wVSmax < w)
                    {
                        wVSmax = w;
                    }
                }
            }
            if (std::isnan(wVSmax) || std::isinf(wVSmax))
            {
                Error("wVS_max : wVSmax is nan or inf, returning 0");
                return 0;
            }
            return wVSmax;
        }

        inline double rhoVS(int n, double Cs, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double b = delta(MIMASSC2, MFMASSC2, mode) - 1;
            double h = b / (n);
            double result = 0;
            for (int i = 0; i < n; i++)
            {
                if (!std::isnan(wVS(1 + h * i, Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode)))
                {

                    double dx1 = 1 + i * h;
                    double dx2 = dx1 + h;
                    result += (h)*wVS(0.5 * (dx1 + dx2), Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
                }
            }
            // std::cout << "b " << b << "\n";
            // std::cout << "res 2 " << wVS(b, Cs, a, MIMASSC2, MFMASSC2, Z, R, mode) << "\n";
            if (std::isnan(result) || std::isinf(result))
            {
                Warning("rhoVS : result is nan or inf, returning 0");
                return 0;
            }
            return result;
        }

        // Equation 5.18 + 5.20
        inline double w0VS(double E2, double Cs, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            return w0(E2, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode) + wVS(E2, Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
        }

        // Calcul du max de w0VS pour Neumann Rejection (spectre)
        inline double w0VS_max(double Cs, double MF, double MGT, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double w0VSmax = 0;
            double n = 5000;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode)) / (n);
            std::vector<double> tableau_E2(n);
            for (int i = 0; i < n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                double w = w0VS(tableau_E2[i], Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);

                if (!std::isnan(w) && !std::isinf(w))
                {
                    if (w0VSmax < w)
                    {
                        w0VSmax = w;
                    }
                }
            }
            if (std::isnan(w0VSmax) || std::isinf(w0VSmax))
            {
                Error("w0VS_max : w0VSmax is nan or inf, returning 0");
                return 0;
            }
            return w0VSmax;
        }

        inline double W0(double E2, double COS, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double res = BETA * E10 * E2 * M0(E2, COS, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);
            if (std::isnan(res) || std::isinf(res))
            {
                Warning("W0 : res is nan or inf, returning 0");
                return 0;
            }
            return res;
        }

        inline double W0_max(double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double W0max = 0;
            double n = 5000;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode) - 1) / (n);
            double intervalle_cos = 2 / (n);
            std::vector<double> tableau_E2(n + 1);
            std::vector<double> tableau_cos(n + 1);
            for (int i = 0; i <= n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                tableau_cos[i] = 1 - i * intervalle_cos;
            }
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n; j++)
                {
                    double w0_value = W0(tableau_E2[i], tableau_cos[j], MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);

                    if (!std::isnan(w0_value) && !std::isinf(w0_value))
                    {
                        if (W0max < w0_value)
                        {
                            W0max = w0_value;
                        }
                    }
                }
            }

            if (std::isnan(W0max) || std::isinf(W0max))
            {
                Warning("W0_max : W0max is nan or inf, returning 0");
                return 0;
            }
            return W0max;
        }

        inline double W0VS(double E2, double C, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            if (E2 == 1.0) return 0; // to avoid division by zero in beta calculation (Kinetic energy of the charged lepton is 0)
            double E10 = delta(MIMASSC2, MFMASSC2, mode) - E2;
            double BETA = std::sqrt((1. - 1. / std::pow(E2, 2)));
            double res = BETA * E10 * E2 * (M0(E2, C, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) + MVS(E2, C, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode));
            return res;
        }

        inline double W0VS_max(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            DecayManager &dm = DecayManager::GetInstance();
            if (dm.configOptions.general.Verbosity >= 2)
                Info("Calculating maximum of W0VS for radiative corrections...", 2);
            double W0VSmax = 0;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode) - 1) / (n);
            double intervalle_cos = 2. / (n);
            std::vector<double> tableau_E2(n + 1);
            std::vector<double> tableau_cos(n + 1);
            for (int i = 0; i <= n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                tableau_cos[i] = 1 - i * intervalle_cos;
            }
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n; j++)
                {
                    double w0vs_value = W0VS(tableau_E2[i], tableau_cos[j], Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);
                    if (!std::isnan(w0vs_value) && !std::isinf(w0vs_value))
                    {
                        if (W0VSmax < w0vs_value)
                        {
                            W0VSmax = w0vs_value;
                        }
                    }
                }
            }
            
            if (std::isnan(W0VSmax) || std::isinf(W0VSmax))
            {
                Error("W0VS_max : W0VSmax is nan or inf, returning 0");
                return 0;
            }

            return W0VSmax;
        }

        inline double W0VS_mean(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            double somme_W0VS = 0;
            int nout = 0;

            double intervalle_E2 = (delta(MIMASSC2, MFMASSC2, mode) - 1) / (n);
            double intervalle_cos = 2. / (n);
            std::vector<double> tableau_E2(n + 1);
            std::vector<double> tableau_cos(n + 1);
            for (int i = 0; i <= n; i++)
            {
                tableau_E2[i] = 1 + i * intervalle_E2;
                tableau_cos[i] = 1 - i * intervalle_cos;
            }
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n; j++)
                {
                    double w0vs_value = W0VS(tableau_E2[i], tableau_cos[j], Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);
                    if (!std::isnan(w0vs_value) && !std::isinf(w0vs_value))
                    {
                        somme_W0VS += w0vs_value;
                        nout += 1;
                    }
                }
            }
            // std::cout << "somme w0vs : " << somme_W0VS << "\n";
            return somme_W0VS / nout;
        }

        inline double PH(double Cs, double MF, double MGT, double a, double MIMASSC2, double MFMASSC2, int Z, double R, int mode)
        {
            DecayManager &dm = DecayManager::GetInstance();
            if (dm.configOptions.general.Verbosity >= 2)
                Info("Caluculation of PH :");
            

            double RHOH = rho_H(1e6, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);
            double RHO0 = rho0(1e6, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
            double RHOVS = rhoVS(1e6, Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, mode);
            double RHO0VS = RHO0 + RHOVS;

            if (dm.configOptions.general.Verbosity >= 2)
            {
                double DELTA_RHOH = delta_rho_H(1e6, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode);
                
                Info(Form("Cs : %.5f", Cs), 1);
                Info(Form("MF : %.5f", MF), 1);
                Info(Form("MGT : %.5f", MGT), 1);
                Info(Form("a : %.5f", a), 1);
                Info(Form("MIMASSC2 : %.5f", MIMASSC2), 1);
                Info(Form("MFMASSC2 : %.5f", MFMASSC2), 1);
                Info(Form("Δ = %.5f", delta(MIMASSC2, MFMASSC2, mode) * EMASSC2), 1);
                Info(Form("Z : %d", Z), 1);
                Info(Form("R : %.5f", R * 1e15), 1);
                Info(Form("Beta : %d", mode), 1);
                Info(Form("rH : %.5f  ±  %.5f", 100 * RHOH / (RHO0VS + RHOH), 100 * RHO0VS / pow(RHO0VS + RHOH, 2) * DELTA_RHOH), 1);
                Info(Form("rρ : %.5f", 100 * (RHOVS + RHOH) / RHO0), 1);
                Info(Form("E0VS: %.5f", 100 * W0VS_mean(1e3, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / W0VS_max(1e3, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode)), 1);
                Info(Form("EH : %.5f", 100 * WH_mean(1e5, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode) / WH_max(1e5, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, mode)), 1);
                std::cout << std::endl;
            }

            // std::cout << "rhoh : " << RHOH << "\n";
            // std::cout << "rho0 : " << RHO0 << "\n";
            // std::cout << "rhoVS : " << RHOVS << "\n";
            // std::cout << "rho0VS : " << RHO0VS << "\n";
            //  std::cout << "r_rho :" << 100*(RHOVS+RHOH)/RHO0 << "\n";
            // std::cout << "pH : " << RHOH/(RHO0VS + RHOH) << "\n";
            // std::cout << "mf : " << MF_2 << "\n";
            // std::cout << "mgt : " << MGT_2 << "\n";
            // std::cout << "wh max : " << WH_max(1000000, Cs, a, MIMASSC2, MFMASSC2, Z, R, mode) << "\n";
            // std::cout << "w0vs max : " << W0VS_max(MF_2, MGT_2, MIMASSC2, MFMASSC2) << "\n";"

            return RHOH / (RHO0VS + RHOH);
        }

    } // End of radiativecorrections namespace
} // End of CRADLE namespace
#endif
