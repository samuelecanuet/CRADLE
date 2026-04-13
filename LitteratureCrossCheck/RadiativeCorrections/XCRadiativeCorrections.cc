#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <iomanip>

#include "CRADLE/Utilities.hh"
#include "CRADLE/Messenger.hh"
#include "CRADLE/RadiativeCorrections.hh"

using namespace CRADLE;
using namespace CRADLE::utilities;
using namespace CRADLE::correlation;
using namespace CRADLE::radiativecorrections;

struct Nucleus
{
    int Z;
    int A;
    double MF;
    double MGT;
    int BetaType; // 0 for pure Fermi, 1 for pure Gamow-Teller, 2 for mixed transition
    int BetaSign; // -1 for beta minus, +1 for beta plus
    double Ex_Daughter; // Excitation energy of the daughter nucleus in keV

    Nucleus(int Z_, int A_, double MF_, double MGT_, int BetaType_, int BetaSign_, double Ex_Daughter_ = 0)
        : Z(Z_), A(A_), MF(MF_), MGT(MGT_), BetaType(BetaType_), BetaSign(BetaSign_), Ex_Daughter(Ex_Daughter_) {}
        
    Nucleus() : Z(0), A(0), MF(0), MGT(0), BetaType(0), BetaSign(0), Ex_Daughter(0) {}       
};

int main(int argc, char* argv[]) {

    // All the RadiativeCorrections based on the 4 body decay are extracted from F.Gluck 1997 (10.1016/S0010-4655(96)00168-3)
    // The result from this paper will be crosschecked with the CRADLE++ RadioactiveCorrections implementation.
    double Cs = 0.001; // Soft Bremstrahlung cut-off
    double n = 1e6; // Number Sample

    // CREATE A LIST OF NUCLEI TO CHECK
    std::map<std::string, Nucleus> List_Of_XC_Nuclei;
    std::vector<std::string> order;
    auto add = [&](const std::string& key, const Nucleus& val)
    {
        List_Of_XC_Nuclei.emplace(key, val);
        order.push_back(key);
    };

    add("n", Nucleus(  0,  1,  sqrt(1),    sqrt(3),    2, -1));
    add("π", Nucleus(     0,  1,  sqrt(2),    0.,         0, -1));
    add("6He", Nucleus(      2,  6,  0,          sqrt(6),    2, -1));
    add("14O", Nucleus(      8,  14, sqrt(2),    0.,         0, 1, 2312.8));
    add("32P", Nucleus(     15,  32, 0.,    sqrt(2),         1, -1));
    add("32Ar", Nucleus(     18,  32, sqrt(2),    0.,         0, 1, 5046.3));
    

    // START CROSS CHECK
    for (const auto& key : order)
    {
        Start(NiceNucleusName(key));

        // Parameters
        double Z = List_Of_XC_Nuclei[key].Z;
        double A = List_Of_XC_Nuclei[key].A;
        double MF = List_Of_XC_Nuclei[key].MF;
        double MGT = List_Of_XC_Nuclei[key].MGT;
        int betaType = List_Of_XC_Nuclei[key].BetaType;
        int mode = List_Of_XC_Nuclei[key].BetaSign;
        double Ex_Daughter = List_Of_XC_Nuclei[key].Ex_Daughter;

        Info("Parameters : ", 1);
        Info(Form("Z : %d", (int)Z), 2);
        Info(Form("A : %d", (int)A), 2);
        Info(Form("MF : %.5f", MF), 2);
        Info(Form("MGT : %.5f", MGT), 2);
        Info(Form("Beta Type : %d", betaType), 2);
        Info(Form("Beta Sign : %d", mode), 2);
        Info(Form("Excitation Energy of Daughter : %.1f keV", Ex_Daughter), 2);
        
        // Calculated parameters
        double R = ApproximateRadius(A);
        double a = CalculateBetaNeutrinoAsymmetry(MF, MGT, CRADLE::utilities::EMASSC2+1, Z, mode);
        double MIMASSC2 = key == "n" ? CRADLE::utilities::NMASSC2 : key == "π" ? 134976.8 : GetAMEMass("../Nuclear_Databases/AMEdata.txt", Z, A);
        double MFMASSC2 = key == "n" ? CRADLE::utilities::PMASSC2+CRADLE::utilities::EMASSC2 : key == "π" ? 134976.8 - 4154 : GetAMEMass("../Nuclear_Databases/AMEdata.txt", Z - mode, A) + Ex_Daughter;
        
        Info("Calculated Parameters : ", 1);
        Info(Form("R : %.5f fm", R*1e15), 2);
        Info(Form("a : %.5f", a), 2);
        Info(Form("MIMASSC2 : %.5f keV", MIMASSC2), 2);
        Info(Form("MFMASSC2 : %.5f keV", MFMASSC2), 2);

        // Observables
        double RHOH = rho_H(n, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode) ;
        double DELTA_RHOH = delta_rho_H(n, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode) ;
        double RHO0 = rho0(n, MF, MGT, MIMASSC2, MFMASSC2, Z, R, betaType, mode);
        double RHOVS = rhoVS(n, Cs, MF, MGT, MIMASSC2, MFMASSC2, Z, R, betaType, mode) ;
        double RHO0VS = RHO0 + RHOVS ;
        Info("Observables : ", 1);
        Info(Form("Δ = %.5f keV", delta(MIMASSC2, MFMASSC2, mode)*CRADLE::utilities::EMASSC2), 2);
        Info(Form("rH : %.5f  ±  %.5f", 100*RHOH/(RHO0VS + RHOH), 100*RHO0VS/pow(RHO0VS+RHOH, 2)*DELTA_RHOH), 2);
        Info(Form("rρ : %.5f", 100*(RHOVS+RHOH)/RHO0), 2);
        double mean = W0VS_mean(1e3, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode);
        double max = W0VS_max(1e3, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode);
        Info(Form("E0VS: %.5f", 100. * mean / max), 2);
        Info(Form("EH : %.5f", 100*WH_mean(1e5, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode) / WH_max(1e5, Cs, MF, MGT, a, MIMASSC2, MFMASSC2, Z, R, betaType, mode)), 2);
        //
    }

    return 0;
}
