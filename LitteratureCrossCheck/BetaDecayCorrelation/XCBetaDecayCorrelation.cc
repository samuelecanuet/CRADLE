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
#include "CRADLE/DecayManager.hh"
#include "CRADLE/DecayChannel.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/Particle.hh"

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace CRADLE;
using namespace CRADLE::utilities;
using namespace CRADLE::correlation;
using namespace CRADLE::radiativecorrections;

inline double F_TF1_COSTHETAEJ(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    ublas::vector<double> ChargedLepton_Dir(3);
    ChargedLepton_Dir(0) = std::sqrt(1 - std::pow(costheta, 2));
    ChargedLepton_Dir(1) = 0;
    ChargedLepton_Dir(2) = costheta;
    ublas::vector<double> NeutralLepton_Dir(3);
    NeutralLepton_Dir(0) = 0;
    NeutralLepton_Dir(1) = 1;
    NeutralLepton_Dir(2) = 0;
    ublas::vector<double> polDir(3);
    polDir(0) = 0;
    polDir(1) = 0;
    polDir(2) = 1;
    double ChargedLepton_Energy = par[7];

    return par[6] * correlation::CalculateAngularCorrelationFactor(a, b, c, A, B, D, ChargedLepton_Energy, ChargedLepton_Dir, NeutralLepton_Dir, polDir);
}

inline double F_TF1_COSTHETAEJ_Ana(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    double ChargedLepton_Energy = par[7];

    double res = ( 1 + b * utilities::EMASSC2 / ChargedLepton_Energy + A * costheta * std::sqrt(1. - utilities::EMASSC2 * utilities::EMASSC2 / ChargedLepton_Energy / ChargedLepton_Energy)) / ( 1. + b * utilities::EMASSC2 / ChargedLepton_Energy) / 2;

    return par[6] * res;
}

inline double F_TF1_COSTHETAENU(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    ublas::vector<double> ChargedLepton_Dir(3);
    ChargedLepton_Dir(0) = std::sqrt(1 - std::pow(costheta, 2));
    ChargedLepton_Dir(1) = 0;
    ChargedLepton_Dir(2) = costheta;
    ublas::vector<double> NeutralLepton_Dir(3);
    NeutralLepton_Dir(0) = 0;
    NeutralLepton_Dir(1) = 0;
    NeutralLepton_Dir(2) = 1;
    ublas::vector<double> polDir(3);
    polDir(0) = 0;
    polDir(1) = 0;
    polDir(2) = 1;
    double ChargedLepton_Energy = par[7];

    return par[6] * correlation::CalculateAngularCorrelationFactor(a, b, c, A, B, D, ChargedLepton_Energy, ChargedLepton_Dir, NeutralLepton_Dir, polDir);
}
inline double F_TF1_COSTHETAENU_Ana(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    double ChargedLepton_Energy = par[7];

    double res = ( 1. + b * utilities::EMASSC2 / ChargedLepton_Energy + (a) * costheta * std::sqrt(1. - utilities::EMASSC2 * utilities::EMASSC2 / ChargedLepton_Energy / ChargedLepton_Energy) ) / ( 1. + b * utilities::EMASSC2 / ChargedLepton_Energy) / 2.; // normalisation factor to have the integral over costheta from -1 to 1 equal to 1 + b * EMASSC2 / E

    return par[6] * res;
}

inline double F_TF1_COSTHETANUJ(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    ublas::vector<double> ChargedLepton_Dir(3);
    ChargedLepton_Dir(0) = 0;
    ChargedLepton_Dir(1) = 1;
    ChargedLepton_Dir(2) = 0;
    ublas::vector<double> NeutralLepton_Dir(3);
    NeutralLepton_Dir(0) = std::sqrt(1 - std::pow(costheta, 2));
    NeutralLepton_Dir(1) = 0;
    NeutralLepton_Dir(2) = costheta;
    ublas::vector<double> polDir(3);
    polDir(0) = 0;
    polDir(1) = 0;
    polDir(2) = 1;
    double ChargedLepton_Energy = par[7];

    return par[6] * correlation::CalculateAngularCorrelationFactor(a, b, c, A, B, D, ChargedLepton_Energy, ChargedLepton_Dir, NeutralLepton_Dir, polDir);
}

inline double F_TF1_COSTHETANUJ_Ana(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double costheta = x[0];
    double ChargedLepton_Energy = par[7];

    double res = ( 1. + b * utilities::EMASSC2 / ChargedLepton_Energy + B * costheta ) / ( 1. + b * utilities::EMASSC2 / ChargedLepton_Energy) / 2.; 

    return par[6] * res;
}

inline double F_TF1_PHI(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double phi = x[0] * utilities::PI / 180.;
    ublas::vector<double> ChargedLepton_Dir(3);
    ChargedLepton_Dir(0) = std::cos(phi);
    ChargedLepton_Dir(1) = std::sin(phi);
    ChargedLepton_Dir(2) = 0;
    ublas::vector<double> NeutralLepton_Dir(3);
    NeutralLepton_Dir(0) = std::sin(phi);
    NeutralLepton_Dir(1) = std::cos(phi);
    NeutralLepton_Dir(2) = 0;
    ublas::vector<double> polDir(3);
    polDir(0) = 0;
    polDir(1) = 0;
    polDir(2) = 1;
    double ChargedLepton_Energy = par[7];

    return par[6] * correlation::CalculateAngularCorrelationFactor(a, b, c, A, B, D, ChargedLepton_Energy, ChargedLepton_Dir, NeutralLepton_Dir, polDir);
}

inline double F_TF1_PHI_Ana(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double A = par[3];
    double B = par[4];
    double D = par[5];

    double phi = x[0] * utilities::PI / 180.;
   
    double res = par[6] * (1. + b * utilities::EMASSC2 / par[7] + M_PI*M_PI/16. * sqrt(1. - pow(utilities::EMASSC2 / par[7], 2)) * ((a + c/3.) * cos(phi) + D*sin(phi))) / (1. + b * utilities::EMASSC2 / par[7]) / 2.* M_PI; // normalisation factor to have the integral over phi from 0 to 360 equal to 1 + b * EMASSC2 / E
    return res;
}

double ComputeChi2Integral(TH1* h, TF1* f) {
    double chi2 = 0.0;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double y = h->GetBinContent(i);
        double err = h->GetBinError(i);
        if (err <= 0) continue;

        double xlow = h->GetBinLowEdge(i);
        double xup  = xlow + h->GetBinWidth(i);

        double fval = f->Integral(xlow, xup) / (xup - xlow);

        chi2 += std::pow(y - fval, 2) / (err * err);
    }

    return chi2 / (h->GetNbinsX());
}

inline void PlotPad(TH1D* H, TF1* f, std::string Sensibility_text_str, TCanvas * c, int pad_number)
{
    c->cd(pad_number);
    double ypad = 0.2;
    TPad *pad_PDF = new TPad("pad_PDF", "pad_PDF", 0, ypad, 1, 1);
    pad_PDF->Draw();
    pad_PDF->cd();
    if (pad_number != 4)
        H->Scale(1. / H->Integral("width"));
    else
        H->Scale(1. / H->Integral("width")*f->Integral(H->GetXaxis()->GetXmin(), H->GetXaxis()->GetXmax()));

    H->GetYaxis()->SetRangeUser(0, -1111);
    H->SetLineColor(kBlack);
    H->Draw("E");
    f->SetLineColor(kRed);
    f->Draw("SAME");
    TLatex *chi2_text = new TLatex(0.2, 0.85, Form("#chi^{2} = %.2f", ComputeChi2Integral(H, f)));
    chi2_text->SetNDC();
    chi2_text->SetTextColor(kRed);
    chi2_text->Draw("SAME");
    TLatex *Sensitivity_text = new TLatex(0.4, 0.85, Form("Sensitivity: %s", Sensibility_text_str.c_str()));
    Sensitivity_text->SetNDC();
    Sensitivity_text->SetTextColor(kBlue);
    Sensitivity_text->Draw("SAME");
    c->cd(pad_number);
    TPad *pad_residuals = new TPad("pad_residuals", "pad_residuals", 0, 0, 1, ypad);
    pad_residuals->Draw();
    pad_residuals->cd();
    TGraphErrors *g_residuals = new TGraphErrors();
    for (int i = 1; i <= H->GetNbinsX(); ++i) {
        double y = H->GetBinContent(i);
        double err = H->GetBinError(i);
        if (err <= 0) continue;
        double xcenter =  H->GetBinCenter(i);

        double fval = f->Eval(xcenter);

        g_residuals->SetPoint(g_residuals->GetN(), xcenter, (y - fval) / fval);
        g_residuals->SetPointError(g_residuals->GetN() - 1, 0, err / fval);
    }
    g_residuals->SetMarkerStyle(20);
    g_residuals->SetMarkerSize(0.8);
    g_residuals->SetMarkerColor(kBlack);
    g_residuals->GetXaxis()->SetLimits(H->GetXaxis()->GetXmin(), H->GetXaxis()->GetXmax());
    g_residuals->Draw("AP");
    TLine *line = new TLine(f->GetXmin(), 0, f->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->Draw("SAME");
    Info(Form("𝜒² = %.2f", ComputeChi2Integral(H, f)), 2);    

    return;
}

struct TEST
{
    std::string name;
    std::string configfilename;
    int Z;
    int A;
    int BetaSign;
    double Qbeta;
    double Ji;
    double Jf;

    TEST(std::string name, std::string configfilename, int Z, int A, int BetaSign, double Qbeta, double Ji, double Jf) : name(name), configfilename(configfilename), Z(Z), A(A), BetaSign(BetaSign), Qbeta(Qbeta), Ji(Ji), Jf(Jf)
    {
        if (BetaSign == 1)
            this->Qbeta = Qbeta - 2 * utilities::EMASSC2;
        else
            this->Qbeta = Qbeta;
    }
    TEST() {}
};

int main(int argc, const char **argv)
{
    std::map<std::string, TEST> Scenarios;

    Scenarios["SM"] = TEST("60Co", "config_60Co_SM.txt", 27, 60, -1, 317, 5, 4);
    Scenarios["ReTensor"] = TEST("60Co", "config_60Co_ReTensor.txt", 27, 60, -1, 317, 5, 4);
    Scenarios["ImTensor"] = TEST("60Co", "config_60Co_ImTensor.txt", 27, 60, -1, 317, 5, 4);
    Scenarios["BSM"] = TEST("39Ca", "config_39Ca_BSM.txt", 20, 39, +1, 6524.488, 1.5, 1.5);

    for (auto &scenario : Scenarios)
    {
        Start("Scenario: " + scenario.first, 0);

        std::string configfilename = "../LitteratureCrossCheck/BetaDecayCorrelation/" + scenario.second.configfilename;
        std::string outputfilename = "output_" + scenario.second.name + "_" + scenario.first + ".root";
        std::string reltivepath = "/mnt/hgfs/shared-2/";

        // RUNNING CRADLE++
        // system(("CRADLE++ nucleus -n " + scenario.second.name + " -A " + std::to_string(scenario.second.A) + " -Z " + std::to_string(scenario.second.Z) + " general -l 10000000 -o " + reltivepath + outputfilename + " -c " + configfilename + " -t 7").c_str());
        // system(("CRADLE++ nucleus -n " + scenario.second.name + " -A " + std::to_string(scenario.second.A) + " -Z " + std::to_string(scenario.second.Z) + " general -l 1 -c " + configfilename + " -t 1 -v 2").c_str());

        // OUTPUT FILE
        std::string readfilename = reltivepath + outputfilename.substr(0, outputfilename.find(".root")) + "_reader.root";
        TFile *file = new TFile(Form("%s", readfilename.c_str()), "RECREATE");
        TH1D *H_CosTheta_e_nu = new TH1D("H_CosTheta_e_nu", "Cosine of angle between electron and neutrino; cos(#theta_{e#nu}); Normalized Counts", 250, -1, 1);
        TH1D *H_CosTheta_e_j = new TH1D("H_CosTheta_e_j", "Cosine of angle between electron and polarisation; cos(#theta_{ej}); Normalized Counts", 250, -1, 1);
        TH1D *H_CosTheta_nu_j = new TH1D("H_CosTheta_nu_j", "Cosine of angle between neutrino and polarisation; cos(#theta_{#nu j}); Normalized Counts", 250, -1, 1);
        TH1D *H_phi = new TH1D("H_phi", "Azimuthal angle between electron and neutrino; #phi (deg); Normalized Counts", 180, -180, 180);

        // INPUT FILE
        TFile *f = new TFile(Form("%s", (reltivepath + outputfilename).c_str()), "READ");
        if (!f->IsOpen())
            Error("Could not open file " + outputfilename);

        // Setting up TTreeReader
        TTree *t = (TTree *)f->Get("ParticleTree");
        if (!t)
            Error("Could not find TTree 'ParticleTree' in file " + outputfilename);
        TTreeReader *Reader = new TTreeReader(t);
        TTreeReaderArray<double> *Time = new TTreeReaderArray<double>(*Reader, "time");
        TTreeReaderArray<int> *Code = new TTreeReaderArray<int>(*Reader, "code");
        TTreeReaderArray<double> *Energy = new TTreeReaderArray<double>(*Reader, "energy");
        TTreeReaderArray<double> *Excitation_energy = new TTreeReaderArray<double>(*Reader, "excitation_energy");
        TTreeReaderArray<double> *p = new TTreeReaderArray<double>(*Reader, "p");
        TTreeReaderArray<double> *Px = new TTreeReaderArray<double>(*Reader, "px");
        TTreeReaderArray<double> *Py = new TTreeReaderArray<double>(*Reader, "py");
        TTreeReaderArray<double> *Pz = new TTreeReaderArray<double>(*Reader, "pz");

        clock_t start = clock();
        int Entries = Reader->GetEntries();
        int step = std::max(1, Entries / 1000);

        TH1D *H_EnergyMean = new TH1D("H_EnergyMean", "Mean energy of particles; Energy (keV); Counts", 10000, 0, 1);
        // Reading TTree
        int Verbosity = 0;
        while (Reader->Next() && Reader->GetCurrentEntry() < Entries)
        {
            if (Verbosity == 0)
                ProgressBar(Reader->GetCurrentEntry(), Entries - 1, start, "", step);

            if (Verbosity > 0)
                Info("#### Event: " + std::to_string(Reader->GetCurrentEntry()) + " #### (N = " + std::to_string(Code->GetSize()) + ")", 0);

            // Lopping on particle in the event
            for (int i = 0; i < Code->GetSize(); i++)
            {
                if (Verbosity > 0)
                    Info(Form("%s \t E = %.1f keV \t Ex = %.1f keV \t px = %.2f \t py = %.2f \t pz = %.2f", PDGtoName((*Code)[i]).c_str(), (*Energy)[i], (*p)[i], (*Px)[i], (*Py)[i], (*Pz)[i]), 1);
            }

            if (abs((*Energy)[1] + (*Energy)[2] - scenario.second.Qbeta) < 1) // 60Co case
            {
                int index_e = 2;
                int index_nu = 1;
                if (abs((*Code)[1]) == 11) // Electron is at index 1
                {
                    index_e = 1;
                    index_nu = 2;
                }

                H_EnergyMean->Fill(std::sqrt(1 - utilities::EMASSC2 * utilities::EMASSC2 / ((*Energy)[index_e] + utilities::EMASSC2) / ((*Energy)[index_e] + utilities::EMASSC2)));

                // cos(theta) between electron and neutrino
                ublas::vector<double> e_dir(3);
                e_dir[0] = (*Px)[index_e];
                e_dir[1] = (*Py)[index_e];
                e_dir[2] = (*Pz)[index_e];
                ublas::vector<double> nu_dir(3);
                nu_dir[0] = (*Px)[index_nu];
                nu_dir[1] = (*Py)[index_nu];
                nu_dir[2] = (*Pz)[index_nu];
                double costheta = inner_prod(e_dir, nu_dir) / utilities::GetNorm(e_dir) / utilities::GetNorm(nu_dir);
                H_CosTheta_e_nu->Fill(costheta);

                // cos(theta) between electron and z
                ublas::vector<double> j_dir(3);
                j_dir[0] = 0;
                j_dir[1] = 0;
                j_dir[2] = 1.;
                costheta = inner_prod(e_dir, j_dir) / utilities::GetNorm(e_dir) / utilities::GetNorm(j_dir);
                H_CosTheta_e_j->Fill(costheta);

                // cos(theta) between neutrino and z
                costheta = inner_prod(nu_dir, j_dir) / utilities::GetNorm(nu_dir) / utilities::GetNorm(j_dir);
                H_CosTheta_nu_j->Fill(costheta);

                // phi between electron and neutrino
                // D * beta_e * inner_prod(polDir, CrossProduct(elDir, enuDir));
                double phi = std::atan2(e_dir[0] * nu_dir[1] - e_dir[1] * nu_dir[0], e_dir[0] * nu_dir[0] + e_dir[1] * nu_dir[1]);
                H_phi->Fill(phi * 180. / M_PI);
            }

            if (Verbosity > 0)
                std::cout << std::endl;
        }

        file->cd();

        // Fake sim to calculate correlation coefficients
        DecayManager &dm = DecayManager::GetInstance();
        int argc = 13;
        const char *argv[13] = {"CRADLE++", "nucleus", "-n", scenario.second.name.c_str(), "-Z", std::to_string(scenario.second.Z).c_str(), "-A", std::to_string(scenario.second.A).c_str(), "general", "-l", "1", "-v", "0"};
        dm.Initialise(configfilename, argc, argv);   
        //

        double mf, mgt, mixing_ratio;
        utilities::FindMatrixElement(scenario.second.Z, scenario.second.A, 0, scenario.second.Z - scenario.second.BetaSign, scenario.second.A, 0, mf, mgt, mixing_ratio);
        Info(Form("Matrix elements: M_F = %.4f, M_GT = %.4f", mf, mgt), 2);
        double Ji = scenario.second.Ji;
        double Jf = scenario.second.Jf;
        int Recoil_Z = scenario.second.Z - scenario.second.BetaSign;
        int BetaSign = scenario.second.BetaSign;
        Info(Form("Initial state: J = %.1f, Z = %d, A = %d", Ji, scenario.second.Z, scenario.second.A), 2);
        Info(Form("Final state: J = %.1f, Z = %d, A = %d", Jf, scenario.second.Z - scenario.second.BetaSign, scenario.second.A), 2);
        double mean_beta = H_EnergyMean->GetMean();
        double mean_Energy = sqrt(utilities::EMASSC2 * utilities::EMASSC2 / (1. - mean_beta * mean_beta));
        Info(Form("Mean energy of the charged lepton: %.2f keV", mean_Energy), 2);
        double a = correlation::CalculateBetaNeutrinoAsymmetry(mf, mgt, mean_Energy, Recoil_Z, -BetaSign);
        double cc = correlation::CalculateAlignmentCorrelation(mf, mgt, Ji, Jf, -BetaSign, Recoil_Z, mean_Energy);
        double b = correlation::CalculateFierz(mf, mgt, Recoil_Z, -BetaSign);
        double A = correlation::CalculateBetaAssymetry(mf, mgt, Ji, Jf, -BetaSign, Recoil_Z, mean_Energy);
        double B = correlation::CalculateNeutrinoAssymetry(mf, mgt, Ji, Jf, -BetaSign, Recoil_Z, mean_Energy);
        double D = correlation::CalculateDTripleCorrelation(mf, mgt, Ji, Jf, -BetaSign, Recoil_Z, mean_Energy);
        Info(Form("Calculated correlation coefficients: a = %.4f, b = %.4f, c = %.4f, A = %.4f, B = %.4f, D = %.4f", a, b, cc, A, B, D), 2);
        TCanvas *c = new TCanvas(Form("c_%s", scenario.first.c_str()), Form("c_%s", scenario.first.c_str()), 1200, 800);

        double y = 0.85;
        double x1 = 0.15;
        double x2 = 0.4;
        c->Divide(2, 2);
        TF1 *f_costhetaej = new TF1("f_costhetaej", F_TF1_COSTHETAEJ_Ana, -1, 1, 8);
        f_costhetaej->SetParameters(a, b, -cc, A, B, D, 1, mean_Energy);
        TF1 *f_costhetaenu = new TF1("f_costhetaenu", F_TF1_COSTHETAENU_Ana, -1, 1, 8);
        f_costhetaenu->SetParameters(a, b, -cc, A, B, D, 1, mean_Energy);
        TF1 *f_costhetanuj  = new TF1("f_costhetanuj", F_TF1_COSTHETANUJ_Ana, -1, 1, 8);
        f_costhetanuj->SetParameters(a, b, -cc, A, B, D, 1, mean_Energy);
        TF1 *f_phi = new TF1("f_phi", F_TF1_PHI_Ana, -180, 180, 8);
        f_phi->SetParameters(a, b, -cc, A, B, D, 1, mean_Energy);
        PlotPad(H_CosTheta_e_j, f_costhetaej, "A", c, 1);
        PlotPad(H_CosTheta_e_nu, f_costhetaenu, "a", c, 2);
        PlotPad(H_CosTheta_nu_j, f_costhetanuj, "B", c, 3);
        PlotPad(H_phi, f_phi, "a, c, D", c, 4);
        c->Write();
        f->Close();

        file->Close();

        // system(("python3 ../LitteratureCrossCheck/BetaDecayCorrelation/PlotCorrelation.py " + readfilename).c_str());

    }
    

    return 0;
}