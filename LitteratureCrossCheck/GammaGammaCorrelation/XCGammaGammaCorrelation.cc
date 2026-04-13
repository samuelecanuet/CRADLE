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

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace CRADLE;
using namespace CRADLE::utilities;
using namespace CRADLE::correlation;
using namespace CRADLE::radiativecorrections;

struct Cascade
{
    std::string Name;
    double Ji;
    double Jm;
    double Jf;

    double delta1 = 0;
    std::pair<int, int> l1;
    double delta2 = 0;
    std::pair<int, int> l2;

    TH1D *H;

    Cascade(std::string name, double ji, double jm, double jl, std::pair<int, int> l1_val, std::pair<int, int> l2_val) : Name(name), Ji(ji), Jm(jm), Jf(jl), l1(l1_val), l2(l2_val)
    {

        // Angle HISTOGRAM
        H = new TH1D(Name.c_str(), Form("%s; #theta (deg); Counts", Name.c_str()), 180*2, 0, 180);

        // l1 = std::make_pair(std::max(abs((int)(Ji - Jm)), abs((int)(Jm - Jf))), std::max(abs((int)(Ji - Jm)), abs((int)(Jm - Jf)))); // default to pure E1
        // l2 = std::make_pair(std::max(abs((int)(Jm - Jf)), abs((int)(Jf - 0))), std::max(abs((int)(Jm - Jf)), abs((int)(Jf - 0)))); // default to pure E1
    }

    Cascade(){}
};

inline double GammaGammaTF1(double *x, double *par)
{
    std::vector<double> parv = {par[0], par[1], par[2]};
    return par[3] * GammaGammaW(cos(x[0] * M_PI / 180.), parv) / GammaGammaW(cos(90. * M_PI / 180.), parv);
}

int main(int argc, char *argv[])
{
    ///////  #### COMPARISING WITH Gamma-Gamma Angular Correlation Measurements With GRIFFIN (Table 1.4) ####  ///////
    Start("Comparing values", 0);
    std::map<std::string, Cascade> DataValues;
    // Paper
    DataValues["0 -> 2 -> 0"] = Cascade("0 -> 2 -> 0", 0., 2., 0., std::make_pair(2, 0), std::make_pair(2, 0));
    DataValues["4 -> 2 -> 0"] = Cascade("4 -> 2 -> 0", 4., 2., 0., std::make_pair(2, 0), std::make_pair(2, 0));
    DataValues["1 -delta> 2 -> 0"] = Cascade("1 -delta> 2 -> 0", 1., 2., 0., std::make_pair(1, 2), std::make_pair(2, 0)); DataValues["1 -delta> 2 -> 0"].delta1 = 5.;
    DataValues["3 -> 2 -> 0"] = Cascade("3 -> 2 -> 0", 3., 2., 0., std::make_pair(1, 0), std::make_pair(2, 0));
    DataValues["2 -> 2 -> 0"] = Cascade("2 -> 2 -> 0", 2., 2., 0., std::make_pair(1, 0), std::make_pair(2, 0));
    DataValues["2 -delta> 2 -> 0"] = Cascade("2 -delta> 2 ->_0", 2., 2., 0., std::make_pair(1, 2), std::make_pair(2, 0)); DataValues["2 -delta> 2 -> 0"].delta1 = -1;

    for (auto &data : DataValues)
    {
        Info(Form("%s : ", data.second.Name.c_str()), 1);
        std::vector<double> ak = correlation::CaluclateGammaCoefficient_a(std::abs(data.second.Ji), std::abs(data.second.Jm), std::abs(data.second.Jf), data.second.l1, data.second.delta1, data.second.l2, data.second.delta2);
        for (size_t i = 0; i < ak.size(); i++)
        {
            Info(Form("a%d : %.4f", 2 * i, ak[i]), 2);
        }
    }

    std::cout << std::endl; 

    ///////  #### COMPARISING WITH Gamma-Gamma Angular Correlation Measurements With GRIFFIN (Figure 1.9) ####  ///////
    // DATA
    Start("Comparing plots", 0);
    std::map<std::string, Cascade> Data;
    Data["Transitions_A"] = Cascade("4 -> 2 -> 0", 4., 2., 0., std::make_pair(2, 0), std::make_pair(2, 0));
    Data["Transitions_B"] = Cascade("0 -> 2 -> 0", 0., 2., 0., std::make_pair(2, 0), std::make_pair(2, 0)); // (ok)
    Data["Transitions_C"] = Cascade("2 -> 1 -> 0", 2., 1., 0., std::make_pair(3, 0), std::make_pair(1, 0)); // (ok)
    Data["Transitions_D"] = Cascade("1 -> 2 -> 0", 1., 2., 0., std::make_pair(2, 0), std::make_pair(2, 0)); // (ok)
    Data["Transitions_E"] = Cascade("3 -> 1 -> 0", 3., 1., 0., std::make_pair(2, 0), std::make_pair(1, 0)); // (ok)
    Data["Transitions_F"] = Cascade("2 -> 2 -> 0", 2., 2., 0., std::make_pair(1, 0), std::make_pair(2, 0)); // (ok)

    
    // OUPUT FILE
    TFile *fout = new TFile("../LitteratureCrossCheck/GammaGammaCorrelation/Cascades.root", "RECREATE");
    TCanvas *cAna = new TCanvas("Analytical", "Analytical", 800, 600);
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9); 
    int N = 1e5;

    int counter = 0;
    for (auto &data : Data)
    {
        Start(data.second.Name);

        // Set Up CRADLE++
        DecayManager &dm = DecayManager::GetInstance();
        dm.configOptions.decay.InFlightDecay = false;        // avoid doppler broadering
        dm.configOptions.decay.GammaGammaCorrelation = true; // enable gamma gamma correlation
        dm.RegisterBasicDecayModes();
        dm.RegisterBasicParticles();

        // State_i (Ji, Ei)
        // γ1 ( E1 = Ei-Em, deltai )
        // State_m (Jm, Em)
        // γ2 ( E2 = Em, deltaf )
        // State_f (Jf, Ef=0)

        // STATES
        double Ei = 3000.+counter;
        double Ji = data.second.Ji;
        double Em = 2000.+counter;
        double Jm = data.second.Jm;
        double Ef = 0.;
        double Jf = data.second.Jf;

        Info("States : ", 1);
        Info(Form("Ji : %.1f", Ji), 2);
        Info(Form("Jm : %.1f", Jm), 2);
        Info(Form("Jf : %.1f", Jf), 2);

        // GAMMA
        std::pair<int, int> L1 = data.second.l1;    
        std::pair<int, int> L2 = data.second.l2;    
        double delta1 = data.second.delta1; // mixing ratio for gamma 1
        double delta2 = data.second.delta2; // mixing ratio for gamma 2
        
        Info("Gammas : ", 1);
        Info(Form("L1 : %d, %d", L1.first, L1.second), 2);
        Info(Form("L2 : %d, %d", L2.first, L2.second), 2);
        Info(Form("delta1 : %.2f", delta1), 2);
        Info(Form("delta2 : %.2f", delta2), 2);

        // ak
        std::vector<double> ak = correlation::CaluclateGammaCoefficient_a(std::abs(Ji), std::abs(Jm), std::abs(Jf), L1, delta1, L2, delta2);
        double W_max = correlation::MaxAnalyticalGammaCorrelation(ak);
        std::vector<std::vector<double>> *dist = new std::vector<std::vector<double>>();
        dist->push_back(ak);
        Info("ak coefficients : ", 1);
        for (size_t i = 0; i < ak.size(); i++)        {
            Info(Form("a%d : %.4f", (int)i, ak[i]), 2);
        }

        Particle *Initstate = new Particle(1000180320, GetApproximateMass(18, 32), 18, 14, Ji, Ei);
        dm.RegisterParticle(Initstate);
    
        std::ostringstream oss_saved;

        clock_t start = clock();
        for (int i = 0; i < N; i++)
        {
            ProgressBar(i, N - 1, start, "", 100);
            // Fictive Nucleus
            Particle *initstate = dm.GetNewParticle(1000180320);
            // decay 1
            std::vector<Particle *> mstate = (&dm.GetDecayMode("Gamma"))->Decay(initstate, Ei - Em, Em);
            Particle *GAMMA1 = mstate.at(1); // first gamma

            std::ostringstream oss;
            oss << "GammaGamma:" << "Z" << initstate->GetCharge() << "A" << initstate->GetCharge() + initstate->GetNeutrons() << "Ei" << Em + GAMMA1->GetKinEnergy() << "Em" << Em << "Ef" << Ef;
            dm.RegisterChannelPropreties(oss.str(), dist, W_max, Ji, Jf, Jm);

            // avoid Doppler Broadering
            ublas::vector<double> m(4, 0);
            m(0) = mstate.at(0)->GetMass();
            mstate.at(0)->SetMomentum(m);
            // decay 2
            std::vector<Particle *> fstate = (&dm.GetDecayMode("Gamma"))->Decay(mstate.at(0), Em - Ef, Ef);
            Particle *GAMMA2 = fstate.at(1); // second gamma

            ublas::vector<double> GAMMA1_dir = GAMMA1->Get3Momentum();
            ublas::vector<double> GAMMA2_dir = GAMMA2->Get3Momentum();

            double costheta = inner_prod(GAMMA1_dir, GAMMA2_dir) / utilities::GetNorm(GAMMA1_dir) / utilities::GetNorm(GAMMA2_dir);
            data.second.H->Fill(acos(costheta) * 180. / M_PI);

            delete initstate;
            delete GAMMA1;
            delete GAMMA2;
        }        

        fout->cd();
        TCanvas *c = new TCanvas(data.second.Name.c_str(), data.second.Name.c_str(), 800, 600);
        TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.9);
        c->cd();
        gStyle->SetOptStat(0);
        // function
        TF1 *fana = new TF1(data.second.Name.c_str(), GammaGammaTF1, 0, 180, 4);
        fana->SetName(Form("%s", data.second.Name.c_str()));
        fana->SetTitle(Form("%s", data.second.Name.c_str()));
        fana->SetParameters(ak[0], ak[1], ak[2], 1.);
        // hist
        data.second.H->Scale(((TH1D *)fana->GetHistogram())->Integral("width") / data.second.H->Integral("width"));
        data.second.H->Draw("E");
        fana->Draw("SAME");
        l->AddEntry(data.second.H, "CRADLE++", "l");
        l->AddEntry(fana, "W(#theta)", "l");
        c->Write();

        cAna->cd(); 
        fana->SetLineColor(cAna->GetListOfPrimitives()->GetEntries() + 1);
        fana->SetMinimum(0);
        fana->SetMaximum(5.);
        if (counter == 0)
            fana->Draw();
        else
            fana->Draw("SAME");
        
        leg->AddEntry(fana, data.second.Name.c_str(), "l");

        counter++;
    }
    cAna->cd();
    leg->Draw();
    cAna->Write();
    fout->Close();

    // run python to create the final plot
    system("cd ../LitteratureCrossCheck/GammaGammaCorrelation/ ; python3 PlotCascade.py ; cd -");
    //

    return 0;
}
