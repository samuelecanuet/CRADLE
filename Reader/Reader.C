/// ## Reader.C ## ///
// Reading Tree from generated ROOT file by CRADLE++

#include "TFile.h"
#include "../include/CRADLE/Messenger.hh"
#include "PDGcode.hh"

map<int, TH1D*> H_time;
map<int, TH1D*> H_E;
map<int, TH1D*> H_Ex;
map<int, TH1D*> H_px;
map<int, TH1D*> H_py;
map<int, TH1D*> H_pz;
map<int, int> ParticleCounter;

struct Correlation
{
    int PDG1;
    int PDG2;
    double E1;
    double E2;

    TH1D *H_Angle = nullptr;

    Correlation(int pdg1, int pdg2, double e1, double e2) : PDG1(pdg1), PDG2(pdg2), E1(e1), E2(e2) {}
    Correlation(string name1, string name2, double e1, double e2) : PDG1(NametoPDG(name1)), PDG2(NametoPDG(name2)), E1(e1), E2(e2) {}

    void Fill(double angle)
    {
        if (!H_Angle)
            H_Angle = new TH1D(Form("H_Angle_%s(%.0fkeV)_%s(%.0fkeV)", PDGtoName(PDG1).c_str(), E1, PDGtoName(PDG2).c_str(), E2), Form("Angle distribution between %s and %s; cos(#theta); Counts", PDGtoName(PDG1).c_str(), PDGtoName(PDG2).c_str()), 1000, -1, 1);
        H_Angle->Fill(angle);
    }
};

void InitHistogram(int pdg)
{
    // Counting
    if (ParticleCounter.count(pdg) == 0)
        ParticleCounter[pdg] = 1;
    else
        ParticleCounter[pdg]++;

    // Init HIstograms for each particle
    if (H_time.count(pdg) != 0)
        return;

    H_time[pdg] = new TH1D(Form("H_time_%s", PDGtoName(pdg).c_str()), Form("Time distribution for %s; Time (s); Counts", PDGtoName(pdg).c_str()), 100, 0, 1);
    H_E[pdg] = new TH1D(Form("H_E_%s", PDGtoName(pdg).c_str()), Form("Energy distribution for %s; Energy (keV); Counts", PDGtoName(pdg).c_str()), 100000, 0, 10000);
    H_Ex[pdg] = new TH1D(Form("H_Ex_%s", PDGtoName(pdg).c_str()), Form("Excitation energy distribution for %s; Excitation Energy (keV); Counts", PDGtoName(pdg).c_str()), 10000, 0, 10000);
    H_px[pdg] = new TH1D(Form("H_px_%s", PDGtoName(pdg).c_str()), Form("px distribution for %s; px; Counts", PDGtoName(pdg).c_str()), 1000, -1, 1);
    H_py[pdg] = new TH1D(Form("H_py_%s", PDGtoName(pdg).c_str()), Form("py distribution for %s; py; Counts", PDGtoName(pdg).c_str()), 1000, -1, 1);
    H_pz[pdg] = new TH1D(Form("H_pz_%s", PDGtoName(pdg).c_str()), Form("pz distribution for %s; pz; Counts", PDGtoName(pdg).c_str()), 1000, -1, 1);
}

int Reader(string filename = "")
{

    // Getting filename
    // string filename = "../build/output.root";
    // if (argc > 1)
    //     filename = argv[1];
    // else
    //     Error("No input file provided. Please provide the path to the ROOT file as an argument.");
    
    // Openning file
    TFile *f = new TFile(filename.c_str(), "READ");
    if (!f->IsOpen())
        Error("Could not open file " + filename);
    
    // Setting up TTreeReader
    TTree *t = (TTree *)f->Get("ParticleTree");
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

    // Creating new file to save histograms 
    string filename_out = filename.substr(0, filename.find_last_of('.')) + "_read.root";
    TFile *fout = new TFile(filename_out.c_str(), "RECREATE");

    // Setup a correlations //
    vector<Correlation> Correlations;
    // Correlation gammagamma_60Co(22, 22, 1173, 1332);
    // Correlations.push_back(gammagamma_60Co);

    InitHistogram(-11); // electron
    InitHistogram(12); // electron
    InitHistogram(-12); // electron
    InitHistogram(11);  // positron
    InitHistogram(22);  // gamma
    InitHistogram(2212); // proton
    InitHistogram(1000170320);
    InitHistogram(1000180320);
    InitHistogram(1000160310);

    // Reading TTree
    int Verbosity = 0;
    while (Reader->Next())
    {
        if (Verbosity == 0)
            ProgressBar(Reader->GetCurrentEntry(), Entries-1, start, "", step);
        
        if (Verbosity > 0)
            Info("#### Event: " + std::to_string(Reader->GetCurrentEntry()) + " #### (N = " + std::to_string(Code->GetSize()) + ")", 0);

        // Lopping on particle in the event
        for (int i = 0; i < Code->GetSize(); i++)
        {
            if (Verbosity > 0)
                Info(Form("%s \t E = %.1f keV \t Ex = %.1f keV \t px = %.2f \t py = %.2f \t pz = %.2f", PDGtoName((*Code)[i]).c_str(), (*Energy)[i], (*p)[i], (*Px)[i], (*Py)[i], (*Pz)[i]), 1);

            // InitHistogram((*Code)[i]);
            
            // Filling histograms
            // H_time[(*Code)[i]]->Fill((*Time)[i]);
            H_E[(*Code)[i]]->Fill((*Energy)[i]);
            H_Ex[(*Code)[i]]->Fill((*Excitation_energy)[i]);
            // H_px[(*Code)[i]]->Fill((*Px)[i]);
            // H_py[(*Code)[i]]->Fill((*Py)[i]);
            // H_pz[(*Code)[i]]->Fill((*Pz)[i]);

            // for (Correlation &correlation : Correlations)
            // {
            //     if (((*Code)[i] == correlation.PDG1 && abs((*Energy)[i] - correlation.E1) < 1))
            //     {
            //         for (int j = i + 1; j < Code->GetSize(); j++)
            //         {
            //             if (((*Code)[j] == correlation.PDG2 && abs((*Energy)[j] - correlation.E2) < 1))
            //             {
            //                 double costheta = ((*Px)[i] * (*Px)[j] + (*Py)[i] * (*Py)[j] + (*Pz)[i] * (*Pz)[j]);
            //                 correlation.Fill(costheta);
            //             }
            //         }
            //     }
            // }
        }
        
        if (Verbosity > 0)
            cout << endl;
    }

    // Writting Histograms
    fout->cd();
    // Particles
    for (map<int, TH1D *>::iterator it = H_time.begin(); it != H_time.end(); ++it)
    {
        fout->mkdir(PDGtoName(it->first).c_str());
        fout->cd(PDGtoName(it->first).c_str());
        H_time[it->first]->Write();
        H_E[it->first]->Write();
        H_Ex[it->first]->Write();
        H_px[it->first]->Write();
        H_py[it->first]->Write();
        H_pz[it->first]->Write();
        fout->cd();
    }
    // Correlations 
    for (Correlation &correlation : Correlations)
    {
        if (fout->GetDirectory(Form("%s_%s", PDGtoName(correlation.PDG1).c_str(), PDGtoName(correlation.PDG2).c_str())) == nullptr)
            fout->mkdir(Form("%s_%s", PDGtoName(correlation.PDG1).c_str(), PDGtoName(correlation.PDG2).c_str()));
        fout->cd(Form("%s_%s", PDGtoName(correlation.PDG1).c_str(), PDGtoName(correlation.PDG2).c_str()));
        if (correlation.H_Angle)
            correlation.H_Angle->Write();
        fout->cd();
    }
    fout->Close();

    // General Parsing 
    Message("", "## General Parsing ##", 0, "cyan");
    Message("", "Total number of events: " + std::to_string(Entries), 1, "blue");
    Message("", "Particle count:", 0, "blue");
    for (map<int, int>::iterator it = ParticleCounter.begin(); it != ParticleCounter.end(); ++it)
    {
        Message("", Form("%s : %d", PDGtoName(it->first).c_str(), it->second), 1, "blue");
    }

    return 0;
}