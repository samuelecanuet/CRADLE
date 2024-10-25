void macro()

{

    TFile *f = new TFile("build/test.root", "READ");

    TTree *t = (TTree *)f->Get("ParticleTree");
    TTreeReader *Reader = new TTreeReader(t);

    TTreeReaderValue<int> *code = new TTreeReaderValue<int>(*Reader, "code");
    TTreeReaderValue<double> *energy = new TTreeReaderValue<double>(*Reader, "energy");

    TH1D *H = new TH1D("h", "h", 10000, 0, 10000);
    while (Reader->Next())
    {
        if (**code == 2212)
        {
            H->Fill(**energy);
        }
    }


    // f = new TFile("../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/32Ar_CS0_CSP0_CV1_CVP1_1.root", "READ");
    // t = (TTree *)f->Get("ParticleTree");
    // Reader = new TTreeReader(t);
    // Reader->Restart();

    // TH1D *H_1 = new TH1D("h1", "h1", 10000, 0, 10000);
    // while (Reader->Next())
    // {
    //     if (**code == 2212)
    //     {
    //         H_1->Fill(**energy);
    //     }
    // }

    TFile *ff = new TFile("testtt.root", "RECREATE");
    ff->cd();
    TCanvas *c = new TCanvas("c", "c", 600, 800);
    c->cd();
    H->Scale(1./H->Integral());
    // H_1->Scale(1./H_1->Integral());
    H->Draw("HIST");
    // H_1->SetLineColor(kRed);
    // H_1->Draw("same");
    c->Write();
    ff->Close();
    f->Close();
}