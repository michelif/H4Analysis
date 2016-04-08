#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis3() {
    Int_t nbins = 800;
    Int_t j;
    char name[20];
    char title[100];
    TH1F *HistoEvent[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"HistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        HistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *PhaseHistoEvent[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"PhaseHistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        PhaseHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 799.9);
    }
    TH1F *NewHistoEvent[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"NewHistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEventFFT[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"NewHistoEventFFT%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEventFFT[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"SignalEventFFT%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1 *SignalEvent[1811];
    for (Int_t z=0;z<1811;z++) {
        sprintf(name,"SignalEvent%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    Float_t mean;
    Float_t rms;
    Float_t rmsprime;
    Int_t count = 0;
    TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
    TTree* h4 = (TTree*) f.Get("h4");
    TFile f1("NormalizedNoiseFFT.root", "read");
    TH1F* NormNoiseFFT = (TH1F*) f1.Get("NormNoiseFFT");
    TFile outputfile("NormalizedSignalNoise.root", "recreate");
    TTree *MyTree = new TTree ("MyTree", "MyTree");
    MyTree->Branch("rms", &rms, "rms/F");
    MyTree->Branch("rmsprime", &rmsprime, "rmsprime/F");
    MyTree->Branch("count", &count, "count/I");
    TString plot;
    TString cut;
    TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
    TH1F* NormSigFFT = new TH1F ("NormSigFFT", "Normalized SignalNoise FFT", nbins, 0, 5);
    TH1F* Phase = new TH1F ("Phase", "Phase", nbins, -0.1, 799.9);
    TStopwatch t;
    t.Start(); //1 hour runtime
    for (j=1;j<1812;j++) {
        TH1 *Throwaway = 0;
        TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        std::vector<Float_t> y;
        Float_t mu = 0;
        plot = "WF_val:WF_time>>TempHisto";
        cut = "WF_ch==2 && spill==1 && event==";
        cut += j;
        h4->Draw(plot, cut, "goff");
        if (TempHisto->GetMaximum() == 0) {
            delete HistoEvent[j];
            delete NewHistoEvent[j];
            delete NewHistoEventFFT[j];
            delete SignalEventFFT[j];
            delete PhaseHistoEvent[j];
            delete SignalEvent[j];
            continue;
        }
        for (Int_t i=0; i<nbins; i++) {
            for (Int_t k=0; k<1000; k++) {
                if (TempHisto->GetBinContent(i+1, k) != 0) {
                    HistoEvent[j]->SetBinContent(i+1,k*0.92-120);
                }
            }
        }
        TempHisto->GetXaxis()->SetRange(0,225);
        mean = TempHisto->GetMean(2);
        rms = TempHisto->GetRMS(2);
        TempHisto->GetXaxis()->SetRange(0,nbins);
        for (Int_t q=0;q<nbins;q++) {
            NewHistoEvent[j]->SetBinContent(q+1, HistoEvent[j]->GetBinContent(q+1)-mean); //NewHistoEvent now becomes the "original"
        }
        NewHistoEvent[j]->FFT(NewHistoEventFFT[j], "MAG"); // Signal+Noise in frequency domain
        NormNoiseFFT->Scale(rms); //Normalize Noise in frequency domain
        for (Int_t q=0;q<nbins;q++) {
            SignalEventFFT[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT->GetBinContent(q+1))); // S+N - N = S in frequency domain
        }
        NormNoiseFFT->Scale(1/rms);
        HistoEvent[j]->FFT(PhaseHistoEvent[j], "PH");
        Double_t *re_full = new Double_t[nbins];
        Double_t *im_full = new Double_t[nbins];
        for (Int_t n=0; n<nbins; n++) {
            (re_full)[n]=(SignalEventFFT[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1))); //magnitude * cos of phase
            (im_full)[n]=(SignalEventFFT[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1))); //magnitude * sin of phase
        }
        invFFT->SetPointsComplex(re_full, im_full);
        invFFT->Transform();
        Throwaway = TH1::TransformHisto(invFFT,Throwaway,"Re");
        for (Int_t p=0; p<nbins; p++) {
            SignalEvent[j]->SetBinContent(p+1, Throwaway->GetBinContent(p+1)/nbins);
        }
        NewHistoEvent[j]->Write();
        SignalEvent[j]->Write();
        for (Int_t q=0; q<225; q++) {
            y.push_back(SignalEvent[j]->GetBinContent(q));
            mu += y.at(q);
        }
        mu = mu/225;
        rmsprime = 0;
        for (Int_t q=0; q<225; q++) {
            rmsprime += (y[q]-mu)*(y[q]-mu);
        }
        rmsprime = sqrt(rmsprime/225);
        y.clear();
        delete HistoEvent[j];
        delete PhaseHistoEvent[j];
        delete NewHistoEventFFT[j];
        delete SignalEventFFT[j];
        cout << "Event " << j << ", Mean = " << mean << ", RMS = " << rms << ", RMS' = " << rmsprime << endl;
        count += 1;
        delete invFFT;
        delete Throwaway;
        MyTree->Fill();
    }
    MyTree->Write();
    //NormSigFFT->Scale(1/count); //scaling by the # of events
    //Phase->Scale(1/count);
    //NormNoiseFFT->Write();
    //TH1F* BetterSignal = new TH1F ("BetterSignal", "Recovered Signal", nbins, -0.1, 159.9);
    //BetterSignal->Write();
    t.Stop();
    t.Print();
}