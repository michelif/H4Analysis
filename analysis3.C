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
    TH1F *HistoEvent[2214];
    for (Int_t z=0;z<2214;z++) {
        sprintf(name,"HistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        HistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *PhaseHistoEvent[2214];
    for (Int_t z=0;z<2214;z++) {
        sprintf(name,"HistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        PhaseHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 799.9);
    }
    TH1F *NewHistoEvent[2214];
    for (Int_t z=0;z<2214;z++) {
        sprintf(name,"NewHistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEventFFT[2214];
    for (Int_t z=0;z<2214;z++) {
        sprintf(name,"NewHistoEventFFT%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEventFFT[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT[2214];
    for (Int_t z=0;z<2214;z++) {
        sprintf(name,"SignalEventFFT%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEventFFT[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    Double_t mean;
    Double_t rms;
    Double_t count = 0;
    TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
    TTree* h4 = (TTree*) f.Get("h4");
    TFile f1("NormalizedNoiseFFT.root", "read");
    TH1F* NormNoiseFFT = (TH1F*) f1.Get("NormNoiseFFT");
    TFile outputfile("NormalizedSignalNoise.root", "recreate");
    TString plot;
    TString cut;
    TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
    TH1F* NormSigFFT = new TH1F ("NormSigFFT", "Normalized SignalNoise FFT", nbins, 0, 5);
    TH1F* Phase = new TH1F ("Phase", "Phase", nbins, -0.1, 799.9);
    TStopwatch t;
    t.Start(); //1 hour runtime
    for (j=10;j<20;j++) {
        plot = "WF_val:WF_time>>TempHisto";
        cut = "WF_ch==2 && event==";
        cut += j;
        cout << plot << " " << cut << endl;
        h4->Draw(plot, cut, "goff");
        if (TempHisto->GetMaximum() == 0) {
            delete HistoEvent[j+1];
            delete NewHistoEvent[j+1];
            delete NewHistoEventFFT[j+1];
            delete SignalEventFFT[j+1];
            continue;
        }
        for (Int_t i=0; i<nbins; i++) {
            for (Int_t k=0; k<1000; k++) {
                if (TempHisto->GetBinContent(i+1, k) != 0) {
                    HistoEvent[j+1]->SetBinContent(i+1,k*0.92-120);
                }
            }
        }
        TempHisto->GetXaxis()->SetRange(0,225);
        mean = TempHisto->GetMean(2);
        rms = TempHisto->GetRMS(2);
        TempHisto->GetXaxis()->SetRange(0,nbins);
        for (Int_t q=0;q<nbins;q++) {
            NewHistoEvent[j+1]->SetBinContent(q+1, HistoEvent[j+1]->GetBinContent(q+1)-mean);
        }
        NewHistoEvent[j+1]->FFT(NewHistoEventFFT[j+1], "MAG"); // Signal+Noise in frequency domain
        //NormNoiseFFT->Scale(rms); //Normalize Noise in frequency domain
        for (Int_t q=0;q<nbins;q++) {
            SignalEventFFT[j+1]->SetBinContent(q+1, (NewHistoEventFFT[j+1]->GetBinContent(q+1)-NormNoiseFFT->GetBinContent(q+1))); // S+N - N = S in frequency domain
        }
        NormSigFFT->Add(NormSigFFT, SignalEventFFT[j+1]); // adding all the signal events in the frequency domain
        HistoEvent[j+1]->FFT(PhaseHistoEvent[j+1], "PH");
        Phase->Add(Phase, PhaseHistoEvent[j+1]);
        NewHistoEvent[j+1]->Write();
        TempHisto->Write();
        NewHistoEventFFT[j+1]->Write();
        SignalEventFFT[j+1]->Write();
        cout << "Event " << j << ", Mean = " << mean << endl;
        count += 1;
    }
    NormSigFFT->Scale(1/count); //scaling by the # of events
    Phase->Scale(1/count);
    NormNoiseFFT->Write();
    Double_t *re_full = new Double_t[nbins];
    Double_t *im_full = new Double_t[nbins];
    for (Int_t n=0; n<nbins; n++) {
        (re_full)[n]=(NormSigFFT->GetBinContent(n+1)*cos(Phase->GetBinContent(n+1))); //magnitude * cos of phase
        (im_full)[n]=(NormSigFFT->GetBinContent(n+1)*sin(Phase->GetBinContent(n+1))); //magnitude * sin of phase
    }
    TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
    invFFT->SetPointsComplex(re_full, im_full);
    invFFT->Transform();
    TH1 *Signal = 0;
    Signal = TH1::TransformHisto(invFFT,Signal,"Re");
    Signal->SetTitle("Recovered Signal 'S'");
    TH1F* BetterSignal = new TH1F ("BetterSignal", "Recovered Signal", nbins, -0.1, 159.9);
    for (Int_t p=0; p<nbins; p++) {
        BetterSignal->SetBinContent(p+1, Signal->GetBinContent(p+1)/nbins);
    }
    BetterSignal->Write();
    t.Stop();
    t.Print();
}