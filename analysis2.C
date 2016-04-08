#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis2() {
    Int_t nbins = 800;
    Int_t j;
    char name[20];
    char title[100];
    TH1F *HistoEvent[2219];
    for (Int_t z=0;z<2219;z++) {
        sprintf(name,"HistoEvent%d",z-1);
        sprintf(title,"Event%d Histo", z-1);
        HistoEvent[z] = new TH1F(name,title,nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEvent[2219];
    for (Int_t z=0;z<2219;z++) {
        sprintf(name,"NewHistoEvent%d",z-1);
        sprintf(title,"Event%d Histo", z-1);
        NewHistoEvent[z] = new TH1F(name,title,nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEventFFT[2219];
    for (Int_t z=0;z<2219;z++) {
        sprintf(name,"NewHistoEventFFT%d",z-1);
        sprintf(title,"Event%d Histo", z-1);
        NewHistoEventFFT[z] = new TH1F(name,title,nbins, 0, 5);
    }
    Double_t mean;
    Double_t rms;
    Double_t count = 0;
    TFile f("/home/marko/H4Analysis/ntuples/analysis_3905.root");
    TTree* h4 = (TTree*) f.Get("h4");
    TFile outputfile("NormalizedNoiseFFT.root", "recreate");
    TString plot;
    TString cut;
    TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", nbins, -0.1, 159.9, 1000, -15, 15); //nanoseconds
    TH1F* NormNoiseFFT = new TH1F ("NormNoiseFFT", "Normalized Noise FFT", nbins, 0, 5);
    TStopwatch t;
    t.Start(); //1 hour runtime
    for (j=0;j<2218;j++) {
                plot = "WF_val:WF_time>>TempHisto";
                cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && event==";
                cut += j;
                h4->Draw(plot, cut, "goff");
                if (TempHisto->GetMaximum() == 0) {
                    delete HistoEvent[j+1];
                    delete NewHistoEvent[j+1];
                    delete NewHistoEventFFT[j+1];
                    continue;
                }
                for (Int_t i=0; i<nbins; i++) {
                    for (Int_t k=0; k<1000; k++) {
                        if (TempHisto->GetBinContent(i+1, k) != 0) {
                            HistoEvent[j+1]->SetBinContent(i+1,k*0.03-15);
                        }
                    }
                }
                mean = TempHisto->GetMean(2);
                rms = TempHisto->GetRMS(2);
                for (Int_t q=0;q<nbins;q++) {
                    NewHistoEvent[j+1]->SetBinContent(q+1, HistoEvent[j+1]->GetBinContent(q+1)-mean);
                }
                NewHistoEvent[j+1]->Scale(1/rms);
                NewHistoEvent[j+1]->FFT(NewHistoEventFFT[j+1], "MAG");
                NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[j+1]);
                //NewHistoEvent[j+1]->Write();
                //NewHistoEventFFT[j+1]->Write();
                cout << "Event " << j << ", Mean = " << mean << ", RMS = " << rms << endl;
                count += 1;
//           }
//        }
    }
    NormNoiseFFT->Scale(1/count);
    NormNoiseFFT->Write();
    t.Stop();
    t.Print();
}