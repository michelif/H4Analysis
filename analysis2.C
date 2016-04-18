#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis2() {
    Int_t nbins = 800;
    Int_t j;
    char name[20];
    char title[100];
    TStopwatch t;
    TFile f1("/home/marko/H4Analysis/ntuples/analysis_3898.root");
    TFile f2("/home/marko/H4Analysis/ntuples/analysis_3902.root");
    TFile f3("/home/marko/H4Analysis/ntuples/analysis_3905.root");
    TTree* h4_1 = (TTree*) f1.Get("h4");
    TTree* h4_2 = (TTree*) f2.Get("h4");
    TTree* h4_3 = (TTree*) f3.Get("h4");
    TFile outputfile("AllPedestalEventFFTs.root", "recreate");
    Int_t nentries1 = h4_1->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t nentries2 = h4_2->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t nentries3 = h4_3->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t entriestotal = nentries1 + nentries2 + nentries3;
    cout << entriestotal << endl;
    Int_t nspill = 15; //how many spills will be analyzed
    TH1F *HistoEvent[entriestotal];
    for (Int_t z=0;z<entriestotal;z++) {
        sprintf(name,"HistoEvent%d",z);
        sprintf(title,"Event%d Histo", z);
        HistoEvent[z] = new TH1F(name,title,nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEvent[entriestotal];
    for (Int_t z=0;z<entriestotal;z++) {
        sprintf(name,"NewHistoEvent%d",z);
        sprintf(title,"Event%d Histo", z);
        NewHistoEvent[z] = new TH1F(name,title,nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEventFFT[entriestotal];
    for (Int_t z=0;z<entriestotal;z++) {
        sprintf(name,"NewHistoEventFFT%d",z);
        sprintf(title,"Event%d Histo", z);
        NewHistoEventFFT[z] = new TH1F(name,title,nbins, 0, 5);
    }
    TH1F* NormNoiseFFT = new TH1F ("NormNoiseFFT", "Normalized Noise FFT", nbins, 0, 5);
    TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", nbins, -0.1, 159.9, 1000, -15, 15); //nanoseconds
    Int_t count = 0;
    t.Start();
    for (Int_t spill=0;spill<nspill;spill++) {
        cout << "Run 3898 Spill " << spill << endl;
        h4_1->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_1->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_1->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_1->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_1->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            plot = "WF_val:WF_time>>TempHisto";
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_1->Draw(plot, cut, "goff");
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            for (Int_t i=0; i<nbins; i++) {
                for (Int_t k=0; k<1000; k++) {
                    if (TempHisto->GetBinContent(i+1, k) != 0) {
                        HistoEvent[count]->SetBinContent(i+1,k*0.03-15);
                    }
                }
            }
            mean = TempHisto->GetMean(2);
            rms = TempHisto->GetRMS(2);
            for (Int_t q=0;q<nbins;q++) {
                NewHistoEvent[count]->SetBinContent(q+1, HistoEvent[count]->GetBinContent(q+1)-mean);
            }
            NewHistoEvent[count]->Scale(1/rms);
            NewHistoEvent[count]->FFT(NewHistoEventFFT[count], "MAG");
            NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[count]);
            NewHistoEventFFT[count]->Write();
            cout << "Event " << count+1 << " out of " << entriestotal << endl;
            count += 1;
        }
    }
    for (Int_t spill=0;spill<nspill;spill++) {
        cout << "Run 3902 Spill " << spill << endl;
        h4_2->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_2->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_2->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_2->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_2->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            plot = "WF_val:WF_time>>TempHisto";
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_2->Draw(plot, cut, "goff");
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            for (Int_t i=0; i<nbins; i++) {
                for (Int_t k=0; k<1000; k++) {
                    if (TempHisto->GetBinContent(i+1, k) != 0) {
                        HistoEvent[count]->SetBinContent(i+1,k*0.03-15);
                    }
                }
            }
            mean = TempHisto->GetMean(2);
            rms = TempHisto->GetRMS(2);
            for (Int_t q=0;q<nbins;q++) {
                NewHistoEvent[count]->SetBinContent(q+1, HistoEvent[count]->GetBinContent(q+1)-mean);
            }
            NewHistoEvent[count]->Scale(1/rms);
            NewHistoEvent[count]->FFT(NewHistoEventFFT[count], "MAG");
            NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[count]);
            NewHistoEventFFT[count]->Write();
            cout << "Event " << count+1 << " out of " << entriestotal << endl;
            count += 1;
        }
    }
    for (Int_t spill=1;spill<2;spill++) {
        cout << "Run 3905 Spill " << spill << endl;
        h4_3->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_3->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_3->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_3->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_3->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            plot = "WF_val:WF_time>>TempHisto";
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_3->Draw(plot, cut, "goff");
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            for (Int_t i=0; i<nbins; i++) {
                for (Int_t k=0; k<1000; k++) {
                    if (TempHisto->GetBinContent(i+1, k) != 0) {
                        HistoEvent[count]->SetBinContent(i+1,k*0.03-15);
                    }
                }
            }
            mean = TempHisto->GetMean(2);
            rms = TempHisto->GetRMS(2);
            for (Int_t q=0;q<nbins;q++) {
                NewHistoEvent[count]->SetBinContent(q+1, HistoEvent[count]->GetBinContent(q+1)-mean);
            }
            NewHistoEvent[count]->Scale(1/rms);
            NewHistoEvent[count]->FFT(NewHistoEventFFT[count], "MAG");
            NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[count]);
            NewHistoEventFFT[count]->Write();
            cout << "Event " << count+1 << " out of " << entriestotal << endl;
            count += 1;
        }
    }
    TFile out("AllNormalizedNoiseFFT.root", "recreate");
    NormNoiseFFT->Scale(1./count);
    NormNoiseFFT->Write();
    t.Stop();
    t.Print();
}