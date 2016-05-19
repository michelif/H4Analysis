#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TH2F.h"
#include <iostream>
#include <fstream>
using namespace std;

TH1F* transform2Dto1D (TH2F *h2) {
    Int_t ybins = h2->GetNbinsY();
    Int_t xbins = h2->GetNbinsX();
    Float_t ymax = h2->GetYaxis()->GetBinUpEdge(h2->GetNbinsY());
    Float_t ymin = h2->GetYaxis()->GetBinLowEdge(1);
    Float_t xmax = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());
    Float_t xmin = h2->GetXaxis()->GetBinLowEdge(1);
    Float_t yrange = ymax - ymin;
    TString name = h2->GetName();
    name += "_Transformed";
    TH1F *h1 = new TH1F (name, "1D Histogram", xbins, xmin, xmax);
    for (Int_t xbin=0; xbin<xbins; xbin++) {
        for (Int_t ybin=0; ybin<ybins; ybin++) {
            if (h2->GetBinContent(xbin+1, ybin) != 0) {
                h1->SetBinContent(xbin+1,ybin*(yrange/ybins)+ymin);
            }
        }
    }
    return h1;
}

int main() {
    Int_t nbins = 800, j;
    char name[20], title[100];
    TStopwatch t;
    TFile f1("/home/marko/Desktop/H4Analysis/ntuples/analysis_3898.root"); //Run 3898 ntuple
    TFile f2("/home/marko/Desktop/H4Analysis/ntuples/analysis_3902.root"); //Run 3902 ntuple
    TFile f3("/home/marko/Desktop/H4Analysis/ntuples/analysis_3905.root"); //Run 3905 ntuple
    TTree* h4_3898 = (TTree*) f1.Get("h4");
    TTree* h4_3902 = (TTree*) f2.Get("h4");
    TTree* h4_3905 = (TTree*) f3.Get("h4");
    TFile outputfile("AllPedestalEventFFTs.root", "recreate"); //individual event noise spectra
    Int_t nentries1 = h4_3898->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t nentries2 = h4_3902->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t nentries3 = h4_3905->GetEntries("WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160");
    Int_t entriestotal = nentries1 + nentries2 + nentries3;
    Int_t nspill = 15; //how many spills will be analyzed (not all spills necessarily have data)
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
    Int_t count = 0;
    t.Start();
    for (Int_t spill=0;spill<nspill;spill++) {
        cout << "Run 3898 Spill " << spill << endl;
        h4_3898->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_3898->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_3898->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_3898->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_3898->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            TString histoname = "Run3898TempHisto_";
            histoname += spill;
            histoname += "_";
            histoname += j;
            TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -15, 15);
            plot = "WF_val:WF_time>>";
            plot += histoname;
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_3898->Draw(plot, cut, "goff");
            TempHisto = (TH2F*) gDirectory->Get(histoname);
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            HistoEvent[count] = transform2Dto1D(TempHisto);
            mean = TempHisto->GetMean(2);
            rms = TempHisto->GetRMS(2);
            for (Int_t q=0;q<nbins;q++) {
                NewHistoEvent[count]->SetBinContent(q+1, HistoEvent[count]->GetBinContent(q+1)-mean); //centering the pedestal at <y> = 0
            }
            NewHistoEvent[count]->Scale(1/rms); //dividing by RMS of pedestal event, later undone when the filter is applied to the wave pulses
            NewHistoEvent[count]->FFT(NewHistoEventFFT[count], "MAG"); //noise power spectrum (frequency domain)
            NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[count]);
            NewHistoEventFFT[count]->Write();
            cout << "Event " << count+1 << " out of " << entriestotal << endl;
            count += 1;
            delete TempHisto;
        }
    }
    for (Int_t spill=0;spill<nspill;spill++) {
        cout << "Run 3902 Spill " << spill << endl;
        h4_3902->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_3902->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_3902->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_3902->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_3902->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            TString histoname = "Run3902TempHisto_";
            histoname += spill;
            histoname += "_";
            histoname += j;
            TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -15, 15);
            plot = "WF_val:WF_time>>";
            plot += histoname;
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_3902->Draw(plot, cut, "goff");
            TempHisto = (TH2F*) gDirectory->Get(histoname);
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            HistoEvent[count] = transform2Dto1D(TempHisto);
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
            delete TempHisto;
        }
    }
    for (Int_t spill=1;spill<2;spill++) {
        cout << "Run 3905 Spill " << spill << endl;
        h4_3905->SetEntryList(0);
        TString listcut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
        listcut += spill;
        TString spillcut = "spill==";
        spillcut += spill;
        h4_3905->Draw(">>myList", listcut, "entrylist");
        TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
        h4_3905->SetEntryList(myList);
        Int_t entriesperspill = myList->GetN();
        h4_3905->Draw("event", spillcut, "goff");
        Double_t *vTemp = h4_3905->GetV1();
        Double_t *vEvent = new Double_t[entriesperspill];
        for (int iEntry = 0; iEntry<entriesperspill; iEntry++){
            vEvent[iEntry] = vTemp[iEntry];
        }
        Double_t mean;
        Double_t rms;
        TString plot;
        TString cut;
        for (j=0;j<entriesperspill;j++) {
            TString histoname = "Run3905TempHisto_";
            histoname += spill;
            histoname += "_";
            histoname += j;
            TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -15, 15);
            plot = "WF_val:WF_time>>";
            plot += histoname;
            cut = "WF_ch==APD1 && amp_max[APD3]<25 && b_rms[APD3]<5. && charge_tot[APD3]<20000 &&  amp_max[APD5]<25 && b_rms[APD5]<5. &&  amp_max[APD6]<25 && b_rms[APD6]<5. &&  amp_max[APD4]<25 && b_rms[APD4]<5. && amp_max[SiPM1]<20 && amp_max[SiPM2]<20 && amp_max[APD1]<40 && amp_max[APD2]<40 && b_rms[APD1]<5. && b_rms[APD2]<5. && WF_time<160 && spill==";
            cut += spill;
            cut += " && event==";
            cut += vEvent[j];
            h4_3905->Draw(plot, cut, "goff");
            TempHisto = (TH2F*) gDirectory->Get(histoname);
            if (TempHisto->GetMaximum() == 0) {
                continue;
            }
            HistoEvent[count] = transform2Dto1D(TempHisto);
            mean = TempHisto->GetMean(2);
            rms = TempHisto->GetRMS(2);
            for (Int_t q=0;q<nbins;q++) {
                NewHistoEvent[count]->SetBinContent(q+1, HistoEvent[count]->GetBinContent(q+1)-mean);
            }
            NewHistoEvent[count]->Scale(1/rms);
            NewHistoEvent[count]->FFT(NewHistoEventFFT[count], "MAG");
            NormNoiseFFT->Add(NormNoiseFFT, NewHistoEventFFT[count]);
            HistoEvent[count]->Write();
            NewHistoEventFFT[count]->Write();
            cout << "Event " << count+1 << " out of " << entriestotal << endl;
            count += 1;
            delete TempHisto;
        }
    }
    TFile out("AllNormalizedNoiseFFT.root", "recreate");
    NormNoiseFFT->Scale(1./count);
    NormNoiseFFT->Write();
    t.Stop();
    t.Print();
}