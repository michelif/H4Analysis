#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TVirtualFFT.h"
#include <iostream>
#include <fstream>
using namespace std;

Float_t pedestalrms1D (TH1 *h1) {
    std::vector<Float_t> y;
    for (Int_t q=0;q<225;q++) {
            y.push_back(h1->GetBinContent(q+1));
    }
    Float_t newrms = TMath::RMS(y.begin(), y.end());
    y.clear();
    return newrms;
}

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

TH1F* inverseFFT (TH1F *h1, TH1F *phase) {
    Int_t nbins = h1->GetNbinsX();
    TString name = h1->GetName();
    name += "_RE";
    TString name2 = h1->GetName();
    name2.Remove(11,3);
    name2 += "_inv";
    Double_t *re_full = new Double_t[nbins];
    Double_t *im_full = new Double_t[nbins];
    TH1 *Throwaway = 0;
    TH1F *invh1 = new TH1F(name2, "Signal Event", nbins, -0.1, 159.9);
    TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
     for (Int_t n=0; n<nbins; n++) {
        (re_full)[n]=(h1->GetBinContent(n+1)*cos(phase->GetBinContent(n+1)));
        (im_full)[n]=(h1->GetBinContent(n+1)*sin(phase->GetBinContent(n+1)));
    }
    invFFT->SetPointsComplex(re_full, im_full);
    invFFT->Transform();
    Throwaway = TH1::TransformHisto(invFFT, Throwaway, "Re");
    Throwaway->SetName(name);
    for (Int_t p=0; p<nbins; p++) {
        invh1->SetBinContent(p+1, Throwaway->GetBinContent(p+1)/nbins);
    }
    return invh1;
}

int main() {
    Int_t nspill = 1;
    Int_t nbins = 800;
    Int_t j;
    TFile f("/home/marko/Desktop/H4Analysis/ntuples/analysis_4443.root");
    TTree* h4 = (TTree*) f.Get("h4");
    h4->SetEntryList(0);
    TString listcut = "WF_ch==2 && spill==1";
    h4->Draw(">>myList", listcut, "entrylist");
    TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    h4->SetEntryList(myList);
    Int_t nevents = myList->GetN();
    h4->Draw("event", "spill==1", "goff");
    Double_t *vTemp = h4->GetV1();
    Double_t *vEvent = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vEvent[iEntry] = vTemp[iEntry];
    }
    char name[20];
    char title[100];
    TH1F *HistoEvent[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"HistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        HistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *PhaseHistoEvent[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"PhaseHistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        PhaseHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 799.9);
    }
    TH1F *NewHistoEvent[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"NewHistoEvent%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEvent[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1F *NewHistoEventFFT[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"NewHistoEventFFT%d",z+1);
        sprintf(title,"Event%d Histo", z+1);
        NewHistoEventFFT[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT100Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT100Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT100Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT250Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT250Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT500Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT500Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT750Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT750Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT750Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1000Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1000Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1000Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1250Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1250Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1500Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1500Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1800Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1800Events%d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1800Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1 *SignalEvent100Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent100Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent100Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent250Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent250Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent500Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent500Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent750Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent750Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent750Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1000Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1000Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1000Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1250Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1250Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1500Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1500Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1800Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1800Events%d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1800Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    Float_t mean;
    Float_t rms;
    Float_t rmsprime100;
    Float_t rmsprime250;
    Float_t rmsprime500;
    Float_t rmsprime750;
    Float_t rmsprime1000;
    Float_t rmsprime1250;
    Float_t rmsprime1500;
    Float_t rmsprime1800;
    Int_t count = 0;
    TFile f1("NormNoiseFFTs.root", "read");
    TFile f2("AllNormalizedNoiseFFT.root", "read");
    TH1F* NormNoiseFFT100Events = (TH1F*) f1.Get("NormNoiseFFT100Events");
    TH1F* NormNoiseFFT250Events = (TH1F*) f1.Get("NormNoiseFFT250Events");
    TH1F* NormNoiseFFT500Events = (TH1F*) f1.Get("NormNoiseFFT500Events");
    TH1F* NormNoiseFFT750Events = (TH1F*) f1.Get("NormNoiseFFT750Events");
    TH1F* NormNoiseFFT1000Events = (TH1F*) f1.Get("NormNoiseFFT1000Events");
    TH1F* NormNoiseFFT1250Events = (TH1F*) f1.Get("NormNoiseFFT1250Events");
    TH1F* NormNoiseFFT1500Events = (TH1F*) f1.Get("NormNoiseFFT1500Events");
    TH1F* NormNoiseFFT1800Events = (TH1F*) f2.Get("NormNoiseFFT");
    TFile out1("NormalizedSignalNoise100Events.root", "recreate");
    TFile out2("NormalizedSignalNoise250Events.root", "recreate");
    TFile out3("NormalizedSignalNoise500Events.root", "recreate");
    TFile out4("NormalizedSignalNoise750Events.root", "recreate");
    TFile out5("NormalizedSignalNoise1000Events.root", "recreate");
    TFile out6("NormalizedSignalNoise1250Events.root", "recreate");
    TFile out7("NormalizedSignalNoise1500Events.root", "recreate");
    TFile out8("NormalizedSignalNoise1800Events.root", "recreate");
    TFile out9("Pre&PostRMS.root", "recreate");
    TTree *MyTree = new TTree ("MyTree", "MyTree");
    MyTree->Branch("rms", &rms, "rms/F");
    MyTree->Branch("rmsprime100", &rmsprime100, "rmsprime100/F");
    MyTree->Branch("rmsprime250", &rmsprime250, "rmsprime250/F");
    MyTree->Branch("rmsprime500", &rmsprime500, "rmsprime500/F");
    MyTree->Branch("rmsprime750", &rmsprime750, "rmsprime750/F");
    MyTree->Branch("rmsprime1000", &rmsprime1000, "rmsprime1000/F");
    MyTree->Branch("rmsprime1250", &rmsprime1250, "rmsprime1250/F");
    MyTree->Branch("rmsprime1500", &rmsprime1500, "rmsprime1500/F");
    MyTree->Branch("rmsprime1800", &rmsprime1800, "rmsprime1800/F");
    MyTree->Branch("count", &count, "count/I");
    TString plot;
    TString cut;
    TH1F* NormSigFFT = new TH1F ("NormSigFFT", "Normalized SignalNoise FFT", nbins, 0, 5);
    TH1F* Phase = new TH1F ("Phase", "Phase", nbins, -0.1, 799.9);
    TStopwatch t;
    t.Start(); //1 hour runtime
    for (j=1;j<nevents;j++) {
    //for (j=1;j<11;j++) {
        TString histoname = "TempHisto_";
        histoname += j;
        TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
        plot = "WF_val:WF_time>>";
        plot += histoname;
        cut = "WF_ch==2 && spill==";
        cut += nspill;
        cut += " && event==";
        cut += j;
        h4->Draw(plot, cut, "goff");
        TempHisto = (TH2F*) gDirectory->Get(histoname);
        if (TempHisto->GetMaximum() == 0) {
            delete HistoEvent[j];
            delete NewHistoEvent[j];
            delete NewHistoEventFFT[j];
            delete SignalEventFFT100Events[j];
            delete SignalEventFFT250Events[j];
            delete SignalEventFFT500Events[j];
            delete SignalEventFFT750Events[j];
            delete SignalEventFFT1000Events[j];
            delete SignalEventFFT1250Events[j];
            delete SignalEventFFT1500Events[j];
            delete SignalEventFFT1800Events[j];
            delete PhaseHistoEvent[j];
            delete SignalEvent100Events[j];
            delete SignalEvent250Events[j];
            delete SignalEvent500Events[j];
            delete SignalEvent750Events[j];
            delete SignalEvent1000Events[j];
            delete SignalEvent1250Events[j];
            delete SignalEvent1500Events[j];
            delete SignalEvent1800Events[j];
            continue;
        }
        HistoEvent[j] = transform2Dto1D(TempHisto);
        TempHisto->GetXaxis()->SetRange(0,225); // sets the range to only view the pedestal
        mean = TempHisto->GetMean(2);
        rms = TempHisto->GetRMS(2);
        TempHisto->GetXaxis()->SetRange(0,nbins);
        for (Int_t q=0;q<nbins;q++) {
            NewHistoEvent[j]->SetBinContent(q+1, HistoEvent[j]->GetBinContent(q+1)-mean); //NewHistoEvent now becomes the "original"
        }
        NewHistoEvent[j]->FFT(NewHistoEventFFT[j], "MAG"); // Signal+Noise in frequency domain
        NormNoiseFFT100Events->Scale(rms); //Normalize Noise in frequency domain
        NormNoiseFFT250Events->Scale(rms);
        NormNoiseFFT500Events->Scale(rms);
        NormNoiseFFT750Events->Scale(rms);
        NormNoiseFFT1000Events->Scale(rms);
        NormNoiseFFT1250Events->Scale(rms);
        NormNoiseFFT1500Events->Scale(rms);
        NormNoiseFFT1800Events->Scale(rms);
        for (Int_t q=0;q<nbins;q++) {
            SignalEventFFT100Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT100Events->GetBinContent(q+1))); // S+N - N = S in frequency domain
            SignalEventFFT250Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT250Events->GetBinContent(q+1)));
            SignalEventFFT500Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT500Events->GetBinContent(q+1)));
            SignalEventFFT750Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT750Events->GetBinContent(q+1)));
            SignalEventFFT1000Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1000Events->GetBinContent(q+1)));
            SignalEventFFT1250Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1250Events->GetBinContent(q+1)));
            SignalEventFFT1500Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1500Events->GetBinContent(q+1)));
            SignalEventFFT1800Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1800Events->GetBinContent(q+1)));
        }
        NormNoiseFFT100Events->Scale(1/rms);
        NormNoiseFFT250Events->Scale(1/rms);
        NormNoiseFFT500Events->Scale(1/rms);
        NormNoiseFFT750Events->Scale(1/rms);
        NormNoiseFFT1000Events->Scale(1/rms);
        NormNoiseFFT1250Events->Scale(1/rms);
        NormNoiseFFT1500Events->Scale(1/rms);
        NormNoiseFFT1800Events->Scale(1/rms);
        HistoEvent[j]->FFT(PhaseHistoEvent[j], "PH");
        SignalEvent100Events[j] = inverseFFT(SignalEventFFT100Events[j], PhaseHistoEvent[j]);
        SignalEvent250Events[j] = inverseFFT(SignalEventFFT250Events[j], PhaseHistoEvent[j]);
        SignalEvent500Events[j] = inverseFFT(SignalEventFFT500Events[j], PhaseHistoEvent[j]);
        SignalEvent750Events[j] = inverseFFT(SignalEventFFT750Events[j], PhaseHistoEvent[j]);
        SignalEvent1000Events[j] = inverseFFT(SignalEventFFT1000Events[j], PhaseHistoEvent[j]);
        SignalEvent1250Events[j] = inverseFFT(SignalEventFFT1250Events[j], PhaseHistoEvent[j]);
        SignalEvent1500Events[j] = inverseFFT(SignalEventFFT1500Events[j], PhaseHistoEvent[j]);
        SignalEvent1800Events[j] = inverseFFT(SignalEventFFT1800Events[j], PhaseHistoEvent[j]);
        out1.cd();
        NewHistoEvent[j]->Write();
        SignalEvent100Events[j]->Write();
        out2.cd();
        NewHistoEvent[j]->Write();
        SignalEvent250Events[j]->Write();
        out3.cd();
        NewHistoEvent[j]->Write();
        SignalEvent500Events[j]->Write();
        out4.cd();
        NewHistoEvent[j]->Write();
        SignalEvent750Events[j]->Write();
        out5.cd();
        NewHistoEvent[j]->Write();
        SignalEvent1000Events[j]->Write();
        out6.cd();
        NewHistoEvent[j]->Write();
        SignalEvent1250Events[j]->Write();
        out7.cd();
        NewHistoEvent[j]->Write();
        SignalEvent1500Events[j]->Write();
        out8.cd();
        NewHistoEvent[j]->Write();
        SignalEvent1800Events[j]->Write();
        rmsprime100 = pedestalrms1D(SignalEvent100Events[j]);
        rmsprime250 = pedestalrms1D(SignalEvent250Events[j]);
        rmsprime500 = pedestalrms1D(SignalEvent500Events[j]);
        rmsprime750 = pedestalrms1D(SignalEvent750Events[j]);
        rmsprime1000 = pedestalrms1D(SignalEvent1000Events[j]);
        rmsprime1250 = pedestalrms1D(SignalEvent1250Events[j]);
        rmsprime1500 = pedestalrms1D(SignalEvent1500Events[j]);
        rmsprime1800 = pedestalrms1D(SignalEvent1800Events[j]);
        delete HistoEvent[j];
        delete PhaseHistoEvent[j];
        delete NewHistoEventFFT[j];
        delete SignalEventFFT100Events[j];
        delete SignalEventFFT250Events[j];
        delete SignalEventFFT500Events[j];
        delete SignalEventFFT750Events[j];
        delete SignalEventFFT1000Events[j];
        delete SignalEventFFT1250Events[j];
        delete SignalEventFFT1500Events[j];
        delete SignalEventFFT1800Events[j];
        cout << "Event " << count+1 << " out of " << nevents << endl;
        count += 1;
        out9.cd();
        MyTree->Fill();
    }
    MyTree->Write();
    t.Stop();
    t.Print();
}