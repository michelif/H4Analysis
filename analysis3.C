#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis3(Int_t nspill = 1) {
    Int_t nbins = 800;
    Int_t j;
    TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
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
    for (int iEntry = 0; iEntry<nevents; iEntry++){
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
        sprintf(name,"SignalEventFFT100Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT100Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT250Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT250Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT500Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT500Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT750Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT750Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT750Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1000Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1000Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1000Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1250Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1250Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1F *SignalEventFFT1500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEventFFT1500Events %d",z+1);
        sprintf(title,"SignalFFT Event%d Histo", z+1);
        SignalEventFFT1500Events[z+1] = new TH1F(name,title, nbins, 0, 5);
    }
    TH1 *SignalEvent100Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent100Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent100Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent250Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent250Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent500Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent500Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent750Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent750Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent750Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1000Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1000Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1000Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1250Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1250Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1250Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
    }
    TH1 *SignalEvent1500Events[nevents];
    for (Int_t z=0;z<nevents;z++) {
        sprintf(name,"SignalEvent1500Events %d",z+1);
        sprintf(title,"Signal Event%d Histo", z+1);
        SignalEvent1500Events[z+1] = new TH1F(name,title, nbins, -0.1, 159.9);
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
    Int_t count = 0;
    TFile f1("NormNoiseFFTs.root", "read");
    TH1F* NormNoiseFFT100Events = (TH1F*) f1.Get("NormNoiseFFT100Events");
    TH1F* NormNoiseFFT250Events = (TH1F*) f1.Get("NormNoiseFFT250Events");
    TH1F* NormNoiseFFT500Events = (TH1F*) f1.Get("NormNoiseFFT500Events");
    TH1F* NormNoiseFFT750Events = (TH1F*) f1.Get("NormNoiseFFT750Events");
    TH1F* NormNoiseFFT1000Events = (TH1F*) f1.Get("NormNoiseFFT1000Events");
    TH1F* NormNoiseFFT1250Events = (TH1F*) f1.Get("NormNoiseFFT1250Events");
    TH1F* NormNoiseFFT1500Events = (TH1F*) f1.Get("NormNoiseFFT1500Events");
    TFile out1("NormalizedSignalNoise100Events.root", "recreate");
    TFile out2("NormalizedSignalNoise250Events.root", "recreate");
    TFile out3("NormalizedSignalNoise500Events.root", "recreate");
    TFile out4("NormalizedSignalNoise750Events.root", "recreate");
    TFile out5("NormalizedSignalNoise1000Events.root", "recreate");
    TFile out6("NormalizedSignalNoise1250Events.root", "recreate");
    TFile out7("NormalizedSignalNoise1500Events.root", "recreate");
    TFile out8("Pre&PostRMS.root", "recreate");
    TTree *MyTree = new TTree ("MyTree", "MyTree");
    MyTree->Branch("rms", &rms, "rms/F");
    MyTree->Branch("rmsprime100", &rmsprime100, "rmsprime100/F");
    MyTree->Branch("rmsprime250", &rmsprime250, "rmsprime250/F");
    MyTree->Branch("rmsprime500", &rmsprime500, "rmsprime500/F");
    MyTree->Branch("rmsprime750", &rmsprime750, "rmsprime750/F");
    MyTree->Branch("rmsprime1000", &rmsprime1000, "rmsprime1000/F");
    MyTree->Branch("rmsprime1250", &rmsprime1250, "rmsprime1250/F");
    MyTree->Branch("rmsprime1500", &rmsprime1500, "rmsprime1500/F");
    MyTree->Branch("count", &count, "count/I");
    TString plot;
    TString cut;
    TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
    TH1F* NormSigFFT = new TH1F ("NormSigFFT", "Normalized SignalNoise FFT", nbins, 0, 5);
    TH1F* Phase = new TH1F ("Phase", "Phase", nbins, -0.1, 799.9);
    TStopwatch t;
    t.Start(); //1 hour runtime
    for (j=1;j<nevents;j++) {
    //for (j=1;j<11;j++) {
        TH1 *Throwaway1 = 0;
        TH1 *Throwaway2 = 0;
        TH1 *Throwaway3 = 0;
        TH1 *Throwaway4 = 0;
        TH1 *Throwaway5 = 0;
        TH1 *Throwaway6 = 0;
        TH1 *Throwaway7 = 0;
        TVirtualFFT *invFFT1 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT2 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT3 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT4 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT5 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT6 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        TVirtualFFT *invFFT7 = TVirtualFFT::FFT(1, &nbins, "C2R M K");
        std::vector<Float_t> y1;
        std::vector<Float_t> y2;
        std::vector<Float_t> y3;
        std::vector<Float_t> y4;
        std::vector<Float_t> y5;
        std::vector<Float_t> y6;
        std::vector<Float_t> y7;
        Float_t mu1 = 0;
        Float_t mu2 = 0;
        Float_t mu3 = 0;
        Float_t mu4 = 0;
        Float_t mu5 = 0;
        Float_t mu6 = 0;
        Float_t mu7 = 0;
        plot = "WF_val:WF_time>>TempHisto";
        cut = "WF_ch==2 && spill==";
        cut += nspill;
        cut += " && event==";
        cut += j;
        h4->Draw(plot, cut, "goff");
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
            delete PhaseHistoEvent[j];
            delete SignalEvent100Events[j];
            delete SignalEvent250Events[j];
            delete SignalEvent500Events[j];
            delete SignalEvent750Events[j];
            delete SignalEvent1000Events[j];
            delete SignalEvent1250Events[j];
            delete SignalEvent1500Events[j];
            continue;
        }
        for (Int_t i=0; i<nbins; i++) {
            for (Int_t k=0; k<1000; k++) {
                if (TempHisto->GetBinContent(i+1, k) != 0) {
                    HistoEvent[j]->SetBinContent(i+1,k*0.92-120); //going from 2D to 1D histogram
                }
            }
        }
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
        for (Int_t q=0;q<nbins;q++) {
            SignalEventFFT100Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT100Events->GetBinContent(q+1))); // S+N - N = S in frequency domain
            SignalEventFFT250Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT250Events->GetBinContent(q+1)));
            SignalEventFFT500Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT500Events->GetBinContent(q+1)));
            SignalEventFFT750Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT750Events->GetBinContent(q+1)));
            SignalEventFFT1000Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1000Events->GetBinContent(q+1)));
            SignalEventFFT1250Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1250Events->GetBinContent(q+1)));
            SignalEventFFT1500Events[j]->SetBinContent(q+1, (NewHistoEventFFT[j]->GetBinContent(q+1)-NormNoiseFFT1500Events->GetBinContent(q+1)));
        }
        NormNoiseFFT100Events->Scale(1/rms);
        NormNoiseFFT250Events->Scale(1/rms);
        NormNoiseFFT500Events->Scale(1/rms);
        NormNoiseFFT750Events->Scale(1/rms);
        NormNoiseFFT1000Events->Scale(1/rms);
        NormNoiseFFT1250Events->Scale(1/rms);
        NormNoiseFFT1500Events->Scale(1/rms);
        HistoEvent[j]->FFT(PhaseHistoEvent[j], "PH");
        Double_t *re_full1 = new Double_t[nbins];
        Double_t *im_full1 = new Double_t[nbins];
        Double_t *re_full2 = new Double_t[nbins];
        Double_t *im_full2 = new Double_t[nbins];
        Double_t *re_full3 = new Double_t[nbins];
        Double_t *im_full3 = new Double_t[nbins];
        Double_t *re_full4 = new Double_t[nbins];
        Double_t *im_full4 = new Double_t[nbins];
        Double_t *re_full5 = new Double_t[nbins];
        Double_t *im_full5 = new Double_t[nbins];
        Double_t *re_full6 = new Double_t[nbins];
        Double_t *im_full6 = new Double_t[nbins];
        Double_t *re_full7 = new Double_t[nbins];
        Double_t *im_full7 = new Double_t[nbins];
        for (Int_t n=0; n<nbins; n++) {
            (re_full1)[n]=(SignalEventFFT100Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1))); //magnitude * cos of phase
            (im_full1)[n]=(SignalEventFFT100Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1))); //magnitude * sin of phase
            (re_full2)[n]=(SignalEventFFT250Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full2)[n]=(SignalEventFFT250Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (re_full3)[n]=(SignalEventFFT500Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full3)[n]=(SignalEventFFT500Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (re_full4)[n]=(SignalEventFFT750Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full4)[n]=(SignalEventFFT750Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (re_full5)[n]=(SignalEventFFT1000Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full5)[n]=(SignalEventFFT1000Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (re_full6)[n]=(SignalEventFFT1250Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full6)[n]=(SignalEventFFT1250Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (re_full7)[n]=(SignalEventFFT1500Events[j]->GetBinContent(n+1)*cos(PhaseHistoEvent[j]->GetBinContent(n+1)));
            (im_full7)[n]=(SignalEventFFT1500Events[j]->GetBinContent(n+1)*sin(PhaseHistoEvent[j]->GetBinContent(n+1)));
        }
        invFFT1->SetPointsComplex(re_full1, im_full1);
        invFFT1->Transform();
        Throwaway1 = TH1::TransformHisto(invFFT1, Throwaway1, "Re");
        Throwaway1->SetName("Re_1");
        invFFT2->SetPointsComplex(re_full1, im_full1);
        invFFT2->Transform();
        Throwaway2 = TH1::TransformHisto(invFFT2, Throwaway2, "Re");
        Throwaway2->SetName("Re_2");
        invFFT3->SetPointsComplex(re_full3, im_full3);
        invFFT3->Transform();
        Throwaway3 = TH1::TransformHisto(invFFT3, Throwaway3, "Re");
        Throwaway3->SetName("Re_3");
        invFFT4->SetPointsComplex(re_full4, im_full4);
        invFFT4->Transform();
        Throwaway4 = TH1::TransformHisto(invFFT4, Throwaway4, "Re");
        Throwaway4->SetName("Re_4");
        invFFT5->SetPointsComplex(re_full5, im_full5);
        invFFT5->Transform();
        Throwaway5 = TH1::TransformHisto(invFFT5, Throwaway5, "Re");
        Throwaway5->SetName("Re_5");
        invFFT6->SetPointsComplex(re_full6, im_full6);
        invFFT6->Transform();
        Throwaway6 = TH1::TransformHisto(invFFT6, Throwaway6, "Re");
        Throwaway6->SetName("Re_6");
        invFFT7->SetPointsComplex(re_full7, im_full7);
        invFFT7->Transform();
        Throwaway7 = TH1::TransformHisto(invFFT7, Throwaway7, "Re");
        Throwaway7->SetName("Re_7");
        for (Int_t p=0; p<nbins; p++) {
            SignalEvent100Events[j]->SetBinContent(p+1, Throwaway1->GetBinContent(p+1)/nbins);
            SignalEvent250Events[j]->SetBinContent(p+1, Throwaway2->GetBinContent(p+1)/nbins);
            SignalEvent500Events[j]->SetBinContent(p+1, Throwaway3->GetBinContent(p+1)/nbins);
            SignalEvent750Events[j]->SetBinContent(p+1, Throwaway4->GetBinContent(p+1)/nbins);
            SignalEvent1000Events[j]->SetBinContent(p+1, Throwaway5->GetBinContent(p+1)/nbins);
            SignalEvent1250Events[j]->SetBinContent(p+1, Throwaway6->GetBinContent(p+1)/nbins);
            SignalEvent1500Events[j]->SetBinContent(p+1, Throwaway7->GetBinContent(p+1)/nbins);
        }
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
        for (Int_t q=0; q<225; q++) {
            y1.push_back(SignalEvent100Events[j]->GetBinContent(q+1));
            mu1 += y1.at(q);
            y2.push_back(SignalEvent250Events[j]->GetBinContent(q+1));
            mu2 += y2.at(q);
            y3.push_back(SignalEvent500Events[j]->GetBinContent(q+1));
            mu3 += y3.at(q);
            y4.push_back(SignalEvent750Events[j]->GetBinContent(q+1));
            mu4 += y4.at(q);
            y5.push_back(SignalEvent1000Events[j]->GetBinContent(q+1));
            mu5 += y5.at(q);
            y6.push_back(SignalEvent1250Events[j]->GetBinContent(q+1));
            mu6 += y6.at(q);
            y7.push_back(SignalEvent1500Events[j]->GetBinContent(q+1));
            mu7 += y7.at(q);
        }
        mu1 = mu1/225;
        mu2 = mu2/225;
        mu3 = mu3/225;
        mu4 = mu4/225;
        mu5 = mu5/225;
        mu6 = mu6/225;
        mu7 = mu7/225;
        rmsprime100 = 0;
        rmsprime250 = 0;
        rmsprime500 = 0;
        rmsprime750 = 0;
        rmsprime1000 = 0;
        rmsprime1250 = 0;
        rmsprime1500 = 0;
        for (Int_t q=0; q<225; q++) {
            rmsprime100 += (y1[q]-mu1)*(y1[q]-mu1);
            rmsprime250 += (y2[q]-mu2)*(y2[q]-mu2);
            rmsprime500 += (y3[q]-mu3)*(y3[q]-mu3);
            rmsprime750 += (y4[q]-mu4)*(y4[q]-mu4);
            rmsprime1000 += (y5[q]-mu5)*(y5[q]-mu5);
            rmsprime1250 += (y6[q]-mu6)*(y6[q]-mu6);
            rmsprime1500 += (y7[q]-mu7)*(y7[q]-mu7);
        }
        rmsprime100 = sqrt(rmsprime100/225);
        rmsprime250 = sqrt(rmsprime250/225);
        rmsprime500 = sqrt(rmsprime500/225);
        rmsprime750 = sqrt(rmsprime750/225);
        rmsprime1000 = sqrt(rmsprime1000/225);
        rmsprime1250 = sqrt(rmsprime1250/225);
        rmsprime1500 = sqrt(rmsprime1500/225);
        y1.clear();
        y2.clear();
        y3.clear();
        y4.clear();
        y5.clear();
        y6.clear();
        y7.clear();
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
        cout << "Event " << count+1 << " out of " << nevents << endl;
        count += 1;
        delete invFFT1;
        delete invFFT2;
        delete invFFT3;
        delete invFFT4;
        delete invFFT5;
        delete invFFT6;
        delete invFFT7;
        delete Throwaway1;
        delete Throwaway2;
        delete Throwaway3;
        delete Throwaway4;
        delete Throwaway5;
        delete Throwaway6;
        delete Throwaway7;
        out8.cd();
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