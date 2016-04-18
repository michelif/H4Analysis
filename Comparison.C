#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;
void Comparison() {
    Int_t nbins = 800;
    Int_t count = 0;
    Float_t muPre100;
    Float_t muPost100;
    Float_t muAveOrig100;
    Float_t muAveSub100;
    Float_t rmsPre100;
    Float_t rmsPost100;
    Float_t rmsAveOrig100 = 0;
    Float_t rmsAveSub100 = 0;
    std::vector<Float_t> vPre100;
    std::vector<Float_t> vPost100;
    std::vector<Float_t> vAveOrig100;
    std::vector<Float_t> vAveSub100;
    Float_t muPre1500;
    Float_t muPost1500;
    Float_t muAveOrig1500;
    Float_t muAveSub1500;
    Float_t rmsPre1500;
    Float_t rmsPost1500;
    Float_t rmsAveOrig1500 = 0;
    Float_t rmsAveSub1500 = 0;
    std::vector<Float_t> vPre1500;
    std::vector<Float_t> vPost1500;
    std::vector<Float_t> vAveOrig1500;
    std::vector<Float_t> vAveSub1500;
    TString HistoName;
    TH1F* original;
    TH1F* subtracted;
    TH1F *originalsum = new TH1F ("originalsum", "sum of original histos", nbins, -0.1, 159.9);
    TH1F *subtractedsum = new TH1F ("subtractedsum", "sum of subtracted histos", nbins, -0.1, 159.9);

    TFile f1("NormalizedSignalNoise100Events.root", "read");
    TFile f2("NormalizedSignalNoise1500Events.root", "read");

    Int_t nevents = f1.GetListOfKeys()->GetEntries();
    nevents = nevents/2; //divide by 2 since there are 2 histograms (before and after subtraction) for each event
    for (Int_t j=1;j<=nevents;j++) {
        if (j==129) continue;
        if (j==670) continue;
        vPre100.clear();
        vPost100.clear();
        muPre100 = 0;
        muPost100 = 0;
        rmsPre100 = 0;
        rmsPost100 = 0;
        HistoName = "NewHistoEvent";
        HistoName += j;
        original = (TH1F*) f1.Get(HistoName);
        HistoName = "SignalEvent100Events ";
        HistoName += j;
        subtracted = (TH1F*) f1.Get(HistoName);
        for (Int_t q=0; q<225; q++) {
            vPre100.push_back(original->GetBinContent(q+1));
            muPre100 += vPre100.at(q);
            vPost100.push_back(subtracted->GetBinContent(q+1));
            muPost100 += vPost100.at(q);
        }
        muPre100 = muPre100/225;
        muPost100 = muPost100/225;
        for (Int_t q=0; q<225; q++) {
            rmsPre100 += (vPre100[q]-muPre100)*(vPre100[q]-muPre100);
            rmsPost100 += (vPost100[q]-muPost100)*(vPost100[q]-muPost100);
        }
        rmsPre100 = sqrt(rmsPre100/225);
        rmsPost100 = sqrt(rmsPost100/225);
        if (rmsPre100 <= 6) {
            originalsum->Add(originalsum, original);
            subtractedsum->Add(subtractedsum, subtracted);
            count += 1;
        }
    }
    originalsum->Scale(1./count);
    subtractedsum->Scale(1./count);
    for (Int_t q=0; q<225; q++) {
        vAveOrig100.push_back(originalsum->GetBinContent(q+1));
        muAveOrig100 += vAveOrig100.at(q);
        vAveSub100.push_back(subtractedsum->GetBinContent(q+1));
        muAveSub100 += vAveSub100.at(q);
    }
    muAveOrig100 = muAveOrig100/225;
    muAveSub100 = muAveSub100/225;
    for (Int_t q=0; q<225; q++) {
            rmsAveOrig100 += (vAveOrig100[q]-muAveOrig100)*(vAveOrig100[q]-muAveOrig100);
            rmsAveSub100 += (vAveSub100[q]-muAveSub100)*(vAveSub100[q]-muAveSub100);
        }
    rmsAveOrig100 = sqrt(rmsAveOrig100/225);
    rmsAveSub100 = sqrt(rmsAveSub100/225);

    //---1500 Events---
    for (Int_t j=1;j<=nevents;j++) {
        if (j==129) continue;
        if (j==670) continue;
        vPre1500.clear();
        vPost1500.clear();
        muPre1500 = 0;
        muPost1500 = 0;
        rmsPre1500 = 0;
        rmsPost1500 = 0;
        HistoName = "NewHistoEvent";
        HistoName += j;
        original = (TH1F*) f2.Get(HistoName);
        HistoName = "SignalEvent1500Events ";
        HistoName += j;
        subtracted = (TH1F*) f2.Get(HistoName);
        for (Int_t q=0; q<225; q++) {
            vPre1500.push_back(original->GetBinContent(q+1));
            muPre1500 += vPre1500.at(q);
            vPost1500.push_back(subtracted->GetBinContent(q+1));
            muPost1500 += vPost1500.at(q);
        }
        muPre1500 = muPre1500/225;
        muPost1500 = muPost1500/225;
        for (Int_t q=0; q<225; q++) {
            rmsPre1500 += (vPre1500[q]-muPre1500)*(vPre1500[q]-muPre1500);
            rmsPost1500 += (vPost1500[q]-muPost1500)*(vPost1500[q]-muPost1500);
        }
        rmsPre1500 = sqrt(rmsPre1500/225);
        rmsPost1500 = sqrt(rmsPost1500/225);
        if (rmsPre1500 <= 6) {
            originalsum->Add(originalsum, original);
            subtractedsum->Add(subtractedsum, subtracted);
            count += 1;
        }
    }
    originalsum->Scale(1./count);
    subtractedsum->Scale(1./count);
    for (Int_t q=0; q<225; q++) {
        vAveOrig1500.push_back(originalsum->GetBinContent(q+1));
        muAveOrig1500 += vAveOrig1500.at(q);
        vAveSub1500.push_back(subtractedsum->GetBinContent(q+1));
        muAveSub1500 += vAveSub1500.at(q);
    }
    muAveOrig1500 = muAveOrig1500/225;
    muAveSub1500 = muAveSub1500/225;
    for (Int_t q=0; q<225; q++) {
            rmsAveOrig1500 += (vAveOrig1500[q]-muAveOrig1500)*(vAveOrig1500[q]-muAveOrig1500);
            rmsAveSub1500 += (vAveSub1500[q]-muAveSub1500)*(vAveSub1500[q]-muAveSub1500);
        }
    rmsAveOrig1500 = sqrt(rmsAveOrig1500/225);
    rmsAveSub1500 = sqrt(rmsAveSub1500/225);
    cout << "100 Event Improvement is " << (rmsAveSub100-rmsAveOrig100)/rmsAveOrig100*-100 << " percent better." << endl;
    cout << "1500 Event Improvement is " << (rmsAveSub1500-rmsAveOrig1500)/rmsAveOrig1500*-100 << " percent better." << endl;
}