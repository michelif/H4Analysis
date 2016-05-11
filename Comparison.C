#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;
void Comparison() {
    Int_t nbins = 800, count = 0;
    Float_t rmsPre100, rmsAveOrig100, rmsAveSub100, rmsPre1800, rmsAveOrig1800, rmsAveSub1800;
    std::vector<Float_t> vPre100, vAveOrig100, vAveSub100, vPre1800, vAveOrig1800, vAveSub1800;
    TString HistoName;
    TH1F *original, *subtracted;
    TH1F *originalsum = new TH1F ("originalsum", "sum of original histos", nbins, -0.1, 159.9);
    TH1F *subtractedsum = new TH1F ("subtractedsum", "sum of subtracted histos", nbins, -0.1, 159.9);

    TFile f1("NormalizedSignalNoise100Events.root", "read");
    TFile f2("NormalizedSignalNoise1800Events.root", "read");

    Int_t nevents = f1.GetListOfKeys()->GetEntries();
    nevents = nevents/2; //divide by 2 since there are 2 histograms (before and after subtraction) for each event

    for (Int_t j=1;j<=nevents;j++) {
        if (j==129) continue;
        if (j==670) continue;
        vPre100.clear();
        HistoName = "NewHistoEvent";
        HistoName += j;
        original = (TH1F*) f1.Get(HistoName);
        HistoName = "SignalEvent100Events";
        HistoName += j;
        HistoName += "_inv";
        subtracted = (TH1F*) f1.Get(HistoName);
        for (Int_t q=0; q<225; q++) {
            vPre100.push_back(original->GetBinContent(q+1));
        }
        rmsPre100 = TMath::RMS(vPre100.begin(), vPre100.end());
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
        vAveSub100.push_back(subtractedsum->GetBinContent(q+1));
    }
    rmsAveOrig100 = TMath::RMS(vAveOrig100.begin(), vAveOrig100.end());
    rmsAveSub100 = TMath::RMS(vAveSub100.begin(), vAveSub100.end());
    //---1800 Events---
    for (Int_t j=1;j<=nevents;j++) {
        if (j==129) continue;
        if (j==670) continue;
        vPre1800.clear();
        HistoName = "NewHistoEvent";
        HistoName += j;
        original = (TH1F*) f2.Get(HistoName);
        HistoName = "SignalEvent1800Events";
        HistoName += j;
        HistoName += "_inv";
        subtracted = (TH1F*) f2.Get(HistoName);
        for (Int_t q=0; q<225; q++) {
            vPre1800.push_back(original->GetBinContent(q+1));
        }
        rmsPre1800 = TMath::RMS(vPre1800.begin(), vPre1800.end());
        if (rmsPre1800 <= 6) {
            originalsum->Add(originalsum, original);
            subtractedsum->Add(subtractedsum, subtracted);
            count += 1;
        }
    }
    originalsum->Scale(1./count);
    subtractedsum->Scale(1./count);
    for (Int_t q=0; q<225; q++) {
        vAveOrig1800.push_back(originalsum->GetBinContent(q+1));
        vAveSub1800.push_back(subtractedsum->GetBinContent(q+1));
    }
    rmsAveOrig1800 = TMath::RMS(vAveOrig1800.begin(), vAveOrig1800.end());
    rmsAveSub1800 = TMath::RMS(vAveSub1800.begin(), vAveSub1800.end());
    cout << "100 Event Improvement is " << (rmsAveSub100-rmsAveOrig100)/rmsAveOrig100*-100 << " percent better." << endl;
    cout << "1800 Event Improvement is " << (rmsAveSub1800-rmsAveOrig1800)/rmsAveOrig1800*-100 << " percent better." << endl;
}