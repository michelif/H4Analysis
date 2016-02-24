#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis() {
	ofstream myfile;
	TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
	TCanvas* canvas = new TCanvas("canvas","Time and Frequency",1200,600);
    TCanvas* canvas2 = new TCanvas("canvas2","Frequency",800,1200);
	canvas->Divide(2,2);
	TTree* h4 = (TTree*) f.Get("h4");
	TH2F* WavePulse = new TH2F ("WavePulse", "Wave Pulse", 1000, -0.1, 199.9, 850, -50, 800);
	TH1F* PulseTime = new TH1F ("PulseTime", "1D Histo", 1000, -0.1, 199.9); //nanoseconds
	h4->Draw("WF_val:WF_time>>WavePulse", "WF_ch==2 && event==1 && spill==1");
	myfile.open ("storage.txt");
	myfile << "x" << "\t \t" << "y" << endl;
	for (Int_t i=0; i<1000; i++) {
		for (Int_t j=0; j<4096; j++) {
			if (WavePulse->GetBinContent(i, j) != 0) {
				myfile << i << "\t \t" << j-50 << endl; //subtract j by distance from lower limit of y to 0 (if negative)
				PulseTime->SetBinContent(i,j-50);
			}
		}
	}
	myfile.close();
	new TFile("storage.txt");
    TGraph *MyGraph = new TGraph("storage.txt");
    //MyGraph->Draw();
    //PulseTime->DrawClone();
    canvas->cd(1);
    WavePulse->GetXaxis()->SetTitle("Time (ns)");
    WavePulse->GetYaxis()->SetTitle("Energy (GeV)");
    WavePulse->DrawClone();
    //h4->Draw("WF_val:WF_time", "WF_ch==2 && event==1 && spill==1");
    canvas->cd(2);
    TH1F* PulseFreq = new TH1F ("PulseFreq", "Pulse FFT", 1000, -0.1, 999.9);
    PulseTime->FFT(PulseFreq, "MAG");
    PulseFreq->SetLineColor(kRed);
    PulseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    PulseFreq->GetYaxis()->SetTitle("Amplitude");
    PulseFreq->DrawClone();
    gPad->SetLogy();
    canvas->cd(3);
    TH1F *NoiseTime = new TH1F ("NoiseTime", "Wave Pulse Noise", 1000, -0.1, 200*0.2-0.1);
    for (Int_t l=0; l<PulseTime->GetNbinsX(); l++) {
    	NoiseTime->SetBinContent(l, PulseTime->GetBinContent(l%200));
    }
    //NoiseTime->GetXaxis()()->SetRangeUser(10,990);
    NoiseTime->GetXaxis()->SetTitle("Time (ns)");
    NoiseTime->GetYaxis()->SetTitle("Energy (GeV)");
    NoiseTime->Draw();
    canvas->cd(4);
    TH1F* NoiseFreq = new TH1F ("NoiseFreq", "Noise FFT", 200, -0.1, 199.9);
    NoiseTime->FFT(NoiseFreq, "MAG");
    NoiseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    NoiseFreq->GetYaxis()->SetTitle("Amplitude");
    NoiseFreq->Draw();
    gPad->SetLogy();
    canvas2->Divide(1,3);
    canvas2->cd(1);
    PulseFreq->DrawClone();
    gPad->SetLogy();
    canvas2->cd(2);
    NoiseFreq->DrawClone();
    gPad->SetLogy();
    canvas2->cd(3);
    Double_t PulseFreqIntegral = PulseFreq->Integral();
    Double_t NoiseFreqIntegral = NoiseFreq->Integral();
    PulseFreq->Scale(1/PulseFreqIntegral);
    NoiseFreq->Scale(1/NoiseFreqIntegral);
    PulseFreq->SetTitle("Normalized Pulse and Noise FFT");
    PulseFreq->Draw();
    NoiseFreq->Draw("same");
    gPad->SetLogy();
}


// Double_t binCenter = xaxis->GetBinCenter(bin); //where bin is a number corresponding to a bin
// Int_t BinSize = htemp->GetNBins();
// TAxis *xaxis = htemp->GetXaxis();
// TAxis *yaxis = htemp->GetYaxis();

// TH2F* WavePulse=(TH2F*)gDirectory->Get("htemp");

// TTree *MyTree = new TTree("MyTree", "MyTree");
// MyTree->ReadFile("storage.txt", "x:y");
// MyTree->Draw("y:x","");
// canvas->cd(1);
// WavePulse->DrawClone("colz");
//NoiseTime->SetTitle("FFT");
