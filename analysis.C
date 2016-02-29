#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void analysis() {
    Int_t nbins = 800;
	ofstream myfile;
	TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
	TCanvas* canvas = new TCanvas("canvas","Time and Frequency",1200,600);
    TCanvas* canvas2 = new TCanvas("canvas2","Frequency",800,1200);
    TCanvas* canvas3 = new TCanvas("canvas3","Pure Signal",800,1200);
	canvas->Divide(2,2);
	TTree* h4 = (TTree*) f.Get("h4");
	TH2F* WavePulse = new TH2F ("WavePulse", "Wave Pulse", nbins, -0.1, 159.9, 850, -50, 800);
	TH1F* PulseTime = new TH1F ("PulseTime", "1D Histo", nbins, -0.1, 159.9); //nanoseconds
	h4->Draw("WF_val:WF_time>>WavePulse", "WF_ch==2 && event==1 && spill==1");
	myfile.open ("storage.txt");
	myfile << "x" << "\t \t" << "y" << endl;
	for (Int_t i=0; i<nbins; i++) {
		for (Int_t j=0; j<4096; j++) {
			if (WavePulse->GetBinContent(i+1, j) != 0) {
				myfile << i << "\t \t" << j-50 << endl; //subtract j by distance from lower limit of y to 0 (if negative)
				PulseTime->SetBinContent(i+1,j-50);
			}
		}
	}
	myfile.close();
	new TFile("storage.txt");
    TGraph *MyGraph = new TGraph("storage.txt");
    //MyGraph->Draw();
    //PulseTime->DrawClone();
    canvas->cd(1);
    //WavePulse->GetXaxis()->SetTitle("Time (ns)");
    //WavePulse->GetYaxis()->SetTitle("Amplitude");
    //WavePulse->DrawClone();
    PulseTime->GetXaxis()->SetTitle("Time (ns)");
    PulseTime->GetYaxis()->SetTitle("Amplitude");
    PulseTime->DrawClone();
    //h4->Draw("WF_val:WF_time", "WF_ch==2 && event==1 && spill==1");
    canvas->cd(2);
    TH1F* PulseFreq = new TH1F ("PulseFreq", "Pulse FFT", nbins, -0.1, 999.9);
    TH1F* PulsePhase = new TH1F ("PulsePhase", "Pulse Phase", nbins, -0.1, 999.9);
    PulseTime->FFT(PulseFreq, "MAG");
    PulseTime->FFT(PulsePhase, "PH");
    PulseFreq->SetLineColor(kRed);
    PulseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    PulseFreq->GetYaxis()->SetTitle("Amplitude");
    PulseFreq->DrawClone();
    gPad->SetLogy();
    canvas->cd(3);
    TH1F *NoiseTime = new TH1F ("NoiseTime", "Wave Pulse Noise", nbins, -0.1, 159.9);
    cout << PulseTime->GetNbinsX() << endl;
    for (Int_t k=0; k<PulseTime->GetNbinsX(); k++) {
        if (k<851) {
            NoiseTime->SetBinContent(k+1, PulseTime->GetBinContent((k+1)%250));
        }
        else  	NoiseTime->SetBinContent(k+1,0);
        //NoiseTime->SetBinContent(k+1, sin(k/10.)/10.);
    }
    //NoiseTime->GetXaxis()()->SetRangeUser(10,990);
    NoiseTime->GetXaxis()->SetTitle("Time (ns)");
    NoiseTime->GetYaxis()->SetTitle("Amplitude");
    NoiseTime->Draw();
    canvas->cd(4);
    TH1F* NoiseFreq = new TH1F ("NoiseFreq", "Noise FFT", nbins, -0.1, 999.9);
    NoiseTime->FFT(NoiseFreq, "MAG");
    NoiseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    NoiseFreq->GetYaxis()->SetTitle("Amplitude");
    NoiseFreq->Draw();
    gPad->SetLogy();
    canvas2->Divide(2,2);
    canvas2->cd(1);
    PulseFreq->DrawClone();
    gPad->SetLogy();
    canvas2->cd(2);
    NoiseFreq->DrawClone();
    gPad->SetLogy();
    canvas2->cd(3);
    //Double_t PulseFreqIntegral = PulseFreq->Integral();
    //Double_t NoiseFreqIntegral = NoiseFreq->Integral();
    //PulseFreq->Scale(1/PulseFreqIntegral);
    //NoiseFreq->Scale(1/NoiseFreqIntegral);
    //PulseFreq->SetTitle("Normalized Pulse and Noise FFT");
    PulseFreq->SetTitle("Pulse and Noise FFT Comparison");
    PulseFreq->Draw();
    NoiseFreq->Draw("same");
    gPad->SetLogy();
    TH1F* UnscaledSignalFreq = new TH1F ("UnscaledSignalFreq", "Unscaled Signal Frequency", nbins, -0.1, 999.9);
    for (Int_t l=0; l<nbins; l++) {
        UnscaledSignalFreq->SetBinContent(l+1, (PulseFreq->GetBinContent(l+1)-NoiseFreq->GetBinContent(l+1))/PulseFreq->GetBinContent(l+1));
    }
    canvas2->cd(4);
    UnscaledSignalFreq->DrawClone();
    TH1F* SignalFreq = new TH1F ("SignalFreq", "Signal Frequency", nbins, -0.1, 999.9);
    for (Int_t m=0; m<nbins; m++) {
        SignalFreq->SetBinContent(m+1, UnscaledSignalFreq->GetBinContent(m+1)*PulseFreq->GetBinContent(m+1));
    }
    Double_t *re_full = new Double_t[nbins];
    Double_t *im_full = new Double_t[nbins];
    for (Int_t n=0; n<nbins; n++) {
        (re_full)[n]=(SignalFreq->GetBinContent(n+1)*cos(PulsePhase->GetBinContent(n+1)));
        (im_full)[n]=(SignalFreq->GetBinContent(n+1)*sin(PulsePhase->GetBinContent(n+1)));
    }
    TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
    invFFT->SetPointsComplex(re_full, im_full);
    invFFT->Transform();
    TH1 *Signal = 0;
    Signal = TH1::TransformHisto(invFFT,Signal,"Re");
    Signal->SetTitle("Recovered Signal 'S'");
    TH1F* fancysignal = new TH1F ("fancysignal", "Recovered Signal", nbins, -0.1, 159.9);
    for (Int_t p=0; p<nbins; p++) {
        fancysignal->SetBinContent(p+1, Signal->GetBinContent(p+1)/nbins);
    }
    canvas3->Divide(1,2);
    canvas3->cd(1);
    PulseTime->DrawClone();
    canvas3->cd(2);
    fancysignal->GetXaxis()->SetTitle("Time (ns)");
    fancysignal->GetYaxis()->SetTitle("Amplitude");
    fancysignal->SetLineColor(kRed);
    fancysignal->Draw();
    PulseTime->DrawClone("same");
    gPad->SetGrid();
    //Double_t N = nbins;
    //TVirtualFFT *SignalTime = TVirtualFFT::FFT(1, &N, "C2R");
    //SignalTime->SetPointsComplex(N,SignalFreq->GetBinContent(N));
    //SignalTime->Transform();

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
