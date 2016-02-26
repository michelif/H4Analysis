#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void testanalysis() {
	//ofstream myfile;
	//TFile f("/home/marko/H4Analysis/ntuples/analysis_4443.root");
	TCanvas* canvas = new TCanvas("canvas","Time & Frequency",1200,600);
    TCanvas* canvas2 = new TCanvas("canvas2","Frequency",800,1200);
    TCanvas* canvas3 = new TCanvas("canvas3","Frequency2", 800, 1200);
	canvas->Divide(2,3);
    Int_t nbins = 100;
    TF1 *func = new TF1("func", "1*TMath::Exp(-((x-TMath::Pi())^2)/(2*0.65^2))", 0, 2*TMath::Pi());
    canvas->cd(1);
    func->SetTitle("Signal");
    func->DrawClone();
    TH1F* funchisto = new TH1F("funchisto", "Signal", nbins, 0, 2*TMath::Pi());
    TH1F* funcFFT = new TH1F("funcFFT", "Signal FFT S", nbins, 0, 2*TMath::Pi());
    for (Int_t i=0; i<=nbins; i++) {
        Double_t x = Double_t(i)/nbins*(2*TMath::Pi());
        funchisto->SetBinContent(i+1, func->Eval(x));
    }
    //func->FFT(funcFFT, "MAG");
    canvas->cd(2);
    funchisto->FFT(funcFFT, "MAG");
    funcFFT->DrawClone();
    gPad->SetLogy();
    TF1* noisy = new TF1("noisy", "1*TMath::Exp(-((x-TMath::Pi())^2)/(2*0.65^2)) + sin(6*x)/10", 0, 2*TMath::Pi());
    canvas->cd(3);
    noisy->SetTitle("Signal & Noise");
    noisy->DrawClone();
    TH1F* noisyhisto = new TH1F("noisyhisto", "Signal & Noise ", nbins, 0, 2*TMath::Pi());
    TH1F* noisyFFT = new TH1F("noisyFFT", "Signal & Noise FFT S+N", nbins, 0, 2*TMath::Pi());
    TH1F* noisyPH = new TH1F("noisyPH", "Noisy Signal Phase", nbins, 0, 2*TMath::Pi());
    for (Int_t j=0; j<=nbins; j++) {
        Double_t x = (Double_t(j)/nbins)*(2*TMath::Pi());
        noisyhisto->SetBinContent(j+1, noisy->Eval(x));
    }
    noisyhisto->FFT(noisyFFT, "MAG");
    noisyhisto->FFT(noisyPH, "PH");
    canvas->cd(4);
    noisyFFT->DrawClone();
    gPad->SetLogy();
    TF1 *noise = new TF1("noise", "sin(6*x)/10", 0, 2*TMath::Pi());
    canvas->cd(5);
    noise->DrawClone();    
    TH1F* noisehisto = new TH1F("noisehisto", "Noise", nbins, 0, 2*TMath::Pi());
    TH1F* noiseFFT = new TH1F("noiseFFT", "Noise FFT", nbins, 0, 2*TMath::Pi());
    for (Int_t k=0; k<=nbins; k++) {
        Double_t x = Double_t(k)/nbins*(2*TMath::Pi());
        noisehisto->SetBinContent(k+1, noise->Eval(x));
    }
    noisehisto->FFT(noiseFFT, "MAG");
    canvas->cd(6);
    noiseFFT->DrawClone();
    gPad->SetLogy();
    canvas2->Divide(1,3);
    canvas2->cd(1);
    noisyFFT->DrawClone();
    gPad->SetLogy();
    canvas2->cd(2);
    noiseFFT->DrawClone();
    gPad->SetLogy();
    //Double_t noisyFFTIntegral = noisyFFT->Integral();
    //Double_t noiseFFTIntegral = noiseFFT->Integral();
    //noisyFFT->Scale(1/noisyFFTIntegral);
    //noiseFFT->Scale(1/noiseFFTIntegral);
    noisyFFT->SetTitle("Normalized Noisy Signal & Noise FFT");
    canvas2->cd(3);
    noisyFFT->Draw();
    noiseFFT->Draw("same");
    gPad->SetLogy();
    TH1F* normsignalFFT = new TH1F ("normsignalFFT", "Normalized Signal FFT S/(S+N)", nbins, 0, 2*TMath::Pi());
    for (Int_t l=0; l<nbins; l++) {
        //normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));
        //normsignalFFT->SetBinContent(l+1, pow((noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1),2));
        normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));

    }
    TH1F* signalFFT = new TH1F ("signalFFT", "Signal FFT S/(S+N)*(S+N)", nbins, 0, 2*TMath::Pi());
    for (Int_t q=0; q<nbins; q++) {
        signalFFT->SetBinContent(q+1, normsignalFFT->GetBinContent(q+1)*noisyFFT->GetBinContent(q+1));
        //signalFFT->SetBinContent(q+1, noisyFFT->GetBinContent(q+1)); //checking if output = input, without taking noise out

    }
    canvas3->Divide(1,3);
    canvas3->cd(1);
    normsignalFFT->Draw();
    canvas3->cd(2);
    signalFFT->Draw();
    canvas3->cd(3);
    Double_t *re_full = new Double_t[nbins];
    Double_t *im_full = new Double_t[nbins];
    for (Int_t m=0; m<nbins; m++) {
        (re_full)[m]=(signalFFT->GetBinContent(m+1)*cos(noisyPH->GetBinContent(m+1)));
        (im_full)[m]=(signalFFT->GetBinContent(m+1)*sin(noisyPH->GetBinContent(m+1)));
    }
    TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
    invFFT->SetPointsComplex(re_full, im_full);
    invFFT->Transform();
    TH1 *Signal = 0;
    Signal = TH1::TransformHisto(invFFT,Signal,"Re");
    Signal->SetTitle("Recovered Signal 'S'");
    //Signal->Draw();
    TH1F* fancysignal = new TH1F ("fancysignal", "Recovered Signal", nbins, 0, 2*TMath::Pi());
    for (Int_t p=0; p<nbins; p++) {
        fancysignal->SetBinContent(p+1, Signal->GetBinContent(p+1)/nbins);
    }
    func->Draw();
    fancysignal->Draw("same");
    



	
    //TTree* h4 = (TTree*) f.Get("h4");
	//TH2F* WavePulse = new TH2F ("WavePulse", "Wave Pulse", nbins0, -0.1, 199.9, 850, -50, 800);
	//TH1F* PulseTime = new TH1F ("PulseTime", "1D Histo", nbins0, -0.1, 199.9);
	//h4->Draw("WF_val:WF_time>>WavePulse", "WF_ch==2 && event==1 && spill==1");
	//myfile.open ("storage.txt");
	//myfile << "x" << "\t \t" << "y" << endl;
	//for (Int_t i=0; i<nbins0; i++) {
	//	for (Int_t j=0; j<4096; j++) {
	//		if (WavePulse->GetBinContent(i, j) != 0) {
	//			myfile << i << "\t \t" << j-50 << endl; //subtract j by distance from lower limit of y to 0 (if negative)
	//			PulseTime->SetBinContent(i,j-50);
	//		}
	//	}
	//}
	//myfile.close();
	//new TFile("storage.txt");
    //TGraph *MyGraph = new TGraph("storage.txt");
    //MyGraph->Draw();
    //PulseTime->DrawClone();
    //canvas->cd(1);
    //func->GetXaxis()->SetTitle("Time (ns)");
    //func->GetYaxis()->SetTitle("Energy (GeV)");
    //func->DrawClone();
    ////h4->Draw("WF_val:WF_time", "WF_ch==2 && event==1 && spill==1");
    //canvas->cd(2);
    //TH1F* PulseFreq = new TH1F ("PulseFreq", "Pulse FFT", nbins0, -0.1, 999.9);
    //PulseTime->FFT(PulseFreq, "MAG");
    //PulseFreq->SetLineColor(kRed);
    //PulseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    //PulseFreq->GetYaxis()->SetTitle("Amplitude");
    //PulseFreq->DrawClone();
    //gPad->SetLogy();
    //canvas->cd(3);
    //TH1F *NoiseTime = new TH1F ("NoiseTime", "Wave Pulse Noise", nbins0, -0.1, 200*0.2-0.1);
    //for (Int_t l=0; l<PulseTime->GetNbinsX(); l++) {
    //	NoiseTime->SetBinContent(l, PulseTime->GetBinContent(l%200));
    //}
    ////NoiseTime->GetXaxis()()->SetRangeUser(10,990);
    //NoiseTime->GetXaxis()->SetTitle("Time (ns)");
    //NoiseTime->GetYaxis()->SetTitle("Energy (GeV)");
    //NoiseTime->Draw();
    //canvas->cd(4);
    //TH1F* NoiseFreq = new TH1F ("NoiseFreq", "Noise FFT", 200, -0.1, 199.9);
    //NoiseTime->FFT(NoiseFreq, "MAG");
    //NoiseFreq->GetXaxis()->SetTitle("Frequency ns^-1");
    //NoiseFreq->GetYaxis()->SetTitle("Amplitude");
    //NoiseFreq->Draw();
    //gPad->SetLogy();
    //canvas2->Divide(1,3);
    //canvas2->cd(1);
    //PulseFreq->DrawClone();
    //gPad->SetLogy();
    //canvas2->cd(2);
    //NoiseFreq->DrawClone();
    //gPad->SetLogy();
    //canvas2->cd(3);
    //Double_t PulseFreqIntegral = PulseFreq->Integral();
    //Double_t NoiseFreqIntegral = NoiseFreq->Integral();
    //PulseFreq->Scale(1/PulseFreqIntegral);
    //NoiseFreq->Scale(1/NoiseFreqIntegral);
    //PulseFreq->SetTitle("Normalized Pulse and Noise FFT");
    //PulseFreq->Draw();
    //NoiseFreq->Draw("same");
    //gPad->SetLogy();
}