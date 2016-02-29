#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

void testanalysis() {
	TCanvas* canvas = new TCanvas("canvas","Time & Frequency",1200,600);
    TCanvas* canvas2 = new TCanvas("canvas2","Frequency",800,1200);
    TCanvas* canvas3 = new TCanvas("canvas3","Frequency2", 800, 1200);
	canvas->Divide(2,3);
    Int_t nbins = 100;
    //TF1 *func = new TF1("func", "800*TMath::Exp(-((x-TMath::Pi())^2)/(2*0.65^2))", 0, 2*TMath::Pi());
    TF1 *func = new TF1("func", "800*(-exp(-3*x) + exp(-1.5*x))", 0, 4);
    canvas->cd(1);
    func->SetTitle("Signal");
    func->DrawClone();
    TH1F* funchisto = new TH1F("funchisto", "Signal", nbins, 0, 2*TMath::Pi());
    TH1F* funcFFT = new TH1F("funcFFT", "Signal FFT S", nbins, 0, 2*TMath::Pi());
    for (Int_t i=0; i<=nbins; i++) {
        Double_t x = Double_t(i)/nbins*(2*TMath::Pi());
        funchisto->SetBinContent(i+1, func->Eval(x));
    }
    canvas->cd(2);
    funchisto->FFT(funcFFT, "MAG");
    funcFFT->DrawClone();
    gPad->SetLogy();
    //TF1* noisy = new TF1("noisy", "800*TMath::Exp(-((x-TMath::Pi())^2)/(2*0.65^2)) + sin(6*x)/10", 0, 2*TMath::Pi());
    TF1* noisy = new TF1("noisy", "800*(-exp(-3*x) + exp(-1.5*x)) + sin(0.01*x)/10 + x/10", 0, 2*TMath::Pi());
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
    TF1 *noise = new TF1("noise", "sin(0.01*x)/10 + x/10", 0, 2*TMath::Pi());
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
    noisyFFT->SetTitle("Normalized Noisy Signal & Noise FFT");
    canvas2->cd(3); 
    noisyFFT->Draw();
    noiseFFT->Draw("same");
    gPad->SetLogy();
    TH1F* normsignalFFT = new TH1F ("normsignalFFT", "Normalized Signal FFT S/(S+N)", nbins, 0, 2*TMath::Pi());
    for (Int_t l=0; l<nbins; l++) {
         normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));
    }
    TH1F* signalFFT = new TH1F ("signalFFT", "Signal FFT S/(S+N)*(S+N)", nbins, 0, 2*TMath::Pi());
    for (Int_t q=0; q<nbins; q++) {
        signalFFT->SetBinContent(q+1, normsignalFFT->GetBinContent(q+1)*noisyFFT->GetBinContent(q+1));
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
    TH1F* fancysignal = new TH1F ("fancysignal", "Recovered Signal", nbins, 0, 2*TMath::Pi());
    for (Int_t p=0; p<nbins; p++) {
        fancysignal->SetBinContent(p+1, Signal->GetBinContent(p+1)/nbins);
    }
    func->Draw();
    fancysignal->Draw("same");
    funchisto->SetLineColor(kGreen);
    funchisto->Draw("same");
    



	

    //Double_t noisyFFTIntegral = noisyFFT->Integral();
    //Double_t noiseFFTIntegral = noiseFFT->Integral();
    //noisyFFT->Scale(1/noisyFFTIntegral);
    //noiseFFT->Scale(1/noiseFFTIntegral);
    //normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));
    //normsignalFFT->SetBinContent(l+1, pow((noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1),2));
    //signalFFT->SetBinContent(q+1, noisyFFT->GetBinContent(q+1)); //checking if output = input, without taking noise out


}