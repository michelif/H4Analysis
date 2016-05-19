{
 //---- copy test file from:
 //       https://github.com/simonepigazzini/H4Analysis
 
 TCanvas* canvasTimeFrequency = new TCanvas("canvasTimeFrequency","canvasTimeFrequency",900,500);
//  TFile f("analysis_4443.root"); 
//  h4->Draw("WF_val:WF_time>>htemp", "WF_ch==2 && event==1 && spill==1");
//  TH1F* htemp=(TH1F*)gDirectory->Get("htemp");
 
 canvasTimeFrequency->Divide(3,1);
 
 canvasTimeFrequency->cd(1);
 TH1F* histogram_time = new TH1F ("histogram_time", "histogram in time domain", 100, 0, 100);

 //---- fill a histogram
//  TF1* funz = new TF1("funz","(x-50)*(x-50)+0.04*x",0,100);
 TF1* funz = new TF1("funz","sin(x/100*6.28)",0,100);
 for (int i = 0; i<1000; i++) {
  histogram_time->Fill(i, funz->Eval(i));
 }
 
 histogram_time->Draw("hist");
 histogram_time->GetXaxis()->SetTitle("time");
 
 canvasTimeFrequency->cd(2);
 
 TVirtualFFT::SetTransform(0);
//  TH1F *histogram_frequency_temp = new TH1F("histogram_frequency_temp", "", 100, 0, 1);
//  histogram_time->FFT(histogram_frequency_temp, "MAG");
 TH1F *histogram_frequency_temp = (TH1F*) histogram_time->FFT(0, "MAG");
 histogram_frequency_temp->Draw();
 std::cout << " n = " << histogram_frequency_temp->GetEntries() << std::endl;
 
 //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
 //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
 
 canvasTimeFrequency->cd(3);
 TH1F* histogram_frequency = new TH1F("histogram_frequency_temp", "histogram in time domain", 100, 0, 1);
 for (int iBin = 0; iBin < histogram_frequency_temp->GetNbinsX(); iBin++ ) {
  histogram_frequency -> SetBinContent (iBin+1, histogram_frequency_temp->GetBinContent(iBin+1));
 }
 histogram_frequency->Draw();
 histogram_frequency->GetXaxis()->SetTitle("frequency");
 
 
}



