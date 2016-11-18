{
  TString run=5880;
  TString outDir = "plots/";

  TFile *_file0 = TFile::Open("ntuples/Run_"+run+".root");
  TFile *outfile = TFile::Open(outDir+"outFile"+run+".root");

  h4->Draw("time[Ch0]-time[Trig]>>h(1000,0,200)","time[Ch0]>18","");
  h->Fit("gaus");
  h->SaveAs(outDir+"TimingResolution.png");
  h->SaveAs(outDir+"TimingResolution.pdf");


}
