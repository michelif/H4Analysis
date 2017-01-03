void fastPlotter(TString run){

  //  TString run = "5580";
  TString outDir = "plots/"+run+"/";
  system("mkdir "+outDir);
  cout<<"Processing Run:"+run<<std::endl;

  TFile *_file0 = TFile::Open("ntuples/Run_"+run+".root");
  TTree* tree = (TTree*) _file0->Get("h4");
  TFile *outfile = TFile::Open(outDir+"outFile"+run+".root","recreate");

  TCanvas*  c1 = new TCanvas();
  //define histos with correct boundaries
  TH1F* h = new TH1F("h","h",10000,0,200);
  TH1F* h2 = new TH1F("h2","h2",1000,0,4000);

  //project on the histo a given variable with a given selection
  tree->Project("h","time[Ch0]-time[Trig]","time[Ch0]>18");
  float mean =  h->GetMean();
  float  rms = h->GetRMS();
  h->GetXaxis()->SetRangeUser(mean-10*rms,mean+10*rms);
  h->Fit("gaus");
  h->Draw();
  c1->SaveAs(outDir+"TimingResolution.png");
  c1->SaveAs(outDir+"TimingResolution.pdf");
  h->Write("TimingResolution");

  tree->Project("h2","amp_max[Ch0]","","");
  h2->Draw();
  c1->SaveAs(outDir+"AmpCh0.png");
  c1->SaveAs(outDir+"AmpCh0.pdf");
  h2->Write("AmpCh0");


  outfile->Write();
  outfile->Close();
}
