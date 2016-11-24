#include "FFTAnalyzer.h"

//----------Utils-------------------------------------------------------------------------
bool FFTAnalyzer::Begin(CfgManager& opts, uint64* index)
{

    //---get all needed information from DigitizerReco
    //   n_channels is fixed by DigiReco since we want to use
    //   the same channel number <-> name mapping
    //   NB: if src is FFTAnalyzer search keep searching for the DigiReco instance    
    vector<string> srcChannels;
    int nChannels;
    float tUnit;
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        cout << ">>> FFTAnalyzer ERROR: no DigitizerReco plugin specified" << endl;
        return false;
    }    
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    if(opts.GetOpt<string>(srcInstance_+".pluginType") != "DigitizerReco")
    {
        string digiInstance = opts.GetOpt<string>(srcInstance_+".srcInstanceName");
        srcChannels = opts.GetOpt<vector<string> >(digiInstance+".channelsNames");
        nSamples_ = (opts.GetOpt<int>(digiInstance+".nSamples"));
        tUnit = (opts.GetOpt<float>(digiInstance+".tUnit"));
    }
    else
    {
        srcChannels = opts.GetOpt<vector<string> >(srcInstance_+".channelsNames");
        nSamples_ = (opts.GetOpt<int>(srcInstance_+".nSamples"));
        tUnit = (opts.GetOpt<float>(srcInstance_+".tUnit"));
    }
    nChannels = srcChannels.size();

    //---register shared FFTs
    //   nSamples is divided by to if FFT is from time to frequency domain
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    fftType_ = opts.OptExist(instanceName_+".FFTType") ?
        opts.GetOpt<string>(instanceName_+".FFTType") : "T2F";
    bool storeFFT = opts.OptExist(instanceName_+".storeFFToutput") ?
        opts.GetOpt<bool>(instanceName_+".storeFFToutput") : false;
    for(auto& channel : channelsNames_)
    {
        if(fftType_ == "T2F")
        {
            FFTs_[channel] = new FFTClass();
            RegisterSharedData(FFTs_[channel], channel, storeFFT);

            if(opts.OptExist(instanceName_+".subtractFFTNoise")){
	      noiseTemplateFile_= TFile::Open(opts.GetOpt<string>(instanceName_+".noiseTemplateFile").c_str());
	      TString noiseTemplateHistoName (opts.GetOpt<string>(instanceName_+".noiseTemplateHisto"));
	      noiseTemplateHistoRe_ = (TH1F*) noiseTemplateFile_->Get(noiseTemplateHistoName+"_Re_tmpl");
	      noiseTemplateHistoIm_ = (TH1F*) noiseTemplateFile_->Get(noiseTemplateHistoName+"_Im_tmpl");
	      if (!noiseTemplateFile_)
		{
		  cout << ">>> FFTAnalyzer ERROR: noiseTemplateFile not open " << endl;
		  return false;
		}

	    }



        }
        else
        {
            WFs_[channel] = new WFClass(1, tUnit);
            RegisterSharedData(WFs_[channel], channel, storeFFT);
            if(opts.OptExist(instanceName_+".subtractFFTNoise")){
	      noiseTemplateFile_= TFile::Open(opts.GetOpt<string>(instanceName_+".noiseTemplateFile").c_str());
	      TString noiseTemplateHistoName (opts.GetOpt<string>(instanceName_+".noiseTemplateHisto"));
	      noiseTemplateHistoRe_ = (TH1F*) noiseTemplateFile_->Get(noiseTemplateHistoName+"_Re_tmpl");
	      noiseTemplateHistoIm_ = (TH1F*) noiseTemplateFile_->Get(noiseTemplateHistoName+"_Im_tmpl");
	      if (!noiseTemplateFile_)
		{
		  cout << ">>> FFTAnalyzer ERROR: noiseTemplateFile not open " << endl;
		  return false;
		}

	    }
            if(opts.OptExist(instanceName_+".wienerFilter")){
	      signalWeinerTemplateFile_= TFile::Open(opts.GetOpt<string>(instanceName_+".signalWeinerTemplateFile").c_str());
	      TString signalWeinerTemplateHistoName (opts.GetOpt<string>(instanceName_+".signalWeinerTemplateHisto"));
	      signalWeinerTemplateHistoAmpl_ = (TH1F*) signalWeinerTemplateFile_->Get(signalWeinerTemplateHistoName+"_Ampl_tmpl");

	      bkgWeinerTemplateFile_= TFile::Open(opts.GetOpt<string>(instanceName_+".bkgWeinerTemplateFile").c_str());
	      TString bkgWeinerTemplateHistoName (opts.GetOpt<string>(instanceName_+".bkgWeinerTemplateHisto"));
	      bkgWeinerTemplateHistoAmpl_ = (TH1F*) bkgWeinerTemplateFile_->Get(bkgWeinerTemplateHistoName+"_Ampl_tmpl");


	      weightHisto_ = new TH1F("dummy","dummy",512,0,512);
	      for(int i=0;i<nSamples_/2;++i){
		float sigPlusBkg = signalWeinerTemplateHistoAmpl_->GetBinContent(i+1);
		float bkg = bkgWeinerTemplateHistoAmpl_->GetBinContent(i+1);
		float sig = sigPlusBkg - bkg;
		float weight = sig*sig/(sig*sig+bkg*bkg);
		if(sig*sig+bkg*bkg>0)		weightHisto_->SetBinContent(i,weight);
	      }

	      weightHisto_->Smooth(4);

	      TCanvas c1;
	      weightHisto_->Draw();
	      c1.SaveAs("weiner.png");


	      //background modeled with a pol2 fit to noise distribution
//	      f_bkg_ = new TF1("f_bkg_","pol2",0.,512.);
//	      f_bkg_->SetParameter(0,249.2);
//	      f_bkg_->SetParameter(1,-0.87);
//	      f_bkg_->SetParameter(2,0.001176);

	      f_filter_=new TF1("f_filter","gaus",0.,512.);//fixme move to config
	      f_filter_->SetParameter(0,1.01292);
	      f_filter_->SetParameter(1,0.0484);
	      f_filter_->SetParameter(2,28.5);


	      if (!signalWeinerTemplateFile_)
		{
		  cout << ">>> FFTAnalyzer ERROR: signalWeinerTemplateFile not open " << endl;
		  return false;
		}

	      if (!bkgWeinerTemplateFile_)
		{
		  cout << ">>> FFTAnalyzer ERROR: bkgWeinerTemplateFile not open " << endl;
		  return false;
		}

	    }


        }
    }

    //---create and register templates istograms
    //   histograms are created with automatic binning alog Y axis
    if(fftType_ == "T2F" && opts.OptExist(instanceName_+".makeTemplates")){
      templatesNames_ =  opts.GetOpt<vector<string> >(instanceName_+".makeTemplates");
    }
    for(auto& channel : channelsNames_)
        for(auto& tmpl : templatesNames_)
        {
            templates2dHistos_[channel+tmpl] = new TH2F((channel+tmpl).c_str(),
                                                      ("Template "+channel+" "+tmpl).c_str(),
                                                      nSamples_/2, 0, nSamples_/2,
                                                      10000, 0, 0);
            templatesHistos_[channel+tmpl] = new TH1F((channel+"_"+tmpl+"_tmpl").c_str(),
                                                      ("Template "+channel+" "+tmpl).c_str(),
                                                      nSamples_/2, 0, nSamples_/2);
            RegisterSharedData(templatesHistos_[channel+tmpl], channel+"_"+tmpl+"_tmpl", true);
        }
    
    //---register output data tree if requested (default true)
    bool storeTree = opts.OptExist(instanceName_+".storeTree") && fftType_ == "T2F" ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : false;
    if(storeTree)
    {
        for(auto& channel : channelsNames_)
        {
            auto point = std::find(srcChannels.begin(), srcChannels.end(), channel);
            if(point != srcChannels.end())
                channelsMap_[channel] = point-srcChannels.begin();            
        }
            
        string fftTreeName = opts.OptExist(instanceName_+".fftTreeName") ?
            opts.GetOpt<string>(instanceName_+".fftTreeName") : "fft";
        RegisterSharedData(new TTree(fftTreeName.c_str(), "fft_tree"), "fft_tree", storeTree);
        //---create tree branches:
        //   array size is determined by DigitizerReco channels
        index_ = index;
        n_tot_ = nChannels*nSamples_/2;
        current_ch_ = new int[n_tot_];
        freqs_ = new float[n_tot_];
        re_ = new float[n_tot_];
        im_ = new float[n_tot_];        
        amplitudes_ = new float[n_tot_];
        phases_ = new float[n_tot_];
        fftTree_ = (TTree*)data_.back().obj;
        fftTree_->Branch("index", index_, "index/l");
        fftTree_->Branch("n_tot", &n_tot_, "n_tot/i");
        fftTree_->Branch("ch", current_ch_, "ch[n_tot]/I");        
        fftTree_->Branch("freq", freqs_, "freq[n_tot]/F");
        fftTree_->Branch("re", re_, "re[n_tot]/F");
        fftTree_->Branch("im", im_, "im[n_tot]/F");
        fftTree_->Branch("ampl", amplitudes_, "ampl[n_tot]/F");
        fftTree_->Branch("phi", phases_, "phi[n_tot]/F");
        //---set default values
        for(int i=0; i<n_tot_; ++i)
        {
	  current_ch_[i]=-1;
            freqs_[i] = (i%nSamples_);
            re_[i] = -10;
            im_[i] = -10;
            amplitudes_[i]=-1;
            phases_[i]=-1;
        }
    }
    else
        fftTree_ = NULL;
    
    return true;
}


bool FFTAnalyzer::ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    for(auto& channel : channelsNames_)
    {
        //---FFT from time to frequency domain /// T2F
        if(fftType_ == "T2F")
        {
            //---get WF from source instance data and reset FFT
            FFTs_[channel]->Reset();
            auto wf = (WFClass*)plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false).at(0).obj;
            auto samples = wf->GetSamples();
            auto samples_norm = *samples;
            if(opts.OptExist(instanceName_+".normalizeInput") && opts.GetOpt<bool>(instanceName_+".normalizeInput"))
            {
                float max = *std::max_element(samples_norm.begin(), samples_norm.end());
                for(auto& sample : samples_norm)
                    sample /= max;
            }
	    int i=0;
            //---build the FFT
            double Re[nSamples_], Im[nSamples_];
	    
	    if(opts.OptExist(instanceName_+".BW4Filter")){//butterworth filter of 4th order
	      auto fftr2c = TVirtualFFT::FFT(1, &nSamples_, "C2C B");//the FFT is done inside the filter so this is C2C B
	      std::pair<Double_t*,Double_t*> realImP=ButterworthFilter4(samples_norm);
	      fftr2c->SetPointsComplex(realImP.first,realImP.second);
	      FFTs_[channel]->SetPointsComplex(nSamples_, realImP.first, realImP.second);
	      delete fftr2c;
	    }else{
	      auto fftr2c = TVirtualFFT::FFT(1, &nSamples_, "R2C");
	      fftr2c->SetPoints(samples_norm.data());
	      fftr2c->Transform();
	      fftr2c->GetPointsComplex(Re, Im);
	      FFTs_[channel]->SetPointsComplex(nSamples_/2, Re, Im);
	      delete fftr2c;
	    }


            map<string, const double*> var_map;
            var_map["Re"] = Re;
            var_map["Im"] = Im;
            var_map["Ampl"] = FFTs_[channel]->GetAmplitudes()->data();
            var_map["Phase"] = FFTs_[channel]->GetPhases()->data();
            if(fftTree_ || templatesNames_.size() != 0)
            {
                for(int k=0; k<nSamples_/2; ++k)
                {
		  for(auto& tmpl : templatesNames_){
                        templates2dHistos_[channel+tmpl]->Fill(k, var_map[tmpl][k]);
		  }
                    if(fftTree_)
                    {
                        int index =  channelsMap_[channel] * nSamples_/2 + k;
                        current_ch_[index] = channelsMap_[channel];
                        re_[index] = Re[k];
                        im_[index] = Im[k];
                        amplitudes_[index] = var_map["Ampl"][index];
                        phases_[index] = var_map["Phase"][index];
                    }
                }
            }

        }
        //---FFT from frequency to time domain /// F2T
        else
        {
            //---get FFT from source instance data and reset old WF
            WFs_[channel]->Reset();
            auto fft = (FFTClass*)plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false).at(0).obj;

            //---build the FFT
            double data[nSamples_];
            auto Re = fft->GetRe();
            auto Im = fft->GetIm();
            auto fftc2r = TVirtualFFT::FFT(1, &nSamples_, "C2R");



	    if(opts.OptExist(instanceName_+".subtractFFTNoise") || opts.OptExist(instanceName_+".wienerFilter")){
	    //---wiener filtering
	    //see http://www.dmf.unisalento.it/~giordano/allow_listing/wiener.pdf
            if(opts.OptExist(instanceName_+".wienerFilter")){

	      float weight=0,newRe=0,newIm=0;
	      for(int i=0;i<nSamples_/2;++i){
		//		float sigPlusBkg = signalWeinerTemplateHistoAmpl_->GetBinContent(i+1);
		//		float bkg = bkgWeinerTemplateHistoAmpl_->GetBinContent(i+1);
		//		bkgWeinerTemplateHistoAmpl_->Print();
		//		std::cout<<i<<" "<<sigPlusBkg<<" "<<bkg<<std::endl;
		//		float bkg = f_bkg_->Eval(i);
		//		float sig = sigPlusBkg - bkg;
		//		weight = sig*sig/(sig*sig+bkg*bkg);
		//		if(sig*sig+bkg*bkg>0)		dummy->SetBinContent(i,weight);
		//		weight = f_filter_->Eval(i);
		//std::cout<<i<<" "<<weight<<std::endl;

		weight = weightHisto_->GetBinContent(i+1);

		newRe = *(Re->data()+i)*weight;
		newIm = *(Im->data()+i)*weight;

		if(!opts.OptExist(instanceName_+".frequencyCut"))fftc2r->SetPoint(i,newRe,newIm);
		else if(opts.OptExist(instanceName_+".frequencyCut") && i<opts.GetOpt<float>(instanceName_+".frequencyCut"))	fftc2r->SetPoint(i,newRe,newIm);
		else if(opts.OptExist(instanceName_+".frequencyCut") && i>opts.GetOpt<float>(instanceName_+".frequencyCut"))    {fftc2r->SetPoint(i,*(Re->data()+i)*TMath::Erfc((i+1-opts.GetOpt<float>(instanceName_+".frequencyCut"))*0.1),*(Im->data()+i)*TMath::Erfc((i+1-opts.GetOpt<float>(instanceName_+".frequencyCut"))*0.1));
		  //		  std::cout<<i<<" "<<TMath::Erfc((i+1-opts.GetOpt<float>(instanceName_+".frequencyCut"))*0.1)<<std::endl;
		}
	      }

//	      TFile* f= TFile::Open("weiner_noCh18Sub.root","recreate"); 
//	      TCanvas c1;
//	      dummy->Draw();
//	      c1.SaveAs("weiner.png");
//	      dummy->Write("weight");
//	      f->Write();
//	      f->Close();


	    }

	    //---subtract FFT of noise from template before going back to time domain
            if(opts.OptExist(instanceName_+".subtractFFTNoise")){
	      float sampleShift=0;
	      if(opts.OptExist(instanceName_+".triggerRefSample")){
		int trigRef=0;
		for(int iSample=nSamples_*8; iSample<nSamples_*9; ++iSample)
		  {
		    if(event.digiSampleValue[iSample] < 1000)
		      {
			trigRef = iSample-nSamples_*8;
			break;
		      }
		  }
		
		float sampleShift=(opts.GetOpt<float>(instanceName_+".triggerRefSample")-trigRef);
		sampleShift = -sampleShift;
	      }
	      double noiseRe=0,noiseIm=0,newRe=0, newIm=0;
	      for(int i=0;i<nSamples_/2;++i){

		noiseRe = noiseTemplateHistoRe_->GetBinContent(i+1);
		noiseIm = noiseTemplateHistoIm_->GetBinContent(i+1);
		//translation in time corresponds to phase shift, check http://dsp.stackexchange.com/questions/509/what-effect-does-a-delay-in-the-time-domain-have-in-the-frequency-domain
		//		float noiseReTranslated=cos(2*3.1415*sampleShift/nSamples_)*noiseRe+sin(2*3.1415*sampleShift/nSamples_)*noiseIm;
		  //		float noiseImTranslated=cos(2*3.1415*sampleShift/nSamples_)*noiseIm-sin(2*3.1415*sampleShift/nSamples_)*noiseRe;

		//		noiseRe = noiseTemplateHistoRe_->GetBinContent(i+1);
		//		noiseIm = noiseTemplateHistoIm_->GetBinContent(i+1);
		newRe = *(Re->data()+i) - noiseRe; 
		newIm = *(Im->data()+i) - noiseIm; 

		if(!opts.OptExist(instanceName_+".frequencyCut"))fftc2r->SetPoint(i,newRe,newIm);
		else if(opts.OptExist(instanceName_+".frequencyCut") && i<opts.GetOpt<float>(instanceName_+".frequencyCut"))	fftc2r->SetPoint(i,newRe,newIm);
		else fftc2r->SetPoint(i,0,0);
		//		if(channel=="xtal11")std::cout<<i<<" "<< *(Re->data()+i)<<" "<<noiseReTranslated<<" "<<newRe<<std::endl;
	      }
	    }
	    }else{
	      fftc2r->SetPointsComplex(Re->data(), Im->data());
	    }
            fftc2r->Transform();
            fftc2r->GetPoints(data);

            //---fill new WF
            for(int iSample=0; iSample<nSamples_; ++iSample)
                WFs_[channel]->AddSample(data[iSample]/nSamples_);

            delete fftc2r;
        }
    }
    //---fill FFT tree
    if(fftTree_)
    {
        fftTree_->Fill();
    }
    
    return true;
}

bool FFTAnalyzer::End(CfgManager& opts)
{
    for(auto& channel : channelsNames_)
        for(auto& tmpl : templatesNames_)
            GetIterativeProfile(templates2dHistos_[channel+tmpl], templatesHistos_[channel+tmpl]);

    return true;
}


std::pair<Double_t*,Double_t*>  FFTAnalyzer::ButterworthFilter4(std::vector<double>& samples_norm){

	      TVirtualFFT *fft_f = TVirtualFFT::FFT(1, &nSamples_, "C2CF M K");
	      TVirtualFFT *fft_r = TVirtualFFT::FFT(1, &nSamples_, "R2C K");

	      int BW4=50;
	      float pi=3.1415;
	      float freq=5.0;
	      float dt=1./freq;
	      float df=freq/nSamples_;
	      // Butterworth filter order 4
	      double o1=0.7654, o2=1.8478;
	      // Cut-off frequency :
	      double oc=2.*pi*BW4*1.e6;
	      double a1=1+o1*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
	      double b1=2.*((oc*dt/2.)*(oc*dt/2.)-1.);
	      double c1=1-o1*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
	      double a2=1+o2*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
	      double b2=2.*((oc*dt/2.)*(oc*dt/2.)-1.);
	      double c2=1-o2*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
	      double ad=a1*a2;
	      double bd=a1*b2+b1*a2;
	      double cd=a1*c2+b1*b2+c1*a2;
	      double dd=b1*c2+c1*b2;
	      double ed=c1*c2;
	      double an=1.;
	      double bn=4.;
	      double cn=6.;
	      double dn=4.;
	      double en=1.;
	      double g=pow(dt*oc,4)/16.;
	      an=g*an/ad;
	      bn=g*bn/ad;
	      cn=g*cn/ad;
	      dn=g*dn/ad;
	      en=g*en/ad;
	      bd=bd/ad;
	      cd=cd/ad;
	      dd=dd/ad;
	      ed=ed/ad;
	      ad=1.;
	      double y4=0.,y3=0.,y2=0.,y1=0.;
	      double x4=0.,x3=0.,x2=0.,x1=0.;

	      double rex_BW4[nSamples_],rey_BW4[nSamples_],rey2_BW4[nSamples_];
	      double imx_BW4[nSamples_],imy_BW4[nSamples_],imy2_BW4[nSamples_];
	      double yref=1.;
	      //step response
	      for(int i=0; i<200; i++)
		{
		  double value=0.;
		  if(i>10)value=1.;
		  double y=x4*en+x3*dn+x2*cn*x1*bn+value*an-y4*ed-y3*dd-y2*cd-y1*bd;
		  x4=x3; x3=x2; x2=x1; x1=value;
		  y4=y3; y3=y2; y2=y1; y1=y;
		  yref=y;
		}


	      y4=0.,y3=0.,y2=0.,y1=0.;
	      x4=0.,x3=0.,x2=0.,x1=0.;


	      for(int is=0; is<nSamples_; is++)
		{
		  double value=*(samples_norm.data()+is);
		  double y=(x4*en+x3*dn+x2*cn*x1*bn+value*an)/yref-y4*ed-y3*dd-y2*cd-y1*bd; 
		  x4=x3; x3=x2; x2=x1; x1=value;
		  y4=y3; y3=y2; y2=y1; y1=y;
		  rex_BW4[is]=y;
		  imx_BW4[is]=0;
		  //		  std::cout<<value<<" "<<y<<" "<<yref<<std::endl;
		}

	      double mod_BW4[nSamples_];
	      fft_f->SetPointsComplex(rex_BW4, imx_BW4);
	      fft_f->Transform();
	      fft_f->GetPointsComplex(rey_BW4, imy_BW4);
	      //	      delete fft_f;

	      fft_r->SetPoints(samples_norm.data());
	      fft_r->Transform();
	      fft_r->GetPointsComplex(rey2_BW4, imy2_BW4);
	      delete fft_r;


	      double xax[nSamples_];
	      float normalization = 0;
	      for(int i=0; i<nSamples_; i++)  {
		xax[i]=i;
		mod_BW4[i]=1/(sqrt(rey_BW4[i]*rey_BW4[i]+imy_BW4[i]*imy_BW4[i])/yref);
		if(i==0)	normalization=mod_BW4[i];//in general the filter has some gain so normalize it
		mod_BW4[i]=mod_BW4[i]/normalization;
		rey2_BW4[i]=mod_BW4[i]*rey2_BW4[i];
		imy2_BW4[i]=mod_BW4[i]*imy2_BW4[i];
	      }
//	      TCanvas c3;
//	      TGraph dummy(512,xax, mod_BW4);
//	      dummy.Draw("AP");
//	      c3.SaveAs("dummy.png");


	      return std::make_pair<Double_t*,Double_t*>(rey2_BW4,imy2_BW4);
}
