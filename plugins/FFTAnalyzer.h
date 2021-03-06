#ifndef __FFT_Analyzer__
#define __FFT_Analyzer__

#include <iostream>
#include <math.h>

#include "TH2F.h"
#include "TCanvas.h"
#include "TVirtualFFT.h"
#include "TGraph.h"

#include "interface/utils.h"
#include "interface/PluginBase.h"
#include "interface/DigiTree.h"
#include "interface/WFTree.h"
#include "interface/WFClass.h"
#include "interface/FFTClass.h"

class FFTAnalyzer: public PluginBase
{
public:
    //---ctors---
    FFTAnalyzer(){};

    //---dtor---
    ~FFTAnalyzer(){};

    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    bool End(CfgManager& opts);

    std::pair<Double_t*,Double_t*> ButterworthFilter4(std::vector<double>&);
          
private:    
    //---internal data
    uint64*                   index_;
    unsigned int              n_tot_;
    int                       nSamples_;
    int*                      current_ch_;
    float*                    freqs_;
    float*                    re_;
    float*                    im_;
    float*                    amplitudes_;
    float*                    phases_;
    string                    fftType_;
    string                    srcInstance_;
    vector<string>            channelsNames_;
    vector<string>            templatesNames_;
    map<string, int>          channelsMap_;
    map<string, TH1F*>        templatesHistos_;
    map<string, TH2F*>        templates2dHistos_;
    map<string, FFTClass*>    FFTs_;
    map<string, WFClass*>     WFs_;
    TTree*                    fftTree_;         
    //---subtract FFT noise
    TFile* noiseTemplateFile_;
    TH1F* noiseTemplateHistoRe_;
    TH1F* noiseTemplateHistoIm_;
    //---wiener filter
    TFile* signalWeinerTemplateFile_;
    TH1F*  signalWeinerTemplateHistoAmpl_;
    TFile* bkgWeinerTemplateFile_;
    TH1F*  bkgWeinerTemplateHistoAmpl_;
    TF1* f_bkg_; 
    TF1* f_filter_; 
    TH1F* weightHisto_;
};

DEFINE_PLUGIN(FFTAnalyzer);

#endif
