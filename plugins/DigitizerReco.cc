#include "DigitizerReco.h"

//----------Utils-------------------------------------------------------------------------
bool DigitizerReco::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    nSamples_ = opts.GetOpt<int>(instanceName_+".nSamples");
    tUnit_ = opts.GetOpt<float>(instanceName_+".tUnit");

    for(auto& channel : channelsNames_)
    {
        if(opts.OptExist(channel+".type") && opts.GetOpt<string>(channel+".type") == "NINO")
            WFs[channel] = new WFClassNINO(opts.GetOpt<int>(channel+".polarity"), tUnit_);
        else
            WFs[channel] = new WFClass(opts.GetOpt<int>(channel+".polarity"), tUnit_);
        RegisterSharedData(WFs[channel], channel, false);
    }
    
    return true;
}

bool DigitizerReco::ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{        
    //---read the digitizer
    //---set time reference from digitized trigger
    int trigRef=0;
    for(int iSample=nSamples_*8; iSample<nSamples_*9; ++iSample)
    {
        if(event.digiSampleValue[iSample] < 1000)
        {
            trigRef = iSample-nSamples_*8;
            break;
        }
    }
    
    //---user channels
    for(auto& channel : channelsNames_)
      {
        //---reset and read new WFs
        WFs[channel]->Reset();
        int digiGr = opts.GetOpt<int>(channel+".digiGroup");
        int digiCh = opts.GetOpt<int>(channel+".digiChannel");
        int offset = event.digiMap.at(make_pair(digiGr, digiCh));
        for(int iSample=offset; iSample<offset+nSamples_; ++iSample)
	  {
            //---H4DAQ bug: sometimes ADC value is out of bound.
            //---skip everything if one channel is bad
	    if(event.digiSampleValue[iSample] > 10000)
	      {
		cout << ">>>DigiReco WARNING: skipped event" <<endl;
		return false;
	      }
            WFs[channel]->AddSample(event.digiSampleValue[iSample]);
	  }
        if(opts.OptExist(channel+".useTrigRef") && opts.GetOpt<bool>(channel+".useTrigRef"))
	  WFs[channel]->SetTrigRef(trigRef);
	
      }
    return true;
}
