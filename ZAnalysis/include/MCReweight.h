#ifndef MCREWEIGHT
#define MCREWEIGHT

#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>

class MCReweight{
  public:

    MCReweight(std::string file);
    ~MCReweight();
    double reweightFactor(float vz);    

  private:
    TFile * f;
    TH1D * vzRatio;

};

MCReweight::MCReweight(std::string file){
  f = TFile::Open(file.c_str(),"read");
  vzRatio = (TH1D*)f->Get("vz_Ratio");
}

MCReweight::~MCReweight(){
  f->Close();
}

double MCReweight::reweightFactor(float vz){
  if(TMath::Abs(vz) > 20){
    std::cout << "Warning, trying to get the vz reweight factor for |vz|>20 cm.  This region is not supported so I am returning 1!!!" << std::endl;
    return 1;
  }

  return vzRatio->GetBinContent( vzRatio->FindBin(vz) );
}


#endif
