#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2mumu(std::vector< std::string > files){
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  
  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
  }  


  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),""); 
    for(unsigned int i = 0; i<v.GetEntries(); i++){
      v.GetEntry(i);

      //check out trigger 
      if( !(v.trigHLT()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12]) ) continue;

      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      if(i%1000==0) std::cout << i << std::endl;

      for(unsigned int j = 0; j<v.candSize(); j++){
        
        if( !(v.pTD1()[j] > 20 )) continue;
        if( !(v.pTD2()[j] > 20 )) continue;
        if( !(v.tightCand(j,"POG"))) continue;//tight Muon 1 && tight Muon 2      
     
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[j] != v.chargeD2()[j];
        if( isOppositeSign){
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin(v.centrality(),k)) massPeakOS[k]->Fill( v.mass()[j] );
          }
        }else{
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin(v.centrality(),k)) massPeakSS[k]->Fill( v.mass()[j] );
          }
        }
      }
    }

  }

  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->SetDirectory(0);
    massPeakSS[i]->SetDirectory(0);
  }
  TFile * output = new TFile("Z2mumu.root","recreate");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
  }
  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_mumu_Channel <fileList>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  doZ2mumu(listOfFiles);
  return 0; 
}
