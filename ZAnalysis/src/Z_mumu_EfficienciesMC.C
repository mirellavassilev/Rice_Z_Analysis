#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

bool isMatched(float eta1, float phi1, float pt1, float eta2, float phi2, float pt2){
  if( TMath::Abs(pt1 - pt2) > 5 ) return false;
  if( TMath::Abs(eta1 - eta2) > 0.1) return false;
  if( TMath::ACos(TMath::Cos(phi1 - phi2)) > 0.1 ) return false;
  return true; 
}

void doZ2mumuMC(std::vector< std::string > files){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoID_pass[nBins];
  TH2D * recoID_net[nBins];
  TH2D * recoTrig_pass[nBins];
  TH2D * recoTrig_net[nBins];
  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
    recoID_pass[k] = new TH2D(Form("recoID_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
    recoID_net[k] = new TH2D(Form("recoID_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
    recoTrig_pass[k] = new TH2D(Form("recoTrig_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
    recoTrig_net[k] = new TH2D(Form("recoTrig_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nMuEffBinsEta,-2.4,2.4,s.nMuEffBinsPt,s.muEffBinsPt);
  }

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),"dimucontana_mc"); 
    for(unsigned int i = 0; i<v.GetEntries(); i++){
    //for(unsigned int i = 0; i<100000; i++){
      v.GetEntry(i);
      
      if(i%1000==0) std::cout << i << std::endl;
      
      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      for(unsigned int j = 0; j<v.candSize_gen(); j++){

        //only look at gen Z's
        if(v.PID_gen()[j] != 23) continue;

        //require both legs to be in acceptance
        if( TMath::Abs( v.EtaD1_gen()[j] ) > 2.4 ) continue;
        if( TMath::Abs( v.EtaD2_gen()[j] ) > 2.4 ) continue;

        //require both muons to be > 20 GeV
        if( v.pTD1_gen()[j] < 20 ) continue;
        if( v.pTD2_gen()[j] < 20 ) continue;

         
        for(int k = 0; k<nBins; k++){ 
          recoEff_net[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j]);
          recoEff_net[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j]);
        }

        bool is1Matched = false;
        bool is2Matched = false;
        bool is1Tight = false;
        bool is2Tight = false;
        bool isDaughter1Trigger = false;
        bool isDaughter2Trigger = false;
        //FIXME, we are assuming we miss reconstruct both legs (which might not be accurate)!
        if(v.RecIdx_gen()[j] < 0 ) continue;
        if( v.chargeD1()[j] == v.chargeD1_gen()[ v.RecIdx_gen()[j] ] ){ //1 matched to 1 and 2 to 2
          is1Matched = isMatched(v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD1_gen()[j], v.EtaD1()[v.RecIdx_gen()[j]], v.PhiD1()[v.RecIdx_gen()[j]], v.pTD1()[v.RecIdx_gen()[j]]);
          is2Matched = isMatched(v.EtaD2_gen()[j], v.PhiD2_gen()[j], v.pTD2_gen()[j], v.EtaD2()[v.RecIdx_gen()[j]], v.PhiD2()[v.RecIdx_gen()[j]], v.pTD2()[v.RecIdx_gen()[j]]); 
        } else { //1 matched to 2 and vice versa
          is1Matched = isMatched(v.EtaD1_gen()[j], v.PhiD1_gen()[j], v.pTD1_gen()[j], v.EtaD2()[v.RecIdx_gen()[j]], v.PhiD2()[v.RecIdx_gen()[j]], v.pTD2()[v.RecIdx_gen()[j]]);
          is2Matched = isMatched(v.EtaD2_gen()[j], v.PhiD2_gen()[j], v.pTD2_gen()[j], v.EtaD1()[v.RecIdx_gen()[j]], v.PhiD1()[v.RecIdx_gen()[j]], v.pTD1()[v.RecIdx_gen()[j]]); 
        }

        for(int k = 0; k<nBins; k++){
          float SF1 = 1;
          float SF2 = 1; 
          if(is1Matched) recoEff_pass[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j], SF1);
          if(is2Matched) recoEff_pass[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j], SF2);
          
          //no SF here
          if(is1Matched) recoID_net[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j]);
          if(is2Matched) recoID_net[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j]);
        }

        if( v.chargeD1()[j] == v.chargeD1_gen()[ v.RecIdx_gen()[j] ] ){ //1 matched to 1 and 2 to 2
          is1Tight = v.tightMuon1()[v.RecIdx_gen()[j]];
          is2Tight = v.tightMuon2()[v.RecIdx_gen()[j]]; 
        } else {
          is1Tight = v.tightMuon2()[v.RecIdx_gen()[j]];
          is2Tight = v.tightMuon1()[v.RecIdx_gen()[j]];
        }
        
        for(int k = 0; k<nBins; k++){
          float SF1 = 1;
          float SF2 = 1; 
          if(is1Matched && is1Tight) recoID_pass[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j], SF1);
          if(is2Matched && is2Tight) recoID_pass[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j], SF2);
          
          if(is1Matched && is1Tight) recoTrig_net[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j]);
          if(is2Matched && is2Tight) recoTrig_net[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j]);
        }
        
        if( v.chargeD1()[j] == v.chargeD1_gen()[ v.RecIdx_gen()[j] ] ){ //1 matched to 1 and 2 to 2
          isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][v.RecIdx_gen()[j]];
          isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][v.RecIdx_gen()[j]];
        } else {
          isDaughter1Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][v.RecIdx_gen()[j]];
          isDaughter2Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][v.RecIdx_gen()[j]];
        }
        
        for(int k = 0; k<nBins; k++){
          float SF1 = 1;
          float SF2 = 1; 
          if(is1Matched && is1Tight && isDaughter1Trigger) recoTrig_pass[k]->Fill(v.EtaD1_gen()[j], v.pTD1_gen()[j], SF1);
          if(is2Matched && is2Tight && isDaughter2Trigger) recoTrig_pass[k]->Fill(v.EtaD2_gen()[j], v.pTD2_gen()[j], SF2);
        }
      }

      //check out trigger 
      //if( !(v.trigHLT()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12]) ) continue;


      /*for(unsigned int j = 0; j<v.candSize(); j++){
        
        if( !(v.pTD1()[j] > 20 )) continue;
        if( !(v.pTD2()[j] > 20 )) continue;
        if( !(v.tightCand(j,"POG"))) continue;//tight Muon 1 && tight Muon 2      
        if( !(v.VtxProb()[j] >0.001)) continue; 
 
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][j];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[j] != v.chargeD2()[j];
        if( isOppositeSign){
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin(v.centrality(),k)){
              massPeakOS[k]->Fill( v.mass()[j] );
            }
          }
        }else{
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin(v.centrality(),k)) massPeakSS[k]->Fill( v.mass()[j] );
          }
        }
      }
      */
    }
  }

  for(int i = 0; i<nBins; i++){
    recoEff_pass[i]->SetDirectory(0);
    recoEff_net[i]->SetDirectory(0);
    recoID_pass[i]->SetDirectory(0);
    recoID_net[i]->SetDirectory(0);
    recoTrig_pass[i]->SetDirectory(0);
    recoTrig_net[i]->SetDirectory(0);
  }

  TFile * output = new TFile("Z2mumu_Efficiencies.root","recreate");
  for(int i = 0; i<nBins; i++){
    recoEff_pass[i]->Write();
    recoEff_net[i]->Write();
    recoID_pass[i]->Write();
    recoID_net[i]->Write();
    recoTrig_pass[i]->Write();
    recoTrig_net[i]->Write();
  }

  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_mumu_EfficienciesMC <fileList>" << std::endl;
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
   
  doZ2mumuMC(listOfFiles);
  return 0; 
}
