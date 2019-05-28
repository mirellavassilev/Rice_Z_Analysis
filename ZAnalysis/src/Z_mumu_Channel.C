#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"

//ROOT stuff
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2mumu(std::vector< std::string > files){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  TProfile * v2Num[nBins];
  TProfile * v2NumVsCent;
  TProfile * v2Denom[nBins];
  TProfile * v2DenomVsCent;
  
  TProfile * v2MuNum[nBins];
  TProfile * v2MuNumVsCent;
  TProfile * v2MuDenom[nBins];
  TProfile * v2MuDenomVsCent;

  TH1D * v2NumVsCentHist;
  TH1D * v2SqrtDenomVsCent;
  TH1D * v2VsCent;
  
  TH1D * v2MuNumVsCentHist;
  TH1D * v2MuSqrtDenomVsCent;
  TH1D * v2MuVsCent;

  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{+}#mu^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{#mu^{#pm}#mu^{#pm}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    
    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
    v2MuNum[i] = new TProfile(Form("v2MuNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2MuDenom[i] = new TProfile(Form("v2MuDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
  }  
  v2NumVsCent = new TProfile("v2NumVsCent","",nBins,0,nBins);
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2SqrtDenomVsCent = new TH1D("v2SqrtDenomVsCent","",nBins,0,nBins);
  v2NumVsCentHist = new TH1D("v2NumVsCentHist","",nBins,0,nBins);
  
  v2MuNumVsCent = new TProfile("v2MuNumVsCent","",nBins,0,nBins);
  v2MuDenomVsCent = new TProfile("v2MuDenomVsCent","",nBins,0,nBins);
  v2MuSqrtDenomVsCent = new TH1D("v2MuSqrtDenomVsCent","",nBins,0,nBins);
  v2MuNumVsCentHist = new TH1D("v2MuNumVsCentHist","",nBins,0,nBins);

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

        //v2 calcuation
        TComplex Q = TComplex(v.ephfpQ()[1], v.ephfpAngle()[1], true) + TComplex(v.ephfmQ()[1], v.ephfmAngle()[1], true);
        TComplex candQ = TComplex(1, 2*v.phi()[j], true);
        TComplex mu1Q = TComplex(1, 2*v.PhiD1()[j], true);
        TComplex mu2Q = TComplex(1, 2*v.PhiD2()[j], true);
        if( isOppositeSign){
          for(int k = 0; k<nBins; k++){
            if(c.isInsideBin(v.centrality(),k)){
              float denom = Q.Rho2();
              float num = (candQ*TComplex::Conjugate(Q)).Re();
              v2Num[k]->Fill(0.5,num);
              v2Denom[k]->Fill(0.5,denom);

              v2NumVsCent->Fill(k,num);
              v2DenomVsCent->Fill(k,denom);

              //muons
              float numMu1 = (mu1Q*TComplex::Conjugate(Q)).Re();
              float numMu2 = (mu2Q*TComplex::Conjugate(Q)).Re();
              v2MuNum[k]->Fill(0.5,numMu1);
              v2MuNum[k]->Fill(0.5,numMu2);
              v2MuDenom[k]->Fill(0.5,denom);

              v2MuNumVsCent->Fill(k,numMu1);
              v2MuNumVsCent->Fill(k,numMu2);
              v2MuDenomVsCent->Fill(k,denom);
            }
          }
        }
      }
    }
  }

  for(int k = 1; k<nBins+1; k++){
    v2NumVsCentHist->SetBinContent(k,v2NumVsCent->GetBinContent(k));
    v2NumVsCentHist->SetBinError(k,v2NumVsCent->GetBinError(k));

    v2SqrtDenomVsCent->SetBinContent(k,TMath::Sqrt(v2DenomVsCent->GetBinContent(k)));
    v2SqrtDenomVsCent->SetBinError(k,v2DenomVsCent->GetBinError(k)/(2*TMath::Sqrt(v2DenomVsCent->GetBinContent(k))));
    
    //muons
    v2MuNumVsCentHist->SetBinContent(k,v2MuNumVsCent->GetBinContent(k));
    v2MuNumVsCentHist->SetBinError(k,v2MuNumVsCent->GetBinError(k));

    v2MuSqrtDenomVsCent->SetBinContent(k,TMath::Sqrt(v2MuDenomVsCent->GetBinContent(k)));
    v2MuSqrtDenomVsCent->SetBinError(k,v2MuDenomVsCent->GetBinError(k)/(2*TMath::Sqrt(v2MuDenomVsCent->GetBinContent(k))));
  }

  v2VsCent = (TH1D*)v2NumVsCentHist->Clone("v2VsCent");
  v2VsCent->Divide(v2SqrtDenomVsCent);
  //muons  
  v2MuVsCent = (TH1D*)v2MuNumVsCentHist->Clone("v2MuVsCent");
  v2MuVsCent->Divide(v2MuSqrtDenomVsCent);

  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->SetDirectory(0);
    massPeakSS[i]->SetDirectory(0);
    v2Num[i]->SetDirectory(0);
    v2Denom[i]->SetDirectory(0);
    v2MuNum[i]->SetDirectory(0);
    v2MuDenom[i]->SetDirectory(0);
  }

  v2DenomVsCent->SetDirectory(0);
  v2SqrtDenomVsCent->SetDirectory(0);
  v2NumVsCentHist->SetDirectory(0);
  v2NumVsCent->SetDirectory(0);
  v2VsCent->SetDirectory(0);
  
  v2MuDenomVsCent->SetDirectory(0);
  v2MuSqrtDenomVsCent->SetDirectory(0);
  v2MuNumVsCentHist->SetDirectory(0);
  v2MuNumVsCent->SetDirectory(0);
  v2MuVsCent->SetDirectory(0);

  TFile * output = new TFile("Z2mumu.root","recreate");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2MuNum[i]->Write();
    v2MuDenom[i]->Write();
  }
  v2NumVsCent->Write();
  v2DenomVsCent->Write();
  v2SqrtDenomVsCent->Write();
  v2NumVsCentHist->Write();
  v2VsCent->Write();
  v2MuNumVsCent->Write();
  v2MuDenomVsCent->Write();
  v2MuSqrtDenomVsCent->Write();
  v2MuNumVsCentHist->Write();
  v2MuVsCent->Write();

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
