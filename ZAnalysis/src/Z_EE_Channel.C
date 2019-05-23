#include "include/electronSelector.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2EE(std::vector< std::string > files){
  ElectronSelector eSel = ElectronSelector();

  TH1D * massPeakOS = new TH1D("massPeakOS","massPeakOS",30,60,120);
  TH1D * massPeakSS = new TH1D("massPeakSS","massPeakSS",30,60,120);

  int nEle;
  int hiBin;
  float vz;
  std::vector< float > * elePt = 0;
  std::vector< float > * eleEta = 0;
  std::vector< float > * elePhi = 0;
  std::vector< float > * eleSigmaIEtaIEta = 0;
  std::vector< int > * eleCharge = 0;
  std::vector< int > * eleMissHits = 0;
  std::vector< float > * eleTrkEta = 0;
  std::vector< float > * eleSCEta = 0;
  std::vector< float > * eleTrkPhi = 0;
  std::vector< float > * eleSCPhi = 0;
  std::vector< float > * eleHoverE = 0;
  std::vector< float > * eleD0 = 0;
  std::vector< float > * eleDz = 0;
  std::vector< float > * eleEoverPInv = 0;

  for(unsigned int f = 0; f<files.size(); f++){
    TFile * in = TFile::Open(files.at(f).c_str(),"read");

    std::cout << f << std::endl;

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    eTree->SetBranchAddress("nEle",&nEle);
    eTree->SetBranchAddress("elePt",&elePt);
    eTree->SetBranchAddress("eleEta",&eleEta);
    eTree->SetBranchAddress("elePhi",&elePhi);
    eTree->SetBranchAddress("eleSigmaIEtaIEta",&eleSigmaIEtaIEta);
    eTree->SetBranchAddress("eleMissHits",&eleMissHits);
    eTree->SetBranchAddress("eleCharge",&eleCharge);
    eTree->SetBranchAddress("eleTrkEta",&eleTrkEta);
    eTree->SetBranchAddress("eleTrkPhi",&eleTrkPhi);
    eTree->SetBranchAddress("eleSCEta",&eleSCEta);
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverE",&eleHoverE);
    eTree->SetBranchAddress("eleD0",&eleD0);
    eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("vz",&vz);

    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      eTree->GetEntry(i);
      if(nEle<2) continue;

      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;
       

      if(hiBin<100) continue;

      std::vector< int > goodElectrons;
   
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        if(elePt->at(j)<20) continue;
        if(TMath::Abs(eleSCEta->at(j)) > 2.4) continue;


        float dEta = TMath::Abs( eleTrkEta->at(j) - eleSCEta->at(j) );
        float dPhi = TMath::Abs( TMath::ACos(TMath::Cos(eleTrkPhi->at(j) - eleSCPhi->at(j))) );
        
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::veto, hiBin, eleEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverE->at(j), eleEoverPInv->at(j), eleD0->at(j), eleDz->at(j))) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      if(goodElectrons.size()>2) std::cout << "Warning: more than 2 good electrons - what do we do here???" << std::endl;
      if(goodElectrons.size()==2){
        TLorentzVector * elec1 = new TLorentzVector();
        elec1->SetPtEtaPhiM(elePt->at(goodElectrons.at(0)), eleEta->at(goodElectrons.at(0)), elePhi->at(goodElectrons.at(0)), 0.000511);

        TLorentzVector * elec2 = new TLorentzVector();
        elec2->SetPtEtaPhiM(elePt->at(goodElectrons.at(1)), eleEta->at(goodElectrons.at(1)), elePhi->at(goodElectrons.at(1)), 0.000511);

        TLorentzVector Zcand = *elec1+*elec2;
        if( eleCharge->at(goodElectrons.at(0)) != eleCharge->at(goodElectrons.at(1)) ){
          massPeakOS->Fill( Zcand.M() );
        }else{
          massPeakSS->Fill( Zcand.M() );
        }
      }

    }

    delete eTree;
    in->Close();
  }

  massPeakOS->SetDirectory(0);
  massPeakSS->SetDirectory(0);
  TFile * output = new TFile("output.root","recreate");
  massPeakOS->Write();
  massPeakSS->Write();
  output->Close();

  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: Z_EE_Channel <fileList>" << std::endl;
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
   
  doZ2EE(listOfFiles);
  return 0; 
}
