#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2EE(std::vector< std::string > files){
  TH1D * massPeak = new TH1D("massPeak","massPeak",30,60,120);

  int nEle;
  std::vector< float > * elePt = 0;
  std::vector< float > * eleEta = 0;
  std::vector< float > * elePhi = 0;


  for(unsigned int f = 0; f<files.size(); f++){
    TFile * in = TFile::Open(files.at(f).c_str(),"read");

    std::cout << f << std::endl;

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    eTree->SetBranchAddress("nEle",&nEle);
    eTree->SetBranchAddress("elePt",&elePt);
    eTree->SetBranchAddress("eleEta",&eleEta);
    eTree->SetBranchAddress("elePhi",&elePhi);

    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      eTree->GetEntry(i);
      if(nEle<2) continue;

      std::vector< int > goodElectrons;
   
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        if(elePt->at(j)<20) continue;

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
        massPeak->Fill( Zcand.M() );
      }

    }

    delete eTree;
    in->Close();
  }

  massPeak->SetDirectory(0);
  TFile * output = new TFile("output.root","recreate");
  massPeak->Write();
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
