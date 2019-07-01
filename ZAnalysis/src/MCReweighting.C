#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doReweighting(std::vector< std::string > filesMC, std::vector< std::string > filesData){
 TH1D * vz_Data = new TH1D("vz_Data",";v_{z};Normalized to unity",50,-25,25); 
 TH1D * vz_MC = new TH1D("vz_MC",";v_{z};Normalized to unity",50,-25,25); 
 TH1D * vz_Ratio; 

 float vz;

  int pprimaryVertexFilter;
  int phfCoincFilter2Th4;
  int pclusterCompatibilityFilter;

  std::cout << "Reading MC Files..." << std::endl;
  for(unsigned int f = 0; f<filesMC.size(); f++){
    std::cout << f << "/" << filesMC.size() << std::endl;
 
    TFile * in = TFile::Open(filesMC.at(f).c_str(),"read");
    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("vz",&vz);
  
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

    for(unsigned int i = 0; i<evtTree->GetEntries(); i++){
      skimTree->GetEntry(i);
      if( !(pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter) ) continue;

      evtTree->GetEntry(i);
      vz_MC->Fill(vz);
    }

    in->Close();
  }

  std::cout << "Reading Data Files..." << std::endl;
  for(unsigned int f = 0; f<filesData.size(); f++){
  //for(unsigned int f = 0; f<40; f++){
    std::cout << f << "/" << filesData.size() << std::endl;
 
    TFile * in = TFile::Open(filesData.at(f).c_str(),"read");
    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("vz",&vz);
  
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
      
    for(unsigned int i = 0; i<evtTree->GetEntries(); i++){
      skimTree->GetEntry(i);
      if( !(pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter) ) continue;

      evtTree->GetEntry(i);
      vz_Data->Fill(vz);
    }

    in->Close();
  }

  float MCEntries = vz_MC->GetEntries();
  float DataEntries = vz_Data->GetEntries();
  
  vz_MC->Scale(1.0/MCEntries);
  vz_Data->Scale(1.0/DataEntries);
  vz_Ratio = (TH1D*) vz_Data->Clone("vz_Ratio");
  vz_Ratio->Divide(vz_MC);

  vz_MC->Print("All");
  vz_Data->Print("All");
  vz_Ratio->Print("All");

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  TCanvas * c1 = new TCanvas("c1","",800,800);
  c1->SetLeftMargin(0.15);
  vz_MC->SetLineColor(kBlack);
  vz_MC->Draw("HIST");
  vz_Data->SetMarkerColor(kRed+1);
  vz_Data->SetLineColor(kRed+1);
  vz_Data->SetMarkerStyle(8);
  vz_Data->Draw("p same");

  TLegend * l = new TLegend(0.2,0.7,0.35,0.8);
  l->SetBorderSize(0);
  l->AddEntry(vz_MC,"MC","l");
  l->AddEntry(vz_Data,"Data","p");
  l->Draw("same");


  c1->SaveAs("plots/MCReweighting/vz.pdf");
  c1->SaveAs("plots/MCReweighting/vz.png");
  c1->SaveAs("plots/MCReweighting/vz.C");

  vz_Ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  vz_Ratio->GetYaxis()->SetTitle("Data/MC");
  vz_Ratio->SetMarkerColor(kRed+1);
  vz_Ratio->SetLineColor(kRed+1);
  vz_Ratio->SetMarkerStyle(8);
  vz_Ratio->Draw();

  c1->SaveAs("plots/MCReweighting/vz_Ratio.pdf");
  c1->SaveAs("plots/MCReweighting/vz_Ratio.png");
  c1->SaveAs("plots/MCReweighting/vz_Ratio.C");

  vz_Ratio->SetDirectory(0);
  TFile * out = TFile::Open("resources/vzReweight.root","recreate");
  vz_Ratio->Write();
  out->Close();

}


int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: MCReweighting.exe <fileListMC> <fileListData>" << std::endl;
    return 1;
  }  

  std::string fListMC = argv[1];
  std::string fListData = argv[2];
  std::string buffer, buffer2;
  std::vector<std::string> listOfFilesMC;
  std::ifstream inFileMC(fListMC.data());

  std::vector<std::string> listOfFilesData;
  std::ifstream inFileData(fListData.data());

  if(!inFileMC.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    while(true)
    {
      inFileMC >> buffer;
      if(inFileMC.eof()) break;
      listOfFilesMC.push_back(buffer);
    }
  }
  
  if(!inFileData.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    while(true)
    {
      inFileData >> buffer2;
      if(inFileData.eof()) break;
      listOfFilesData.push_back(buffer2);
    }
  }
   
  doReweighting(listOfFilesMC, listOfFilesData);
  return 0; 
}
