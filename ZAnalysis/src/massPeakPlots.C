#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void plotMassPeaks(std::string Zee, std::string Zmumu){
  //Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS_EE[nBins]; 
  TH1D * massPeakSS_EE[nBins]; 
  TH1D * massPeakOS_MuMu[nBins]; 
  TH1D * massPeakSS_MuMu[nBins]; 
  
  TFile * ZeeFile = TFile::Open(Zee.c_str(),"read");
  TFile * ZmumuFile = TFile::Open(Zmumu.c_str(),"read");
  
  for(int i = 0; i<nBins; i++){
    massPeakOS_EE[i] = (TH1D*)ZeeFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_EE[i] = (TH1D*)ZeeFile->Get(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakOS_MuMu[i] = (TH1D*)ZmumuFile->Get(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    massPeakSS_MuMu[i] = (TH1D*)ZmumuFile->Get(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));

    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.2);
    c1->SetBottomMargin(0.2);
    massPeakOS_MuMu[i]->GetYaxis()->SetRangeUser(0, massPeakOS_MuMu[i]->GetMaximum()*1.3);  
    massPeakOS_MuMu[i]->GetXaxis()->SetTitle("m_{ll} (GeV)");
    massPeakOS_MuMu[i]->GetYaxis()->SetTitle("Counts");
    massPeakOS_MuMu[i]->GetXaxis()->CenterTitle();
    massPeakOS_MuMu[i]->GetYaxis()->CenterTitle();

    massPeakOS_MuMu[i]->SetStats(0);
    
    massPeakOS_MuMu[i]->SetLineColor(kBlack);
    massPeakOS_MuMu[i]->SetFillColor(kOrange-2);
    massPeakOS_MuMu[i]->Draw("HIST");
    massPeakOS_EE[i]->SetLineColor(kBlack);
    massPeakOS_EE[i]->SetFillColor(kBlack);
    massPeakOS_EE[i]->SetFillStyle(3444);
    massPeakOS_EE[i]->Draw("HIST same");

    massPeakSS_MuMu[i]->SetMarkerStyle(8);
    massPeakSS_MuMu[i]->SetMarkerColor(kRed+1);

    massPeakSS_EE[i]->SetMarkerStyle(25);
    massPeakSS_EE[i]->SetMarkerColor(kBlack);
    
    massPeakSS_MuMu[i]->Draw("p same");
    massPeakSS_EE[i]->Draw("p same");

    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.png",c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.pdf",c.getCentBinLow(i),c.getCentBinHigh(i)));
    c1->SaveAs(Form("plots/massPeaks/massPeak_%d_%d.C",c.getCentBinLow(i),c.getCentBinHigh(i)));

    delete c1;
  }


  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: massPeakPlots <Z2EE file> <Z2mumu file>" << std::endl;
    return 1;
  }  

  std::string Zee = argv[1];
  std::string Zmumu = argv[2];
   
  plotMassPeaks(Zee, Zmumu);
  return 0; 
}
