#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void v2Plots(std::string Zmumu){
  //Settings s = Settings();

  //CentralityTool c = CentralityTool();

  TH1D * v2; 
  TH1D * v2Mu; 
  
  TFile * ZmumuFile = TFile::Open(Zmumu.c_str(),"read");
 
  v2 = (TH1D*)ZmumuFile->Get("v2VsCent");  
  v2Mu = (TH1D*)ZmumuFile->Get("v2MuVsCent");  

  TH1D * v2Plot = new TH1D("v2Plot","",7,0,7);
  TH1D * v2MuPlot = new TH1D("v2MuPlot","",7,0,7);
  TH1D * ATLAS = new TH1D("ATLAS","",7,0,7);

  v2Plot->SetBinContent(2,v2->GetBinContent(v2->FindBin(21)));
  v2Plot->SetBinError(2,v2->GetBinError(v2->FindBin(21)));
  v2Plot->SetBinContent(3,v2->GetBinContent(v2->FindBin(11)));
  v2Plot->SetBinError(3,v2->GetBinError(v2->FindBin(11)));
  v2Plot->SetBinContent(4,v2->GetBinContent(v2->FindBin(12)));
  v2Plot->SetBinError(4,v2->GetBinError(v2->FindBin(12)));
  v2Plot->SetBinContent(5,v2->GetBinContent(v2->FindBin(13)));
  v2Plot->SetBinError(5,v2->GetBinError(v2->FindBin(13)));
  v2Plot->SetBinContent(6,v2->GetBinContent(v2->FindBin(14)));
  v2Plot->SetBinError(6,v2->GetBinError(v2->FindBin(14)));
  v2Plot->SetBinContent(7,v2->GetBinContent(v2->FindBin(20)));
  v2Plot->SetBinError(7,v2->GetBinError(v2->FindBin(20)));

  v2MuPlot->SetBinContent(2,v2Mu->GetBinContent(v2Mu->FindBin(21)));
  v2MuPlot->SetBinError(2,v2Mu->GetBinError(v2Mu->FindBin(21)));
  v2MuPlot->SetBinContent(3,v2Mu->GetBinContent(v2Mu->FindBin(11)));
  v2MuPlot->SetBinError(3,v2Mu->GetBinError(v2Mu->FindBin(11)));
  v2MuPlot->SetBinContent(4,v2Mu->GetBinContent(v2Mu->FindBin(12)));
  v2MuPlot->SetBinError(4,v2Mu->GetBinError(v2Mu->FindBin(12)));
  v2MuPlot->SetBinContent(5,v2Mu->GetBinContent(v2Mu->FindBin(13)));
  v2MuPlot->SetBinError(5,v2Mu->GetBinError(v2Mu->FindBin(13)));
  v2MuPlot->SetBinContent(6,v2Mu->GetBinContent(v2Mu->FindBin(14)));
  v2MuPlot->SetBinError(6,v2Mu->GetBinError(v2Mu->FindBin(14)));
  v2MuPlot->SetBinContent(7,v2Mu->GetBinContent(v2Mu->FindBin(20)));
  v2MuPlot->SetBinError(7,v2Mu->GetBinError(v2Mu->FindBin(20)));
   
  ATLAS->SetBinContent(1,-0.015);
  ATLAS->SetBinError(1,-0.018);

  const char * labels[7] = {"0-80%","0-80%","0-100%","0-10%", "10-30%", "30-50%", "50-100%"};

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);

  gStyle->SetErrorX(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  v2Plot->SetStats(0);
  v2Plot->GetYaxis()->SetRangeUser(-0.05,0.1);

  v2Plot->SetMarkerStyle(8);
  v2Plot->SetMarkerSize(1.3);
  v2Plot->SetMarkerColor(kBlack);
  v2Plot->SetLineColor(kBlack);
  v2Plot->SetLineWidth(2);
  v2Plot->GetXaxis()->CenterTitle();
  v2Plot->GetYaxis()->CenterTitle();

  v2MuPlot->SetMarkerStyle(21);
  v2MuPlot->SetMarkerColor(kRed+1);
  v2MuPlot->SetMarkerSize(1.5);
  v2MuPlot->SetLineColor(kRed+1);
  v2MuPlot->SetLineWidth(4);

  v2Plot->GetYaxis()->SetTitle("v_{2}");
  v2Plot->GetXaxis()->SetTitle("Centrality");

  ATLAS->SetMarkerColor(kBlue);
  ATLAS->SetLineColor(kBlue);
  ATLAS->SetMarkerStyle(21);
  ATLAS->SetLineWidth(3);
  ATLAS->SetMarkerSize(1.5);

  for(int i = 1; i<8; i++){
    v2Plot->GetXaxis()->SetBinLabel(i, labels[i-1]);
    v2Plot->GetXaxis()->ChangeLabel(i,45);
  }
 
  TLine * l = new TLine(0,0,7,0);
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);

  
  TLegend * leg = new TLegend(0.28,0.62,0.83,0.87);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0,"2018 PbPb","");
  leg->AddEntry((TObject*)0,"Scalar Product Method","");
  leg->AddEntry((TObject*)0,"p_{T}^{#mu} > 20 GeV","");
  leg->AddEntry(v2Plot,"Z candidates (#mu^{+}#mu^{-} channel)","p");
  leg->AddEntry(v2MuPlot,"Daughter muons","p");
  leg->AddEntry(ATLAS,"2012 ATLAS (#mu#mu + ee, stat error only)","p");

  v2Plot->Draw();
  l->Draw("same");
  leg->Draw("same");
  ATLAS->Draw("same");
  v2MuPlot->Draw("same");
  v2Plot->Draw("same");
 
  c1->SaveAs("plots/v2/v2Summary.png");
  c1->SaveAs("plots/v2/v2Summary.pdf");
  c1->SaveAs("plots/v2/v2Summary.C");
 

  return;
}

int main(int argc, const char* argv[])
{
  if(argc != 2)
  {
    std::cout << "Usage: massPeakPlots <Z2mumu file>" << std::endl;
    return 1;
  }  

  std::string Zmumu = argv[1];
   
  v2Plots(Zmumu);
  return 0; 
}
