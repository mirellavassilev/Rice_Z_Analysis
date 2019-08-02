#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLegend.h>
#include <iostream>
using std::cout;
using std::endl;

#include <TCanvas.h>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include "RooUnfold/RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfold/src/RooUnfoldBinByBin.h"
#include "RooUnfold/RooUnfold/src/RooUnfoldInvert.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#endif

void unfold()
{
#ifdef __CINT__
 gSystem->Load("include/libRooUnfold");
 #endif


//get gen/reco/Response
TFile * f=new TFile("resources/Z2ee_EfficiencyMC_0.root");
TH1D* genn=(TH1D*)f->Get("genn");
TH1D* rico=(TH1D*)f->Get("rico");
TH2D* Response=(TH2D*)f->Get("Response");

genn->Draw("COLZ");
rico->Draw("COLZ");
Response->Draw("COLZ");

//get measured data
TFile* r=new TFile("backgroundSubtraction_ee_isMu0.root");
TH1D *pTOS_minusAll_0_90= (TH1D*)r->Get("pTOS_minusAll_0_90");
TH1D* data=(TH1D*) pTOS_minusAll_0_90->Clone("data");
int nBins=14;
Float_t v[]={};
for( int i=1; i<nBins; i++){
v[i]=pTOS_minusAll_0_90->GetXaxis()->GetBinWidth(i);
data->SetBinContent(i,pTOS_minusAll_0_90->GetBinContent(i)*v[i]);
}
data->Draw("COLZ");  

//Sow pT ditsribution genn v. rico v. measured data 
TCanvas *plots=new TCanvas("plots","");
genn->Draw();
genn->SetLineColor(kGreen);
rico->Draw("SAME");
rico->SetLineColor(kRed);
data->Draw("SAME");
plots->SaveAs("plots/plotse.pdf");

//Bayes Unfolding 
RooUnfoldResponse response (0,0,Response,"","");
RooUnfoldBayes    unfold (&response,data,4);
TH1D* hReco= (TH1D*) unfold.Hreco();
 
 TCanvas *bu=new TCanvas("bu","Bayes Unfold");
 hReco->SetTitle("pT Spectrum Bayes Reconstructed");
  unfold.PrintTable (cout, genn);
  hReco->Draw();
  hReco->SetLineColor(kRed);
  data->Draw("SAME");
  genn->SetLineColor(8);
  genn->Draw("SAME");
  rico->Draw("SAME");
  rico->SetLineColor(kBlack);
auto legend=new TLegend(.3,.75,.7,.9);
legend->AddEntry(hReco,"Bayes Reconstructed Measurement","l");
legend->AddEntry(data,"Measurement","l");
legend->AddEntry(rico,"reco(MC)","l");
legend->AddEntry(genn,"gen(MC)","l");
legend->Draw();


bu->SaveAs("plots/bu1e.pdf");
bu->SaveAs("plots/bu1e.png");


//Bin by Bin Plot//
RooUnfoldBinByBin unfold2 (&response,data);
TH1D* hReco2= (TH1D*) unfold2.Hreco();
TCanvas *bbu=new TCanvas("bbu","BinbyBin Unfold");
hReco2->SetTitle("pT Spectrum Bin-by-Bin Reconstructed");
unfold2.PrintTable (cout, genn);
hReco2->Draw();
hReco2->SetLineColor(kRed);
data->Draw("SAME");
 genn->SetLineColor(8);
 genn->Draw("SAME");
 rico->SetLineColor(kBlack);
 rico->Draw("SAME");
 auto legend2=new TLegend(.3,.75,.7,.9);
 legend2->AddEntry(hReco2,"Bin-by-Bin Reconstructed Measurement","l");
 legend2->AddEntry(data,"Measurement","l");
 legend2->AddEntry(rico,"reco(MC)","l");
 legend2->AddEntry(genn,"gen(MC)","l");
 legend2->Draw();

 bbu->SaveAs("plots/bbu1e.pdf");
 bbu->SaveAs("plots/bbu1e.png");
//INvert 
RooUnfoldInvert unfold3 (&response,data);
TH1D* hReco3 =(TH1D*) unfold3.Hreco();
 TCanvas * invert=new TCanvas("invert","Invert Unfold");
hReco3->SetTitle("pT Spectrum Matrix Inversion Reconstructed");
  unfold3.PrintTable (cout, genn);
  hReco3->Draw();
  hReco3->SetLineColor(kRed);
  data->Draw("SAME");
  genn->SetLineColor(8);
  genn->Draw("SAME");
  rico->SetLineColor(kBlack);
 rico->Draw("SAME");
auto legend4=new TLegend(.3,.75,.7,.9);
legend4->AddEntry(hReco2,"Matrix Inversion Reconstructed Measurement","l");
legend4->AddEntry(data,"Measurement","l");
legend4->AddEntry(rico,"reco(MC)","l");
legend4->AddEntry(genn,"gen(MC)","l");
legend4->Draw();
invert->SaveAs("plots/inversione.pdf");
invert->SaveAs("plots/inversione.png");


//Plots Bayes
 TCanvas *bbur=new TCanvas("bbur","Bayes Unfold");

 TH1D* recratiob=(TH1D*) hReco->Clone("recratiob");
 recratiob->Divide(data);
 recratiob->SetLineColor(kBlue);
 recratiob->SetTitle("Bayes v. Bin-by-Bin Reconstruction");
 recratiob->Draw();
 TH1D* recratiobb=(TH1D*) hReco2->Clone("recratiobb");
 recratiobb->Divide(data);
 recratiobb->SetLineColor(kGreen);
 recratiobb->Draw("SAME");
 TH1D* ratiob=(TH1D*) genn->Clone("ratiob");
 ratiob->Divide(rico);
 ratiob->SetLineColor(kRed);
 ratiob->Draw("SAME");
 auto legend3=new TLegend(.3,.75,.7,.9);
 legend3->AddEntry(recratiob,"Bayes Reconstructed Data/Data","l");
 legend3->AddEntry(recratiobb,"Bin by Bin Reconstructed Data/Data","l");
 legend3->AddEntry(ratiob,"gen/reco","l");
 legend3->Draw();
 bbur->SaveAs("plots/bbure.pdf");

TCanvas *bbb=new TCanvas("bbb","Bayes v BB");
TH1D* bayesbb=(TH1D*) hReco2->Clone("bayesbb");
bayesbb->Divide(hReco);
bayesbb->SetTitle("Bin-by-Bin pT Reconstructrion/Bayes Reconstruction");
bayesbb->Draw();
bbb->SaveAs("plots/bbbe.pdf");

//PLot binbybin vs inversion 

 TCanvas *bbinv=new TCanvas("bbinv","Bayes Unfold");

 TH1D* invrat=(TH1D*) hReco3->Clone("recratiob");
 invrat->Divide(data);
 invrat->SetLineColor(kBlue);
 invrat->SetTitle("MatrixInversion  v. Bin-by-Bin Reconstruction");
 invrat->Draw();
  recratiobb->Draw("SAME");
  ratiob->Draw("SAME");
  auto legend5=new TLegend(.3,.75,.7,.9);
  legend5->AddEntry(invrat,"MatrixInversion Reconstructed Data/Data","l");
  legend5->AddEntry(recratiobb,"Bin by Bin Reconstructed Data/Data","l");
  legend5->AddEntry(ratiob,"gen/reco","l");
  legend5->Draw();
  bbinv->SaveAs("plots/bbinve.pdf");

 TCanvas *bbbinv=new TCanvas("bbbinv","Inversion v BB");
  TH1D* invbb=(TH1D*) hReco2->Clone("bayesbb");
  invbb->Divide(hReco3);
  invbb->SetTitle("Bin-by-Bin pT Reconstructrion/MatrixInversion Reconstruction");
  invbb->Draw();
  bbbinv->SaveAs("plots/bbbinve.pdf");

}


#ifndef __CINT__
int main()
{unfold();
  return 0;
}
#endif     
