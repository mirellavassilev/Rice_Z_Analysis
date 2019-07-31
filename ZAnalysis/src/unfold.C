#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLegend.h>
#include <iostream>
using std::cout;
using std::endl;
#include <TMatrixTBase.h>
#include "TVectorT.h"
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
#include "RooUnfold/RooUnfold/src/RooUnfold.h"
#endif

void unfold()
{
#ifdef __CINT__
  gSystem->Load("include/libRooUnfold");
#endif



//Get Response Histogram
TFile* r=new TFile("backgroundSubtraction_24_isMu1.root");
TH1D *pTOS_minusAll_0_90= (TH1D*)r->Get("pTOS_minusAll_0_90");

TH1D* data=(TH1D*) pTOS_minusAll_0_90->Clone("data");
int nBins=14;
Float_t v[]={};
for( int i=1; i<nBins; i++){
v[i]=pTOS_minusAll_0_90->GetXaxis()->GetBinWidth(i);
data->SetBinContent(i,pTOS_minusAll_0_90->GetBinContent(i)*v[i]);
}


//// Get genn/rico
TFile* f=new TFile("resources/Z2mumu_Efficiencies.root");
TH1D* genn=(TH1D*)f->Get("genn");
TH1D* rico=(TH1D*)f->Get("rico");
TH2D* Response=(TH2D*)f->Get("Response");

genn->Draw("COLZ");
rico->Draw("COLZ");
Response->Draw("COLZ");
data->Draw("COLZ");

TCanvas *plots=new TCanvas("plots","");
genn->Draw();
genn->SetLineColor(kGreen);
rico->Draw("SAME");
rico->SetLineColor(kRed);
data->Draw("SAME");
plots->SaveAs("plots/plots.pdf");
//Bayes Plot///
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


bu->SaveAs("plots/bu1.pdf");
bu->SaveAs("plots/bu1.png");

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

hReco2->Print("All");
bbu->SaveAs("plots/bbu1.pdf");
bbu->SaveAs("plots/bbu1.png");

//Error Plots Bin-by Bin//
 TVectorD e0=(TVectorD) unfold2.ErecoV(RooUnfold::ErrorTreatment::kErrors);
TVectorD e1=(TVectorD) unfold2.ErecoV( RooUnfold::ErrorTreatment::kCovariance);
TVectorD e2=(TVectorD) unfold2.ErecoV(RooUnfold::ErrorTreatment::kCovToy);
//TVectorD* e3=(TVectorD*) unfold2.ErecoV( witherror=3);
TCanvas *error=new TCanvas("error","error compare");
e0.Draw();
e1.Draw("SAME");
e2.Draw("SAME");
//e3->Draw("SAME");

//auto legend6=new TLegend(.3,.75,.7,.9);
//legend6->AddEntry(e0,"e0","l");
//legend6->AddEntry(e1,"e1","l");
//legend6->AddEntry(e2,"e2","l");
//legend6->AddEntry(e3,"e3","l");
//legend6->Draw();
error->SaveAs("plots/error.pdf");

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
invert->SaveAs("plots/inversion.pdf");
invert->SaveAs("plots/inversion.png");

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
bbur->SaveAs("plots/bbur.pdf");

TCanvas *bbb=new TCanvas("bbb","Bayes v BB");
TH1D* bayesbb=(TH1D*) hReco2->Clone("bayesbb");
bayesbb->Divide(hReco);
bayesbb->SetTitle("Bin-by-Bin pT Reconstructrion/Bayes Reconstruction");
bayesbb->Draw();
bbb->SaveAs("plots/bbb.pdf");

//Plots Binbybin vs inversion
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
  bbinv->SaveAs("plots/bbinv.pdf");

 TCanvas *bbbinv=new TCanvas("bbbinv","Inversion v BB");
  TH1D* invbb=(TH1D*) hReco2->Clone("bayesbb");
  invbb->Divide(hReco3);
  invbb->SetTitle("Bin-by-Bin pT Reconstructrion/MatrixInversion Reconstruction");
  invbb->Draw();
  bbbinv->SaveAs("plots/bbbinv.pdf");


////////////////////////////////////////Systematics BInbybin//////////////////////////
TFile* sys=new TFile("systematics_mu24_isMu210.root");
TH1D *pT_totalError_0_90= (TH1D*)sys->Get("pT_totalError_0_90");
pT_totalError_0_90->Draw("COLZ");

TH1D* mod=(TH1D*) data->Clone("mod");
TH1D*dataplus=(TH1D*) data->Clone("dataplus");
TH1D*dataminus=(TH1D*) data->Clone("dataminus");

Float_t w[]={};
for( int i=1; i<nBins; i++){
w[i]=pT_totalError_0_90->GetBinContent(i);
mod->SetBinContent(i,data->GetBinContent(i)*w[i]);
dataplus->SetBinContent(i,data->GetBinContent(i)+mod->GetBinContent(i));
dataminus->SetBinContent(i,data->GetBinContent(i)-mod->GetBinContent(i));
}
//Unfold BIn-by-Bin
RooUnfoldBinByBin unfoldsys (&response,dataplus);
TH1D* dataplusu= (TH1D*) unfoldsys.Hreco();
RooUnfoldBinByBin unfoldsys2 (&response,dataminus);
TH1D* dataminusu= (TH1D*) unfoldsys2.Hreco();
TCanvas *sysefferr=new TCanvas ("sysefferr","");
dataplusu->Draw("");
dataminusu->SetLineColor(kRed);
dataminusu->Draw("SAME");
hReco2->SetLineColor(kBlack);
hReco2->Draw("SAME");
sysefferr->SaveAs("plots/systoterr.pdf");
//Ratio
TCanvas *totsysbb=new TCanvas ("totsysrat");
 TH1D* lower=(TH1D*) dataminusu->Clone("lower");
 lower->Divide(hReco2);
 TH1D* upper=(TH1D*) dataplusu->Clone("upper");
 upper->Divide(hReco2);
lower->SetTitle("Bin-by-Bin Propagated Systematic Error");
lower->GetYaxis()->SetRangeUser(0.7,1.2);
lower->SetLineColor(kBlack);
lower->Draw();
upper->SetLineColor(kRed);
upper->Draw("SAME");
  auto legend7=new TLegend(.3,.75,.7,.9);
legend7->AddEntry(upper,"upper","l");
legend7->AddEntry(lower,"lower","l");
legend7->Draw();
totsysbb->SaveAs("plots/totsysbb.pdf");
//Unfold Invert
RooUnfoldInvert unfoldsysi (&response,dataplus);
TH1D* dataplusui= (TH1D*) unfoldsysi.Hreco();
RooUnfoldBinByBin unfoldsysi2 (&response,dataminus);
TH1D* dataminusui= (TH1D*) unfoldsys2.Hreco();
TCanvas *sysefferri=new TCanvas ("sysefferri","");
dataplusui->Draw("");
dataminusui->SetLineColor(kRed);
dataminusui->Draw("SAME");
hReco3->SetLineColor(kBlack);
hReco3->Draw("SAME");
sysefferri->SaveAs("plots/systoterri.pdf");
//Ratio
TCanvas *totsysi=new TCanvas ("totsysi");
 TH1D* loweri=(TH1D*) dataminusui->Clone("loweri");
 loweri->Divide(hReco3);
 TH1D* upperi=(TH1D*) dataplusui->Clone("upperi");
 upperi->Divide(hReco3);
loweri->SetTitle("Matrix Inversion Propagated Systematic Error");
loweri->SetLineColor(kBlack);
loweri->Draw();
upperi->SetLineColor(kRed);
upperi->Draw("SAME");
  auto legend8=new TLegend(.3,.75,.7,.9);
legend8->AddEntry(upperi,"upper","l");
legend8->AddEntry(loweri,"lower","l");
legend8->Draw();
totsysi->SaveAs("plots/totsysi.pdf");

}

#ifndef __CINT__
int main()
{unfold();
  return 0;
}
#endif
