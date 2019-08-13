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


/////////////////////////Systematics/////////////////////////////
TFile* sys=new TFile("systematics_ee_isMu210.root");
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
////////////Unfold Bin-by-Bin//////////////////
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
sysefferr->SaveAs("plots/systoterre.pdf");
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
totsysbb->SaveAs("plots/totsysbbe.pdf");
//Unfold Matrix Inversion//
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
sysefferri->SaveAs("plots/systoterrie.pdf");
//Ratio//
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
totsysi->SaveAs("plots/totsysie.pdf");
//Smeared Response Matrix
TH2D* SmearResponse=(TH2D*)f->Get("SmearResponse");
SmearResponse->Draw("COLZ");
RooUnfoldResponse response2 (0,0,SmearResponse,"","");

RooUnfoldBinByBin unfoldsbb (&response2,data);
TH1D* smearbb= (TH1D*) unfoldsbb.Hreco();
TCanvas *bbsmear=new TCanvas ("bbsmear");
TH1D*bbsmearratio=(TH1D*) smearbb->Clone("bbsmearratio");
bbsmearratio->Divide(hReco2);
RooUnfoldInvert unfoldsmi (&response2,data);
TH1D* smearmi =(TH1D*) unfoldsmi.Hreco();
TH1D*mismearratio=(TH1D*) smearmi->Clone("mismearratio");
mismearratio->Divide(hReco3);
bbsmearratio->SetTitle("Smeared Bin-by-Bin Unfolding/Bin-by-Bin Unfolding");
bbsmearratio->SetLineColor(kRed);
bbsmearratio->Draw();
mismearratio->SetLineColor(kBlue);
mismearratio->Draw("SAME");

auto legend9=new TLegend(.3,.75,.7,.9);
legend9->AddEntry(bbsmearratio,"Bin-by-Bin","l");
legend9->AddEntry(mismearratio,"Matrix Inversion","l");
legend9->Draw();
bbsmear->SaveAs("plots/smeare.pdf");
}


#ifndef __CINT__
int main()
{unfold();
  return 0;
}
#endif     
