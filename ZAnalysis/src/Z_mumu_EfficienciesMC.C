#include "include/VertexCompositeNtuple.h"
#include "include/PbPb_5TeV_2018_eventUtils.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "include/MCReweight.h"

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
#include "TRandom3.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitter.h"
#include "TStyle.h"
#include "TUnfold.h"
//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


void doZ2mumuMC(std::vector< std::string > files){
  TH1::SetDefaultSumw2();
  Settings s = Settings();

  MCReweight vzRW = MCReweight("resources/vzReweight.root");

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TRandom3 r = TRandom3();
  TRandom3 * r1 = new TRandom3();
////////////////////////////////////////Create Histograms//////////////////////////////////////////////////
  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoEff[nBins];
  TEfficiency * eff[nBins];

  
//Resolution
  TH1D *res= new TH1D("res","Resolution",30,-.2,.2);
  TH1D *yreso= new TH1D("yreso","Rapidity Resolution",35,-.2,.2);
   

  TH1D *rc1= new TH1D("rc1","Centrality 0-20",30,-.2,.2);
  TF1 *g1  = new TF1("g1","gaus",-.05,.05);
  TH1D *rc2= new TH1D("rc2","Centrality 20-60",30,-.2,.2);
  TF1 *g2  = new TF1("g2","gaus",-.05,.05);
  TH1D *rc3= new TH1D("rc3","Centrality 60-100",30,-.2,.2);
  TF1 *g3  = new TF1("g3","gaus",-.05,.05);
  TH1D *rc4= new TH1D("rc4","Centrality 100-200",30,-.2,.2);
  TF1 *g4  = new TF1("g4","gaus",-.05,.05);
 
  TH1D *ry1= new TH1D("ry1","y -2.4/-2",30,-.2,.2);
  TF1 *y1  = new TF1("y1","gaus",-.05,.05);
  TH1D *ry2= new TH1D("ry2","y -2/-1.6",30,-.2,.2); 
  TF1 *y2  = new TF1("y2","gaus",-.05,.05);
  TH1D *ry3= new TH1D("ry3","y -1.6/-1.2",30,-.2,.2);
  TF1 *y3  = new TF1("y3","gaus",-.05,.05);
  TH1D *ry4= new TH1D("ry4","y -1.2/-.8",30,-.2,.2);
  TF1 *y4  = new TF1("y4","gaus",-.05,.05);
  TH1D *ry5= new TH1D("ry5","y -.8/-.4",30,-.2,.2);
  TF1 *y5  = new TF1("y5","gaus",-.05,.05);
  TH1D *ry6= new TH1D("ry6","y -.4/0",30,-.2,.2);
  TF1 *y6  = new TF1("y6","gaus",-.05,.05);
  TH1D *ry7= new TH1D("ry7","y 0/.4",30,-.2,.2);
  TF1 *y7  = new TF1("y7","gaus",-.05,.05);
  TH1D *ry8= new TH1D("ry8","y .4/.8",30,-.2,.2);
  TF1 *y8  = new TF1("y8","gaus",-.05,.05);
  TH1D *ry9= new TH1D("ry9","y .8/1.2",30,-.2,.2);
  TF1 *y9  = new TF1("y9","gaus",-.05,.05); 
  TH1D *ry10= new TH1D("ry10","y 1.2/1.6",30,-.2,.2);
  TF1 *y10  = new TF1("y10","gaus",-.05,.05);
  TH1D *ry11= new TH1D("ry11","y 1.6/2",30,-.2,.2);
  TF1 *y11  = new TF1("y11","gaus",-.05,.05);
  TH1D *ry12= new TH1D("ry12","y 2/2.4",30,-.2,.2);
  TF1 *y12  = new TF1("y12","gaus",-.05,.05);
//Rapidity 
  TH1D *rc1y= new TH1D("rc1y","Centrality 0-20",30,-.2,.2);
  TF1 *g1y  = new TF1("g1y","gaus",-.1,.1);
  TH1D *rc2y= new TH1D("rc2y","Centrality 20-60",30,-.2,.2);
  TF1 *g2y  = new TF1("g2y","gaus",-.1,.1);
  TH1D *rc3y= new TH1D("rc3y","Centrality 60-100",30,-.2,.2);
  TF1 *g3y  = new TF1("g3y","gaus",-.1,.1);
  TH1D *rc4y= new TH1D("rc4y","Centrality 100-200",30,-.2,.2);
  TF1 *g4y  = new TF1("g4y","gaus",-.1,.1);


// pT 
Float_t Bins[]={0,1.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,90.0,120.0,150.0,200.0};
Int_t  nptBins = 13;
TH1D *genn=new TH1D("genn","",nptBins,Bins);
TH1D *rico=new TH1D("rico","",nptBins,Bins);
TH1D *genny=new TH1D("genny","",13,-2.1,2.1);
TH1D *ricoy=new TH1D("ricoy","",13,-2.1,2.1);

Float_t centBins[]={0,10,20,40,60,80,100,140,180};
TH1D *acccorr= new TH1D("acccorr","",8,centBins);
TH1D *acccorr2= new TH1D("acccorr2","",8,centBins);
//Respose Matrix
TH2D *Response= new TH2D("Response","",nptBins,Bins,nptBins,Bins);
TH2D * SmearResponse=new TH2D ("SmearResponse","",nptBins,Bins,nptBins,Bins);
  TH2D *yResponse= new TH2D ("Rapidity Response Matrix","",13,-2.1,2.1,13,-2.1,2.1);
//////////////////////////////////////////Select Events////////////////////////////////////////////
    for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap,s.nZPtBins-1,s.zPtBins);

   }

  //starting looping over the file
  for(unsigned int f = 0; f<files.size(); f++){
    VertexCompositeNtuple v = VertexCompositeNtuple();
    v.GetTree(files.at(f),"dimucontana_mc"); 
    //for(unsigned int i = 0; i<v.GetEntries(); i++){
    for(unsigned int i = 0; i<20000; i++){
      v.GetEntry(i);
      
      if(i%1000==0) std::cout << i << std::endl;
      
      //event selection
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::hfCoincFilter2Th4 ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::primaryVertexFilter ])) continue;
      if( !(v.evtSel()[ PbPb::R5TeV::Y2018::clusterCompatibilityFilter ])) continue;
      if( TMath::Abs(v.bestvtxZ()) > 15 ) continue;

      double eventWeight = vzRW.reweightFactor( v.bestvtxZ() ) * c.findNcoll( v.centrality() );

      for(unsigned int j = 0; j<v.candSize_gen(); j++){

        //only look at gen Z's
        if(v.PID_gen()[j] != 23) continue;
 
//Acceptance Correction//
	if( TMath::Abs( v.y_gen()[j])<2.4 && v.mass()[j]>20 && v.mass()[j]<120){
        acccorr->Fill(v.centrality());
	}
       if( TMath::Abs( v.y_gen()[j])<2.4 && v.mass()[j]>20 && v.mass()[j]<120 && v.pTD1_gen()[j]>20 && v.pTD2_gen()[j]>20 &&TMath::Abs( v.EtaD1_gen()[j] ) < 2.4 && TMath::Abs( v.EtaD2_gen()[j] ) < 2.4){
        acccorr2->Fill(v.centrality());
        }
        
        //require both legs to be in acceptance
        if( TMath::Abs( v.EtaD1_gen()[j] ) > 2.4 ) continue;
        if( TMath::Abs( v.EtaD2_gen()[j] ) > 2.4 ) continue;

        //require both muons to be > 20 GeV
        if( v.pTD1_gen()[j] < 20 ) continue;
        if( v.pTD2_gen()[j] < 20 ) continue;

        //Fill denominator 
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin(v.centrality(),k)){
            recoEff_net[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight );
          }
        }

        //see if we have a matched candidate
        if(v.RecIdx_gen()[j] < 0 ) continue;
        int indx = v.RecIdx_gen()[j];

        //see if it passes all our selections
        if( v.mass()[indx] < s.zMassRange[0] || v.mass()[indx] > s.zMassRange[1]) continue; 
        if( !(v.pTD1()[indx] > s.minMuonPt )) continue;
        if( !(v.pTD2()[indx] > s.minMuonPt )) continue;
        if( !(v.tightCand(indx,"POG"))) continue;//tight Muon 1 && tight Muon 2      
        if( !(v.VtxProb()[indx] >0.001)) continue; 
 
        //make sure that one of the daughters was the trigger muon
        bool isDaughter1Trigger = v.trigMuon1()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        bool isDaughter2Trigger = v.trigMuon2()[PbPb::R5TeV::Y2018::HLT_HIL3Mu12][indx];
        if( !(isDaughter1Trigger || isDaughter2Trigger) ) continue;
 
        bool isOppositeSign =  v.chargeD1()[indx] != v.chargeD2()[indx];
        if( !isOppositeSign ) continue;

        //Looks like this is a good candidate match! Let's get the scale factor
        double scaleFactor = 1.0; //FIXME getScaleFactor();

        //Fill numerator (and apply the scale factor here!)
        for(int k = 0; k<nBins; k++){ 
          if(c.isInsideBin(v.centrality(),k)){
            recoEff_pass[k]->Fill( v.y_gen()[j], v.pT_gen()[j], eventWeight * scaleFactor );
          }
        }
	
	//Fill Response Histogram
	Response->Fill(v.pT()[v.RecIdx_gen()[j]],v.pT_gen()[j]);
	SmearResponse->Fill(v.pT()[v.RecIdx_gen()[j]] * r1->Gaus(1,0.03),v.pT_gen()[j]);
          yResponse->Fill(v.y()[v.RecIdx_gen()[j]],v.y_gen()[j]);
	//Fill gen and rico
	genn->Fill(v.pT_gen()[j]);
	rico->Fill(v.pT()[v.RecIdx_gen()[j]]);	   
        genny->Fill(v.y_gen()[j]);
        ricoy->Fill(v.y()[v.RecIdx_gen()[j]]);   

   //Fill Resolution
       res->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
       yreso->Fill(v.y_gen()[j]-v.y()[v.RecIdx_gen()[j]]);
    
//Fill Resolution by Centrality and Rapidity pT
	if (v.centrality()>0&&v.centrality()<20) { 
	rc1->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);	
	rc1->SetLineColor(kGreen);
	}
	else if (v.centrality()>20&&v.centrality()<60) {
        rc2->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);  
	rc2->SetLineColor(kBlue);
	}
	else if (v.centrality()>60&&v.centrality()<100) {
        rc3->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        rc3->SetLineColor(kRed);
	}
	else if (v.centrality()>100&&v.centrality()<200) {
        rc4->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        rc4->SetLineColor(kYellow);
	}



	if(v.y_gen()[j]>-2.4&&v.y_gen()[j]<-2.0){
	ry1->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
	}
	else if(v.y_gen()[j]>-2.0&&v.y_gen()[j]<-1.6){
        ry2->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }	
	else if(v.y_gen()[j]>-1.6&&v.y_gen()[j]<-1.2){
        ry3->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>-1.2&&v.y_gen()[j]<-.8){
        ry4->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>-.8&&v.y_gen()[j]<-.4){
        ry5->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>-.4&&v.y_gen()[j]<0.0){
        ry6->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>0.0&&v.y_gen()[j]<.4){
        ry7->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>.4&&v.y_gen()[j]<.8){
        ry8->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>.8&&v.y_gen()[j]<1.2){
        ry9->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>1.2&&v.y_gen()[j]<1.6){
        ry10->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>1.6&&v.y_gen()[j]<2.0){
        ry11->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
        }
	else if(v.y_gen()[j]>2.0&&v.y_gen()[j]<2.4){
        ry12->Fill((v.pT_gen()[j]-v.pT()[v.RecIdx_gen()[j]])/v.pT_gen()[j]);
//Fill Resolution by centrality Rapidity
 
          if (v.centrality()>0&&v.centrality()<20) {
        rc1y->Fill(v.y_gen()[j]-v.y()[v.RecIdx_gen()[j]]);
        rc1y->SetLineColor(kGreen);
        }
        else if (v.centrality()>20&&v.centrality()<60) {
        rc2y->Fill(v.y_gen()[j]-v.y()[v.RecIdx_gen()[j]]);
        rc2y->SetLineColor(kBlue);
        }
        else if (v.centrality()>60&&v.centrality()<100) {
        rc3y->Fill(v.y_gen()[j]-v.y()[v.RecIdx_gen()[j]]);
        rc3y->SetLineColor(kRed);
        }
        else if (v.centrality()>100&&v.centrality()<200) {
        rc4y->Fill(v.y_gen()[j]-v.y()[v.RecIdx_gen()[j]]);
        rc4y->SetLineColor(kYellow);
        }   
     }

      }
    }
  }
rc1->Fit(g1,"","",-.05,.05);
rc2->Fit(g2,"","",-.05,.05);
rc3->Fit(g3,"","",-.05,.05);
rc4->Fit(g4,"","",-.05,.05);

rc1y->Fit(g1y,"","",-.05,.05);
rc2y->Fit(g2y,"","",-.05,.05);
rc3y->Fit(g3y,"","",-.05,.05);
rc4y->Fit(g4y,"","",-.05,.05);

ry1->Fit(y1,"","",-.05,.05);
ry2->Fit(y2,"","",-.05,.05);
ry3->Fit(y3,"","",-.05,.05);
ry4->Fit(y4,"","",-.05,.05);
ry5->Fit(y5,"","",-.05,.05);
ry6->Fit(y6,"","",-.05,.05);
ry7->Fit(y7,"","",-.05,.05);
ry8->Fit(y8,"","",-.05,.05);
ry9->Fit(y9,"","",-.05,.05);
ry10->Fit(y10,"","",-.05,.05);
ry11->Fit(y11,"","",-.05,.05);
ry12->Fit(y12,"","",-.05,.05);


//Histogram sig v y bin entries/error 
Float_t y1b =y1 ->GetParameter(2);
Float_t y2b =y2->GetParameter(2);
Float_t y3b =y3->GetParameter(2);
Float_t y4b =y4->GetParameter(2);
Float_t y5b =y5->GetParameter(2);
Float_t y6b =y6->GetParameter(2);
Float_t y7b =y7->GetParameter(2);
Float_t y8b =y8->GetParameter(2);
Float_t y9b =y9->GetParameter(2);
Float_t y10b =y10->GetParameter(2);
Float_t y11b =y11->GetParameter(2);
Float_t y12b =y12->GetParameter(2);
Float_t y1e=y1->GetParError(2);
Float_t y2e=y2->GetParError(2); 
Float_t y3e=y3->GetParError(2);
Float_t y4e=y4->GetParError(2);
Float_t y5e=y5->GetParError(2);
Float_t y6e=y6->GetParError(2);
Float_t y7e=y7->GetParError(2);
Float_t y8e=y8->GetParError(2);
Float_t y9e=y9->GetParError(2);
Float_t y10e=y10->GetParError(2);
Float_t y11e=y11->GetParError(2);
Float_t y12e=y12->GetParError(2);

//histogram sig v. y 
Float_t biny[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0,2.4};
Int_t binnumy=12;
TH1D*sigy=new TH1D("sigy","Sigma vs. y",binnumy,biny);
sigy->SetBinContent(1,y1b);
sigy->SetBinContent(2,y2b);
sigy->SetBinContent(3,y3b);
sigy->SetBinContent(4,y4b);
sigy->SetBinContent(5,y5b);
sigy->SetBinContent(6,y6b);
sigy->SetBinContent(7,y7b);
sigy->SetBinContent(8,y8b);
sigy->SetBinContent(9,y9b);
sigy->SetBinContent(10,y10b);
sigy->SetBinContent(11,y11b);
sigy->SetBinContent(12,y12b);
sigy->SetBinError(1,y1e);
sigy->SetBinError(2,y2e);
sigy->SetBinError(3,y3e);
sigy->SetBinError(4,y4e);
sigy->SetBinError(5,y5e);
sigy->SetBinError(6,y6e);
sigy->SetBinError(7,y7e);
sigy->SetBinError(8,y8e);
sigy->SetBinError(9,y9e);
sigy->SetBinError(10,y10e);
sigy->SetBinError(11,y11e);
sigy->SetBinError(12,y12e);
sigy->SetYTitle("sigma");
sigy->SetXTitle("y");

TCanvas *sigyh = new TCanvas("sigyh","Sigma vs. y Resolution");
sigy->Draw();
 sigyh->SaveAs("plots/sigyh.pdf");
 sigyh->SaveAs("plots/sigyh.png");


//Graph sig v. centrality

Float_t c1 =g1->GetParameter(2);
Float_t c1e=g1->GetParError(2);
Float_t c2 =g2->GetParameter(2);
Float_t c2e=g2->GetParError(2);
Float_t c3 =g3->GetParameter(2);
Float_t c3e=g3->GetParError(2);
Float_t c4 =g4->GetParameter(2);
Float_t c4e=g4->GetParError(2);

//histogram sig v. centrality 
Float_t bins[]={0,20,60,100,200};
Int_t binnum=4;
TH1D*sigc= new TH1D("sigc","Sigma v. Centrality",binnum,bins);
sigc->SetBinContent(1,c1);
sigc->SetBinContent(2,c2);
sigc->SetBinContent(3,c3);
sigc->SetBinContent(4,c4);
sigc->SetBinError(1,c1e);
sigc->SetBinError(2,c2e);
sigc->SetBinError(3,c3e);
sigc->SetBinError(4,c4e);

sigy->SetYTitle("sigma");
sigy->SetXTitle("Centrality");

TCanvas *sigcen = new TCanvas("sigcen","Sigma vs. Centrality Resolution");
sigc->Draw();
 sigcen->SaveAs("plots/sigcenhis.pdf");
 sigcen->SaveAs("plots/sigcenhis.png");

TCanvas *acceptancecent=new TCanvas("acceptancecent", "Acceptance Correction");
TH1D*acccent=(TH1D*) acccorr->Clone("acccent");
acccent->Divide(acccorr2);
acccent->Draw();
acceptancecent->SaveAs("plots/acceptancecent.pdf");


    for(int i = 0; i<nBins; i++){
    recoEff[i] = (TH2D*)recoEff_pass[i]->Clone(Form("recoEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
     recoEff[i]->Divide(recoEff_net[i]);
     recoEff[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pass[i]), *(recoEff_net[i]),"w") ){
      eff[i] = new TEfficiency(*(recoEff_pass[i]), *(recoEff_net[i]));
      eff[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff[i]->SetName(Form("eff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff[i]->SetDirectory(0);
    }
    recoEff_pass[i]->SetDirectory(0);
    recoEff_net[i]->SetDirectory(0);
}
    res->SetDirectory(0);
    yreso->SetDirectory(0);
    rc1->SetDirectory(0);
    rc2->SetDirectory(0);
    rc3->SetDirectory(0);
    rc4->SetDirectory(0);
    ry1->SetDirectory(0);
    ry2->SetDirectory(0);
    ry3->SetDirectory(0);
    ry4->SetDirectory(0);
    ry5->SetDirectory(0);
    ry6->SetDirectory(0);
    ry7->SetDirectory(0);
    ry8->SetDirectory(0);
    ry8->SetDirectory(0);
    ry9->SetDirectory(0);
    ry10->SetDirectory(0);
    ry11->SetDirectory(0);
    ry12->SetDirectory(0);
    Response->SetDirectory(0);
    genn->SetDirectory(0);
    rico->SetDirectory(0);
    SmearResponse->SetDirectory(0); 

THStack *hs = new THStack("hs","Resolution");
hs->Add(res);
hs->Add(rc1);
hs->Add(rc2);
hs->Add(rc3);
hs->Add(rc4);

TCanvas *Yrat=new TCanvas("Yrat","");
TH1D *yrat=(TH1D*) genny->Clone("yrat");
yrat->Divide(ricoy);
yrat->GetXaxis()->SetTitle("y");
yrat->GetYaxis()->SetTitle("y_gen/y");
yrat->Draw();
Yrat->SaveAs ("plots/Yrat.pdf");

  TCanvas *ccent=new TCanvas("ccent","Resolution centrality dependance");
  ccent->Divide(2,2);
  ccent->cd(1);
 rc1y->GetXaxis()->SetTitle("y_gen-y");
 rc1y->Draw();
  g1y->Draw("same");
  ccent->cd(2);
rc2y->GetXaxis()->SetTitle("y_gen-y");
  rc2y->Draw();
  g2y->Draw("same");
  ccent->cd(3);
rc3y->GetXaxis()->SetTitle("y_gen-y");  
rc3y->Draw();
  g3y->Draw("same");
  ccent->cd(4);
rc4y->GetXaxis()->SetTitle("y_gen-y");  
rc4y->Draw();
  g4y->Draw("same");

 ccent->SaveAs("plots/ycentres.pdf");
 ccent->SaveAs("plots/ycentres.png");
 
 TCanvas * yResponsec= new TCanvas("yResponsec","");
yResponse->GetYaxis()->SetTitle("y");
yResponse->GetXaxis()->SetTitle("y_gen");
yResponsec->SetLogz();
yResponse->Draw("COLZ");
yResponsec->SaveAs("plots/yResponse.pdf");

TCanvas*yresoc=new TCanvas("yresoc","");
yreso->GetXaxis()->SetTitle("y_gen-y");
yreso->Draw();
yresoc->SaveAs("plots/yresolution.pdf");


  TCanvas *ccenty=new TCanvas("ccenty","Rapidity Resolution centrality dependance");
  ccenty->Divide(2,2);
  ccenty->cd(1);
  rc1y->Draw();
  g1y->Draw("same");
  ccenty->cd(2);
  rc2y->Draw();
  g2y->Draw("same");
  ccenty->cd(3);
  rc3y->Draw();
  g3y->Draw("same");
  ccenty->cd(4);
  rc4y->Draw();
  g4y->Draw("same");

 ccenty->SaveAs("plots/centresy.pdf");
 ccenty->SaveAs("plots/centresy.png");

  TCanvas *cy=new TCanvas("cy","Resolution y dependance");
  cy->Divide(3,4);
  cy->cd(1);
  ry1->Draw();
  y1->Draw("same");
  cy->cd(2);
  ry2->Draw();
  y2->Draw("same");
  cy->cd(3);
  ry3->Draw();
  y3->Draw("same");
  cy->cd(4);
  ry4->Draw();
  y4->Draw("same");
  cy->cd(5);
  ry5->Draw();
  y5->Draw("same");
  cy->cd(6);
  ry6->Draw();
  y6->Draw("same");
  cy->cd(7);
  ry7->Draw();
  y7->Draw("same");
  cy->cd(8);
  ry8->Draw();
  y8->Draw("same");
  cy->cd(9);
  ry9->Draw();
  y9->Draw("same");
  cy->cd(10);
  ry10->Draw();
  y10->Draw("same");
  cy->cd(11);
  ry11->Draw();
  y11->Draw("same");
  cy->cd(12);
  ry12->Draw();
  y12->Draw("same");
  

 cy->SaveAs("plots/yres.pdf");
 cy->SaveAs("plots/yres.png");


  TFile * output = new TFile("resources/Z2mumu_Efficiencies.root","recreate");
  for(int i = 0; i<nBins; i++){
   recoEff[i]->Write();
   recoEff_pass[i]->Write();
   recoEff_net[i]->Write();
   eff[i]->Write();
  }
TH1D* rationum=(TH1D*) rico->Clone("rationum");
rationum->Divide(genn);
    
     res->Write();
     yreso->Write();

  Response->Write();
  genn->Write();
  rico->Write();
  rationum->Write();
  SmearResponse->Write();
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
