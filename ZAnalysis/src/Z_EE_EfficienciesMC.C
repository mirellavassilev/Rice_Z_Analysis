#include "include/centralityBine.h"
#include "include/electronEnergyScalee.h"
#include "include/forceConsistencye.h"
#include "include/electronSelectore.h"
#include "include/electronTriggerMatchinge.h"
#include "include/centralityToole.h"
#include "include/Settingse.h"
#include "include/Timere.h"
#include "include/MCReweighte.h"
#include "include/ElectronTnPe.h"

#include "TRandom3.h"
#include <TLegend.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TF1.h>
#include <TCanvas.h>


void doZ2EE(std::vector< std::string > files, int jobNumber){
  Timer timer = Timer();
  timer.Start();
  timer.StartSplit("Start Up");

  TH1::SetDefaultSumw2();
  ElectronEnergyScale energyScale = ElectronEnergyScale("MC");
  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  ElectronTnP eTnP = ElectronTnP();  

  MCReweight vzRW = MCReweight("resources/vzReweighte.root");
  
  Settings s = Settings();

  CentralityBin cb = CentralityBin();
  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();
  
  TRandom3 * r1 = new TRandom3();

  TH1D * recoEff_pt_pass[nBins];
  TH1D * recoEff_pt_net[nBins];
  TH1D * recoEff_pt[nBins];  
  TEfficiency * eff_pt[nBins];
  
  TH1D * recoEff_y_pass[nBins];
  TH1D * recoEff_y_net[nBins];
  TH1D * recoEff_y[nBins];  
  TEfficiency * eff_y[nBins];
  
  TH1D * recoEff_phi_pass[nBins];
  TH1D * recoEff_phi_net[nBins];
  TH1D * recoEff_phi[nBins];  
  TEfficiency * eff_phi[nBins];
  
  TH1D * recoEff_cent_pass[nBins];
  TH1D * recoEff_cent_net[nBins];
  TH1D * recoEff_cent[nBins];  
  TEfficiency * eff_cent[nBins];

  TH2D * recoEff_pass[nBins];
  TH2D * recoEff_net[nBins];
  TH2D * recoEff[nBins];
  TEfficiency * eff[nBins];
  
  TH2D * recoEff_noSF_pass[nBins];
  TH2D * recoEff_noSF_net[nBins];
  TH2D * recoEff_noSF[nBins];
  TEfficiency * eff_noSF[nBins];
  ///////Histograms Resolution///////////////////////////
  TH1D *res= new TH1D("res","Resolution",30,-.2,.2);
    TH1D *yreso= new TH1D("yreso","Rapidity Resolution",35,-.2,.2);

  // pT 
  Float_t Bins[]={0,1.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,90.0,120.0,150.0,200.0};
  Int_t  nptBins = 13; 
  TH1D *genn=new TH1D("genn","",nptBins,Bins);
  TH1D *rico=new TH1D("rico","",nptBins,Bins);
    TH1D *genny=new TH1D("genny","",13,-2.1,2.1);
  TH1D *ricoy=new TH1D("ricoy","",13,-2.1,2.1);
 //Respose Matrix
  TH2D *Response= new TH2D("Response","",nptBins,Bins,nptBins,Bins);
  TH2D * SmearResponse=new TH2D ("SmearResponse","",nptBins,Bins,nptBins,Bins);
  TH2D *yResponse= new TH2D ("Rapidity Response Matrix","",13,-2.4,2.4,13,-2.4,2.4);
  TH1D *ybb=new TH1D ("Rapidity Resolution","",12,-2.4,2.4);
 
 TH1D *rc1= new TH1D("rc1","Centrality 0-20",30,-.2,.2);
  TF1 *g1  = new TF1("g1","gaus",-.1,.1);
  TH1D *rc2= new TH1D("rc2","Centrality 20-60",30,-.2,.2);
  TF1 *g2  = new TF1("g2","gaus",-.1,.1);
  TH1D *rc3= new TH1D("rc3","Centrality 60-100",30,-.2,.2);
  TF1 *g3  = new TF1("g3","gaus",-.1,.1);
  TH1D *rc4= new TH1D("rc4","Centrality 100-200",30,-.2,.2);
  TF1 *g4  = new TF1("g4","gaus",-.1,.1);

////////////////////////////////////////////////////////////
  for(int k = 0; k<nBins; k++){
    recoEff_pass[k] = new TH2D(Form("recoEff_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_net[k] = new TH2D(Form("recoEff_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_noSF_pass[k] = new TH2D(Form("recoEff_noSF_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    recoEff_noSF_net[k] = new TH2D(Form("recoEff_noSF_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBinsEle,-s.maxZRapEle,s.maxZRapEle,s.nZPtBins-1,s.zPtBins);
    
    recoEff_pt_pass[k] = new TH1D(Form("recoEff_pt_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_pt_net[k] = new TH1D(Form("recoEff_pt_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZPtBins-1,s.zPtBins);
    recoEff_y_pass[k] = new TH1D(Form("recoEff_y_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_y_net[k] = new TH1D(Form("recoEff_y_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",s.nZRapBins,-s.maxZRap,s.maxZRap);
    recoEff_phi_pass[k] = new TH1D(Form("recoEff_phi_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_phi_net[k] = new TH1D(Form("recoEff_phi_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",30,-TMath::Pi(),TMath::Pi());
    recoEff_cent_pass[k] = new TH1D(Form("recoEff_cent_pass_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
    recoEff_cent_net[k] = new TH1D(Form("recoEff_cent_net_%d_%d",c.getCentBinLow(k),c.getCentBinHigh(k)),"",20,0,100);
  }

  int nEle;
  int hiBin;
  float hiHF;
  float vz;

  int pprimaryVertexFilter;
  int phfCoincFilter2Th4;
  int pclusterCompatibilityFilter;

  int HLT_DoubleEle10;
  int HLT_SingleEle20;

  std::vector< float > * elePt = 0;
  std::vector< float > * eleEta = 0;
  std::vector< float > * elePhi = 0;
  std::vector< float > * eleSigmaIEtaIEta = 0;
  std::vector< int > * eleCharge = 0;
  std::vector< int > * eleMissHits = 0;
  std::vector< float > * eledEtaAtVtx = 0;
  std::vector< float > * eledPhiAtVtx = 0;
  std::vector< float > * eleSCEta = 0;
  std::vector< float > * eleSCPhi = 0;
  std::vector< float > * eleHoverEBc = 0;
  //std::vector< float > * eleD0 = 0;
  //std::vector< float > * eleDz = 0;
  std::vector< float > * eleIP3D = 0;
  std::vector< float > * eleEoverPInv = 0;

  int nMC = 0;
  //std::vector< int > * mcStatus = 0;
  std::vector< int > * mcPID = 0;
  std::vector< float > * mcPt = 0;
  std::vector< float > * mcEta = 0;
  //std::vector< float > * mcPhi = 0;
  std::vector< int > * mcMomPID = 0;
  std::vector< int > * mcGMomPID = 0;
  std::vector< float > * mcMomPt = 0;
  std::vector< float > * mcMomEta = 0;
  std::vector< float > * mcMomPhi = 0; 
  std::vector< float > * mcMomMass = 0;

  for(unsigned int f = 0; f<files.size(); f++){
    timer.StartSplit("Opening Files");

    TFile * in = TFile::Open(files.at(f).c_str(),"read");
    if(f%5 == 0)  std::cout << f << "/" << files.size() << std::endl;

    TTree * hltTree = (TTree*)in->Get("hltanalysis/HltTree");
    hltTree->SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v1",&HLT_DoubleEle10); 
    hltTree->SetBranchAddress("HLT_HIEle20Gsf_v1",&HLT_SingleEle20); 

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    TTree * eTreeMC = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
    eTreeMC->SetBranchAddress("nMC",&nMC);
    //eTreeMC->SetBranchAddress("mcStatus",&mcStatus);
    eTreeMC->SetBranchAddress("mcPID",&mcPID);
    eTreeMC->SetBranchAddress("mcMomPID",&mcMomPID);
    eTreeMC->SetBranchAddress("mcGMomPID",&mcGMomPID);
   
    eTreeMC->SetBranchAddress("mcPt",&mcPt);
    eTreeMC->SetBranchAddress("mcEta",&mcEta);
    //eTreeMC->SetBranchAddress("mcPhi",&mcPhi);
    eTreeMC->SetBranchAddress("mcMomPt",&mcMomPt);
    eTreeMC->SetBranchAddress("mcMomEta",&mcMomEta);
    eTreeMC->SetBranchAddress("mcMomPhi",&mcMomPhi);
    eTreeMC->SetBranchAddress("mcMomMass",&mcMomMass);

    eTree->SetBranchAddress("nEle",&nEle);
    eTree->SetBranchAddress("elePt",&elePt);
    eTree->SetBranchAddress("eleEta",&eleEta);
    eTree->SetBranchAddress("elePhi",&elePhi);
    eTree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta);
    eTree->SetBranchAddress("eleMissHits",&eleMissHits);
    eTree->SetBranchAddress("eleCharge",&eleCharge);
    eTree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
    eTree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
    eTree->SetBranchAddress("eleSCEta",&eleSCEta);
    eTree->SetBranchAddress("eleSCPhi",&eleSCPhi);
    eTree->SetBranchAddress("eleHoverEBc",&eleHoverEBc);
    //eTree->SetBranchAddress("eleD0",&eleD0);
    //eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleIP3D",&eleIP3D);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    //evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("hiHF",&hiHF);
    evtTree->SetBranchAddress("vz",&vz);
  
    TTree * skimTree = (TTree*)in->Get("skimanalysis/HltTree");
    skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

    TTree * L1Tree = (TTree*)in->Get("l1object/L1UpgradeFlatTree");
    L1Tree->SetBranchAddress("nEGs",&(eTrig.L1nEGs));
    L1Tree->SetBranchAddress("egEta", &(eTrig.L1egEta));
    L1Tree->SetBranchAddress("egPhi", &(eTrig.L1egPhi));
    L1Tree->SetBranchAddress("egEt", &(eTrig.L1egEt));

    TTree * HLTObjTree;
    HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIEle20Gsf_v");
    HLTObjTree->SetBranchAddress("eta",&(eTrig.HLTEta));
    HLTObjTree->SetBranchAddress("phi",&(eTrig.HLTPhi));
    HLTObjTree->SetBranchAddress("pt",&(eTrig.HLTPt));

    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      if(i%5000 == 0) std::cout << i << "/" << eTree->GetEntries() << std::endl;
      timer.StartSplit("Checking Evt Selections");
      //event selections
      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;
      hiBin = cb.getHiBinFromhiHF(hiHF,3);

      skimTree->GetEntry(i);
      if(! (pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter)) continue;
      
      double eventWeight = vzRW.reweightFactor( vz ) * c.findNcoll( hiBin );
      
      timer.StartSplit("Loading GEN electron tree");
      eTreeMC->GetEntry(i);

      timer.StartSplit("Gen Loop");
      int nGenElectronsFound = 0;
      
      TLorentzVector mom = TLorentzVector();
      bool foundGen = false;
      for(int j = 0; j<nMC; j++){
        //break out if you find a Z->tautau
        if( mcGMomPID->at(j) == 23 && TMath::Abs(mcMomPID->at(j)) == 15) break;
        
        //only look for electrons coming directly from a Z
        if( TMath::Abs(mcMomPID->at(j)) != 23) continue;
        
        //if it decays to muons, break out
        if( TMath::Abs(mcPID->at(j)) == 13) break;       
   
        //make sure it's an electron
        if( TMath::Abs(mcPID->at(j)) != 11) continue;       
        nGenElectronsFound++;

        //if they are not in our acceptance, break out
        if( mcPt->at(j) < s.minElectronPt ) break; 
        if( TMath::Abs(mcEta->at(j)) > s.maxZRapEle ) break; 

        //we found both daughters and they are in our acceptance, lets fill our histogram
        if( nGenElectronsFound == 2 ){
          //Fill denominator
          foundGen = true; 
          mom.SetPtEtaPhiM( mcMomPt->at(j), mcMomEta->at(j), mcMomPhi->at(j), mcMomMass->at(j));
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              recoEff_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_noSF_net[k]->Fill( mom.Rapidity(), mcMomPt->at(j), eventWeight );
              recoEff_pt_net[k]->Fill( mom.Pt(), eventWeight);
              recoEff_y_net[k]->Fill( mom.Rapidity(), eventWeight);
              recoEff_cent_net[k]->Fill( hiBin/2.0, eventWeight);
              recoEff_phi_net[k]->Fill( mom.Phi(), eventWeight);

              //for the last few pt bins, halve the y binning so we have better stats
              if(mcMomPt->at(j) > s.zPtBins[s.nZPtBins - s.nPtBinsToRebinRapEff]){
                int bin = recoEff_net[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                if( bin%2 ==1){
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin+1), mcMomPt->at(j), eventWeight );
                } else {
                  recoEff_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
                  recoEff_noSF_net[k]->Fill( recoEff_net[k]->GetXaxis()->GetBinCenter(bin-1), mcMomPt->at(j), eventWeight );
                }
              }
            }
          }
          break; 
        }//end of if statement
      }//end of gen loop      
      if( !foundGen ) continue;
      
      timer.StartSplit("Loading RECO electron tree");
      eTree->GetEntry(i);

      timer.StartSplit("Checking Number of electrons");
      if(nEle<2) continue;
      
      timer.StartSplit("Checking HLT Selections");
      //check for the trigger we want (double ele 10 or single ele 20)
      hltTree->GetEntry(i);
      if( !HLT_SingleEle20) continue;

      timer.StartSplit("Electron Cuts");
      //make a list of electrons passing our cuts
      std::vector< int > goodElectrons;
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        //correct Pt
        elePt->at(j) = energyScale.correctPt( elePt->at(j), eleSCEta->at(j), hiBin);
        
        if(elePt->at(j)< s.minElectronPt) continue;
        if(TMath::Abs(eleSCEta->at(j)) > s.maxZRapEle) continue;
        //veto on transition region
        if(TMath::Abs(eleSCEta->at(j)) > 1.442 && TMath::Abs(eleSCEta->at(j)) < 1.556 ) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        float dEta = TMath::Abs( eledEtaAtVtx->at(j) );
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverEBc->at(j), eleEoverPInv->at(j), eleIP3D->at(j) )) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      bool moreThan2 = false;
      if(goodElectrons.size()>2) moreThan2 = true;

      timer.StartSplit("Loading HLT/L1 Object stuff");
      //get trigger matching stuff
      L1Tree->GetEntry(i);
      HLTObjTree->GetEntry(i);
    
      //make Z candidates 
      timer.StartSplit("Z candidates");
      TLorentzVector * elec1 = new TLorentzVector();
      TLorentzVector * elec2 = new TLorentzVector();
      for(unsigned int j = 0; j<goodElectrons.size(); j++){
        elec1->SetPtEtaPhiM(elePt->at(goodElectrons.at(j)), eleEta->at(goodElectrons.at(j)), elePhi->at(goodElectrons.at(j)), 0.000511);
        for(unsigned int j2 = j+1; j2<goodElectrons.size(); j2++){
          elec2->SetPtEtaPhiM(elePt->at(goodElectrons.at(j2)), eleEta->at(goodElectrons.at(j2)), elePhi->at(goodElectrons.at(j2)), 0.000511);
          TLorentzVector Zcand = *elec1+*elec2;
          if(Zcand.M() < s.zMassRange[0] || Zcand.M() > s.zMassRange[1]) continue;      

          //L1 trigger matching (1 L1 EG > 15 GeV)
          bool isFirstElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 15.0);
          bool isSecondElectronL1Matched =  matcher.isL1Matched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 15.0);
          if(! (isFirstElectronL1Matched || isSecondElectronL1Matched)) continue;

          //HLT trigger matching (1 HLT match > 20 GeV)
          bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
          bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
          if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
  

          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " " << Zcand.Eta() << " " << Zcand.Phi() << " " << mom.Pt() << " " << mom.Eta() << " " << mom.Phi() <<  " isOS? " << (int)isOppositeSign << std::endl;
          if(!isOppositeSign) continue;
          if(moreThan2 && TMath::ACos(TMath::Cos(Zcand.Phi() - mom.Phi())) > 0.1 ) continue;
          //Looks like this is a good candidate match! Let's get the scale factor
          float scaleFactor = eTnP.getZSF(hiBin, elePt->at(goodElectrons.at(j)), eleSCEta->at(goodElectrons.at(j)), elePt->at(goodElectrons.at(j2)), eleSCEta->at(goodElectrons.at(j2)), 0) ;

//////////////////////Fill Resolution Histograms////////////////////////
  //Fill Response Histogram
  Response->Fill(Zcand.Pt(),mom.Pt());
  SmearResponse->Fill(Zcand.Pt() *r1->Gaus(1,0.05),mom.Pt());  
  yResponse->Fill(Zcand.Rapidity(),mom.Rapidity());
//Fill gen and rico
  genn->Fill(mom.Pt());
  rico->Fill(Zcand.Pt());
  genny->Fill(mom.Rapidity());
  ricoy->Fill(Zcand.Rapidity());
 //Fill Resolution
  res->Fill((mom.Pt()-Zcand.Pt())/mom.Pt());
  yreso->Fill(mom.Rapidity()-Zcand.Rapidity());
ybb->Fill(Zcand.Rapidity()/mom.Rapidity());


        if (hiBin>0&&hiBin<20) {
        rc1->Fill(mom.Rapidity()-Zcand.Rapidity());
        rc1->SetLineColor(kGreen);
        }
        else if (hiBin>20&&hiBin<60) {
        rc2->Fill(mom.Rapidity()-Zcand.Rapidity());
        rc2->SetLineColor(kBlue);
        }
        else if (hiBin>60&&hiBin<100) {
        rc3->Fill(mom.Rapidity()-Zcand.Rapidity());
        rc3->SetLineColor(kRed);
        }
        else if (hiBin>100&&hiBin<200) {
        rc4->Fill(mom.Rapidity()-Zcand.Rapidity());
        rc4->SetLineColor(kYellow);
        }


///////////////////////////////////////////////////////////////////////




          //Fill numerator (and apply the scale factor here!)
          for(int k = 0; k<nBins; k++){ 
            if(c.isInsideBin(hiBin,k)){
              //make sure this is in our fiducial histogram range otherwise CheckConsistency can freak out
              if( mom.Pt() < s.zPtBins[ s.nZPtBins-1 ] && TMath::Abs( mom.Rapidity() ) < s.maxZRap ){
                recoEff_pass[k]->Fill( mom.Rapidity(), mom.Pt(), eventWeight * scaleFactor );
                recoEff_pt_pass[k]->Fill( mom.Pt(), eventWeight * scaleFactor);
                recoEff_y_pass[k]->Fill( mom.Rapidity(), eventWeight * scaleFactor);
                recoEff_cent_pass[k]->Fill( hiBin/2.0, eventWeight * scaleFactor);
                recoEff_phi_pass[k]->Fill( mom.Phi(), eventWeight * scaleFactor);
          
                //for the last few pt bins, halve the y binning so we have better stats
                if(mom.Pt()> s.zPtBins[s.nZPtBins- s.nPtBinsToRebinRapEff]){
                  int bin = recoEff_pass[k]->GetXaxis()->FindBin( mom.Rapidity() ); 
                  if( bin%2 ==1){
                    recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight * scaleFactor );
                    recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin+1), mom.Pt(), eventWeight );
                  } else { 
                    recoEff_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight * scaleFactor );
                    recoEff_noSF_pass[k]->Fill( recoEff_pass[k]->GetXaxis()->GetBinCenter(bin-1), mom.Pt(), eventWeight );
                  }
                }
              }
            }
          }
        }
      }
      delete elec1;
      delete elec2;
    }

    delete eTree;
    in->Close();
  }

  timer.StartSplit("End of analysis");
  std::vector< bool > isConsistent;
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_pass[i], recoEff_net[i]);
    recoEff[i] = (TH2D*)recoEff_pass[i]->Clone(Form("recoEff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff[i]->Divide(recoEff_net[i]);
    recoEff[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pass[i]), *(recoEff_net[i]),"w") ){
      eff[i] = new TEfficiency(*(recoEff_pass[i]), *(recoEff_net[i]));
      eff[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff[i]->SetName(Form("eff_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff[i]->SetDirectory(0);
      isConsistent.push_back(true);
    }
    else{
      isConsistent.push_back(false);
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }
    recoEff_net[i]->SetDirectory(0);
    recoEff_pass[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    recoEff_noSF[i] = (TH2D*)recoEff_noSF_pass[i]->Clone(Form("recoEff_noSF_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_noSF[i]->Divide(recoEff_noSF_net[i]);
    recoEff_noSF[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_noSF_pass[i]), *(recoEff_noSF_net[i]),"w") ){
      eff_noSF[i] = new TEfficiency(*(recoEff_noSF_pass[i]), *(recoEff_noSF_net[i]));
      eff_noSF[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_noSF[i]->SetName(Form("eff_noSF_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_noSF[i]->SetDirectory(0);
    }
    recoEff_noSF_net[i]->SetDirectory(0);
    recoEff_noSF_pass[i]->SetDirectory(0);
  }
  
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_pt_pass[i], recoEff_pt_net[i]);
    recoEff_pt[i] = (TH1D*)recoEff_pt_pass[i]->Clone(Form("recoEff_pt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_pt[i]->Divide(recoEff_pt_net[i]);
    recoEff_pt[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_pt_pass[i]), *(recoEff_pt_net[i]),"w") ){
      eff_pt[i] = new TEfficiency(*(recoEff_pt_pass[i]), *(recoEff_pt_net[i]));
      eff_pt[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_pt[i]->SetName(Form("eff_pt_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_pt[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_pt_pass[i]->SetDirectory(0);
    recoEff_pt_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_y_pass[i], recoEff_y_net[i]);
    recoEff_y[i] = (TH1D*)recoEff_y_pass[i]->Clone(Form("recoEff_y_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_y[i]->Divide(recoEff_y_net[i]);
    recoEff_y[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_y_pass[i]), *(recoEff_y_net[i]),"w") ){
      eff_y[i] = new TEfficiency(*(recoEff_y_pass[i]), *(recoEff_y_net[i]));
      eff_y[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_y[i]->SetName(Form("eff_y_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_y[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_y_pass[i]->SetDirectory(0);
    recoEff_y_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_phi_pass[i], recoEff_phi_net[i]);
    recoEff_phi[i] = (TH1D*)recoEff_phi_pass[i]->Clone(Form("recoEff_phi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_phi[i]->Divide(recoEff_phi_net[i]);
    recoEff_phi[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_phi_pass[i]), *(recoEff_phi_net[i]),"w") ){
      eff_phi[i] = new TEfficiency(*(recoEff_phi_pass[i]), *(recoEff_phi_net[i]));
      eff_phi[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_phi[i]->SetName(Form("eff_phi_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_phi[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_phi_pass[i]->SetDirectory(0);
    recoEff_phi_net[i]->SetDirectory(0);
  }
  for(int i = 0; i<nBins; i++){
    forceConsistency(recoEff_cent_pass[i], recoEff_cent_net[i]);
    recoEff_cent[i] = (TH1D*)recoEff_cent_pass[i]->Clone(Form("recoEff_cent_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
    recoEff_cent[i]->Divide(recoEff_cent_net[i]);
    recoEff_cent[i]->SetDirectory(0);

    if( TEfficiency::CheckConsistency(*(recoEff_cent_pass[i]), *(recoEff_cent_net[i]),"w") ){
      eff_cent[i] = new TEfficiency(*(recoEff_cent_pass[i]), *(recoEff_cent_net[i]));
      eff_cent[i]->SetStatisticOption(TEfficiency::kBJeffrey);
      eff_cent[i]->SetName(Form("eff_cent_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)));
      eff_cent[i]->SetDirectory(0);
    }
    else{
      std::cout << "Warning, these histograms are not consistent!" << std::endl;
    }

    recoEff_cent_pass[i]->SetDirectory(0);
    recoEff_cent_net[i]->SetDirectory(0);
  }
    
rc1->Fit(g1,"","",-.05,.05);
rc2->Fit(g2,"","",-.05,.05);
rc3->Fit(g3,"","",-.05,.05);
rc4->Fit(g4,"","",-.05,.05);

 
res->SetDirectory(0);
yreso->SetDirectory(0);
genn->SetDirectory(0);
rico->SetDirectory(0);
Response->SetDirectory(0);
ricoy->SetDirectory(0);
genny->SetDirectory(0);
SmearResponse->SetDirectory(0);
yResponse->SetDirectory(0);
ybb->SetDirectory(0);
    rc1->SetDirectory(0);
    rc2->SetDirectory(0);
    rc3->SetDirectory(0);
    rc4->SetDirectory(0);


  TCanvas *ccent=new TCanvas("ccent","Resolution centrality dependance");
  ccent->Divide(2,2);
  ccent->cd(1);
  rc1->GetXaxis()->SetTitle("mom.Rapidity-Zcand.Rapidity");
  rc1->Draw();
  g1->Draw("same");
  ccent->cd(2);
  rc2->GetXaxis()->SetTitle("mom.Rapidity-Zcand.Rapidity");
  rc2->Draw();
  g2->Draw("same");
  ccent->cd(3);
  rc3->GetXaxis()->SetTitle("mom.Rapidity-Zcand.Rapidity");
  rc3->Draw();
  g3->Draw("same");
  ccent->cd(4);
  rc4->GetXaxis()->SetTitle("mom.Rapidity-Zcand.Rapidity");
  rc4->Draw();
  g4->Draw("same");

 ccent->SaveAs("plots/ycentrese.pdf");
 ccent->SaveAs("plots/ycentrese.png");

TCanvas*yresoc=new TCanvas("yresoc","");
yreso->GetXaxis()->SetTitle("mom.Rapidity-Zcand.Rapidity");
yreso->Draw();
yresoc->SaveAs("plots/yresolutione.pdf");

  TFile * output = new TFile(Form("resources/Z2ee_EfficiencyMC_%d.root",jobNumber),"recreate");
  for(int i = 0; i<nBins; i++){
    recoEff_net[i]->Write();
    recoEff_pass[i]->Write();
    recoEff[i]->Write();
    eff[i]->Write();
    
    recoEff_noSF_net[i]->Write();
    recoEff_noSF_pass[i]->Write();
    recoEff_noSF[i]->Write();
    eff_noSF[i]->Write();
    
    recoEff_pt[i]->Write();
    recoEff_pt_pass[i]->Write();
    recoEff_pt_net[i]->Write();
    eff_pt[i]->Write();
    
    recoEff_y[i]->Write();
    recoEff_y_pass[i]->Write();
    recoEff_y_net[i]->Write();
    eff_y[i]->Write();
    
    recoEff_phi[i]->Write();
    recoEff_phi_pass[i]->Write();
    recoEff_phi_net[i]->Write();
    eff_phi[i]->Write();
    
    recoEff_cent[i]->Write();
    recoEff_cent_pass[i]->Write();
    recoEff_cent_net[i]->Write();
    eff_cent[i]->Write();
  }
  
 res->Write();
 yreso->Write(); 
genn->Write();
 rico->Write();
genny->Write();
ricoy->Write();
Response->Write();
  SmearResponse->Write();   
yResponse->Write();
ybb->Write();
  output->Close();

  timer.Stop();
  timer.Report();


////////////////////////Make Canvas////////////////////////////////////

TCanvas *Res = new TCanvas("Res","");
res->Draw();
Res->SaveAs("plots/rese.pdf");

TCanvas *response = new TCanvas("responsee","");
Response->Draw();
response->SaveAs("plots/responsee.pdf");

 TCanvas * yResponsec= new TCanvas("yResponsec","");
yResponse->GetYaxis()->SetTitle("Zcand.Rapidity");
yResponse->GetXaxis()->SetTitle("mom.Rapidity");
yResponsec->SetLogz();
yResponse->Draw("COLZ");
yResponsec->SaveAs("plots/yResponsee.pdf");

TCanvas *Yrate=new TCanvas("Yrate","");
TH1D *yrat=(TH1D*) genny->Clone("yrat");
yrat->Divide(ricoy);
yrat->GetXaxis()->SetTitle("Rapidity");
yrat->GetYaxis()->SetTitle("mom.Rapidity/Zcand.Rapidity");
yrat->Draw();
Yrate->SaveAs ("plots/Yrate.pdf");
/////////////////////////////////////////////////////////////////////////


  return;
}


int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: Z_EE_EfficiencyMC <fileList> <job #> <total number of jobs>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);

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
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  doZ2EE(listOfFiles, job);
  return 0; 
}
