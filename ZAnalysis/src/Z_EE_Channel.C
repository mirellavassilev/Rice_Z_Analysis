#include "include/electronSelector.h"
#include "include/electronTriggerMatching.h"
#include "include/centralityTool.h"
#include "include/Settings.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void doZ2EE(std::vector< std::string > files){
  //switch between single electron 20 and double electron 10
  bool doSingleEle20 = true;

  ElectronSelector eSel = ElectronSelector();
  ElectronTriggerMatcher matcher = ElectronTriggerMatcher();
  ElecTrigObject eTrig = ElecTrigObject();
  Settings s = Settings();

  CentralityTool c = CentralityTool();
  const int nBins = c.getNCentBins();

  TH1D * massPeakOS[nBins]; 
  TH1D * massPeakSS[nBins]; 
  
  TProfile * v2Num[nBins];
  TProfile * v2NumVsCent;
  TProfile * v2Denom[nBins];
  TProfile * v2DenomVsCent;
  
  TProfile * v2EleNum[nBins];
  TProfile * v2EleNumVsCent;
  TProfile * v2EleDenom[nBins];
  TProfile * v2EleDenomVsCent;

  TH1D * v2NumVsCentHist;
  TH1D * v2SqrtDenomVsCent;
  TH1D * v2VsCent;
  
  TH1D * v2EleNumVsCentHist;
  TH1D * v2EleSqrtDenomVsCent;
  TH1D * v2EleVsCent;
  
  for(int i = 0; i<nBins; i++){
    massPeakOS[i] = new TH1D(Form("massPeakOS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{+}e^{-}};counts",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    massPeakSS[i] = new TH1D(Form("massPeakSS_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),";m_{e^{#pm}e^{#pm}}",s.nZMassBins,s.zMassRange[0],s.zMassRange[1]);
    v2Num[i] = new TProfile(Form("v2Num_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2Denom[i] = new TProfile(Form("v2Denom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    
    v2EleNum[i] = new TProfile(Form("v2EleNum_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
    v2EleDenom[i] = new TProfile(Form("v2EleDenom_%d_%d",c.getCentBinLow(i),c.getCentBinHigh(i)),"",1,0,1);
  }  
  v2NumVsCent = new TProfile("v2NumVsCent","",nBins,0,nBins);
  v2DenomVsCent = new TProfile("v2DenomVsCent","",nBins,0,nBins);
  v2SqrtDenomVsCent = new TH1D("v2SqrtDenomVsCent","",nBins,0,nBins);
  v2NumVsCentHist = new TH1D("v2NumVsCentHist","",nBins,0,nBins);
  
  v2EleNumVsCent = new TProfile("v2EleNumVsCent","",nBins,0,nBins);
  v2EleDenomVsCent = new TProfile("v2EleDenomVsCent","",nBins,0,nBins);
  v2EleSqrtDenomVsCent = new TH1D("v2EleSqrtDenomVsCent","",nBins,0,nBins);
  v2EleNumVsCentHist = new TH1D("v2EleNumVsCentHist","",nBins,0,nBins);

  int nEle;
  int hiBin;
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
  std::vector< float > * eleHoverE = 0;
  std::vector< float > * eleD0 = 0;
  std::vector< float > * eleDz = 0;
  std::vector< float > * eleEoverPInv = 0;

  int hiNevtPlane;
  float hiQVecMag[200];
  float hiQVecAngle[200];

  for(unsigned int f = 0; f<files.size(); f++){
    TFile * in = TFile::Open(files.at(f).c_str(),"read");
    if(f%10 == 0)  std::cout << f << std::endl;

    TTree * hltTree = (TTree*)in->Get("hltanalysis/HltTree");
    hltTree->SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v1",&HLT_DoubleEle10); 
    hltTree->SetBranchAddress("HLT_HIEle20Gsf_v1",&HLT_SingleEle20); 

    TTree * eTree = (TTree*)in->Get("ggHiNtuplizerGED/EventTree");
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
    eTree->SetBranchAddress("eleHoverE",&eleHoverE);
    eTree->SetBranchAddress("eleD0",&eleD0);
    eTree->SetBranchAddress("eleDz",&eleDz);
    eTree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);

    TTree * evtTree = (TTree*)in->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("vz",&vz);
    evtTree->SetBranchAddress("hiNevtPlane",&hiNevtPlane);  
    evtTree->SetBranchAddress("hiQVecMag",hiQVecMag);  
    evtTree->SetBranchAddress("hiQVecAngle",hiQVecAngle);  
  
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
    if(!doSingleEle20){
      HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIDoubleEle10GsfMass50_v");
    } else {
      HLTObjTree = (TTree*)in->Get("hltobject/HLT_HIEle20Gsf_v");
    }
    HLTObjTree->SetBranchAddress("eta",&(eTrig.HLTEta));
    HLTObjTree->SetBranchAddress("phi",&(eTrig.HLTPhi));
    HLTObjTree->SetBranchAddress("pt",&(eTrig.HLTPt));

    for(unsigned int i = 0; i < eTree->GetEntries(); i++){
      hltTree->GetEntry(i);

      //check for the trigger we want (double ele 10 or single ele 20)
      if((!doSingleEle20) && (!HLT_DoubleEle10)) continue;
      if( doSingleEle20 && (!HLT_SingleEle20)) continue;

      eTree->GetEntry(i);
      if(nEle<2) continue;

      skimTree->GetEntry(i);
      if(! (pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter)) continue;

      evtTree->GetEntry(i);
      if(TMath::Abs(vz)>15) continue;

      L1Tree->GetEntry(i);
      HLTObjTree->GetEntry(i);

      std::vector< int > goodElectrons;
   
      for(unsigned int j = 0; j < (unsigned int) nEle; j++){
        if(elePt->at(j)<20) continue;
        if(TMath::Abs(eleSCEta->at(j)) > 2.4) continue;
        //veto on dead endcap region
        if(eleSCEta->at(j) < -1.39 && eleSCPhi->at(j) < -0.9 && eleSCPhi->at(j) > -1.6) continue;

        //check electron qualty variables
        float dEta = TMath::Abs( eledEtaAtVtx->at(j) );
        float dPhi = TMath::Abs( eledPhiAtVtx->at(j) );
        if(!eSel.isGoodElectron(ElectronSelector::WorkingPoint::loose, hiBin, eleSCEta->at(j), eleSigmaIEtaIEta->at(j), dEta, dPhi, eleMissHits->at(j), eleHoverE->at(j), eleEoverPInv->at(j), eleD0->at(j), eleDz->at(j))) continue;

        goodElectrons.push_back(j);
      }

      if(goodElectrons.size()<2) continue;
      bool moreThan2 = false;
      if(goodElectrons.size()>2) moreThan2 = true;
     
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
          if(doSingleEle20){
            bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 20.0);
            bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 20.0);
            if(! (isFirstElectronHLTMatched || isSecondElectronHLTMatched)) continue;
          } else {  //(2 HLT matches > 10 GeV)
            bool isFirstElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j)), eleSCPhi->at(goodElectrons.at(j)), eTrig, 10.0);
            bool isSecondElectronHLTMatched = matcher.isHLTMatched(eleSCEta->at(goodElectrons.at(j2)), eleSCPhi->at(goodElectrons.at(j2)), eTrig, 10.0);
            if(! (isFirstElectronHLTMatched && isSecondElectronHLTMatched)) continue;
          }
   
          bool isOppositeSign =  eleCharge->at(goodElectrons.at(j)) != eleCharge->at(goodElectrons.at(j2));
          if(moreThan2) std::cout << j << " " << j2 << " " << Zcand.M() <<" " << Zcand.Pt() << " isOS? " << (int)isOppositeSign << std::endl;
          if( isOppositeSign){
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBin,k)) massPeakOS[k]->Fill( Zcand.M() );
            }
          }else{
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBin,k)) massPeakSS[k]->Fill( Zcand.M() );
            }
          }

          //v2 calcuation
          if( isOppositeSign){
            TComplex Qp = TComplex(hiQVecMag[7], hiQVecAngle[7], true);
            TComplex Qn = TComplex(hiQVecMag[6], hiQVecAngle[6], true);

            TComplex candQ = TComplex(1, 2*Zcand.Phi(), true);
            TComplex ele1Q = TComplex(1, 2*elePhi->at(goodElectrons.at(j)), true);
            TComplex ele2Q = TComplex(1, 2*elePhi->at(goodElectrons.at(j2)), true);
            for(int k = 0; k<nBins; k++){
              if(c.isInsideBin(hiBin,k)){
                TComplex Q = Qp;
                if(Zcand.Eta()>0) Q = Qn;

                float denom = Q.Rho2();
                float num = (candQ*TComplex::Conjugate(Q)).Re();
                v2Num[k]->Fill(0.5,num);
                v2Denom[k]->Fill(0.5,denom);

                v2NumVsCent->Fill(k,num);
                v2DenomVsCent->Fill(k,denom);

                //electrons
                TComplex Q1 = Qp;
                if(eleEta->at(goodElectrons.at(j))>0) Q1 = Qn; 
                float numEle1 = (ele1Q*TComplex::Conjugate(Q1)).Re();
                v2EleNum[k]->Fill(0.5,numEle1);
                v2EleDenom[k]->Fill(0.5,Q1.Rho2());
              
                TComplex Q2 = Qp;
                if(eleEta->at(goodElectrons.at(j2))>0) Q2 = Qn;
                float numEle2 = (ele2Q*TComplex::Conjugate(Q2)).Re();
                v2EleNum[k]->Fill(0.5,numEle2);
                v2EleDenom[k]->Fill(0.5,Q2.Rho2());

                v2EleNumVsCent->Fill(k,numEle1);
                v2EleNumVsCent->Fill(k,numEle2);
                v2EleDenomVsCent->Fill(k,Q1.Rho2());
                v2EleDenomVsCent->Fill(k,Q2.Rho2());
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
  
  for(int k = 1; k<nBins+1; k++){
    v2NumVsCentHist->SetBinContent(k,v2NumVsCent->GetBinContent(k));
    v2NumVsCentHist->SetBinError(k,v2NumVsCent->GetBinError(k));

    v2SqrtDenomVsCent->SetBinContent(k,TMath::Sqrt(v2DenomVsCent->GetBinContent(k)));
    v2SqrtDenomVsCent->SetBinError(k,v2DenomVsCent->GetBinError(k)/(2*TMath::Sqrt(v2DenomVsCent->GetBinContent(k))));
    
    //electrons
    v2EleNumVsCentHist->SetBinContent(k,v2EleNumVsCent->GetBinContent(k));
    v2EleNumVsCentHist->SetBinError(k,v2EleNumVsCent->GetBinError(k));

    v2EleSqrtDenomVsCent->SetBinContent(k,TMath::Sqrt(v2EleDenomVsCent->GetBinContent(k)));
    v2EleSqrtDenomVsCent->SetBinError(k,v2EleDenomVsCent->GetBinError(k)/(2*TMath::Sqrt(v2EleDenomVsCent->GetBinContent(k))));
  }

  v2VsCent = (TH1D*)v2NumVsCentHist->Clone("v2VsCent");
  v2VsCent->Divide(v2SqrtDenomVsCent);
  //electrons  
  v2EleVsCent = (TH1D*)v2EleNumVsCentHist->Clone("v2EleVsCent");
  v2EleVsCent->Divide(v2EleSqrtDenomVsCent);

  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->SetDirectory(0);
    massPeakSS[i]->SetDirectory(0);
    v2Num[i]->SetDirectory(0);
    v2Denom[i]->SetDirectory(0);
    v2EleNum[i]->SetDirectory(0);
    v2EleDenom[i]->SetDirectory(0);
  }

  v2DenomVsCent->SetDirectory(0);
  v2SqrtDenomVsCent->SetDirectory(0);
  v2NumVsCentHist->SetDirectory(0);
  v2NumVsCent->SetDirectory(0);
  v2VsCent->SetDirectory(0);
  
  v2EleDenomVsCent->SetDirectory(0);
  v2EleSqrtDenomVsCent->SetDirectory(0);
  v2EleNumVsCentHist->SetDirectory(0);
  v2EleNumVsCent->SetDirectory(0);
  v2EleVsCent->SetDirectory(0);
  
  TFile * output = new TFile("Z2ee.root","recreate");
  for(int i = 0; i<nBins; i++){
    massPeakOS[i]->Write();
    massPeakSS[i]->Write();
    v2Num[i]->Write();
    v2Denom[i]->Write();
    v2EleNum[i]->Write();
    v2EleDenom[i]->Write();
  }
  v2NumVsCent->Write();
  v2DenomVsCent->Write();
  v2SqrtDenomVsCent->Write();
  v2NumVsCentHist->Write();
  v2VsCent->Write();
  v2EleNumVsCent->Write();
  v2EleDenomVsCent->Write();
  v2EleSqrtDenomVsCent->Write();
  v2EleNumVsCentHist->Write();
  v2EleVsCent->Write();
  
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
