#include "TF1.h"
#include "TH1D.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"

using namespace std;

//**********************************************************
// Trigger-Combine the data in order to unfold properly later
//**********************************************************

//[0] = Jet20, [1] = Jet40, [2] = Jet60, [3] = Jet80
double trigComb(bool *trg, double *pscl, double pt, int combinationMethod){
  //int combinationMethod = 4;
    double weight=0;

    //HIN-12-017 (charged part RpA) combination method - solid but loses events that slip through the pt bins
    if(combinationMethod==1){
      if(trg[0] && pt<40) weight = pscl[0];
    if(trg[1] && pt>40 && pt<60) weight = pscl[1];
    if(trg[2] && pt>60 && pt<75) weight = pscl[2];
    if(trg[3] && pt>75 && pt<95) weight = pscl[3];
    if(trg[4] && pt>95 && pt<120) weight = pscl[4];
    if(trg[5] && pt>120) weight = pscl[5];
    }

    //Experimental method from STAR - requires trigger decision pt, though!
    if(combinationMethod==2){
      
        if(trg[0] && pt<20) weight = 1./(1./pscl[0]);
        if((trg[0] || trg[1]) && pt>=20 && pt<40) weight = 1./(1./pscl[0] + 1./pscl[1] - 1./(pscl[0]*pscl[1]));
        if((trg[0] || trg[1] || trg[2]) && pt>=40 && pt<60) weight = 1./(1./pscl[0] + 1./pscl[1] + 1./pscl[2] - 1./(pscl[0]*pscl[1]) - 1./(pscl[0]*pscl[2]) - 1./(pscl[1]*pscl[2]) + 1./(pscl[0]*pscl[1]*pscl[2]));
        if((trg[0] || trg[1] || trg[2] || trg[3]) && pt>=60 && pt<80) weight = 1./(1./pscl[0] + 1./pscl[1] + 1./pscl[2] + 1./pscl[3] - 1./(pscl[0]*pscl[1]) - 1./(pscl[0]*pscl[2]) - 1./(pscl[0]*pscl[3]) - 1./(pscl[1]*pscl[2]) - 1./(pscl[1]*pscl[3]) - 1./(pscl[2]*pscl[3]) + 1./(pscl[0]*pscl[1]*pscl[2]) + 1./(pscl[0]*pscl[1]*pscl[3]) + 1./(pscl[0]*pscl[2]*pscl[3]) + 1./(pscl[1]*pscl[2]*pscl[3]) - 1./(pscl[0]*pscl[1]*pscl[2]*pscl[3]));
        if((trg[0] || trg[1] || trg[2] || trg[3] || trg[4]) && pt>=80 && pt<100) weight = 1./(1./pscl[0] + 1./pscl[1] + 1./pscl[2] + 1./pscl[3] + 1./pscl[4] - 1./(pscl[0]*pscl[1]) - 1./(pscl[0]*pscl[2]) - 1./(pscl[0]*pscl[3]) - 1./(pscl[0]*pscl[4]) - 1./(pscl[1]*pscl[2]) - 1./(pscl[1]*pscl[3]) - 1./(pscl[1]*pscl[4])- 1./(pscl[2]*pscl[3]) - 1./(pscl[2]*pscl[4]) - 1./(pscl[3]*pscl[4]) + 1./(pscl[0]*pscl[1]*pscl[2]) + 1./(pscl[0]*pscl[1]*pscl[3]) + 1./(pscl[0]*pscl[1]*pscl[4]) + 1./(pscl[0]*pscl[2]*pscl[3]) + 1./(pscl[0]*pscl[2]*pscl[4]) + 1./(pscl[0]*pscl[3]*pscl[4])+ 1./(pscl[1]*pscl[2]*pscl[3]) + 1./(pscl[1]*pscl[2]*pscl[4]) + 1./(pscl[1]*pscl[3]*pscl[4]) + 1./(pscl[2]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[2]*pscl[3]) - 1./(pscl[0]*pscl[2]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[2]*pscl[4]) - 1./(pscl[1]*pscl[2]*pscl[3]*pscl[4]) + 1./(pscl[0]*pscl[1]*pscl[2]*pscl[3]*pscl[4]));
	if((trg[0] || trg[1] || trg[2] || trg[3] || trg[4] || trg[5]) && pt>=100) weight = 1.;
    }

    //Workaround method for HIN-12-003 (bjets).  Only works with 3 triggers and tends to lose events and assumes Jet80 is unprescaled (false for pA).  For cross-check only.
    if(combinationMethod==3){
    // if(triggerDecision[0] && !triggerDecision[1] && !triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[0]); //Removing finnicky Jet20 sample
        if(trg[1] && !trg[2] && !trg[3]) weight = pscl[1]; //1./(1./pscl[1]);
        if(trg[2] && !trg[3]) weight = 1.; //1./(1./pscl[1] + 1./pscl[2] - (1./(pscl[1]*pscl[2])));
        if(trg[3]) weight = 1.;
    }

    //Totally new experimental method
    if(combinationMethod==4){
      if(trg[5]) weight = 1.;
      if(trg[4] && !trg[5]) weight = 1.;
      if(trg[3] && !trg[4] && !trg[5]) weight = 1.;
      if(trg[2] && !trg[3] && !trg[4] && !trg[5]) weight = 1.;
      if(trg[1] && !trg[2] && !trg[3] && !trg[4] && !trg[5]) weight = 1.;
      if(trg[0] && !trg[1] && !trg[2] && !trg[3] && !trg[4] && !trg[5]) weight = 1./(1./pscl[0]);
    }

    if(combinationMethod==5){
      if((trg[0] || trg[1] || trg[2] || trg[3] || trg[4])) weight = 1./(1./pscl[0] + 1./pscl[1] + 1./pscl[2] + 1./pscl[3] + 1./pscl[4] - 1./(pscl[0]*pscl[1]) - 1./(pscl[0]*pscl[2]) - 1./(pscl[0]*pscl[3]) - 1./(pscl[0]*pscl[4]) - 1./(pscl[1]*pscl[2]) - 1./(pscl[1]*pscl[3]) - 1./(pscl[1]*pscl[4])- 1./(pscl[2]*pscl[3]) - 1./(pscl[2]*pscl[4]) - 1./(pscl[3]*pscl[4]) + 1./(pscl[0]*pscl[1]*pscl[2]) + 1./(pscl[0]*pscl[1]*pscl[3]) + 1./(pscl[0]*pscl[1]*pscl[4]) + 1./(pscl[0]*pscl[2]*pscl[3]) + 1./(pscl[0]*pscl[2]*pscl[4]) + 1./(pscl[0]*pscl[3]*pscl[4])+ 1./(pscl[1]*pscl[2]*pscl[3]) + 1./(pscl[1]*pscl[2]*pscl[4]) + 1./(pscl[1]*pscl[3]*pscl[4]) + 1./(pscl[2]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[2]*pscl[3]) - 1./(pscl[0]*pscl[2]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[3]*pscl[4]) - 1./(pscl[0]*pscl[1]*pscl[2]*pscl[4]) - 1./(pscl[1]*pscl[2]*pscl[3]*pscl[4]) + 1./(pscl[0]*pscl[1]*pscl[2]*pscl[3]*pscl[4]));
    }
    
    if(combinationMethod==6){
      if(trg[5] && pt>=100) weight = pscl[5];
      if(!trg[5] && trg[4] && pt>=80 && pt<100) weight = pscl[4];
      if(!trg[5] && !trg[4] && trg[3] && pt>=60 && pt<80) weight = pscl[3];
      if(!trg[5] && !trg[4] && !trg[3] && trg[2] && pt>=40 && pt<60) weight = pscl[2];
      if(!trg[5] && !trg[4] && !trg[3] && !trg[2] && trg[1] && pt>=20 && pt<40) weight = pscl[1];
      if(!trg[5] && !trg[4] && !trg[3] && !trg[2] && !trg[1] && trg[0] && pt>=0 && pt<20) weight = pscl[0];
    }

    if(combinationMethod==7){
      if(trg[5] && pt>=100) weight = pscl[5];
      if(trg[4] && pt>=80 && pt<100) weight = pscl[4];
      if(trg[3] && pt>=60 && pt<80) weight = pscl[3];
      if(trg[2] && pt>=40 && pt<60) weight = pscl[2];
      if(trg[1] && pt>=20 && pt<40) weight = pscl[1];
      if(trg[0] && pt>=0 && pt<20) weight = pscl[0];
    }

    return weight;
}


void trgCombToy(int nEvents=50000000, bool testRTbenefit=true){
  
  //needed to create different random no's for different runs!
  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  TStopwatch t1;
  TStopwatch t2_1;
  TStopwatch t2_2;
  TStopwatch t2;
  TStopwatch t3;
  TStopwatch t4;
  TStopwatch tpt;

  const int nTrgs = 6;

  double trgs[nTrgs] = {0,20,40,60,80,100};
  double pscls[nTrgs] = {100.,16.,8.,4.,2.,1.};
  TRandom3* r1 = new TRandom3(0);
  r1->SetSeed(0);
  TF1 *f1 = new TF1("f1","x*TMath::Power((TMath::Exp(-0.001*x-0.001*x*x)+x/3),-3.1)",0,300);
  TH1D *totalData = new TH1D("totalData","",500,0,300);
  TH1D *method1 = new TH1D("method1","",500,0,300); method1->Sumw2();
  TH1D *method2 = new TH1D("method2","",500,0,300); method2->Sumw2();
  TH1D *method3 = new TH1D("method3","",500,0,300); method3->Sumw2();
  TH1D *method4 = new TH1D("method4","",500,0,300); method4->Sumw2();
  TH1D *method5 = new TH1D("method5","",500,0,300); method5->Sumw2();

  TH1D *hTrg40 = new TH1D("hTrg40","",500,0,300); hTrg40->Sumw2();
  totalData->Sumw2();
  
  bool trgDecisions[nTrgs] = {0,0,0,0,0,0};
  int trgCounter[nTrgs] = {0,0,0,0,0,0};

  double w1=0, w2=0, w3=0, w4=0, w5=0, w6=0;

  //Function for realistic number of events as a function of pt
  TF1* fNjetDistr = new TF1("fNjetDistr","-0.022 + 0.685*TMath::Log(x)",0,300);
  TF1* NjetSigma = new TF1("NjetSigma","gaus(0)",-4,4);
  TF1 *NjetSigmaErr = new TF1("NjetSigmaErr","-0.00357142 + 0.00167857*x",0,300);
  NjetSigma->SetParameter(0,1);
  NjetSigma->SetParameter(1,0);
  TF1 *trigObjPtSmear = new TF1("trigObjPtSmear","gaus(0)",-1,1);
  trigObjPtSmear->SetParameter(0,1);
  trigObjPtSmear->SetParameter(1,0);
  trigObjPtSmear->SetParameter(2,0.2); //15 GeV/c smear between trigger pt and actual pt is sort of conservative...

  //first = reco jet pt.... second = trigger object pt (smeared by realistic gaussian)
  pair<double,double> jetPts[10];
  t2.Reset();
  t1.Start();


  for(int i=0; i<nEvents; i++){
    if(i%(nEvents/10)==0 && i>0) cout << (float)i/(float)nEvents*100. << "% complete..." << endl;
    //jetPts = {0,0,0,0,0,0,0,0,0,0};
    tpt.Start(0);
    double pt= f1->GetRandom(3,300);
    tpt.Stop();
    jetPts[0].first = pt;
    jetPts[0].second = pt*(1+trigObjPtSmear->GetRandom());
    unsigned int nJets;
    t2.Start(0);
    if(testRTbenefit){
      //Smear number of jets in each event by a realistic distribution based on pt
      NjetSigma->SetParameter(2, TMath::Power(NjetSigmaErr->Eval(pt),2));
      t2_2.Start(0);
      double tmp=0;
      if(pt>150) tmp = NjetSigma->GetRandom(); //Warning! Tuning less than 150 makes this perform very poorly!!
      t2_2.Stop();
      nJets = TMath::Nint(fNjetDistr->Eval(pt)+tmp);
      if(nJets<1) nJets = 1;
    }
    else nJets=1;
    t2.Stop();

    //start from 1, because we've already pushed back one jet
    for(unsigned int ijet=1; ijet<nJets; ijet++){
      double ptsub = f1->GetRandom(0.1,pt); //0,pt
      jetPts[ijet].first = ptsub;
      jetPts[ijet].second = ptsub*(1+trigObjPtSmear->GetRandom());
    }
    
    for(unsigned int ijet=0; ijet<nJets; ijet++){
      totalData->Fill(jetPts[ijet].second);
    }

    t3.Start(0);
    //Now check all the trigger decisions for the various jets in the collection.
    
    //THIS IS THE CORRECT METHOD FOR RANDOMIZED TRIGGERS!
    for(int j=0; j<nTrgs; j++){
      double rnum = r1->Rndm();
      for(unsigned int ijet=0; ijet<nJets; ijet++){
	if(jetPts[ijet].first >= trgs[j] && rnum <= 1./pscls[j]){
	  trgDecisions[j] = 1;
	}
      }
    }

    //THIS IS FOR UNRANDOMIZED TRIGGERS!!
    /*bool trgChecker[nTrgs] = {0,0,0,0,0,0};
    for(int j=0; j<nTrgs; j++){
      for(unsigned int ijet=0; ijet<nJets; ijet++){
	if(jetPts[ijet].first >= trgs[j]){
	  //cout << "jet pt: "<< jetPts[ijet].first << " passes trg: " << trgs[j] << endl;
	  if(trgChecker[j]==0){
	    trgChecker[j] = 1;
	    trgCounter[j]++;
	  }
	}
      }
    }
    for(int j=0; j<nTrgs; j++){
      if(trgCounter[j]%(int)pscls[j]==0 && trgCounter[j]!=0){
	trgDecisions[j]=1;
	trgCounter[j]=0;
      }
      }*/
    //END TRIGGER LOGIC!

    t3.Stop();
    
    if(!trgDecisions[0] && !trgDecisions[1] && !trgDecisions[2] && !trgDecisions[3] && !trgDecisions[4] && !trgDecisions[5]){
      continue;
    }

    t4.Start(0);
    //do the weighting considering all the trigger decisions
    //  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for(unsigned int ijet=0; ijet<nJets; ijet++){
      double trigpt = jetPts[0].first;
      w1 = trigComb(trgDecisions, pscls, trigpt, 1);
      w2 = trigComb(trgDecisions, pscls, trigpt, 2);
      w3 = trigComb(trgDecisions, pscls, trigpt, 3);
      w4 = trigComb(trgDecisions, pscls, trigpt, 4);
      w6 = trigComb(trgDecisions, pscls, trigpt, 7);
 
      //  cout << "jet " << ijet << " recoPt " << jetPts[ijet].first << ", triggerPt: "<< jetPts[ijet].second << endl;
      //  cout << "weight1: "<< w1 << " weight2: "<< w2 << " weight3: "<< w3 << " weight4: "<< w4<< endl;

      //Fill all the jets using the event weight
      method1->Fill(jetPts[ijet].second,w1);
      method2->Fill(jetPts[ijet].second,w2);
      method3->Fill(jetPts[ijet].second,w3);
      method4->Fill(jetPts[ijet].second,w4);
      method5->Fill(jetPts[ijet].second,w6);

      if(trgDecisions[1] && jetPts[ijet].first>40) hTrg40->Fill(jetPts[ijet].second);
    }
    t4.Stop();
    
    // for(int itrig=0; itrig<nTrgs; itrig++){
    //  cout << "trigger["<<itrig<<"] decision: " << trgDecisions[itrig] << endl;
    //}
    // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    //reset prescale factors for next event!
    for(int j=0; j<nTrgs; j++){
      trgDecisions[j] = 0;
    }
  }
  t1.Stop();

  //Diagnostics
  cout << "total time per event: "<< t1.RealTime()/(float)nEvents << endl;
  cout << "to initialize RT benefit: "<< t2.RealTime()/(float)nEvents << endl;
  cout << "RT benefit step 1: "<< t2_1.RealTime()/(float)nEvents << endl;
  cout << "RT benefit step 2: "<< t2_2.RealTime()/(float)nEvents << endl;
  cout << "to check trg decisions: "<< t3.RealTime()/(float)nEvents << endl;
  cout << "to weight triggers: " << t4.RealTime()/(float)nEvents << endl;
  cout << "to f1->GetRandom(): "<< tpt.RealTime()/(float)nEvents << endl;

  TFile *fout = new TFile("trgCombTest.root","recreate");
  fout->cd();
  totalData->Write();
  method1->Write();
  method2->Write();
  method3->Write();
  method4->Write();
  method5->Write();
  
  hTrg40->Write();
}
