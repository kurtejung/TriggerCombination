void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.525,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.53);
  c->GetPad(2)->SetBottomMargin(0.2);
  c->GetPad(2)->SetGridy(1);
}
void formatRatioHist(TH1D *h){
  h->GetYaxis()->SetRangeUser(0.92,1.08);
  h->GetXaxis()->SetRangeUser(0,140);
  //h->GetYaxis()->SetNdivisions(505);
  h->SetYTitle("Orig/Combined");
  h->SetXTitle("p_{T} (GeV/c)");
  //h->SetTitle("asdfasdf");
  h->GetYaxis()->CenterTitle(1);
  h->GetYaxis()->SetLabelSize(18);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(14);
  h->GetXaxis()->SetLabelFont(43);
}
void toyModelPlot(){
  TH1D *origpt;
  TH1D *trgComb[5];

  char* histoname = new char[50];
  TFile *f1 = new TFile("trgCombTest.root");
  origpt = (TH1D*)(f1->Get("totalData"))->Clone("totalData");
  for(int i=1; i<6; i++){
  sprintf(histoname,"%s%d","method",i);
    trgComb[i-1] = (TH1D*)(f1->Get(histoname))->Clone(histoname);
    trgComb[i-1]->SetMarkerColor(1);
    trgComb[i-1]->SetMarkerSize(1);
    trgComb[i-1]->SetLineColor(1);
    trgComb[i-1]->SetMarkerStyle(20);
    trgComb[i-1]->SetTitle("");
  }

  TLegend *leg1 = new TLegend(0.453,0.657,0.824,0.855);
  leg1->SetFillColor(0);
  leg1->AddEntry(origpt,"Original Spectrum");
  leg1->AddEntry(trgComb[0],"Combined Spectrum [HIN-12-017]");
  //leg1->AddEntry(trgComb[1],"Combined Spectrum [STAR]");
  leg1->AddEntry(trgComb[4],"Combined Spectrum [New Method]");
  TCanvas *c1 = new TCanvas("c1","",800,600);
  formatCanvas(c1);
  c1->cd(1);
  origpt->SetYTitle("Counts");
  origpt->SetXTitle("p_{T} (GeV/c)");
  origpt->Draw();
  trgComb[0]->SetMarkerColor(2);
  trgComb[0]->SetLineColor(2);
  trgComb[0]->Draw("same");
  trgComb[1]->SetMarkerColor(6);
  trgComb[1]->SetLineColor(6);
  //trgComb[1]->Draw("same");
  trgComb[4]->SetMarkerColor(4);
  trgComb[4]->SetLineColor(4);
  trgComb[4]->Draw("same");
  //cout << "original mean: " << origpt[0]->GetMean() << endl;
  //cout << "scaled mean: "<< smearpt[0]->GetMean() << endl;
  leg1->Draw();
  c1->cd(2);
  TH1D *orgcln = (TH1D*)origpt->Clone("orgcln");
  TH1D *orgcln2 = (TH1D*)origpt->Clone("orgcln2");
  TH1D *orgcln3 = (TH1D*)origpt->Clone("orgcln3");
  orgcln->Divide(orgcln,trgComb[0]);
  //orgcln2->Divide(orgcln2,trgComb[2]);
  orgcln3->Divide(orgcln3,trgComb[4]);
  formatRatioHist(orgcln);
  //formatRatioHist(orgcln2);
  formatRatioHist(orgcln3);
  orgcln->SetMarkerColor(2);
  orgcln2->SetMarkerColor(6);
  orgcln3->SetMarkerColor(4);
  orgcln3->Draw();
  orgcln->Draw("same");
  //orgcln2->Draw("same");
}
