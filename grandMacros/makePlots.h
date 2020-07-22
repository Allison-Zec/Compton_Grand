#ifndef makePlots_h
#define makePlots_h

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "../vars.h"

TString experimentCode(Int_t prexOrCrex){
  TString expt("prex");
  if(prexOrCrex == 2){
    expt = "crex";
  }
  else if(prexOrCrex != 1 && prexOrCrex != 2){
    printf("Please enter correct parameter for experiment:\n");
    printf("    PREX == 1\n");
    printf("    CREX == 2\n\n");
    exit(1);
  }
  return expt;
}

Int_t getMarkerStyle(Int_t index){
  if(index == 0){return 23;}
  else if(index == 1){return 22;}
  else if(index == 2){return 20;}
  else if(index == 3){return 21;}
  else{
    printf("Invalid graph index (marker)!\n");
    return 0;
  }
}

Int_t getColor(Int_t index){
  if(index == 0){return kViolet;}
  else if(index == 1){return kOrange;}
  else if(index == 2){return kRed;}
  else if(index == 3){return kBlue;}
  else{
    printf("Invalid graph index (color)!\n");
    return 0;
  }
}

TString getLegendEntry(Int_t index){
  if(index == 0){return "Left Out";}
  else if(index == 1){return "Left In";}
  else if(index == 2){return "Right Out";}
  else if(index == 3){return "Right In";}
  else{
    printf("Invalid graph index (legend)!\n");
    return "";
  }
}

void adjustMinMaxs(vector<Float_t> ymins, vector<Float_t> ymaxs, Float_t val, Int_t msmt){
  if(val < ymins[msmt]){ymins[msmt] = val;}
  if(val > ymaxs[msmt]){ymaxs[msmt] = val;}
}

Bool_t isCloseTo(Float_t num, Float_t ref){
  Bool_t val = TMath::Abs(num)>=0.99*TMath::Abs(ref) && TMath::Abs(num)<=1.01*TMath::Abs(ref);
  return val;
}

Bool_t isFlipLeft(Int_t num, Float_t HWienAngle, Float_t VWienAngle, Float_t PhiFG, Bool_t isSnail=true){
  if((num < 99 && isSnail) || (num < 4900 && !isSnail)){
    if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, -89.9008) && isCloseTo(PhiFG, 86.9014))
      return true;
    else if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 91.1902))
      return false;
  }
  else{
    if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, -90.5996) && isCloseTo(PhiFG, 88.0277))
      return true;
    else if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 89.9558))
      return false;
  }
  return false;
}

Int_t getGraphInd(Int_t snailNum, Float_t hWien, Float_t vWien, Float_t solWien, Float_t ihwp, Bool_t isSnail=true){
  if(isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp <= 0.5){return 0;}
  else if(isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp > 0.5){return 1;}
  else if(!isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp <= 0.5){return 2;}
  else if(!isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp > 0.5){return 3;}
  else{printf("Invalid graph index (graph)!\n"); return 0;}
}

void adjustGraphLimits(vector<TGraphErrors *> graphs, Float_t ymin, Float_t ymax, Float_t smallFac, Float_t largeFac, Float_t xmin, Float_t xmax){
  Float_t yminFac = smallFac; Float_t ymaxFac = largeFac;
  if(ymin < 0){yminFac = largeFac;}
  if(ymax < 0){ymaxFac = smallFac;}
  for(Int_t i = 0; i < 4; i++){
    graphs[i]->GetXaxis()->SetLimits(xmin, xmax);
    graphs[i]->GetYaxis()->SetRangeUser(yminFac*ymin, ymaxFac*ymax);
  }
}

void drawAndPrintGraphs(TString msmt, Int_t msmtNum, TCanvas *can, vector<TGraphErrors *> graphs){
  can->cd();
  can->SetGridx(); can->SetGridy();
  TLegend *leg = new TLegend(0.9, 0.75, 0.98, 0.9, "", "NDC");
  for(Int_t i = 0; i < 4; i++){
    leg->AddEntry(graphs[i], getLegendEntry(i).Data());
    if(i == 0){graphs[i]->Draw("ap");}
    else{graphs[i]->Draw("p && same");}
  }
  leg->Draw("same");
  can->Print(Form("%s/plots/msmt%04i_%s.pdf", getenv("COMPMON_GRAND"), msmtNum, msmt.Data()), "pdf");
}

void plotPolSnl(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=true, bool chi2=false){
  FitPolVar var;
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  TFile *grand = new TFile(Form("%s/%s", getenv("COMPMON_GRAND"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("snl");
  TString chiID(""); TString signID("");
  if(chi2) chiID = "chi2";
  if(signCorr) signID = "sign";
  TCanvas *can = new TCanvas(Form("can%s%s_%s", signID.Data(), chiID.Data(), msmt.Data()), "Some Title", 1200, 400);
  TString titleAdd("");
  if(signCorr) titleAdd = " (Sign Corrected)";
  TString chiAdd("");
  if(chi2) chiAdd = " (Chi2 / NDF) ";

  Int_t snailNum, sign;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("snailNum", &snailNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress(msmt.Data(), &var);
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s%s%s vs Snail", msmt.Data(), chiAdd.Data(), titleAdd.Data()));
    g->GetXaxis()->SetTitle("snailNum"); g->GetYaxis()->SetTitle(msmt.Data());
    if(chi2){g->GetYaxis()->SetTitle(Form("%s (Chi2 / NDF)", msmt.Data()));}
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(snailNum, hWien, vWien, solWien, ihwp, true);
    Float_t polVar = var.mean; Float_t polErr = var.meanErr;
    if(chi2){polVar = var.Chi2*1.0/var.NDF; polErr = 0.0;}
    if(signCorr){polVar = var.mean*sign; polErr = var.meanErr*sign;}
    graphs[ind]->SetPoint(graphCounts[ind], snailNum, polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, TMath::Abs(polErr));
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    if(snailNum - 1 < xmin){xmin = snailNum - 1;}
    if(snailNum + 1 > xmax){xmax = snailNum + 1;}
    //printf("<tr class=\"myRow\"><td class=\"myCell\">%i</td><td class=\"myCell\">%.4f +/ %.4f</td></tr>\n", snailNum, 1000*polVar, 1000*polErr);
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotStandardSnl(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool isFloat=true){
  Int_t snailNum, sign, intVar;
  Float_t hWien, vWien, solWien, ihwp, fltVar;

  TFile *grand = new TFile(Form("%s/%s", getenv("COMPMON_GRAND"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("snl");
  TCanvas *can = new TCanvas(Form("can_%s", msmt.Data()), "", 1200, 400); 
  vector<TGraphErrors *> graphs; vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  tree->SetBranchAddress("snailNum", &snailNum);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("sign", &sign);
  if(isFloat)
    tree->SetBranchAddress(msmt.Data(), &fltVar);
  else
    tree->SetBranchAddress(msmt.Data(), &intVar);
  
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs Snail", msmt.Data()));
    g->GetXaxis()->SetTitle("snailNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(snailNum, hWien, vWien, solWien, ihwp, true);
    if(isFloat){
      graphs[ind]->SetPoint(graphCounts[ind], snailNum, fltVar);
      if(fltVar < ymin){ymin = fltVar;}
      if(fltVar > ymax){ymax = fltVar;}
    }
    else{
      graphs[ind]->SetPoint(graphCounts[ind], snailNum, intVar);
      if(((Float_t)intVar) < ymin){ymin = (Float_t)intVar;}
      if(((Float_t)intVar) > ymax){ymax = (Float_t)intVar;}
    }
    if(snailNum - 1 < xmin){xmin = snailNum - 1;}
    if(snailNum + 1 > xmax){xmax = snailNum + 1;}
    graphCounts[ind] += 1;
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotPolRun(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=true, bool chi2=false){
  FitPolVar var;
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  TFile *grand = new TFile(Form("%s/%s", getenv("COMPMON_GRAND"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("run");
  TString chiID(""); TString signID("");
  if(chi2) chiID = "chi2";
  if(signCorr) signID = "sign";
  TString canName = Form("can%s%s_%s", signID.Data(), chiID.Data(), msmt.Data());
  TCanvas *can = new TCanvas(canName.Data(), "Some Title", 1200, 400);
  TString titleAdd(""); TString chiAdd("");
  if(signCorr) titleAdd = " (Sign Corrected)";
  if(chi2) chiAdd = " (Chi2 / NDF) ";

  Int_t runNum, sign;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress(msmt.Data(), &var);
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s%s%s vs Run", msmt.Data(), chiAdd.Data(), titleAdd.Data()));
    g->GetXaxis()->SetTitle("runNum"); g->GetYaxis()->SetTitle(msmt.Data());
    if(chi2) g->GetYaxis()->SetTitle(Form("%s (Chi2 / NDF)", msmt.Data()));
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    Float_t polVar = var.mean; Float_t polErr = var.meanErr;
    if(chi2){polVar = var.Chi2*1.0/var.NDF; polErr = 0.0;}
    if(signCorr){polVar = var.mean*sign; polErr = var.meanErr*sign;}
    if(TMath::IsNaN(polVar)) continue;
    graphs[ind]->SetPoint(graphCounts[ind], runNum, polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, TMath::Abs(polErr));
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    if(runNum - 5 < xmin){xmin = runNum - 5;}
    if(runNum + 5 > xmax){xmax = runNum + 5;}
    //printf("<tr class=\"myRow\"><td class=\"myCell\">%i</td><td class=\"myCell\">%.4f +/ %.4f</td></tr>\n", snailNum, 1000*polVar, 1000*polErr);
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotStandardRun(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool isFloat=true){
  Int_t runNum, sign, intVar;
  Float_t hWien, vWien, solWien, ihwp, fltVar;

  TFile *grand = new TFile(Form("%s/%s", getenv("COMPMON_GRAND"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("run");
  TCanvas *can = new TCanvas(Form("can_%s", msmt.Data()), "", 1200, 400); 
  vector<TGraphErrors *> graphs; vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("sign", &sign);
  if(isFloat)
    tree->SetBranchAddress(msmt.Data(), &fltVar);
  else
    tree->SetBranchAddress(msmt.Data(), &intVar);
  
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs Run", msmt.Data()));
    g->GetXaxis()->SetTitle("runNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    if(isFloat){
      graphs[ind]->SetPoint(graphCounts[ind], runNum, fltVar);
      if(fltVar < ymin){ymin = fltVar;}
      if(fltVar > ymax){ymax = fltVar;}
    }
    else{
      graphs[ind]->SetPoint(graphCounts[ind], runNum, intVar);
      if(((Float_t)intVar) < ymin){ymin = (Float_t)intVar;}
      if(((Float_t)intVar) > ymax){ymax = (Float_t)intVar;}
    }
    if(runNum - 5 < xmin){xmin = runNum - 5;}
    if(runNum + 5 > xmax){xmax = runNum + 5;}
    graphCounts[ind] += 1;
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

#endif
