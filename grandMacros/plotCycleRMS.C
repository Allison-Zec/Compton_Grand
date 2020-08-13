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

using namespace std;

void plotCycleRMS(Int_t prexOrCrex){
  TString exptStr("");
  TString exptName("");
  if(prexOrCrex == 1){
    exptStr = "prex";
    exptName = "PREX-II";
  }
  else if(prexOrCrex == 2){
    exptStr = "crex";
    exptName = "CREX";
  }
  else{
    printf("Invalid Experiment Code\n");
    exit(1);
  }
  
  TFile *f = TFile::Open(Form("%s/%sComptonGrand.root", getenv("COMPMON_GRAND"), exptStr.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t runNum, cycNum;
  DataVar fOff, lOff;
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Acc0LasOff1", &fOff);
  cyc->SetBranchAddress("Acc0LasOff2", %lOff);
  
  TGraphErrors *hF = new TH1F("hF", Form("%s Acc0 RMS by Cycle", exptName.Data()),
                              (Int_t)cyc->GetEntries(), 0, (Int_t)cyc->GetEntries());
  TGraphErrors *hL = new TH1F("hL", Form("%s Acc0 RMS by Cycle", exptName.Data()),
                              (Int_t)cyc->GetEntries(), 0, (Int_t)cyc->GetEntries());

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    gF->SetPoint(i, i+1, fOff.rms);
    gF->SetPointError(i, 0, fOff.rmsErr);
    gL->SetPoint(i, i+1, lOff.rms);
    gL->SetPointError(i, 0, lOff.rmsErr);
  }

  TCanvas *cRMS = new TCanvas("cRMS", "RMS Graph", 1200, 800);
  cRMS->cd();
  hF->GetXaxis()->SetTitle("cycle num");
  hL->GetXaxis()->SetTitle("cycle num");
}
