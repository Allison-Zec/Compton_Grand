#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "../buildGrandRootfile.h"
#include "makePlots.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

using namespace std;

void snailwisePlots(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TString fname = Form("%sGrandCompton.root", expt.Data());
  const Int_t nPolPlots = 6;
  const Int_t nStdPlots = 6;
  Int_t nMsmt = 0;
  TString polVar[nPolPlots] = {"Pol0", "Asym0", "Asym0LasOff", "Pol0", "Asym0", "Asym0LasOff"};
  TString stdVar[nStdPlots] = {"numRuns", "numCycles", "snailTime", "qw1", "hw1", "qw2"};
  Float_t polYmin[nPolPlots] = {0.99, 0.99, 0.99, 0.9, 0.9, 0.9};
  Float_t polYmax[nPolPlots] = {1.01, 1.01, 1.01, 1.1, 1.1, 1.1};
  Float_t polFactors[nPolPlots] = {100., 1000., 1000., 100., 1000., 1000.};
  Float_t stdYmin[nStdPlots] = {0.0, 0.0, 0.0, 0.9, 0.9, 0.9};
  Float_t stdYmax[nStdPlots] = {1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
  Bool_t polSign[nPolPlots] = {true, true, true, false, false, false};
  Bool_t stdFloat[nStdPlots] = {false, false, true, true, true, true};

  for(Int_t i = 0; i < nPolPlots; i++){
    plotPolSnl(fname.Data(), polVar[i].Data(), nMsmt++, polYmin[i], polYmax[i], polSign[i], false, polFactors[i]);
  }
  plotPolSnl(fname.Data(), "Pol0", nMsmt++, 0.0, 1.1, false, true);
  plotPolSnl(fname.Data(), "Asym0LasOff", nMsmt++, 0.0, 1.1, false, true);
  for(Int_t i = 0; i < nStdPlots; i++){
    plotStandardSnl(fname.Data(), stdVar[i].Data(), nMsmt++, stdYmin[i], stdYmax[i], stdFloat[i]);
  }

  gSystem->Exec(Form("pdfunite %s/plots/msmt*.pdf %s/aggregates/%sGrandSnailwise.pdf", 
                getenv("COMPMON_GRAND"), getenv("COMPMON_WEB"), expt.Data()));
  gSystem->Exec(Form("rm -f %s/plots/msmt*.pdf", getenv("COMPMON_GRAND")));
}
