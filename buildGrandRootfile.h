#ifndef buildGrandRootfile_h
#define buildGrandRootfile_h

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "vars.h"

// cyc tree
// basics
Int_t snailNum, runNum, cycleNum;
Int_t firstOffStartMPS, firstOffEndMPS, onStartMPS, onEndMPS, lastOffStartMPS, lastOffEndMPS;
// mpswise
//vector<vector<TString>> cycMPSVars;
vector<DataVar> cycMPSData;
// quartetwise
//vector<vector<TString>> cycQrtVars;
vector<DataVar> cycQrtData;
//triggerwise
//Float_t meanPedF, meanErrPedF, rmsPedF, rmsErrPedF;
//Float_t meanPedL, meanErrPedL, rmsPedL, rmsErrPedL;
DataVar firstOffPedestal, lastOffPedestal;
// calculations
Float_t cycleTime;
vector<TString> cycQrtPols;
vector<PolVar> cycQrtCalc;
Float_t meanDiff0LasOn, meanErrDiff0LasOn, meanSum0LasOn, meanErrSum0LasOn;
Float_t meanDiff0LasOff, meanErrDiff0LasOff, meanDiff0LasOff1, meanErrDiff0LasOff1, meanDiff0LasOff2, meanErrDiff0LasOff2;
Float_t meanSum0LasOff, meanErrSum0LasOff, meanSum0LasOff1, meanErrSum0LasOff1, meanSum0LasOff2, meanErrSum0LasOff2;
Float_t meanDiff4LasOn, meanErrDiff4LasOn, meanSum4LasOn, meanErrSum4LasOn;
Float_t meanDiff4LasOff, meanErrDiff4LasOff, meanDiff4LasOff1, meanErrDiff4LasOff1, meanDiff4LasOff2, meanErrDiff4LasOff2;
Float_t meanSum4LasOff, meanErrSum4LasOff, meanSum4LasOff1, meanErrSum4LasOff1, meanSum4LasOff2, meanErrSum4LasOff2;
Float_t anPow;
PolVar bkSubSum0, bkSubSum4;
PolVar bkSubAsym0LasOn, bkSubAsym0LasOff, bkSubAsym0LasOff1, bkSubAsym0LasOff2, asym0;
PolVar bkSubAsym4LasOn, bkSubAsym4LasOff, bkSubAsym4LasOff1, bkSubAsym4LasOff2, asym4;
PolVar asym0LasOn, asym0LasOff, asym0LasOff1, asym0LasOff2;
PolVar asym4LasOn, asym4LasOff, asym4LasOff1, asym4LasOff2;
PolVar pol0, pol4;

// run tree
// basics
Int_t numRunCycles;
// mpswise
//vector<vector<TString>> runMPSVars;
vector<DataVar> runMPSData;
// epicswise
//vector<vector<TString>> runEpcVars;
vector<StdVar> runEpcData;
// calculations
Float_t runTime;
Int_t runSign;
vector<Float_t> runPol0Avg, runPol0Err;
vector<Float_t> runAsym0Avg, runAsym0Err;
vector<Float_t> runAsym0OffAvg, runAsym0OffErr;
FitPolVar runAsym0, runAsym0Off, runPol0;

// snl tree
// basics
Int_t numSnlCycles;
Int_t numSnlRuns;
//calculations
Float_t snailTime;
Int_t snailSign;
//vector<Float_t> snlPol0Avg, snlPol0Err, snlPol4Avg, snlPol4Err;
//vector<Float_t> snlAsym0Avg, snlAsym0Err, snlAsym4Avg, snlAsym4Err;
//FitPolVar snlAsym0, snlAsym4, snlPol0, snlPol4;
vector<Float_t> snlPol0Avg, snlPol0Err;
vector<Float_t> snlAsym0Avg, snlAsym0Err;
vector<Float_t> snlAsym0OffAvg, snlAsym0OffErr;
FitPolVar snlAsym0, snlAsym0Off, snlPol0;
vector<Float_t> qw1Lst, hw1Lst, qw2Lst, ihwpLst, VWienAngleLst, HWienAngleLst, PhiFGLst;
Float_t qw1, hw1, qw2, ihwp, VWienAngle, HWienAngle, PhiFG;

Int_t helicityFreq(Int_t runNum){
  if(runNum > 4242 && runNum < 4621)
    return 240;
  else
    return 120;
}

Double_t getAnalyzingPower(int run_num){
  if(run_num < 3800){return 0.1105;}
  else if(run_num >= 4232 && run_num <= 4929){return 0.01655915;}
  else if(run_num >= 4930){return 0.036017934;}
  else{
    printf("Run doesn't have a defined analyzing power.\n");
    return 0.0;
  }
}

vector<vector<TString>> readVarsFile(TString outTree, TString inTree){
  ifstream infile(Form("%s/grandConfigs/%s_%s_vars.csv", getenv("COMPMON_GRAND"), outTree.Data(), inTree.Data()));
  string readStr;
  vector<vector<TString>> vars;
  while(getline(infile, readStr)){
    vector<TString> oneVar; stringstream ss(readStr);
    string token;
    while(getline(ss, token, ',')){
      oneVar.push_back(token.c_str());
    }
    vars.push_back(oneVar);
  }
  return vars;
}

vector<TString> readPolsFile(TString outTree, TString inTree){
  ifstream infile(Form("%s/grandConfigs/%s_%s_pol.csv", getenv("COMPMON_GRAND"), outTree.Data(), inTree.Data()));
  string readStr;
  vector<TString> vars;
  while(getline(infile, readStr)){
    //printf("Read from file: %s\n", readStr.c_str());
    vars.push_back(readStr.c_str());
  }
  return vars;
}

vector<vector<int>> findCycles(int runNum){
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), runNum));
  string readStr;
  int cyclesInThisRun = 0;
  vector<vector<int>> run;
  while(getline(infile, readStr)){
    vector<int> limits; stringstream ss(readStr);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    cyclesInThisRun++;
    vector<int> cycle;
    for(int i = 0; i < limits.size(); i++)
      cycle.push_back(limits[i]);
    run.push_back(cycle);
  }
  printf("Looked in %s/cycles_%i.dat\n", getenv("COMPMON_MINIRUNS"), runNum);
  printf("Found %i cycles\n", (Int_t)run.size());
  return run;
}

#endif
