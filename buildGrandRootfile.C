#include "buildGrandRootfile.h"
#include "plot.h"
#include "../CompMon/runs.h"

using namespace std;

void readKeysFile(TString expt){
  ifstream keysfile(Form("%s/%s_rmsCut.key", getenv("COMPMON_MAPS"), expt.Data()));
  if(keysfile.is_open()){
    string line;
    while(getline(keysfile, line)){
      vector<Float_t> keyRange; stringstream ss(line);
      for(Float_t i; ss >> i;){
        keyRange.push_back(i);
        if(ss.peek() == ',')
          ss.ignore();
      }
      keys.push_back(keyRange);
    }
  }
  keysfile.close();
}

Bool_t acceptCycle(Int_t runNum){
  if(keys.size() == 0){
    acc0OnMode = 0.0; acc0OnOffset = 0.0;
    acc0OffMode = 0.0; acc0OffOffset = 0.0;
    return true;
  }
  else{
    for(Int_t i = 0; i < keys.size(); i++){
      if(runNum >= (Int_t)keys[i][0] && runNum <= (Int_t)keys[i][1]){
        acc0OnMode = keys[i][3]; acc0OnOffset = keys[i][4];
        acc0OffMode = keys[i][2]; acc0OffOffset = keys[i][4];
        Float_t offLimit = keys[i][2] + keys[i][4];
        Float_t onLimit = keys[i][3] + keys[i][4];
        return cycMPSData[1].rms < offLimit && cycMPSData[2].rms < offLimit && cycMPSData[0].rms < onLimit;
      }
    }
    return true;
  }
}

void calcCyclePol(Int_t runNum){
  setAnalyzingPower(runNum);

  Int_t diff0On = 0; Int_t sum0On = 1; Int_t diff0Off1 = 2; Int_t sum0Off1 = 3; Int_t diff0Off2 = 4; Int_t sum0Off2 = 5;
  Int_t diff4On = 6; Int_t sum4On = 7; Int_t diff4Off1 = 8; Int_t sum4Off1 = 9; Int_t diff4Off2 = 10; Int_t sum4Off2 = 11;

  meanDiff0LasOn = cycQrtCalc[diff0On].mean; meanErrDiff0LasOn = cycQrtCalc[diff0On].meanErr; 
  meanSum0LasOn = cycQrtCalc[sum0On].mean; meanErrSum0LasOn = cycQrtCalc[sum0On].meanErr;
  meanDiff0LasOff1 = cycQrtCalc[diff0Off1].mean; meanErrDiff0LasOff1 = cycQrtCalc[diff0Off1].meanErr; 
  meanSum0LasOff1 = cycQrtCalc[sum0Off1].mean; meanErrSum0LasOff1 = cycQrtCalc[sum0Off1].meanErr; 
  meanDiff0LasOff2 = cycQrtCalc[diff0Off2].mean; meanErrDiff0LasOff2 = cycQrtCalc[diff0Off2].meanErr;
  meanSum0LasOff2 = cycQrtCalc[sum0Off2].mean; meanErrSum0LasOff2 = cycQrtCalc[sum0Off2].meanErr;
  meanDiff0LasOff = (meanDiff0LasOff1 + meanDiff0LasOff2)/2.0;
  meanErrDiff0LasOff = TMath::Sqrt(TMath::Power(meanErrDiff0LasOff1, 2) + TMath::Power(meanErrDiff0LasOff2, 2))/2.0;
  meanSum0LasOff = (meanSum0LasOff1 + meanSum0LasOff2)/2.;
  meanErrSum0LasOff = TMath::Sqrt(TMath::Power(meanErrSum0LasOff1, 2) + TMath::Power(meanErrSum0LasOff2, 2))/2.0;
  bkSubSum0.mean = meanSum0LasOn - meanSum0LasOff;
  bkSubSum0.meanErr = TMath::Sqrt(TMath::Power(meanErrSum0LasOn, 2) + TMath::Power(meanErrSum0LasOff, 2));
  bkSubAsym0LasOn.mean = meanDiff0LasOn/bkSubSum0.mean;
  bkSubAsym0LasOn.meanErr = TMath::Abs(bkSubAsym0LasOn.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOn/meanDiff0LasOn, 2) + TMath::Power(bkSubSum0.meanErr/bkSubSum0.mean, 2));
  bkSubAsym0LasOff.mean = meanDiff0LasOff/bkSubSum0.mean;
  bkSubAsym0LasOff.meanErr = TMath::Abs(bkSubAsym0LasOff.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff/meanDiff0LasOff, 2) + TMath::Power(bkSubSum0.meanErr/bkSubSum0.mean, 2));
  bkSubAsym0LasOff1.mean = meanDiff0LasOff1/bkSubSum0.mean;
  bkSubAsym0LasOff1.meanErr = TMath::Abs(bkSubAsym0LasOff1.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff1/meanDiff0LasOff1, 2) + TMath::Power(bkSubSum0.meanErr/bkSubSum0.mean, 2));
  bkSubAsym0LasOff2.mean = meanDiff0LasOff2/bkSubSum0.mean;
  bkSubAsym0LasOff2.meanErr = TMath::Abs(bkSubAsym0LasOff2.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff2/meanDiff0LasOff2, 2) + TMath::Power(bkSubSum0.meanErr/bkSubSum0.mean, 2));
  //asym0.mean = asym0LasOn.mean - asym0LasOff.mean;
  //asym0.meanErr = TMath::Sqrt(TMath::Power(asym0LasOn.meanErr, 2) + TMath::Power(asym0LasOff.meanErr, 2));
  asym0LasOff1.mean = bkSubAsym0LasOff1.mean; asym0LasOff1.meanErr = bkSubAsym0LasOff1.meanErr;
  asym0LasOff2.mean = bkSubAsym0LasOff2.mean; asym0LasOff2.meanErr = bkSubAsym0LasOff2.meanErr;
  asym0LasOff.mean = (asym0LasOff1.mean + asym0LasOff2.mean)/2.0;
  asym0LasOff.meanErr = TMath::Sqrt(TMath::Power(asym0LasOff1.meanErr, 2) + TMath::Power(asym0LasOff2.meanErr, 2))/2.0;
  asym0.mean = asym0LasOn.mean - asym0LasOff.mean; 
  asym0.meanErr = TMath::Sqrt(TMath::Power(asym0LasOn.meanErr, 2) + TMath::Power(asym0LasOff.meanErr, 2));
  pol0.mean = asym0.mean/anPow.mean;
  pol0.meanErr = TMath::Abs(pol0.mean)*TMath::Sqrt(TMath::Power(asym0.meanErr/asym0.mean, 2) + TMath::Power(anPow.meanErr/anPow.mean, 2));

  meanDiff4LasOn = cycQrtCalc[diff4On].mean; meanErrDiff4LasOn = cycQrtCalc[diff4On].meanErr; 
  meanSum4LasOn = cycQrtCalc[sum4On].mean; meanErrSum4LasOn = cycQrtCalc[sum4On].meanErr;
  meanDiff4LasOff1 = cycQrtCalc[diff4Off1].mean; meanErrDiff4LasOff = cycQrtCalc[diff4Off1].meanErr; 
  meanSum4LasOff1 = cycQrtCalc[sum4Off1].mean; meanErrSum4LasOff1 = cycQrtCalc[sum4Off1].meanErr;
  meanDiff4LasOff2 = cycQrtCalc[diff4Off2].mean; meanErrDiff4LasOff = cycQrtCalc[diff4Off2].meanErr;
  meanSum4LasOff2 = cycQrtCalc[sum4Off2].mean; meanErrSum4LasOff2 = cycQrtCalc[sum4Off2].meanErr;
  meanDiff4LasOff = (meanDiff4LasOff1 + meanDiff4LasOff2)/2.0;
  meanErrDiff4LasOff = TMath::Sqrt(TMath::Power(meanErrDiff4LasOff1, 2) + TMath::Power(meanErrDiff4LasOff2, 2))/2.0;
  meanSum4LasOff = (meanSum4LasOff1 + meanSum4LasOff2)/2.;
  meanErrSum4LasOff = TMath::Sqrt(TMath::Power(meanErrSum4LasOff1, 2) + TMath::Power(meanErrSum4LasOff2, 2))/2.0;
  bkSubSum4.mean = meanSum4LasOn - meanSum4LasOff;
  bkSubSum4.meanErr = TMath::Sqrt(TMath::Power(meanErrSum4LasOn, 2) + TMath::Power(meanErrSum4LasOff, 2));
  bkSubAsym4LasOn.mean = meanDiff4LasOn/bkSubSum4.mean;
  bkSubAsym4LasOn.meanErr = TMath::Abs(bkSubAsym4LasOn.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOn/meanDiff4LasOn, 2) + TMath::Power(bkSubSum4.meanErr/bkSubSum4.mean, 2));
  bkSubAsym4LasOff.mean = meanDiff4LasOff/bkSubSum4.mean;
  bkSubAsym4LasOff.meanErr = TMath::Abs(bkSubAsym4LasOff.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff/meanDiff4LasOff, 2) + TMath::Power(bkSubSum4.meanErr/bkSubSum4.mean, 2));
  bkSubAsym4LasOff1.mean = meanDiff4LasOff1/bkSubSum4.mean;
  bkSubAsym4LasOff1.meanErr = TMath::Abs(bkSubAsym4LasOff1.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff1/meanDiff4LasOff1, 2) + TMath::Power(bkSubSum4.meanErr/bkSubSum4.mean, 2));
  bkSubAsym4LasOff2.mean = meanDiff4LasOff2/bkSubSum4.mean;
  bkSubAsym4LasOff2.meanErr = TMath::Abs(bkSubAsym4LasOff2.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff2/meanDiff4LasOff2, 2) + TMath::Power(bkSubSum4.meanErr/bkSubSum4.mean, 2));
  //asym4.mean = asym4LasOn.mean - asym4LasOff.mean;
  //asym4.meanErr = TMath::Sqrt(TMath::Power(asym4LasOn.meanErr, 2) + TMath::Power(asym4LasOff.meanErr, 2));
  asym4LasOff1.mean = bkSubAsym4LasOff1.mean; asym4LasOff1.meanErr = bkSubAsym4LasOff1.meanErr;
  asym4LasOff2.mean = bkSubAsym4LasOff2.mean; asym4LasOff2.meanErr = bkSubAsym4LasOff2.meanErr;
  asym4LasOff.mean = (asym4LasOff1.mean + asym4LasOff2.mean)/2.0;
  asym4LasOff.meanErr = TMath::Sqrt(TMath::Power(asym4LasOff1.meanErr, 2) + TMath::Power(asym4LasOff2.meanErr, 2))/2.0;
  asym4.mean = asym4LasOn.mean - asym4LasOff.mean; 
  asym4.meanErr = TMath::Sqrt(TMath::Power(asym4LasOn.meanErr, 2) + TMath::Power(asym4LasOff.meanErr, 2));
  pol4.mean = asym4.mean/anPow.mean;
  pol4.meanErr = TMath::Abs(pol4.mean)*TMath::Sqrt(TMath::Power(asym4.meanErr/asym4.mean, 2) + TMath::Power(anPow.meanErr/anPow.mean, 2));

  if(acceptCycle(runNum)){
    runAsym0Avg.push_back(asym0.mean); runAsym0Err.push_back(asym0.meanErr);
    runAsym0OffAvg.push_back(bkSubAsym0LasOff.mean); runAsym0OffErr.push_back(bkSubAsym0LasOff.meanErr);
    runPol0Avg.push_back(pol0.mean); runPol0Err.push_back(pol0.meanErr);
    snlAsym0Avg.push_back(asym0.mean); snlAsym0Err.push_back(asym0.meanErr);
    snlAsym0OffAvg.push_back(bkSubAsym0LasOff.mean); snlAsym0OffErr.push_back(bkSubAsym0LasOff.meanErr);
    //snlAsym4Avg.push_back(asym4.mean); snlAsym4Err.push_back(asym4.meanErr);
    snlPol0Avg.push_back(pol0.mean); snlPol0Err.push_back(pol0.meanErr);
    //snlPol4Avg.push_back(pol4.mean); snlPol4Err.push_back(pol4.meanErr);
    cycleCut = 0;
  }
  else{
    printf("      Rejecting cycle for polarization!\n");
    cycleCut = 1;
  }
}

void calcCyclePedestals(TFile *plotfile){
  TString hNameF = Form("h%i.%i_pedF", runNum, cycleNum);
  TString hNameL = Form("h%i.%i_pedL", runNum, cycleNum);
  TString fNameF1 = Form("h%i.%i_pedF1", runNum, cycleNum);
  TString fNameF2 = Form("h%i.%i_pedF2", runNum, cycleNum);
  TString fNameL1 = Form("h%i.%i_pedL1", runNum, cycleNum);
  TString fNameL2 = Form("h%i.%i_pedL2", runNum, cycleNum);

  TH1F *hF = (TH1F *)plotfile->Get(hNameF.Data());
  TH1F *hL = (TH1F *)plotfile->Get(hNameL.Data());

  TF1 *fF1 = new TF1(fNameF1.Data(), "gaus");
  TF1 *fF2 = new TF1(fNameF2.Data(), "gaus");
  TF1 *fL1 = new TF1(fNameL1.Data(), "gaus");
  TF1 *fL2 = new TF1(fNameL2.Data(), "gaus");

  Float_t meanF1 = hF->GetMean(); Float_t rmsF1 = hF->GetRMS();
  Float_t meanL1 = hL->GetMean(); Float_t rmsL1 = hL->GetRMS();
  hF->Fit(fNameF1.Data(), "Q", "goff", meanF1 - rmsF1, meanF1 + rmsF1);
  hL->Fit(fNameL1.Data(), "Q", "goff", meanL1 - rmsL1, meanL1 + rmsL1);
  Float_t meanF2 = fF1->GetParameter(1); Float_t rmsF2 = fF1->GetParameter(2);
  Float_t meanL2 = fL1->GetParameter(1); Float_t rmsL2 = fF1->GetParameter(2);
  hF->Fit(fNameF2.Data(), "Q", "goff", meanF2 - rmsF2, meanF2 + rmsF2);
  hL->Fit(fNameL2.Data(), "Q", "goff", meanL2 - rmsL2, meanL2 + rmsL2);
  firstOffPedestal.mean = fF2->GetParameter(1); firstOffPedestal.meanErr = fF2->GetParError(1);
  firstOffPedestal.rms  = fF2->GetParameter(2); firstOffPedestal.rmsErr  = fF2->GetParError(2);
  lastOffPedestal.mean  = fL2->GetParameter(1); lastOffPedestal.meanErr  = fL2->GetParError(1);
  lastOffPedestal.rms   = fL2->GetParameter(2); lastOffPedestal.rmsErr   = fL2->GetParError(2);
}

Bool_t isCloseTo(Float_t num, Float_t ref){
  Bool_t val = TMath::Abs(num)>=0.99*TMath::Abs(ref) && TMath::Abs(num)<=1.01*TMath::Abs(ref);
  //printf("Close to checks: Num: %f, Ref: %f, Result: %s\n", num, ref, val ? "true" : "false");
  return val;
}

void calcSnailSign(Int_t snailNum){
  Int_t ihwpSign = 1 - 2*ihwp;
  Int_t wienSign = 0;
  if(snailNum < 99){
    if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, -89.9008) && isCloseTo(PhiFG, 86.9014))
      wienSign = 1;
    else if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 91.1902))
      wienSign = -1;
  }
  else{
    if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, -90.5996) && isCloseTo(PhiFG, 88.0277))
      wienSign = 1;
    else if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 89.9558))
      wienSign = -1;
  }
  snailSign = ihwpSign*wienSign;
  //printf("Snail sign: %i\n", snailSign);
}

void calcRunSign(Int_t runNum, Float_t runHWien, Float_t runVWien, Float_t runSolWien, Int_t runIHWP){
  Int_t ihwpSign = 1 - 2*runIHWP;
  Int_t wienSign = 0;
  if(runNum < 4800){
    if(isCloseTo(runHWien, -13.0) && isCloseTo(runVWien, -89.9008) && isCloseTo(runSolWien, 86.9014))
      wienSign = 1;
    else if(isCloseTo(runHWien, -13.0) && isCloseTo(runVWien, 88.0008) && isCloseTo(runSolWien, 91.1902))
      wienSign = -1;
  }
  else{
    if(isCloseTo(runHWien, -29.6402) && isCloseTo(runVWien, -90.5996) && isCloseTo(runSolWien, 88.0277))
      wienSign = 1;
    else if(isCloseTo(runHWien, -29.6402) && isCloseTo(runVWien, 88.0008) && isCloseTo(runSolWien, 89.9558))
      wienSign = -1;
  }
  runSign = ihwpSign*wienSign;
  //printf("Snail sign: %i\n", snailSign);
}

void initCycleTree(TTree *cyc){
  cyc->Branch("snailNum", &snailNum, "snailNum/I");
  cyc->Branch("runNum", &runNum, "runNum/I");
  cyc->Branch("cycleNum", &cycleNum, "cycleNum/I");
  cyc->Branch("firstOffStartMPS", &firstOffStartMPS, "firstOffStartMPS/I");
  cyc->Branch("firstOffEndMPS", &firstOffEndMPS, "firstOffEndMPS/I");
  cyc->Branch("onStartMPS", &onStartMPS, "onStartMPS/I");
  cyc->Branch("onEndMPS", &onEndMPS, "onEndMPS/I");
  cyc->Branch("lastOffStartMPS", &lastOffStartMPS, "lastOffStartMPS/I");
  cyc->Branch("lastOffEndMPS", &lastOffEndMPS, "lastOffEndMPS/I");
  cyc->Branch("cycleTime", &cycleTime, "cycleTime/F");
  cyc->Branch("sign", &cycSign, "sign/I");
  cyc->Branch("ihwp", &cycIHWP, "ihwp/F");
  cyc->Branch("VWienAngle", &cycVWien, "VWienAngle/F");
  cyc->Branch("HWienAngle", &cycHWien, "HWienAngle/F");
  cyc->Branch("PhiFG", &cycSolWien, "PhiFG/F");
  for(Int_t i = 0; i < cycMPSVars; i++){DataVar data; cycMPSData.push_back(data);}
  for(Int_t i = 0; i < cycMPSVars; i++){cyc->Branch(cycMPSTitles[i].Data(), &cycMPSData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  //cyc->Branch("PedestalMeanFirstOff", &firstOffPedestal, "mean/F:meanErr/F:rms/F:rmsErr/F");
  //cyc->Branch("PedestalMeanLastOff", &lastOffPedestal, "mean/F:meanErr/F:rms/F:rmsErr/F");
  for(Int_t i = 0; i < cycQrtVars; i++){DataVar data; cycQrtData.push_back(data);}
  for(Int_t i = 0; i < cycQrtVars; i++){cyc->Branch(cycQrtTitles[i].Data(), &cycQrtData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  cyc->Branch("AnalyzingPower", &anPow, "mean/F:meanErr/F");
  cyc->Branch("BackSubSum0", &bkSubSum0, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOn", &asym0LasOn, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff", &asym0LasOff, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff1", &asym0LasOff1, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff2", &asym0LasOff2, "mean/F:meanErr/F");
  cyc->Branch("Asym0", &asym0, "mean/F:meanErr/F");
  cyc->Branch("Pol0", &pol0, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym0LasOn", &bkSubAsym0LasOn, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym0LasOff1", &bkSubAsym0LasOff1, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym0LasOff2", &bkSubAsym0LasOff2, "mean/F:meanErr/F");
  cyc->Branch("BackSubSum4", &bkSubSum4, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOn", &asym4LasOn, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff", &asym4LasOff, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff1", &asym4LasOff1, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff2", &asym4LasOff2, "mean/F:meanErr/F");
  cyc->Branch("Asym4", &asym4, "mean/F:meanErr/F");
  cyc->Branch("Pol4", &pol4, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym4LasOn", &bkSubAsym4LasOn, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym4LasOff1", &bkSubAsym4LasOff1, "mean/F:meanErr/F");
  cyc->Branch("BackSubAsym4LasOff2", &bkSubAsym4LasOff2, "mean/F:meanErr/F");
  cyc->Branch("CycleCut", &cycleCut, "CycleCut/I");
  cyc->Branch("Acc0OnMode", &acc0OnMode, "Acc0OnMode/F");
  cyc->Branch("Acc0OnOffset", &acc0OnOffset, "Acc0OnOffset/F");
  cyc->Branch("Acc0OffMode", &acc0OffMode, "Acc0OffMode/F");
  cyc->Branch("Acc0OffOffset", &acc0OffOffset, "Acc0OffOffset/F");
}

void initRunTree(TTree *run){
  run->Branch("snailNum", &snailNum, "snailNum/I");
  run->Branch("runNum", &runNum, "runNum/I");
  run->Branch("numCycles", &numRunCycles, "numCycles/I");
  run->Branch("runTime", &runTime);
  for(Int_t i = 0; i < runMPSVars; i++){DataVar data; runMPSData.push_back(data);}
  for(Int_t i = 0; i < runMPSVars; i++){run->Branch(runMPSTitles[i].Data(), &runMPSData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  for(Int_t i = 0; i < runEpcVars; i++){StdVar data; runEpcData.push_back(data);}
  for(Int_t i = 0; i < runEpcVars; i++){run->Branch(runEpcTitles[i].Data(), &runEpcData[i], "mean/F");}
  run->Branch("sign", &runSign, "sign/I");
  run->Branch("Asym0", &runAsym0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  run->Branch("Asym0LasOff", &runAsym0Off, "mean/F:meanErr/F:Chi2/F:NDF/I");
  run->Branch("Pol0", &runPol0, "mean/F:meanErr/F:Chi2/F:NDF/I");
}

void initSnailTree(TTree *snl){
  snl->Branch("snailNum", &snailNum, "snailNum/I");
  snl->Branch("numRuns", &numSnlRuns, "numRuns/I");
  snl->Branch("numCycles", &numSnlCycles, "numCycles/I");
  snl->Branch("snailTime", &snailTime);
  snl->Branch("Asym0", &snlAsym0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  snl->Branch("Asym0LasOff", &snlAsym0Off, "mean/F:meanErr/F:Chi2/F:NDF/I");
  snl->Branch("Pol0", &snlPol0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Asym4", &snlAsym4, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Pol4", &snlPol4, "mean/F:meanErr/F:Chi2/F:NDF/I");
  snl->Branch("qw1", &qw1, "qw1/F");
  snl->Branch("hw1", &hw1, "hw1/F");
  snl->Branch("qw2", &qw2, "qw2/F");
  snl->Branch("ihwp", &ihwp, "ihwp/F");
  snl->Branch("VWienAngle", &VWienAngle, "VWienAngle/F");
  snl->Branch("HWienAngle", &HWienAngle, "HWienAngle/F");
  snl->Branch("PhiFG", &PhiFG, "PhiFG/F");
  snl->Branch("sign", &snailSign, "sign/I");
}

/**
  Called before run and cycle loops
**/
void snailIterSet(Int_t base, Int_t snlInd){
  printf("Analyzing snail %i...\n", base + snlInd);
  snailNum = base + snlInd;
  numSnlCycles = 0;
  numSnlRuns = 0;
  snailTime = 0;
  snailSign = 0;
  snlAsym0Avg.clear(); snlAsym0Err.clear();
  snlAsym0OffAvg.clear(); snlAsym0OffErr.clear();  
  snlPol0Avg.clear(); snlPol0Err.clear();
  //snlAsym4Avg.clear(); snlAsym4Err.clear();
  //snlPol4Err.clear(); snlPol4Err.clear();
  qw1 = 0; hw1 = 0; qw2 = 0; ihwp = 0; HWienAngle = 0; VWienAngle = 0; PhiFG = 0;
}

void snailIterAfterSet(Int_t base, Int_t snlInd){
  printf("Finalizing snail %i...\n", base + snlInd);
  TH1F *hAsym0 = new TH1F(Form("hAsym0_%i", base + snlInd), "Asym0", numSnlCycles, 0, numSnlCycles);
  TH1F *hAsym0Off = new TH1F(Form("hAsym0Off_%i", base + snlInd), "Asym0 Off", numSnlCycles, 0, numSnlCycles);
  TH1F *hPol0 = new TH1F(Form("hPol0_%i", base + snlInd), "Pol0", numSnlCycles, 0, numSnlCycles);
  //TH1F *hAsym4 = new TH1F(Form("hAsym4_%i", base + snlInd), "Asym4", numSnlCycles, 0, numSnlCycles);
  //TH1F *hPol4 = new TH1F(Form("hPol4_%i", base + snlInd), "Pol4", numSnlCycles, 0, numSnlCycles);
  TF1 *fAsym0 = new TF1(Form("fAsym0_%i", base + snlInd), "pol0");
  TF1 *fAsym0Off = new TF1(Form("fAsym0_%i", base + snlInd), "pol0");
  TF1 *fPol0 = new TF1(Form("fPol0_%i", base + snlInd), "pol0");
  //TF1 *fAsym4 = new TF1(Form("fAsym4_%i", base + snlInd), "pol0");
  //TF1 *fPol4 = new TF1(Form("fPol4_%i", base + snlInd), "pol0");
  for(Int_t i = 0; i < snlAsym0Avg.size(); i++){
    hAsym0->SetBinContent(i + 1, snlAsym0Avg[i]); hAsym0->SetBinError(i + 1, snlAsym0Err[i]);
    hAsym0Off->SetBinContent(i + 1, snlAsym0OffAvg[i]); hAsym0Off->SetBinError(i + 1, snlAsym0OffErr[i]);
    hPol0->SetBinContent(i + 1, snlPol0Avg[i]); hPol0->SetBinError(i + 1, snlPol0Err[i]);
    //hAsym4->SetBinContent(i + 1, snlAsym4Avg[i]); hAsym4->SetBinError(i + 1, snlAsym4Err[i]);
    //hPol4->SetBinContent(i + 1, snlPol4Avg[i]); hPol4->SetBinError(i + 1, snlPol4Err[i]);
    //printf("  Asyms, Cycle %i: %f +/- %f, %f +/- %f\n", i + 1, snlAsym0Avg[i], snlAsym0Err[i], snlAsym4Avg[i], snlAsym4Err[i]);
    //printf("  Pol, Cycle %i: %f +/- %f, %f +/- %f\n", i + 1, snlPol0Avg[i], snlPol0Err[i], snlPol4Avg[i], snlPol4Err[i]);
    //printf("  Asyms, Cycle %i: %f +/- %f\n", i + 1, snlAsym0Avg[i], snlAsym0Err[i]);
    //printf("  Pol, Cycle %i: %f +/- %f\n", i + 1, snlPol0Avg[i], snlPol0Err[i]);
  }
  hAsym0->Fit(fAsym0, "Q", "", 0, numSnlCycles);
  hAsym0Off->Fit(fAsym0Off, "Q", "", 0, numSnlCycles);
  hPol0->Fit(fPol0, "Q", "", 0, numSnlCycles);
  snlAsym0.mean = fAsym0->GetParameter(0); snlAsym0.meanErr = fAsym0->GetParError(0);
  snlAsym0.Chi2 = fAsym0->GetChisquare(); snlAsym0.NDF = fAsym0->GetNDF();
  snlAsym0Off.mean = fAsym0Off->GetParameter(0); snlAsym0Off.meanErr = fAsym0Off->GetParError(0);
  snlAsym0Off.Chi2 = fAsym0Off->GetChisquare(); snlAsym0Off.NDF = fAsym0Off->GetNDF();
  snlPol0.mean = fPol0->GetParameter(0); snlPol0.meanErr = fPol0->GetParError(0);
  snlPol0.Chi2 = fPol0->GetChisquare(); snlPol0.NDF = fPol0->GetNDF();
  //hAsym4->Fit(fAsym4, "Q", "", 0, numSnlCycles);
  //hPol4->Fit(fPol4, "Q", "", 0, numSnlCycles);
  //snlAsym4.mean = fAsym4->GetParameter(0); snlAsym4.meanErr = fAsym4->GetParError(0);
  //snlAsym4.Chi2 = fAsym4->GetChisquare(); snlAsym4.NDF = fAsym4->GetNDF();
  //snlPol4.mean = fPol4->GetParameter(0); snlPol4.meanErr = fPol4->GetParError(0);
  //snlPol4.Chi2 = fPol4->GetChisquare(); snlPol4.NDF = fPol4->GetNDF();

  //printf("Sign Vars snail %i: IHWP %f, VWien: %f, HWien: %f, PhiFG: %f; Num Runs: %i\n", base + snlInd, ihwp, VWienAngle, HWienAngle, PhiFG, numSnlRuns);

  qw1 /= numSnlRuns; hw1 /= numSnlRuns; qw2 /= numSnlRuns; ihwp /= numSnlRuns;
  HWienAngle /= numSnlRuns; VWienAngle /= numSnlRuns; PhiFG /= numSnlRuns;
  calcSnailSign(base + snlInd);
}

/**
  Called before cycle loop
**/
void runIterSet(Int_t runNumber, Int_t numCycles, TFile *plotFile){
  printf("  Analyzing run %i...\n", runNumber);
  runAsym0Avg.clear(); runAsym0Err.clear();
  runAsym0OffAvg.clear(); runAsym0OffErr.clear();
  runPol0Avg.clear(); runPol0Err.clear();
  runNum = runNumber;
  numRunCycles = numCycles;
  numSnlCycles += numCycles;
  numSnlRuns += 1;
  runSign = 0;
  //runTime = 1.0*mpswise->GetEntries()/helicityFreq(runNumber);
  for(int i = 0; i < runMPSVars; i++){
    TString hName = Form("h%i_%s", runNumber, runMPSTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    runMPSData[i].mean = h->GetMean(); runMPSData[i].meanErr = h->GetMeanError();
    runMPSData[i].rms = h->GetRMS(); runMPSData[i].rmsErr = h->GetRMSError();
  }
  for(int i = 0; i < runEpcVars; i++){
    TString hName = Form("h%i_%s", runNumber, runEpcTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    runEpcData[i].mean = h->GetMean();
  }
  runTime = 1.0*((TH1F *)plotFile->Get(Form("h%i_mps", runNumber)))->GetEntries()/helicityFreq(runNumber);
  snailTime += runTime;
  Float_t ihwpAvg = ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "ihwp")))->GetMean();
  Int_t ihwpVar;
  if(ihwpAvg > 0.5) ihwpVar = 1;
  else ihwpVar = 0;
  qw1 += ((TH1F *)plotFile->Get(Form("h%i_qw1", runNumber)))->GetMean();
  hw1 += ((TH1F *)plotFile->Get(Form("h%i_hw1", runNumber)))->GetMean();
  qw2 += ((TH1F *)plotFile->Get(Form("h%i_qw2", runNumber)))->GetMean();
  ihwp += ihwpVar;
  Float_t hWienMean = ((TH1F *)plotFile->Get(Form("h%i_HWienAngle", runNumber)))->GetMean();
  Float_t vWienMean = ((TH1F *)plotFile->Get(Form("h%i_VWienAngle", runNumber)))->GetMean();
  Float_t solWienMean = ((TH1F *)plotFile->Get(Form("h%i_PhiFG", runNumber)))->GetMean();
  HWienAngle += hWienMean;
  VWienAngle += vWienMean;
  PhiFG += solWienMean;
  //printf("  Sign settings: IHWP: %i, VWien: %f, HWien: %f, Sol: %f\n", 
  //       ihwpVar, ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "VWienAngle")))->GetMean(), 
  //       ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "HWienAngle")))->GetMean(), ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "PhiFG")))->GetMean());
  calcRunSign(runNumber, hWienMean, vWienMean, solWienMean, ihwpVar);
}

void runIterAfterSet(Int_t runNumber){
  printf("  Finalizing run %i...\n", runNumber);
  TH1F *hAsym0 = new TH1F(Form("hAsym0_%i", runNumber), "Asym0", numRunCycles, 0, numRunCycles);
  TH1F *hAsym0Off = new TH1F(Form("hAsym0Off_%i", runNumber), "Asym0 Off", numRunCycles, 0, numRunCycles);
  TH1F *hPol0 = new TH1F(Form("hPol0_%i", runNumber), "Pol0", numRunCycles, 0, numRunCycles);
  TF1 *fAsym0 = new TF1(Form("fAsym0_%i", runNumber), "pol0");
  TF1 *fAsym0Off = new TF1(Form("fAsym0Off_%i", runNumber), "pol0");
  TF1 *fPol0 = new TF1(Form("fPol0_%i", runNumber), "pol0");
  for(Int_t i = 0; i < runAsym0Avg.size(); i++){
    hAsym0->SetBinContent(i + 1, runAsym0Avg[i]); hAsym0->SetBinError(i + 1, runAsym0Err[i]);
    hAsym0Off->SetBinContent(i + 1, runAsym0OffAvg[i]); hAsym0Off->SetBinError(i + 1, runAsym0OffErr[i]);
    hPol0->SetBinContent(i + 1, runPol0Avg[i]); hPol0->SetBinError(i + 1, runPol0Err[i]);
  }
  hAsym0->Fit(fAsym0, "Q", "", 0, numRunCycles);
  hAsym0Off->Fit(fAsym0Off, "Q", "", 0, numRunCycles);
  hPol0->Fit(fPol0, "Q", "", 0, numRunCycles);
  runAsym0.mean = fAsym0->GetParameter(0); runAsym0.meanErr = fAsym0->GetParError(0);
  runAsym0.Chi2 = fAsym0->GetChisquare(); runAsym0.NDF = fAsym0->GetNDF();
  runAsym0Off.mean = fAsym0Off->GetParameter(0); runAsym0Off.meanErr = fAsym0Off->GetParError(0);
  runAsym0Off.Chi2 = fAsym0Off->GetChisquare(); runAsym0Off.NDF = fAsym0Off->GetNDF();
  runPol0.mean = fPol0->GetParameter(0); runPol0.meanErr = fPol0->GetParError(0);
  runPol0.Chi2 = fPol0->GetChisquare(); runPol0.NDF = fPol0->GetNDF();
}

/**
  Main body of cycle loop
**/
void cycleIterSet(vector<vector<int>> cycles, Int_t cycInd, Int_t runNumber, TFile *plotFile){
  printf("    Analyzing cycle %i...\n", cycInd + 1);
  cycleNum = cycInd + 1;
  firstOffStartMPS = cycles[cycInd][0];
  firstOffEndMPS = cycles[cycInd][1];
  onStartMPS = cycles[cycInd][2];
  onEndMPS = cycles[cycInd][3];
  lastOffStartMPS = cycles[cycInd][4];
  lastOffEndMPS = cycles[cycInd][5];
  Int_t totMPS = (firstOffEndMPS - firstOffStartMPS) + (onEndMPS - onStartMPS) + (lastOffEndMPS - lastOffStartMPS);
  cycleTime = 1.0*totMPS/helicityFreq(runNumber);
  cycSign = runSign;
  cycIHWP = runEpcData[5].mean;
  cycHWien = runEpcData[8].mean;
  cycVWien = runEpcData[7].mean;
  cycSolWien = runEpcData[9].mean;
  cycQrtCalc.clear();
  for(Int_t i = 0; i < cycMPSVars; i++){
    TString hName = Form("h%i.%i_%s", runNumber, cycleNum, cycMPSTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    cycMPSData[i].mean = h->GetMean(); cycMPSData[i].meanErr = h->GetMeanError();
    cycMPSData[i].rms = h->GetRMS(); cycMPSData[i].rmsErr = h->GetRMSError();
  }
  for(Int_t i = 0; i < cycQrtVars; i++){
    TString hName = Form("h%i.%i_%s", runNumber, cycleNum, cycQrtTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    cycQrtData[i].mean = h->GetMean(); cycQrtData[i].meanErr = h->GetMeanError();
    cycQrtData[i].rms = h->GetRMS(); cycQrtData[i].rmsErr = h->GetRMSError();
    for(Int_t j = 0; j < cycPolVars; j++){
      if(cycQrtTitles[i].CompareTo(cycPolTitles[j]) == 0){
        PolVar pol; pol.mean = h->GetMean(); pol.meanErr = h->GetMeanError();
        cycQrtCalc.push_back(pol);
      }
    }
  }
  TString hAsym0Name = Form("h%i.%i_AsymAcc0LasOn", runNumber, cycleNum);
  TH1F *h0 = (TH1F *)plotFile->Get(hAsym0Name.Data());
  asym0LasOn.mean = h0->GetMean(); asym0LasOn.meanErr = h0->GetMeanError();
  TString hAsym4Name = Form("h%i.%i_AsymAcc4LasOn", runNumber, cycleNum);
  TH1F *h4 = (TH1F *)plotFile->Get(hAsym4Name.Data());
  asym4LasOn.mean = h4->GetMean(); asym4LasOn.meanErr = h4->GetMeanError();
  calcCyclePol(runNumber);
  //calcCyclePedestals(plotFile);
}

/**
  Main function
**/
void buildGrandRootfile(Int_t prexOrCrex){
  TString expt("prex");
  if(prexOrCrex == 2){
    expt = "crex";
  }
  else if(prexOrCrex != 1 && prexOrCrex != 2){
    printf("Please enter correct parameter for experiment:\n");
    printf("    PREX == 1\n");
    printf("    CREX == 2\n\n");
    printf("Get it right this time, %s.\n", randInsult());
  }
  readKeysFile(expt);

  vector<vector<int>> runList = productionRunList(prexOrCrex);
  TFile *out = new TFile(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()), "RECREATE");
  TTree *cyc = new TTree("cyc", "cycle testing tree");
  TTree *run = new TTree("run", "run number testing tree");
  TTree *snl = new TTree("snl", "snail number testing tree");

  //initVars();
  initCycleTree(cyc);
  initRunTree(run);
  initSnailTree(snl);

  Int_t base = 100*(prexOrCrex - 1) + 1;
  for(Int_t s = 0; s < runList.size(); s++){
    //if(s < 29 || s > 39) continue;
    snailIterSet(base, s);
    for(Int_t r = 0; r < runList[s].size(); r++){
      TFile *plotFile = new TFile(Form("%s/Run%i_Plots.root", getenv("COMPMON_RUNPLOTS"), runList[s][r]), "READ");
      vector<vector<int>> cycles = findCycles(runList[s][r]);
      runIterSet(runList[s][r], (Int_t)cycles.size(), plotFile);
      for(Int_t c = 0; c < cycles.size(); c++){
        cycleIterSet(cycles, c, runList[s][r], plotFile);
        cyc->Fill();
      }
      runIterAfterSet(runList[s][r]);
      run->Fill();
      plotFile->Close();
    }
    snailIterAfterSet(base, s);
    snl->Fill();
  }

  out->cd();
  //printf("Alright to here!\n");
  cyc->Write();
  //printf("Still okay by me!\n");
  run->Write();
  //printf("No problems detected!\n");
  snl->Write();
  //printf("We're still going!\n");
  out->Close();
  //printf("Made it all the way!\n");
}
