#ifndef plot_h
#define plot_h

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

TString B1L1("beamState==1 && (laserState==0 || laserState==1)");
TString B1L0("beamState==1 && (laserState==2 || laserState==3)");
TString B1("beamState==1");
TString B0("beamState==0");
TString L1("(laserState==0 || laserState==1)");
TString pos0("PosHelAcc0/PosHelNSamples0");
TString neg0("NegHelAcc0/NegHelNSamples0");
TString pos4("PosHelAcc4");
TString neg4("NegHelAcc4");
const Int_t cycMPSVars = 24;
TString cycMPSTitles[cycMPSVars] = {"Acc0LasOn", "Acc0LasOff1", "Acc0LasOff2", "Acc0BeamOff", "Acc4LasOn", "Acc4LasOff1",
                                    "Acc4LasOff2", "Acc4BeamOff", "LaserPower", "BeamCurrent", "bpmAx", "bpmAy", 
                                    "bpmBx", "bpmBy", "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2",
                                    "DSbg1", "DSbg2", "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff"};
TString cycMPSNames[cycMPSVars] = {"Acc0/NAcc0", "Acc0/NAcc0", "Acc0/NAcc0", "Acc0/NAcc0", "Acc4", "Acc4",
                                   "Acc4", "Acc4", "cavPowerCalibrated", "bcm", "bpmAx", "bpmAy",
                                   "bpmBx", "bpmBy", "scaler_ip13", "scaler_ip13", "scaler_run4", "scaler_run13",
                                   "scaler_run0", "scaler_run1", "scaler_run10", "scaler_run10", "scaler_run11", "scaler_run11"};

const Int_t cycQrtVars = 36;
TString cycQrtTitles[cycQrtVars] = {"PosAcc0LasOn",   "NegAcc0LasOn",   "DiffAcc0LasOn",   "SumAcc0LasOn",
                                    "PosAcc0LasOff1", "NegAcc0LasOff1", "DiffAcc0LasOff1", "SumAcc0LasOff1",
                                    "PosAcc0LasOff2", "NegAcc0LasOff2", "DiffAcc0LasOff2", "SumAcc0LasOff2",
                                    "PosAcc0BeamOff", "NegAcc0BeamOff", "DiffAcc0BeamOff", "SumAcc0BeamOff",
                                    "PosAcc4LasOn",   "NegAcc4LasOn",   "DiffAcc4LasOn",   "SumAcc4LasOn",
                                    "PosAcc4LasOff1", "NegAcc4LasOff1", "DiffAcc4LasOff1", "SumAcc4LasOff1",
                                    "PosAcc4LasOff2", "NegAcc4LasOff2", "DiffAcc4LasOff2", "SumAcc4LasOff2",
                                    "PosAcc4BeamOff", "NegAcc4BeamOff", "DiffAcc4BeamOff", "SumAcc4BeamOff",
                                    "diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy"};
TString cycQrtNames[cycQrtVars] =  {pos0, neg0, Form("%s - %s", pos0.Data(), neg0.Data()), Form("%s + %s", pos0.Data(), neg0.Data()),
                                    pos0, neg0, Form("%s - %s", pos0.Data(), neg0.Data()), Form("%s + %s", pos0.Data(), neg0.Data()),
                                    pos0, neg0, Form("%s - %s", pos0.Data(), neg0.Data()), Form("%s + %s", pos0.Data(), neg0.Data()),
                                    pos0, neg0, Form("%s - %s", pos0.Data(), neg0.Data()), Form("%s + %s", pos0.Data(), neg0.Data()),
                                    pos4, neg4, Form("%s - %s", pos4.Data(), neg4.Data()), Form("%s + %s", pos4.Data(), neg4.Data()),
                                    pos4, neg4, Form("%s - %s", pos4.Data(), neg4.Data()), Form("%s + %s", pos4.Data(), neg4.Data()),
                                    pos4, neg4, Form("%s - %s", pos4.Data(), neg4.Data()), Form("%s + %s", pos4.Data(), neg4.Data()),
                                    pos4, neg4, Form("%s - %s", pos4.Data(), neg4.Data()), Form("%s + %s", pos4.Data(), neg4.Data()),
                                    "PosHelBPMAx - NegHelBPMAx", "PosHelBPMAy - NegHelBPMAy", "PosHelBPMBx - NegHelBPMBx", "PosHelBPMBy - NegHelBPMBy"};

const Int_t runMPSVars = 23;
TString runMPSTitles[runMPSVars] =  {"Acc0LasOn", "Acc0LasOff", "Acc0BeamOff", "Acc4LasOn", "Acc4LasOff", "Acc4BeamOff",
                                     "LaserPower", "BeamCurrent", "bpmAx", "bpmAy", "bpmBx", "bpmBy",
                                     "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2", "DSbg1", "DSbg2",
                                     "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff", "mps"};
TString runMPSNames[runMPSVars] = {"Acc0/NAcc0", "Acc0/NAcc0", "Acc0/NAcc0", "Acc4", "Acc4", "Acc4",
                                   "cavPowerCalibrated", "bcm", "bpmAx", "bpmAy", "bpmBx", "bpmBy",
                                   "scaler_ip13", "scaler_ip13", "scaler_run4", "scaler_run13", "scaler_run0", "scaler_run1",
                                   "scaler_run10", "scaler_run10", "scaler_run11", "scaler_run11", "mpsCoda"};

const Int_t runEpcVars = 10;
TString runEpcTitles[runEpcVars] = {"tablePosX", "tablePosY", "qw1", "hw1", "qw2", "ihwp", 
                                    "targetPos", "VWienAngle", "HWienAngle", "PhiFG"};
TString runEpcNames[runEpcVars]  = {"epics_tablePosX", "epics_tablePosY", "epics_qw1", "epics_hw1", "epics_qw2", "epics_ihwp_in", 
                                    "epics_targetPos", "epics_VWienAngle", "epics_HWienAngle", "epics_PhiFG"};
const Int_t runSignVars = 4;
Int_t runSignInds[runSignVars] = {8, 7, 9, 5};

Int_t cycPolVars = 12;
TString cycPolTitles[12] =  {"DiffAcc0LasOn", "SumAcc0LasOn", "DiffAcc0LasOff1", "SumAcc0LasOff1", "DiffAcc0LasOff2", "SumAcc0LasOff2",
                             "DiffAcc4LasOn", "SumAcc4LasOn", "DiffAcc4LasOff1", "SumAcc4LasOff1", "DiffAcc4LasOff2", "SumAcc0LasOff2"};

int translate_month(TString month){
  if(month.EqualTo("Jan")){return 1;}
  else if(month.EqualTo("Feb")){return 2;}
  else if(month.EqualTo("Mar")){return 3;}
  else if(month.EqualTo("Apr")){return 4;}
  else if(month.EqualTo("May")){return 5;}
  else if(month.EqualTo("Jun")){return 6;}
  else if(month.EqualTo("Jul")){return 7;}
  else if(month.EqualTo("Aug")){return 8;}
  else if(month.EqualTo("Sep")){return 9;}
  else if(month.EqualTo("Oct")){return 10;}
  else if(month.EqualTo("Nov")){return 11;}
  else if(month.EqualTo("Dec")){return 12;}
  else{return 0;}
}

vector<vector<int>> readEventCuts(int runNum){
  vector<vector<int>> eventCuts;
  ifstream mapfile(Form("%s/compmon_event_cuts_%i.map", getenv("COMPMON_MAPS"), runNum));
  if(!mapfile.good()){return eventCuts;}
  string readStr;
  while(getline(mapfile, readStr)){
    vector<int> cutRegion; stringstream ss(readStr);
    for(int i; ss >> i;){
      cutRegion.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    eventCuts.push_back(cutRegion);
  }
  mapfile.close();
  printf("Applying eventcuts file for run %i\n", runNum);
  return eventCuts;
}

TString evtCutStr(Int_t runNum){
  TString mapCuts("");
  vector<vector<int>> evtCuts = readEventCuts(runNum);
  for(Int_t i = 0; i < evtCuts.size(); i++){
    if(i == 0)
      mapCuts = Form("!(mpsCoda>%i && mpsCoda<%i)", evtCuts[i][0], evtCuts[i][1]);
    else
      mapCuts = mapCuts + Form(" && !(mpsCoda>%i && mpsCoda<%i)", evtCuts[i][0], evtCuts[i][1]);
  }
  return mapCuts;
}

void cycMPSPlots(TTree *mpswise, Int_t runNum, Int_t cycNum, vector<int> cyc, TFile *runOut){
  Int_t fStart = cyc[0]; Int_t fEnd = cyc[1];
  Int_t oStart = cyc[2]; Int_t oEnd = cyc[3];
  Int_t lStart = cyc[4]; Int_t lEnd = cyc[5];
  TString per1 = Form("mpsCoda>=%i && mpsCoda<=%i", fStart, fEnd);
  TString per2 = Form("mpsCoda>=%i && mpsCoda<=%i", oStart, oEnd);
  TString per3 = Form("mpsCoda>=%i && mpsCoda<=%i", lStart, lEnd);
  TString cycCut = Form("(mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i)", fStart, fEnd, oStart, oEnd, lStart, lEnd);
  
  TString varCuts[cycMPSVars] =  {Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per3.Data()), Form("(%s) && (%s)", B0.Data(), cycCut.Data()),
                                  Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per3.Data()), Form("(%s) && (%s)", B0.Data(), cycCut.Data()),
                                  L1, B1, B1, B1, B1, B1, B1L1, B1L0, B1, B1, B1, B1, B1L1, B1L0, B1L1, B1L0};

  for(Int_t i = 0; i < cycMPSVars; i++){
    TString hName = Form("h%i.%i_%s", runNum, cycNum, cycMPSTitles[i].Data());
    mpswise->Draw(Form("%s>>%s", cycMPSNames[i].Data(), hName.Data()), varCuts[i].Data(), "goff");
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
  }
}

void cycQrtPlots(TTree *quartetwise, Int_t runNum, Int_t cycNum, vector<int> cyc, TFile *runOut){
  Int_t fStart = cyc[0]; Int_t fEnd = cyc[1];
  Int_t oStart = cyc[2]; Int_t oEnd = cyc[3];
  Int_t lStart = cyc[4]; Int_t lEnd = cyc[5];
  TString per1 = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", fStart, fEnd);
  TString per2 = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", oStart, oEnd);
  TString per3 = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", lStart, lEnd);
  TString cycCut = Form("(firstMPSnumber>=%i && firstMPSnumber<=%i) || (firstMPSnumber>=%i && firstMPSnumber<=%i) || (firstMPSnumber>=%i && firstMPSnumber<=%i)", fStart, fEnd, oStart, oEnd, lStart, lEnd);

  TString varCuts[cycQrtVars] =  {Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L1.Data(), per2.Data()), 
                                  Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L1.Data(), per2.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per1.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per1.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per3.Data()), Form("(%s) && (%s)", B1L0.Data(), per3.Data()),
                                  Form("(%s) && (%s)", B1L0.Data(), per3.Data()), Form("(%s) && (%s)", B1L0.Data(), per3.Data()),
                                  Form("(%s) && (%s)", B0.Data(), cycCut.Data()), Form("(%s) && (%s)", B0.Data(), cycCut.Data()),
                                  Form("(%s) && (%s)", B0.Data(), cycCut.Data()), Form("(%s) && (%s)", B0.Data(), cycCut.Data()),
                                  B1, B1, B1, B1};
  Int_t sum0OnInd = 3;
  Int_t sum0Off1Ind = 7;
  Int_t sum0Off2Ind = 11;
  Int_t sum4OnInd = 19;
  Int_t sum4Off1Ind = 23;
  Int_t sum4Off2Ind = 27;

  Float_t sumOnAvg0 = 0.0; Float_t sumOff1Avg0 = 0.0; Float_t sumOff2Avg0 = 0.0;
  Float_t sumOnAvg4 = 0.0; Float_t sumOff1Avg4 = 0.0; Float_t sumOff2Avg4 = 0.0;
  for(Int_t i = 0; i < cycQrtVars; i++){
    TString hName = Form("h%i.%i_%s", runNum, cycNum, cycQrtTitles[i].Data());
    quartetwise->Draw(Form("%s>>%s", cycQrtNames[i].Data(), hName.Data()), varCuts[i].Data(), "goff");
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
    if(i==3){sumOnAvg0 = h->GetMean();}
    if(i==7){sumOff1Avg0 = h->GetMean();}
    if(i==11){sumOff2Avg0 = h->GetMean();}
    if(i==19){sumOnAvg4 = h->GetMean();}
    if(i==23){sumOff1Avg4 = h->GetMean();}
    if(i==27){sumOff2Avg4 = h->GetMean();}
  }

  Float_t sumOffAvg0 = (sumOff1Avg0 + sumOff2Avg0)/2.;
  Float_t sumOffAvg4 = (sumOff1Avg4 + sumOff2Avg4)/2.;
  Int_t nAsyms = 6;
  TString asymTitles[6] =  {"AsymAcc0LasOn", "AsymAcc0LasOff1", "AsymAcc0LasOff2",
                            "AsymAcc4LasOn", "AsymAcc4LasOff1", "AsymAcc4LasOff2"};
  TString asymNames[6] = {Form("(%s - %s)/(%s + %s - %f)", pos0.Data(), neg0.Data(), pos0.Data(), neg0.Data(), sumOffAvg0),
                          Form("(%s - %s)/(%f - %f)", pos0.Data(), neg0.Data(), sumOnAvg0, sumOffAvg0),
                          Form("(%s - %s)/(%f - %f)", pos0.Data(), neg0.Data(), sumOnAvg0, sumOffAvg0),
                          Form("(%s - %s)/(%s + %s - %f)", pos4.Data(), neg4.Data(), pos4.Data(), neg4.Data(), sumOffAvg4),
                          Form("(%s - %s)/(%f - %f)", pos4.Data(), neg4.Data(), sumOnAvg0, sumOffAvg4),
                          Form("(%s - %s)/(%f - %f)", pos4.Data(), neg4.Data(), sumOnAvg0, sumOffAvg4)};
  TString asymCuts[6] =  {Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()), Form("(%s) && (%s)", B1L0.Data(), per3.Data()),
                          Form("(%s) && (%s)", B1L1.Data(), per2.Data()), Form("(%s) && (%s)", B1L0.Data(), per1.Data()), Form("(%s) && (%s)", B1L0.Data(), per3.Data())};
  
  for(Int_t i = 0; i < nAsyms; i++){
    TString hName = Form("h%i.%i_%s", runNum, cycNum, asymTitles[i].Data());
    quartetwise->Draw(Form("%s>>%s", asymNames[i].Data(), hName.Data()), asymCuts[i].Data(), "goff");
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
  }
}

void runMPSPlots(TTree *mpswise, Int_t runNum, TFile *runOut){
  TString varCuts[23] =  {B1L1, B1L0, B0, B1L1, B1L0, B0,
                          L1, B1, B1, B1, B1, B1,
                          B1L1, B1L0, B1, B1, B1, B1,
                          B1L1, B1L0, B1L1, B1L0, "Entry$>=0"};

  TString mapCuts = evtCutStr(runNum);
  
  for(Int_t i = 0; i < runMPSVars; i++){
    TString hName = Form("h%i_%s", runNum, runMPSTitles[i].Data());
    TString cuts = varCuts[i];
    if(mapCuts.CompareTo("") != 0){
      cuts = Form("%s && %s", varCuts[i].Data(), mapCuts.Data());
    }
    mpswise->Draw(Form("%s>>%s", runMPSNames[i].Data(), hName.Data()), cuts.Data(), "goff");
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
  }
}

void runEpicsPlots(TTree *epicswise, TTree* mpswise, Int_t runNum, TFile *runOut){
  string *date_ptr = 0;
  epicswise->SetBranchAddress("epics_datestring", &date_ptr);
  epicswise->GetEntry(0);
  
  TString date_str(*date_ptr); TObjArray *date_split = date_str.Tokenize(" ");
  TString day(        ( (TObjString *)(date_split->At(0)) )->String());
  TString month(      ( (TObjString *)(date_split->At(1)) )->String());
  TString date_num_str( ( (TObjString *)(date_split->At(2)) )->String());
  TString time_str    ( ( (TObjString *)(date_split->At(3)) )->String());
  TString year_str    ( ( (TObjString *)(date_split->At(5)) )->String());

  TObjArray *time_ptr = time_str.Tokenize(":");
  TString hour_str(   ( (TObjString *)(time_ptr->At(0)) )->String());
  TString minute_str( ( (TObjString *)(time_ptr->At(1)) )->String());
  TString second_str( ( (TObjString *)(time_ptr->At(2)) )->String());

  Int_t date_num = atoi(date_num_str.Data());
  Int_t year     = atoi(year_str.Data());
  Int_t hour     = atoi(hour_str.Data());
  Int_t minute   = atoi(minute_str.Data());
  Int_t second   = atoi(second_str.Data());

  TH1F *h_month    = new TH1F(Form("h%i_month", runNum),    "Month of Year", 12,    1,   13); h_month->Fill(translate_month(month));
  TH1F *h_date_num = new TH1F(Form("h%i_date_num", runNum), "Day Number",    31,    1,   32); h_date_num->Fill(date_num);
  TH1F *h_year     = new TH1F(Form("h%i_year", runNum),     "Year",          50, 2001, 2051); h_year->Fill(year);
  TH1F *h_hour     = new TH1F(Form("h%i_hour", runNum),     "Hour",          24,    0,   24); h_hour->Fill(hour);
  TH1F *h_minute   = new TH1F(Form("h%i_minute", runNum),   "Minute",        60,    0,   60); h_minute->Fill(minute);
  TH1F *h_second   = new TH1F(Form("h%i_second", runNum),   "Second",        60,    0,   60); h_second->Fill(second);

  runOut->cd(); h_year->Write(); h_month->Write(); h_date_num->Write();
  h_hour->Write(); h_minute->Write(); h_second->Write();

  TString varCuts[10] =  {"Entry$>0", "Entry$>0", "Entry$>0", "Entry$>0", "Entry$>0", "Entry$>0", 
                          "Entry$>0", "epics_HWienAngle!=0", "epics_VWienAngle!=0", "epics_PhiFG!=0"};

  TString mapCuts = evtCutStr(runNum);

  for(Int_t i = 0; i < runEpcVars; i++){
    TString hName = Form("h%i_%s", runNum, runEpcTitles[i].Data());
    TString cuts = varCuts[i];
    if(mapCuts.CompareTo("") != 0){
      cuts = Form("%s && %s", varCuts[i].Data(), mapCuts.Data());
    }
    mpswise->Project(hName.Data(), runEpcNames[i].Data(), cuts.Data());
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
  }
}

#endif
