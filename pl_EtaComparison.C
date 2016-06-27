#include "anaConst14.h"
int readSysFile(TString,float*,float*,float*,float*,float*,float*);
int numActiveBins[4];

int pl_EtaComparison(){
  // Set Flags for Code/Plots
  Bool_t DEBUG = kFALSE;
  int isLogY = 1;

  // Access outside required info
  const Int_t numPtBins = anaConst::nPtBins;
  Float_t lowpt[numPtBins],highpt[numPtBins];
  for(Int_t c=0; c< numPtBins; c++){
    lowpt[c] = anaConst::lpt[c];
    highpt[c] = anaConst::hpt[c];
  }
  Float_t hptCut=anaConst::hptCut;
  const Int_t numTrigs = 4;
  Float_t hptMax=anaConst::hptMax; // Set max above range to allow overflow
  Float_t lowPhi=anaConst::lowPhi, highPhi=anaConst::highPhi;
  const Int_t numEtaBins = anaConst::nEtaBins;
  Float_t loweta[numEtaBins],higheta[numEtaBins];
  for(Int_t c=0; c< numEtaBins; c++){
    loweta[c] = anaConst::etaLow[c];
    higheta[c] = anaConst::etaHigh[c];
  }

  // create generic structures for use in all plots
  TPaveText* lbl[numEtaBins][numPtBins];
  TPaveText* lblE[numEtaBins];
  TPaveText* runInfo[numEtaBins];
  TPaveText* runInfo2[numEtaBins];
  char textLabel[100]; 


  // Load all trigger type files
  const char* baseName = "outputs/wConstraint/purityHists_June9";
  char BEMCName[4][100],SMDLName[4][100],SMDTName[4][100],BEMCHFTName[4][100],SMDLHFTName[4][100],SMDTHFTName[4][100];
  char trigName[4][100] = {"MB","BHT1","BHT2","BHT3"};
  TFile* BEMC[4];
  TFile* SMDL[4];
  TFile* SMDT[4];
  TFile* BEMCh[4];
  TFile* SMDLh[4];
  TFile* SMDTh[4];
  for(int q=0; q<4; q++)
  {
    //without HFT
    sprintf(BEMCName[q],"%s_Eta_%s_BEMC_processed.root",baseName,trigName[q]);
    cout<<BEMCName[q]<<endl;
    BEMC[q] = new TFile(BEMCName[q],"READ");
    sprintf(SMDLName[q],"%s_Eta_%s_SMD_processed.root",baseName,trigName[q]);
    SMDL[q] = new TFile(SMDLName[q],"READ");
    sprintf(SMDTName[q],"%s_Eta_%s_SMD2_processed.root",baseName,trigName[q]);
    SMDT[q] = new TFile(SMDTName[q],"READ");

    /* // with HFT
       sprintf(BEMCName[q],"%s_Eta_%s_BEMC_wHFT_processed.root",baseName,trigName[q]);
       BEMCh[q] = new TFile(BEMCName[q],"READ");
       sprintf(SMDLName[q],"%s_Eta_%s_SMD_wHFT_processed.root",baseName,trigName[q]);
       SMDLh[q] = new TFile(SMDLName[q],"READ");
       sprintf(SMDTName[q],"%s_Eta_%s_SMD2_wHFT_processed.root",baseName,trigName[q]);
       SMDTh[q] = new TFile(SMDTName[q],"READ");*/

  }

  for(int i=0; i<4; i++)
  {
    if(BEMC[i]->IsOpen() && SMDL[i]->IsOpen() && SMDT[i]->IsOpen())
      cout << "Trig " << i << " Files Open." << endl;
    else
    {
      cout << "Files Not Opened Properly." << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Access histograms
  TGraphErrors* nSigBEMC[4][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDL[4][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDT[4][numEtaBins][numPtBins];
  TGraphErrors* PurityBEMC[4][numEtaBins];
  TGraphErrors* PuritySMDL[4][numEtaBins];
  TGraphErrors* PuritySMDT[4][numEtaBins];
  TGraphErrors* dNdpTBEMC[4][numEtaBins];
  TGraphErrors* dNdpTSMDL[4][numEtaBins];
  TGraphErrors* dNdpTSMDT[4][numEtaBins];
  TGraphErrors* Purity2BEMC[4][numEtaBins];
  TGraphErrors* Purity2SMDL[4][numEtaBins];
  TGraphErrors* Purity2SMDT[4][numEtaBins];
  TGraphErrors* dNdpT2BEMC[4][numEtaBins];
  TGraphErrors* dNdpT2SMDL[4][numEtaBins];
  TGraphErrors* dNdpT2SMDT[4][numEtaBins];
  TGraphErrors* nSigBEMCh[4][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDLh[4][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDTh[4][numEtaBins][numPtBins];
  TGraphErrors* PurityBEMCh[4][numEtaBins];
  TGraphErrors* PuritySMDLh[4][numEtaBins];
  TGraphErrors* PuritySMDTh[4][numEtaBins];
  TGraphErrors* dNdpTBEMCh[4][numEtaBins];
  TGraphErrors* dNdpTSMDLh[4][numEtaBins];
  TGraphErrors* dNdpTSMDTh[4][numEtaBins];
  TGraphAsymmErrors* purBSys[4][numEtaBins];
  TGraphAsymmErrors* purSLSys[4][numEtaBins];
  TGraphAsymmErrors* purSTSys[4][numEtaBins];

  float pT[numEtaBins][numPtBins], pTErr[numEtaBins][numPtBins];
  float purFitB[numEtaBins][numPtBins],purFitSL[numEtaBins][numPtBins], purFitST[numEtaBins][numPtBins];
  float purFitBn[numEtaBins][numPtBins],purFitSLn[numEtaBins][numPtBins], purFitSTn[numEtaBins][numPtBins];
  float purFitErrB[numEtaBins][numPtBins],purFitErrSL[numEtaBins][numPtBins], purFitErrST[numEtaBins][numPtBins];
  float purFitErrBn[numEtaBins][numPtBins],purFitErrSLn[numEtaBins][numPtBins], purFitErrSTn[numEtaBins][numPtBins];
  float purCntB[numEtaBins][numPtBins],purCntSL[numEtaBins][numPtBins], purCntST[numEtaBins][numPtBins];
  float purCntBn[numEtaBins][numPtBins],purCntSLn[numEtaBins][numPtBins], purCntSTn[numEtaBins][numPtBins];
  float purCntErrB[numEtaBins][numPtBins],purCntErrSL[numEtaBins][numPtBins], purCntErrST[numEtaBins][numPtBins];
  float purCntErrBn[numEtaBins][numPtBins],purCntErrSLn[numEtaBins][numPtBins], purCntErrSTn[numEtaBins][numPtBins];
  float fitCountSysB[numEtaBins][numPtBins],wwoConstSysB[numEtaBins][numPtBins], purBSysUp[numEtaBins][numPtBins],purBSysDown[numEtaBins][numPtBins];
  float fitCountSysSL[numEtaBins][numPtBins],wwoConstSysSL[numEtaBins][numPtBins],purSLSysUp[numEtaBins][numPtBins],purSLSysDown[numEtaBins][numPtBins];;
  float fitCountSysST[numEtaBins][numPtBins],wwoConstSysST[numEtaBins][numPtBins],purSTSysUp[numEtaBins][numPtBins],purSTSysDown[numEtaBins][numPtBins];;

  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int etabin=0; etabin<numEtaBins; etabin++)
    {
      for(int ptbin=0; ptbin<numPtBins; ptbin++)
      {
        nSigBEMC[w][etabin][ptbin] = (TGraphErrors*)BEMC[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDL[w][etabin][ptbin] = (TGraphErrors*)SMDL[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDT[w][etabin][ptbin] = (TGraphErrors*)SMDT[w]->Get(Form("drawnSigmaE_%i%i",etabin,ptbin)); 
        /* nSigBEMCh[w][etabin][ptbin] = (TGraphErrors*)BEMCh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
           nSigSMDLh[w][etabin][ptbin] = (TGraphErrors*)SMDLh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
           nSigSMDTh[w][etabin][ptbin] = (TGraphErrors*)SMDTh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));*/
      }

      PurityBEMC[w][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawPurityFit_%i",etabin));
      PuritySMDL[w][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawPurityFit_%i",etabin));
      PuritySMDT[w][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawPurityFit_%i",etabin));
      dNdpTBEMC[w][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpTSMDL[w][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpTSMDT[w][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawdNdpT_%i",etabin));
      Purity2BEMC[w][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawPurityFit_%i",etabin));
      Purity2SMDL[w][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawPurityFit_%i",etabin));
      Purity2SMDT[w][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawPurityFit_%i",etabin));
      dNdpT2BEMC[w][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpT2SMDL[w][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpT2SMDT[w][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawdNdpT_%i",etabin));
      /* PurityBEMCh[w][etabin] = (TH1F*)BEMCh[w]->Get(Form("drawPurity_%i",etabin));
         PuritySMDLh[w][etabin] = (TH1F*)SMDLh[w]->Get(Form("drawPurity_%i",etabin));
         PuritySMDTh[w][etabin] = (TH1F*)SMDTh[w]->Get(Form("drawPurity_%i",etabin));*/

      // Read in information for systematics
      TString wConstB   = Form("outputs/wConstraint/systematicInformation_%s_BEMC_EtaBin%i.txt",trigName[w],etabin);
      TString wConstS   = Form("outputs/wConstraint/systematicInformation_%s_SMD_EtaBin%i.txt",trigName[w],etabin);
      TString wConstS2  = Form("outputs/wConstraint/systematicInformation_%s_SMD2_EtaBin%i.txt",trigName[w],etabin);
      TString noConstB  = Form("outputs/noConstraint/systematicInformation_%s_BEMC_EtaBin%i.txt",trigName[w],etabin);
      TString noConstS  = Form("outputs/noConstraint/systematicInformation_%s_SMD_EtaBin%i.txt",trigName[w],etabin);
      TString noConstS2 = Form("outputs/noConstraint/systematicInformation_%s_SMD2_EtaBin%i.txt",trigName[w],etabin);

      numActiveBins[w] = readSysFile(wConstB,pT[etabin],pTErr[etabin],purFitB[etabin],purFitErrB[etabin],purCntB[etabin],purCntErrB[etabin]);
      numActiveBins[w] = readSysFile(wConstS,pT[etabin],pTErr[etabin],purFitSL[etabin],purFitErrSL[etabin],purCntSL[etabin],purCntErrSL[etabin]);
      numActiveBins[w] = readSysFile(wConstS2,pT[etabin],pTErr[etabin],purFitST[etabin],purFitErrST[etabin],purCntST[etabin],purCntErrST[etabin]);
      numActiveBins[w] = readSysFile(noConstB,pT[etabin],pTErr[etabin],purFitBn[etabin],purFitErrBn[etabin],purCntBn[etabin],purCntErrBn[etabin]);
      numActiveBins[w] = readSysFile(noConstS,pT[etabin],pTErr[etabin],purFitSLn[etabin],purFitErrSLn[etabin],purCntSLn[etabin],purCntErrSLn[etabin]);
      numActiveBins[w] = readSysFile(noConstS2,pT[etabin],pTErr[etabin],purFitSTn[etabin],purFitErrSTn[etabin],purCntSTn[etabin],purCntErrSTn[etabin]);
      for(int ptbin =0; ptbin<numPtBins; ptbin++)
      {
        if(DEBUG) cout << pT[etabin][ptbin] << " " << purFitB[etabin][ptbin] <<" " << purFitBn[etabin][ptbin] << endl;
        fitCountSysB[etabin][ptbin] = purFitB[etabin][ptbin] - purCntB[etabin][ptbin];
        if(fitCountSysB[etabin][ptbin] > 0) purBSysDown[etabin][ptbin] += fitCountSysB[etabin][ptbin];
        wwoConstSysB[etabin][ptbin] = purFitB[etabin][ptbin] - purFitBn[etabin][ptbin];
        if(wwoConstSysB[etabin][ptbin] > 0) purBSysDown[etabin][ptbin] += wwoConstSysB[etabin][ptbin];
        fitCountSysSL[etabin][ptbin] = purFitSL[etabin][ptbin] - purCntSL[etabin][ptbin];
        if(fitCountSysSL[etabin][ptbin] > 0) purSLSysDown[etabin][ptbin] += fitCountSysSL[etabin][ptbin];
        wwoConstSysSL[etabin][ptbin] = purFitSL[etabin][ptbin] - purFitSLn[etabin][ptbin];
        if(wwoConstSysSL[etabin][ptbin] > 0) purSLSysDown[etabin][ptbin] += wwoConstSysSL[etabin][ptbin];
        fitCountSysST[etabin][ptbin] = purFitST[etabin][ptbin] - purCntST[etabin][ptbin];
        if(fitCountSysST[etabin][ptbin] > 0) purSTSysDown[etabin][ptbin] += fitCountSysST[etabin][ptbin];
        wwoConstSysST[etabin][ptbin] = purFitST[etabin][ptbin] - purFitSTn[etabin][ptbin];
        if(wwoConstSysST[etabin][ptbin] > 0) purSTSysDown[etabin][ptbin] += wwoConstSysST[etabin][ptbin];
      }
     purBSys[w][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[etabin],purFitB[etabin],0,0,purBSysDown[etabin],purBSysUp[etabin]); 
     purSLSys[w][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[etabin],purFitSL[etabin],0,0,purSLSysDown[etabin],purSLSysUp[etabin]); 
     purSTSys[w][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[etabin],purFitST[etabin],0,0,purSTSysDown[etabin],purSTSysUp[etabin]); 
    } 
  }
  if(DEBUG)
    cout << "Found Hists." << endl;


  // Make Canvas
  TCanvas* purityOL[4];
  TCanvas* purityTrigOL[4][numEtaBins];
  TCanvas* purityTrigOL2[4][numEtaBins];
  for(int etabin=0; etabin<numEtaBins; etabin++)
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      lbl[etabin][ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
      if(ptbin==0)lblE[etabin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
      lbl[etabin][ptbin]->AddText(textLabel);
      if(ptbin==0)lblE[etabin]->AddText(textLabel);
      sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
      lbl[etabin][ptbin]->AddText(textLabel);
      lbl[etabin][ptbin]->SetFillColor(kWhite); 
      if(ptbin==0)lblE[etabin]->SetFillColor(kWhite); 
    }
    runInfo[etabin] = new TPaveText(.55,.77,.85,.87,Form("NB NDC%i",etabin));
    sprintf(textLabel,"Run 14 Au+Au 200 GeV, 0-80%% Centrality");
    runInfo[etabin]->AddText(textLabel);
    sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
    runInfo[etabin]->AddText(textLabel);
    runInfo[etabin]->SetFillColor(kWhite);
    runInfo2[etabin] = new TPaveText(.55,.81,.85,.87,Form("NB NDC%i",etabin));
    sprintf(textLabel,"Run 14 Au+Au 200 GeV, 0-80%% Centrality");
    runInfo2[etabin]->AddText(textLabel);
    sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
    runInfo2[etabin]->AddText(textLabel);
    runInfo2[etabin]->SetFillColor(kWhite);
  }
  if(DEBUG)
    cout << "Canvas Made." << endl;

  // Plot Settings
  for(int r=0; r<4; r++)
  {
    if(DEBUG) cout << "plot settings trig: " << r << endl;
    purityOL[r]= new TCanvas(Form("purityOL_%i",r),"Purity Overlays",50,50,1050,1050);
    purityOL[r]->Divide(2,2);

    for(int etabin=0;etabin<numEtaBins;etabin++)
    {
      if(DEBUG) cout << "plot settings eta: " << etabin << endl;

      purityTrigOL[r][etabin]= new TCanvas(Form("purityTrigOL_%i_%i",r,etabin),"Purity Overlays",50,50,1050,1050);
      purityTrigOL2[r][etabin]= new TCanvas(Form("purityTrigOL2_%i_%i",r,etabin),"Purity Overlays",50,50,1050,1050);
      purityTrigOL[r][etabin]->Divide(1,2);
      purityTrigOL2[r][etabin]->Divide(1,2);
      PurityBEMC[r][etabin]->SetLineColor(1+etabin);
      PurityBEMC[r][etabin]->SetMarkerColor(1+etabin);
      PurityBEMC[r][etabin]->SetMarkerStyle(4+etabin);
      PurityBEMC[r][etabin]->SetTitle(Form("%s BEMC Purity",trigName[r]));
      if(DEBUG) cout << "purityBEMC" << endl;

      PuritySMDL[r][etabin]->SetLineColor(1+etabin);
      PuritySMDL[r][etabin]->SetMarkerColor(1+etabin);
      PuritySMDL[r][etabin]->SetMarkerStyle(4+etabin);
      PuritySMDL[r][etabin]->SetTitle(Form("%s SMD Purity",trigName[r]));
      if(DEBUG) cout << "puritySMDL" << endl;

      PuritySMDT[r][etabin]->SetLineColor(1+etabin);
      PuritySMDT[r][etabin]->SetMarkerColor(1+etabin);
      PuritySMDT[r][etabin]->SetMarkerStyle(4+etabin);
      PuritySMDT[r][etabin]->SetTitle(Form("%s SMD2 Purity",trigName[r]));
      PurityBEMC[r][etabin]->SetMarkerSize(1.0);
      PuritySMDL[r][etabin]->SetMarkerSize(1.0);
      PuritySMDT[r][etabin]->SetMarkerSize(1.0);
      if(DEBUG) cout << "puritySMDT" << endl;

      Purity2BEMC[r][etabin]->SetLineColor(kBlack);
      Purity2SMDL[r][etabin]->SetLineColor(kRed);
      Purity2SMDT[r][etabin]->SetLineColor(kBlue);
      Purity2BEMC[r][etabin]->SetMarkerColor(kBlack);
      Purity2SMDL[r][etabin]->SetMarkerColor(kRed);
      Purity2SMDT[r][etabin]->SetMarkerColor(kBlue);
      Purity2BEMC[r][etabin]->SetMarkerStyle(20);
      Purity2SMDL[r][etabin]->SetMarkerStyle(21);
      Purity2SMDT[r][etabin]->SetMarkerStyle(22);
      Purity2BEMC[r][etabin]->SetMarkerSize(1.0);
      Purity2SMDL[r][etabin]->SetMarkerSize(1.0);
      Purity2SMDT[r][etabin]->SetMarkerSize(1.0);
      Purity2BEMC[r][etabin]->SetTitle(Form("%s Purity Trigger Compare",trigName[r]));
      if(DEBUG) cout << "purity2" << endl;

      dNdpT2BEMC[r][etabin]->GetYaxis()->SetRangeUser(1e1,1e6);
      dNdpT2BEMC[r][etabin]->SetLineColor(kBlack);
      dNdpT2SMDL[r][etabin]->SetLineColor(kRed);
      dNdpT2SMDT[r][etabin]->SetLineColor(kBlue);
      dNdpT2BEMC[r][etabin]->SetMarkerColor(kBlack);
      dNdpT2SMDL[r][etabin]->SetMarkerColor(kRed);
      dNdpT2SMDT[r][etabin]->SetMarkerColor(kBlue);
      dNdpT2BEMC[r][etabin]->SetMarkerStyle(20);
      dNdpT2SMDL[r][etabin]->SetMarkerStyle(21);
      dNdpT2SMDT[r][etabin]->SetMarkerStyle(22);
      dNdpT2BEMC[r][etabin]->SetMarkerSize(1.0);
      dNdpT2SMDL[r][etabin]->SetMarkerSize(1.0);
      dNdpT2SMDT[r][etabin]->SetMarkerSize(1.0);
      dNdpT2BEMC[r][etabin]->SetTitle(Form("%s dN/dpT Trigger Compare",trigName[r]));
      if(DEBUG) cout << "dndpt2" << endl;

      purBSys[r][etabin]->SetLineColor(kBlack);
      purBSys[r][etabin]->SetLineWidth(1);
      purSLSys[r][etabin]->SetLineColor(kRed);
      purSLSys[r][etabin]->SetLineWidth(1);
      purSTSys[r][etabin]->SetLineColor(kBlue);
      purSTSys[r][etabin]->SetLineWidth(1);
    }
  }

  if(DEBUG)
    cout << "Settings Assigned." << endl;

  // Actually Draw

  // Use input for only one ptbin in legend. They're all the same, no need to loop.
  TLegend* leg = new TLegend(0.15,0.15,0.45,0.45);
  for(int etabin=0;etabin<numEtaBins;etabin++)
  {
    sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
    leg->AddEntry(PurityBEMC[1][etabin],textLabel,"lpe");
  }

  TLegend* leg2 = new TLegend(0.15,0.15,0.30,0.45);
  leg2->AddEntry(Purity2BEMC[1][0],"BEMC","lpe");
  leg2->AddEntry(Purity2SMDL[1][0],"SMD","lpe");
  leg2->AddEntry(Purity2SMDT[1][0],"SMD2","lpe");
  TLegend* leg3 = new TLegend(0.15,0.15,0.30,0.30);
  leg3->AddEntry(Purity2BEMC[1][0],"BEMC","lpe");
  leg3->AddEntry(Purity2SMDL[1][0],"SMD","lpe");
  leg3->AddEntry(Purity2SMDT[1][0],"SMD2","lpe");
  TLegend* leg4 = new TLegend(0.65,0.65,0.85,0.85);
  leg4->AddEntry(Purity2BEMC[1][0],"BEMC","lpe");
  leg4->AddEntry(Purity2SMDL[1][0],"SMD","lpe");
  leg4->AddEntry(Purity2SMDT[1][0],"SMD2","lpe");
  TLegend* leg5 = new TLegend(0.65,0.45,0.80,0.75);
  leg5->AddEntry(Purity2BEMC[1][0],"BEMC","lpe");
  leg5->AddEntry(Purity2SMDL[1][0],"SMD","lpe");
  leg5->AddEntry(Purity2SMDT[1][0],"SMD2","lpe");
  if(DEBUG)
    cout << "Legend Made." << endl;

  for(int t=0; t<4; t++)
  {
    if(DEBUG) cout << "Trig: " << t << endl;
    for(int etabin=0;etabin<numEtaBins;etabin++)
    {
      purityOL[t]->cd(1);
      if(etabin==0) PurityBEMC[t][etabin]->Draw("Ape"); 
      else PurityBEMC[t][etabin]->Draw("same pe");
      if(etabin==(numEtaBins-1))leg->Draw("same");

      purityOL[t]->cd(2);
      if(etabin==0)PuritySMDL[t][etabin]->Draw("Ape");
      else PuritySMDL[t][etabin]->Draw("same pe");
      if(etabin==(numEtaBins-1))leg->Draw("same");

      purityOL[t]->cd(3);
      if(etabin==0)PuritySMDT[t][etabin]->Draw("Ape");
      else PuritySMDT[t][etabin]->Draw("same pe");
      if(etabin==(numEtaBins-1))leg->Draw("same");

      purityTrigOL[t][etabin]->cd(1);
      Purity2BEMC[t][etabin]->Draw("Ape");
      purBSys[t][etabin]->Draw("[]");
      Purity2SMDL[t][etabin]->Draw("same pe");
      purSLSys[t][etabin]->Draw("[]");
      Purity2SMDT[t][etabin]->Draw("same pe");
      purSTSys[t][etabin]->Draw("[]");
      leg2->Draw("same");
      runInfo[etabin]->Draw("same");
      purityTrigOL[t][etabin]->cd(2);
      gPad->SetLogy(1);
      dNdpT2BEMC[t][etabin]->Draw("Ape");
      dNdpT2SMDL[t][etabin]->Draw("same pe");
      dNdpT2SMDT[t][etabin]->Draw("same pe");
      leg2->Draw("same");
      //lblE[etabin]->Draw("same");
      runInfo[etabin]->Draw("same");

      purityTrigOL2[t][etabin]->cd(1);
      gPad->SetLogy(1);
      leg2->Draw("same");
      runInfo[etabin]->Draw("same");

      purityTrigOL2[t][etabin]->cd(2);
      gPad->SetLogy(1);
      leg5->Draw("same");
      runInfo[etabin]->Draw("same");
    }
  }
  if(DEBUG) cout << "Drawing Finished. " << endl;

  // Place canvases in order on PDF
  //Set front page
  TCanvas* fp = new TCanvas("fp","Front Page",100,0,1000,900);
  fp->cd();
  TBox *bLabel = new TBox(0.01, 0.88, 0.99, 0.99);
  bLabel->SetFillColor(38);
  bLabel->Draw();
  TLatex tl;
  tl.SetNDC();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.033);
  char tlName[100];
  char tlName2[100];

  TString titlename = baseName;
  int found = titlename.Last('/');
  if(found >= 0){
    titlename.Replace(0, found+1, "");
  } 
  sprintf(tlName, "RUN 14 AuAu 200 GeV NPE Purity");
  tl.SetTextSize(0.05);
  tl.SetTextColor(kWhite);
  tl.DrawLatex(0.05, 0.92,tlName);

  TBox *bFoot = new TBox(0.01, 0.01, 0.99, 0.12);
  bFoot->SetFillColor(38);
  bFoot->Draw();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.05);
  tl.DrawLatex(0.05, 0.05, (new TDatime())->AsString());
  tl.SetTextColor(kBlack);
  tl.SetTextSize(0.03);
  tl.DrawLatex(0.1, 0.14, titlename);
  sprintf(tlName,"Study of Various Cut Sets on Inclusive Electron Purity.");
  tl.DrawLatex(0.1, 0.75,tlName);
  sprintf(tlName,"Z.W. Miller - UIC");
  tl.DrawLatex(0.1, 0.67,tlName);


  TCanvas* temp = new TCanvas();
  char name[100];
  sprintf(name, "%s_purityOverlays_CompareEta.pdf[", baseName);
  temp->Print(name);
  sprintf(name, "%s_purityOverlays_CompareEta.pdf", baseName);
  temp = fp; // print front page
  temp->Print(name);

  for(int q=0; q<4; q++)
  {
    if(DEBUG) cout << "Place Canvas in PDF trig: " << q << endl;
    temp = purityOL[q]; // print data canvases
    //temp->Print(name);
    for(int e=0; e<numEtaBins; e++)
    {  
      if(DEBUG) cout << "Place Canvas in PDF eta: " << e << endl;
      temp = purityTrigOL[q][e]; // print data canvases
      temp->Print(name);
      temp = purityTrigOL2[q][e]; // print data canvases
      //temp->Print(name);
    }
  }

  sprintf(name, "%s_purityOverlays_CompareEta.pdf]", baseName);
  temp->Print(name);


  return 1;
}


int readSysFile(TString fn, float* pT, float* ptErr, float* purFit, float* purFitErr, float* purCnt, float* purCntErr)
{
  char line[100];
  int k=0;
  ifstream f(fn,ios::in);
  if(f.is_open())
  {
    while(!f.eof())
    {
      f.getline(line,100);
      sscanf(line,"%f %f %f %f %f %f",&pT[k],&ptErr[k],&purFit[k],&purFitErr[k],&purCnt[k],&purCntErr[k]);
//      cout << "[" << line << "]" << endl;
      if(pT[k] > 0.5) k++;
    }
    f.close();
  }
  else
    cout << "a text file didn't open." << endl;
  return --k;
}
