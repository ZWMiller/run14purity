#include "anaConst14.h"
int readSysFile(TString,float*,float*,float*,float*,float*,float*);
int getLowCentralityLabel(int);
int getHighCentralityLabel(int);
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
  const Int_t numCentBins = anaConst::nCentBins;
  Float_t loweta[numEtaBins],higheta[numEtaBins];
  Float_t lowcent[numCentBins],highcent[numCentBins];
  for(Int_t c=0; c< numEtaBins; c++){
    loweta[c] = anaConst::etaLow[c];
    higheta[c] = anaConst::etaHigh[c];
  }
  for(Int_t c=0; c< numCentBins; c++){
    lowcent[c] = anaConst::centLow[c];
    highcent[c] = anaConst::centHigh[c];
  }

  // create generic structures for use in all plots
  TPaveText* lbl[numCentBins][numEtaBins][numPtBins];
  TPaveText* lblE[numCentBins][numEtaBins];
  TPaveText* runInfo[numCentBins][numEtaBins];
  TPaveText* runInfo2[numCentBins][numEtaBins];
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
  TGraphErrors* nSigBEMC[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDL[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDT[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* PurityBEMC[4][numCentBins][numEtaBins];
  TGraphErrors* PuritySMDL[4][numCentBins][numEtaBins];
  TGraphErrors* PuritySMDT[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTBEMC[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTSMDL[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTSMDT[4][numCentBins][numEtaBins];
  TGraphErrors* Purity2BEMC[4][numCentBins][numEtaBins];
  TGraphErrors* Purity2SMDL[4][numCentBins][numEtaBins];
  TGraphErrors* Purity2SMDT[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpT2BEMC[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpT2SMDL[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpT2SMDT[4][numCentBins][numEtaBins];
  TGraphErrors* nSigBEMCh[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDLh[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* nSigSMDTh[4][numCentBins][numEtaBins][numPtBins];
  TGraphErrors* PurityBEMCh[4][numCentBins][numEtaBins];
  TGraphErrors* PuritySMDLh[4][numCentBins][numEtaBins];
  TGraphErrors* PuritySMDTh[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTBEMCh[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTSMDLh[4][numCentBins][numEtaBins];
  TGraphErrors* dNdpTSMDTh[4][numCentBins][numEtaBins];
  TGraphAsymmErrors* purBSys[4][numCentBins][numEtaBins];
  TGraphAsymmErrors* purSLSys[4][numCentBins][numEtaBins];
  TGraphAsymmErrors* purSTSys[4][numCentBins][numEtaBins];

  float pT[numCentBins][numEtaBins][numPtBins], pTErr[numCentBins][numEtaBins][numPtBins];
  float purFitB[numCentBins][numEtaBins][numPtBins],purFitSL[numCentBins][numEtaBins][numPtBins], purFitST[numCentBins][numEtaBins][numPtBins];
  float purFitBn[numCentBins][numEtaBins][numPtBins],purFitSLn[numCentBins][numEtaBins][numPtBins], purFitSTn[numCentBins][numEtaBins][numPtBins];
  float purFitErrB[numCentBins][numEtaBins][numPtBins],purFitErrSL[numCentBins][numEtaBins][numPtBins], purFitErrST[numCentBins][numEtaBins][numPtBins];
  float purFitErrBn[numCentBins][numEtaBins][numPtBins],purFitErrSLn[numCentBins][numEtaBins][numPtBins], purFitErrSTn[numCentBins][numEtaBins][numPtBins];
  float purCntB[numCentBins][numEtaBins][numPtBins],purCntSL[numCentBins][numEtaBins][numPtBins], purCntST[numCentBins][numEtaBins][numPtBins];
  float purCntBn[numCentBins][numEtaBins][numPtBins],purCntSLn[numCentBins][numEtaBins][numPtBins], purCntSTn[numCentBins][numEtaBins][numPtBins];
  float purCntErrB[numCentBins][numEtaBins][numPtBins],purCntErrSL[numCentBins][numEtaBins][numPtBins], purCntErrST[numCentBins][numEtaBins][numPtBins];
  float purCntErrBn[numCentBins][numEtaBins][numPtBins],purCntErrSLn[numCentBins][numEtaBins][numPtBins], purCntErrSTn[numCentBins][numEtaBins][numPtBins];
  float fitCountSysB[numCentBins][numEtaBins][numPtBins],wwoConstSysB[numCentBins][numEtaBins][numPtBins], purBSysUp[numCentBins][numEtaBins][numPtBins],purBSysDown[numCentBins][numEtaBins][numPtBins];
  float fitCountSysSL[numCentBins][numEtaBins][numPtBins],wwoConstSysSL[numCentBins][numEtaBins][numPtBins],purSLSysUp[numCentBins][numEtaBins][numPtBins],purSLSysDown[numCentBins][numEtaBins][numPtBins];;
  float fitCountSysST[numCentBins][numEtaBins][numPtBins],wwoConstSysST[numCentBins][numEtaBins][numPtBins],purSTSysUp[numCentBins][numEtaBins][numPtBins],purSTSysDown[numCentBins][numEtaBins][numPtBins];;

  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int centbin = 0; centbin<numCentBins; centbin++){
      for(int etabin=0; etabin<numEtaBins; etabin++)
      {
        for(int ptbin=0; ptbin<numPtBins; ptbin++)
        {
          nSigBEMC[w][centbin][etabin][ptbin] = (TGraphErrors*)BEMC[w]->Get(Form("drawnSigmaE_%i_%i_%i",centbin,etabin,ptbin));
          nSigSMDL[w][centbin][etabin][ptbin] = (TGraphErrors*)SMDL[w]->Get(Form("drawnSigmaE_%i_%i_%i",centbin,etabin,ptbin));
          nSigSMDT[w][centbin][etabin][ptbin] = (TGraphErrors*)SMDT[w]->Get(Form("drawnSigmaE_%i_%i%i",centbin,etabin,ptbin)); 
          /* nSigBEMCh[w][etabin][ptbin] = (TGraphErrors*)BEMCh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
             nSigSMDLh[w][etabin][ptbin] = (TGraphErrors*)SMDLh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
             nSigSMDTh[w][etabin][ptbin] = (TGraphErrors*)SMDTh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));*/
        }

        PurityBEMC[w][centbin][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        PuritySMDL[w][centbin][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        PuritySMDT[w][centbin][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        dNdpTBEMC[w][centbin][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        dNdpTSMDL[w][centbin][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        dNdpTSMDT[w][centbin][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        Purity2BEMC[w][centbin][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        Purity2SMDL[w][centbin][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        Purity2SMDT[w][centbin][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawPurityFit_%i_%i",centbin,etabin));
        dNdpT2BEMC[w][centbin][etabin] = (TGraphErrors*)BEMC[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        dNdpT2SMDL[w][centbin][etabin] = (TGraphErrors*)SMDL[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        dNdpT2SMDT[w][centbin][etabin] = (TGraphErrors*)SMDT[w]->Get(Form("drawdNdpT_%i_%i",centbin,etabin));
        /* PurityBEMCh[w][centbin][etabin] = (TH1F*)BEMCh[w]->Get(Form("drawPurity_%i",etabin));
           PuritySMDLh[w][centbin][etabin] = (TH1F*)SMDLh[w]->Get(Form("drawPurity_%i",etabin));
           PuritySMDTh[w][centbin][etabin] = (TH1F*)SMDTh[w]->Get(Form("drawPurity_%i",etabin));*/

        // Read in information for systematics
        TString wConstB   = Form("outputs/wConstraint/systematicInformation_%s_BEMC_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);
        TString wConstS   = Form("outputs/wConstraint/systematicInformation_%s_SMD_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);
        TString wConstS2  = Form("outputs/wConstraint/systematicInformation_%s_SMD2_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);
        TString noConstB  = Form("outputs/noConstraint/systematicInformation_%s_BEMC_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);
        TString noConstS  = Form("outputs/noConstraint/systematicInformation_%s_SMD_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);
        TString noConstS2 = Form("outputs/noConstraint/systematicInformation_%s_SMD2_CentBin%i_EtaBin%i.txt",trigName[w],centbin,etabin);

        numActiveBins[w] = readSysFile(wConstB,pT[centbin][etabin],pTErr[centbin][etabin],purFitB[centbin][etabin],purFitErrB[centbin][etabin],purCntB[centbin][etabin],purCntErrB[centbin][etabin]);
        numActiveBins[w] = readSysFile(wConstS,pT[centbin][etabin],pTErr[centbin][etabin],purFitSL[centbin][etabin],purFitErrSL[centbin][etabin],purCntSL[centbin][etabin],purCntErrSL[centbin][etabin]);
        numActiveBins[w] = readSysFile(wConstS2,pT[centbin][etabin],pTErr[centbin][etabin],purFitST[centbin][etabin],purFitErrST[centbin][etabin],purCntST[centbin][etabin],purCntErrST[centbin][etabin]);
        numActiveBins[w] = readSysFile(noConstB,pT[centbin][etabin],pTErr[centbin][etabin],purFitBn[centbin][etabin],purFitErrBn[centbin][etabin],purCntBn[centbin][etabin],purCntErrBn[centbin][etabin]);
        numActiveBins[w] = readSysFile(noConstS,pT[centbin][etabin],pTErr[centbin][etabin],purFitSLn[centbin][etabin],purFitErrSLn[centbin][etabin],purCntSLn[centbin][etabin],purCntErrSLn[centbin][etabin]);
        numActiveBins[w] = readSysFile(noConstS2,pT[centbin][etabin],pTErr[centbin][etabin],purFitSTn[centbin][etabin],purFitErrSTn[centbin][etabin],purCntSTn[centbin][etabin],purCntErrSTn[centbin][etabin]);
        for(int ptbin =0; ptbin<numPtBins; ptbin++)
        {
          if(DEBUG) cout << pT[centbin][etabin][ptbin] << " " << purFitB[centbin][etabin][ptbin] <<" " << purFitBn[centbin][etabin][ptbin] << endl;
          fitCountSysB[centbin][etabin][ptbin] = purFitB[centbin][etabin][ptbin] - purCntB[centbin][etabin][ptbin];
          if(fitCountSysB[centbin][etabin][ptbin] > 0) purBSysDown[centbin][etabin][ptbin] += fitCountSysB[centbin][etabin][ptbin];
          wwoConstSysB[centbin][etabin][ptbin] = purFitB[centbin][etabin][ptbin] - purFitBn[centbin][etabin][ptbin];
          if(wwoConstSysB[centbin][etabin][ptbin] > 0) purBSysDown[centbin][etabin][ptbin] += wwoConstSysB[centbin][etabin][ptbin];
          fitCountSysSL[centbin][etabin][ptbin] = purFitSL[centbin][etabin][ptbin] - purCntSL[centbin][etabin][ptbin];
          if(fitCountSysSL[centbin][etabin][ptbin] > 0) purSLSysDown[centbin][etabin][ptbin] += fitCountSysSL[centbin][etabin][ptbin];
          wwoConstSysSL[centbin][etabin][ptbin] = purFitSL[centbin][etabin][ptbin] - purFitSLn[centbin][etabin][ptbin];
          if(wwoConstSysSL[centbin][etabin][ptbin] > 0) purSLSysDown[centbin][etabin][ptbin] += wwoConstSysSL[centbin][etabin][ptbin];
          fitCountSysST[centbin][etabin][ptbin] = purFitST[centbin][etabin][ptbin] - purCntST[centbin][etabin][ptbin];
          if(fitCountSysST[centbin][etabin][ptbin] > 0) purSTSysDown[centbin][etabin][ptbin] += fitCountSysST[centbin][etabin][ptbin];
          wwoConstSysST[centbin][etabin][ptbin] = purFitST[centbin][etabin][ptbin] - purFitSTn[centbin][etabin][ptbin];
          if(wwoConstSysST[centbin][etabin][ptbin] > 0) purSTSysDown[centbin][etabin][ptbin] += wwoConstSysST[centbin][etabin][ptbin];
        }
        purBSys[w][centbin][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[centbin][etabin],purFitB[centbin][etabin],0,0,purBSysDown[centbin][etabin],purBSysUp[centbin][etabin]); 
        purSLSys[w][centbin][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[centbin][etabin],purFitSL[centbin][etabin],0,0,purSLSysDown[centbin][etabin],purSLSysUp[centbin][etabin]); 
        purSTSys[w][centbin][etabin] = new TGraphAsymmErrors(numActiveBins[w],pT[centbin][etabin],purFitST[centbin][etabin],0,0,purSTSysDown[centbin][etabin],purSTSysUp[centbin][etabin]); 
      }
    } 
  }
  if(DEBUG)
    cout << "Found Hists." << endl;


  // Make Canvas
  TCanvas* purityOL[4];
  TCanvas* purityTrigOL[4][numCentBins][numEtaBins];
  TCanvas* purityTrigOL2[4][numCentBins][numEtaBins];
  for(int centbin = 0; centbin<numCentBins; centbin++)
  {
    for(int etabin=0; etabin<numEtaBins; etabin++)
    {
      for(int ptbin=0; ptbin<numPtBins; ptbin++)
      {
        lbl[centbin][etabin][ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
        if(ptbin==0)lblE[centbin][etabin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
        sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
        lbl[centbin][etabin][ptbin]->AddText(textLabel);
        if(ptbin==0)lblE[centbin][etabin]->AddText(textLabel);
        sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
        lbl[centbin][etabin][ptbin]->AddText(textLabel);
        lbl[centbin][etabin][ptbin]->SetFillColor(kWhite); 
        if(ptbin==0)lblE[centbin][etabin]->SetFillColor(kWhite); 
      }
      runInfo[centbin][etabin] = new TPaveText(.55,.77,.85,.87,Form("NB NDC%i",etabin));
      sprintf(textLabel,"Run 14 Au+Au 200 GeV, %i-%i%% Centrality",getLowCentralityLabel(lowcent[centbin]),getHighCentralityLabel(highcent[centbin]));
      runInfo[centbin][etabin]->AddText(textLabel);
      sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
      runInfo[centbin][etabin]->AddText(textLabel);
      runInfo[centbin][etabin]->SetFillColor(kWhite);
      runInfo2[centbin][etabin] = new TPaveText(.55,.81,.85,.87,Form("NB NDC%i",etabin));
      sprintf(textLabel,"Run 14 Au+Au 200 GeV, 0-80%% Centrality");
      runInfo2[centbin][etabin]->AddText(textLabel);
      sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
      runInfo2[centbin][etabin]->AddText(textLabel);
      runInfo2[centbin][etabin]->SetFillColor(kWhite);
    }
  }
  if(DEBUG)
    cout << "Canvas Made." << endl;

  // Plot Settings
  for(int r=0; r<4; r++)
  {
    if(DEBUG) cout << "plot settings trig: " << r << endl;
    purityOL[r]= new TCanvas(Form("purityOL_%i",r),"Purity Overlays",50,50,1050,1050);
    purityOL[r]->Divide(2,2);

    for(int centbin = 0; centbin<numCentBins; centbin++)
    {
      for(int etabin=0;etabin<numEtaBins;etabin++)
      {
        if(DEBUG) cout << "plot settings eta: " << etabin << endl;

        purityTrigOL[r][centbin][etabin]= new TCanvas(Form("purityTrigOL_%i_%i_%i",r,centbin,etabin),"Purity Overlays",50,50,1050,1050);
        purityTrigOL2[r][centbin][etabin]= new TCanvas(Form("purityTrigOL2_%i_%i_%i",r,centbin,etabin),"Purity Overlays",50,50,1050,1050);
        purityTrigOL[r][centbin][etabin]->Divide(1,2);
        purityTrigOL2[r][centbin][etabin]->Divide(1,2);
        PurityBEMC[r][centbin][etabin]->SetLineColor(1+etabin);
        PurityBEMC[r][centbin][etabin]->SetMarkerColor(1+etabin);
        PurityBEMC[r][centbin][etabin]->SetMarkerStyle(4+etabin);
        PurityBEMC[r][centbin][etabin]->SetTitle(Form("%s BEMC Purity",trigName[r]));
        if(DEBUG) cout << "purityBEMC" << endl;

        PuritySMDL[r][centbin][etabin]->SetLineColor(1+etabin);
        PuritySMDL[r][centbin][etabin]->SetMarkerColor(1+etabin);
        PuritySMDL[r][centbin][etabin]->SetMarkerStyle(4+etabin);
        PuritySMDL[r][centbin][etabin]->SetTitle(Form("%s SMD Purity",trigName[r]));
        if(DEBUG) cout << "puritySMDL" << endl;

        PuritySMDT[r][centbin][etabin]->SetLineColor(1+etabin);
        PuritySMDT[r][centbin][etabin]->SetMarkerColor(1+etabin);
        PuritySMDT[r][centbin][etabin]->SetMarkerStyle(4+etabin);
        PuritySMDT[r][centbin][etabin]->SetTitle(Form("%s SMD2 Purity",trigName[r]));
        PurityBEMC[r][centbin][etabin]->SetMarkerSize(1.0);
        PuritySMDL[r][centbin][etabin]->SetMarkerSize(1.0);
        PuritySMDT[r][centbin][etabin]->SetMarkerSize(1.0);
        if(DEBUG) cout << "puritySMDT" << endl;

        Purity2BEMC[r][centbin][etabin]->SetLineColor(kBlack);
        Purity2SMDL[r][centbin][etabin]->SetLineColor(kRed);
        Purity2SMDT[r][centbin][etabin]->SetLineColor(kBlue);
        Purity2BEMC[r][centbin][etabin]->SetMarkerColor(kBlack);
        Purity2SMDL[r][centbin][etabin]->SetMarkerColor(kRed);
        Purity2SMDT[r][centbin][etabin]->SetMarkerColor(kBlue);
        Purity2BEMC[r][centbin][etabin]->SetMarkerStyle(20);
        Purity2SMDL[r][centbin][etabin]->SetMarkerStyle(21);
        Purity2SMDT[r][centbin][etabin]->SetMarkerStyle(22);
        Purity2BEMC[r][centbin][etabin]->SetMarkerSize(1.0);
        Purity2SMDL[r][centbin][etabin]->SetMarkerSize(1.0);
        Purity2SMDT[r][centbin][etabin]->SetMarkerSize(1.0);
        Purity2BEMC[r][centbin][etabin]->SetTitle(Form("%s Purity Trigger Compare",trigName[r]));
        if(DEBUG) cout << "purity2" << endl;

        dNdpT2BEMC[r][centbin][etabin]->GetYaxis()->SetRangeUser(1e1,1e6);
        dNdpT2BEMC[r][centbin][etabin]->SetLineColor(kBlack);
        dNdpT2SMDL[r][centbin][etabin]->SetLineColor(kRed);
        dNdpT2SMDT[r][centbin][etabin]->SetLineColor(kBlue);
        dNdpT2BEMC[r][centbin][etabin]->SetMarkerColor(kBlack);
        dNdpT2SMDL[r][centbin][etabin]->SetMarkerColor(kRed);
        dNdpT2SMDT[r][centbin][etabin]->SetMarkerColor(kBlue);
        dNdpT2BEMC[r][centbin][etabin]->SetMarkerStyle(20);
        dNdpT2SMDL[r][centbin][etabin]->SetMarkerStyle(21);
        dNdpT2SMDT[r][centbin][etabin]->SetMarkerStyle(22);
        dNdpT2BEMC[r][centbin][etabin]->SetMarkerSize(1.0);
        dNdpT2SMDL[r][centbin][etabin]->SetMarkerSize(1.0);
        dNdpT2SMDT[r][centbin][etabin]->SetMarkerSize(1.0);
        dNdpT2BEMC[r][centbin][etabin]->SetTitle(Form("%s dN/dpT Trigger Compare",trigName[r]));
        if(DEBUG) cout << "dndpt2" << endl;

        purBSys[r][centbin][etabin]->SetLineColor(kBlack);
        purBSys[r][centbin][etabin]->SetLineWidth(1);
        purSLSys[r][centbin][etabin]->SetLineColor(kRed);
        purSLSys[r][centbin][etabin]->SetLineWidth(1);
        purSTSys[r][centbin][etabin]->SetLineColor(kBlue);
        purSTSys[r][centbin][etabin]->SetLineWidth(1);
      }
    }
  }

  if(DEBUG)
    cout << "Settings Assigned." << endl;

  // Actually Draw

  // Use input for only one ptbin in legend. They're all the same, no need to loop.
  TLegend* leg = new TLegend(0.15,0.15,0.45,0.45);
  for(int centbin=0;centbin<numCentBins;centbin++)
  {
    for(int etabin=0;etabin<numEtaBins;etabin++)
    {
      sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
      leg->AddEntry(PurityBEMC[1][centbin][etabin],textLabel,"lpe");
    }
  }

  TLegend* leg2 = new TLegend(0.15,0.15,0.30,0.45);
  leg2->AddEntry(Purity2BEMC[1][0][0],"BEMC","lpe");
  leg2->AddEntry(Purity2SMDL[1][0][0],"SMD","lpe");
  leg2->AddEntry(Purity2SMDT[1][0][0],"SMD2","lpe");
  TLegend* leg3 = new TLegend(0.15,0.15,0.30,0.30);
  leg3->AddEntry(Purity2BEMC[1][0][0],"BEMC","lpe");
  leg3->AddEntry(Purity2SMDL[1][0][0],"SMD","lpe");
  leg3->AddEntry(Purity2SMDT[1][0][0],"SMD2","lpe");
  TLegend* leg4 = new TLegend(0.65,0.65,0.85,0.85);
  leg4->AddEntry(Purity2BEMC[1][0][0],"BEMC","lpe");
  leg4->AddEntry(Purity2SMDL[1][0][0],"SMD","lpe");
  leg4->AddEntry(Purity2SMDT[1][0][0],"SMD2","lpe");
  TLegend* leg5 = new TLegend(0.65,0.45,0.80,0.75);
  leg5->AddEntry(Purity2BEMC[1][0][0],"BEMC","lpe");
  leg5->AddEntry(Purity2SMDL[1][0][0],"SMD","lpe");
  leg5->AddEntry(Purity2SMDT[1][0][0],"SMD2","lpe");
  if(DEBUG)
    cout << "Legend Made." << endl;

  for(int t=0; t<4; t++)
  {
    if(DEBUG) cout << "Trig: " << t << endl;
    for(int centbin = 0; centbin<numCentBins; centbin++)
    {
      for(int etabin=0;etabin<numEtaBins;etabin++)
      {
        purityOL[t]->cd(1);
        if(etabin==0) PurityBEMC[t][centbin][etabin]->Draw("Ape"); 
        else PurityBEMC[t][centbin][etabin]->Draw("same pe");
        if(etabin==(numEtaBins-1))leg->Draw("same");

        purityOL[t]->cd(2);
        if(etabin==0)PuritySMDL[t][centbin][etabin]->Draw("Ape");
        else PuritySMDL[t][centbin][etabin]->Draw("same pe");
        if(etabin==(numEtaBins-1))leg->Draw("same");

        purityOL[t]->cd(3);
        if(etabin==0)PuritySMDT[t][centbin][etabin]->Draw("Ape");
        else PuritySMDT[t][centbin][etabin]->Draw("same pe");
        if(etabin==(numEtaBins-1))leg->Draw("same");

        purityTrigOL[t][centbin][etabin]->cd(1);
        Purity2BEMC[t][centbin][etabin]->Draw("Ape");
        purBSys[t][centbin][etabin]->Draw("[]");
        Purity2SMDL[t][centbin][etabin]->Draw("same pe");
        purSLSys[t][centbin][etabin]->Draw("[]");
        Purity2SMDT[t][centbin][etabin]->Draw("same pe");
        purSTSys[t][centbin][etabin]->Draw("[]");
        leg2->Draw("same");
        runInfo[centbin][etabin]->Draw("same");
        purityTrigOL[t][centbin][etabin]->cd(2);
        gPad->SetLogy(1);
        dNdpT2BEMC[t][centbin][etabin]->Draw("Ape");
        dNdpT2SMDL[t][centbin][etabin]->Draw("same pe");
        dNdpT2SMDT[t][centbin][etabin]->Draw("same pe");
        leg2->Draw("same");
        //lblE[centbin][etabin]->Draw("same");
        runInfo[centbin][etabin]->Draw("same");

        purityTrigOL2[t][centbin][etabin]->cd(1);
        gPad->SetLogy(1);
        leg2->Draw("same");
        runInfo[centbin][etabin]->Draw("same");

        purityTrigOL2[t][centbin][etabin]->cd(2);
        gPad->SetLogy(1);
        leg5->Draw("same");
        runInfo[centbin][etabin]->Draw("same");
      }
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
    for(int c=0; c<numCentBins; c++)
    for(int e=0; e<numEtaBins; e++)
    {  
      if(DEBUG) cout << "Place Canvas in PDF eta: " << e << endl;
      temp = purityTrigOL[q][c][e]; // print data canvases
      temp->Print(name);
      temp = purityTrigOL2[q][c][e]; // print data canvases
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

int getLowCentralityLabel(int x)
{
  return 75-x*5;
}
int getHighCentralityLabel(int x)
{
  return 80-x*5;
}
