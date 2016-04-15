
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////
//// Offline Plots - Z. Miller Jan 6, 2016
//// --------------------
//// Updated 3/3/16  to allow for all cut sets in a single call - ZWM
//// Updated 3/8/16  to produce eta dependent binning in production - ZWM
//// Updated 3/16/16 to use THnSparse input, project to TH3, then fit - ZWM
//// --------------------
////
//// root -l
//// .L offline.C
//// offline("FILENAME", triggerType, Cutset, isHFT) 
//// ## File name without .root Extension
//// 
//// Trigger type: 0 = MB, 1=BHT1, 2=BHT2, 3=BHT3, 4=All
//// Cut Set Allowed Values: "BEMC", "SMD", "SMD2", "ALL"
//// isHFT either True or False. If true, only process tracks with HFT Match
////
//// Each cut set and trigger type generates it's own file. Using 4 and "ALL"
//// will create all possible cut sets for each trigger. 
////
//// Example: offline("Feb6_PuritySample",4,"ALL") - This will generate all 
//// cutset-trigger combinations for the Feb6_PuritySample. Depending on runtime
//// options selected (you will be prompted) it will generate a PDF for each set
//// with the relevant plots and a .root file containing all histograms created.
////
//// For general use, should just need to change the path where it looks for the 
//// main input files. Search for "NPEPurity" to find the path.
////
//// anaConst14.h contains most of the constants - pT binning, eta Binning, etc
////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "anaConst14.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "TString.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TF1.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include <sstream>
#include "TLatex.h"

using namespace std;

// Declare functions
void makeHist(const char*,Int_t, const char*);
void offlineEtaBin_ErrorUpdate(const char*, Int_t, const char*, Bool_t);
void checkBatchMode();
Bool_t checkMakePDF();
Bool_t checkMakeRoot();
Bool_t makePDF,makeROOT;
int isLogY = 1;
Bool_t drawAll = kFALSE;
Bool_t withMergedPion = kFALSE;
Bool_t DEBUG = kFALSE;
Bool_t wHFT = kFALSE;

void offline(const char* FileName="test", Int_t trig=4, const char* Cuts="BEMC",Bool_t ishft=kFALSE)
{
  offlineEtaBin_ErrorUpdate(FileName,trig,Cuts,ishft);
}

void offlineEtaBin_ErrorUpdate(const char* FileName="test", Int_t trig=4, const char* Cuts="BEMC",Bool_t ishft=kFALSE) //0=MB,1=HT1,2=HT2,3=HT3,4=ALL, BEMC, SMD, TOF
{
  wHFT = ishft;
  const char* cutTypes[3] = {"BEMC","SMD","SMD2"};
  if(wHFT) cout << "---Running with HFT Match Requirement!---" << endl;
  if(trig == 4)
  {
    checkBatchMode();
    makePDF = checkMakePDF();
    makeROOT= checkMakeRoot();

    if(!strncmp(Cuts,"ALL",3))
    {
      for(int q=0;q<3;q++)
      {
        for(Int_t i=1;i<trig;i++)
          makeHist(FileName,i,cutTypes[q]);
      }
    }
    else
    {
      for(Int_t i=1;i<trig;i++)
        makeHist(FileName,i,Cuts);
    }
  }
  else if(trig > 4)
    cout << "Bad Trigger Input, 0-4 only" << endl;
  else
  {
    checkBatchMode();
    makePDF = checkMakePDF();
    makeROOT= checkMakeRoot();
    if(DEBUG) cout << "Make Root: " << makeROOT << endl;
    if(!strncmp(Cuts,"ALL",3))
    {
      for(int q=0;q<3;q++)
      {
        makeHist(FileName,trig,cutTypes[q]);
      }
    }
    makeHist(FileName,trig,Cuts);
  }
}

void makeHist(const char* FileName="test", Int_t trig=4,const char* Cut="BEMC")
{
  TH1F::SetDefaultSumw2();
  TH3D::SetDefaultSumw2();
  // Set Style parameters for this macro
  gStyle->SetOptTitle(1); // Show Title (off by default for cleanliness)
  gErrorIgnoreLevel = kError; // Set Verbosity Level (kPrint shows all)
  char Cuts[100];
  sprintf(Cuts,"%s",Cut);
  if(strncmp(Cut,"BEMC",4) && strncmp(Cut,"TOF",3) && strncmp(Cut,"SMD",3) && strncmp(Cut,"SMD2",4))
  {
    cout << "Wrong Input for cut type: default to BEMC" << endl;
    sprintf(Cuts, "BEMC");
  }

  char FileLabel[100];
  if(trig == 0)
    sprintf(FileLabel, "MB");
  if(trig == 1)
    sprintf(FileLabel, "BHT1");
  if(trig == 2)
    sprintf(FileLabel, "BHT2");
  if(trig == 3)
    sprintf(FileLabel, "BHT3");

  Int_t number;
  // Open ROOT File
  char name[1000];
  sprintf(name,"/Users/zach/Research/rootFiles/run14NPEpurity/%s.root",FileName);
  TFile *f = new TFile(name,"READ");
  if (f->IsOpen()==kFALSE)
  { std::cout << "!!! File Not Found !!!" << std::endl;
    exit(1); }
  if (f->IsOpen())
  {
    cout << "Opened " << name << endl;
  }

  // f->ls(); // - DEBUG by printing all objects in ROOT file

  char fname[100];
  TFile* file;
  if(DEBUG) cout << "Make Root: " << makeROOT << endl;
  if(makeROOT){
    sprintf(fname,"/Users/zach/Research/rootFiles/run14NPEpurity/%s_Eta_%s_%s_processed.root",FileName,FileLabel,Cuts);
    if(wHFT)
      sprintf(fname,"/Users/zach/Research/rootFiles/run14NPEpurity/%s_Eta_%s_%s_wHFT_processed.root",FileName,FileLabel,Cuts);
    file = new TFile(fname,"RECREATE");
    if (file->IsOpen()==kFALSE)
    {
      std::cout << "!!! Outfile Not Opened !!!" << std::endl;
      makeROOT = kFALSE;
    }
  }

  const Int_t numPtBins  = anaConst::nPtBins;
  const Int_t numEtaBins = anaConst::nEtaBins;
  const Int_t numCanvas = numPtBins/9 + 1;// Each eta gets its own pt bins
  Float_t lowpt[numPtBins],highpt[numPtBins];
  Float_t loweta[numEtaBins],higheta[numEtaBins];
  for(Int_t c=0; c< numPtBins; c++){
    lowpt[c] = anaConst::lpt[c];
    highpt[c] = anaConst::hpt[c];
  }
  for(Int_t c=0; c< numEtaBins; c++){
    loweta[c] = anaConst::etaLow[c];
    higheta[c] = anaConst::etaHigh[c];
  }
  Float_t hptCut=anaConst::hptCut;
  const Int_t numTrigs = 4;
  Float_t hptMax=anaConst::hptMax; // Set max above range to allow overflow
  Float_t lowPhi=anaConst::lowPhi, highPhi=anaConst::highPhi;

  TPaveText* lbl[numEtaBins][numPtBins];
  char textLabel[100];

  TH1D* projnSigmaE[numEtaBins][numPtBins];
  TH1D* drawnSigmaE[numEtaBins][numPtBins];
  TF1 *fitPi[numEtaBins][numPtBins];
  TF1 *fitKP[numEtaBins][numPtBins];
  TF1 *fitmPi[numEtaBins][numPtBins];
  TF1 *fitE[numEtaBins][numPtBins];
  TF1 *fitCom[numEtaBins][numPtBins];
  TF1 *fitPiD[numEtaBins][numPtBins];
  TF1 *fitmPiD[numEtaBins][numPtBins];
  TF1 *fitKPD[numEtaBins][numPtBins];
  TF1 *fitED[numEtaBins][numPtBins];
  TF1 *fitComD[numEtaBins][numPtBins];
  double eInte[numEtaBins][numPtBins];
  double piInte[numEtaBins][numPtBins];
  double mpiInte[numEtaBins][numPtBins];
  double kpInte[numEtaBins][numPtBins];
  double sum[numEtaBins][numPtBins], purity[numEtaBins][numPtBins],
         pT[numEtaBins][numPtBins], dNdpT[numEtaBins][numPtBins],
         XnDOF[numEtaBins][numPtBins], aYield[numEtaBins][numPtBins],
         aPur[numEtaBins][numPtBins];
  double parPlot[12][numEtaBins][numPtBins], 
         errPlot[12][numEtaBins][numPtBins];
  double dx[numEtaBins][numPtBins],dy[numEtaBins][numPtBins],dy2[numEtaBins][numPtBins];

  // Make Canvas
  TCanvas* nSigE[numEtaBins][numCanvas]; 
  TCanvas* purC[numEtaBins];
  TCanvas* xSquare[numEtaBins];
  TCanvas* purityGaussians[numEtaBins];
  TCanvas* parameterGaussians[numEtaBins];
  TCanvas* parameterFitCanvas[numEtaBins];
  TCanvas* junk = new TCanvas("junk","junk",50,50,1050,1050);

  for(Int_t etaBin=0; etaBin < numEtaBins; etaBin++){

    purC[etaBin] = new TCanvas(Form("purC_%i",etaBin),"Electron Purity",50,50,1050,1050);
    parameterFitCanvas[etaBin] = new TCanvas(Form("parFitCanvas_%i",etaBin),"Fit Params vs pT",50,50,1050,1050);

    if(withMergedPion)
      parameterFitCanvas[etaBin]->Divide(4,3);
    else
      parameterFitCanvas[etaBin]->Divide(3,3);
    purC[etaBin]->Divide(1,2);
    xSquare[etaBin] = new TCanvas(Form("xSquare_%i",etaBin),"Chi Square Check",50,50,1050,1050);
    purityGaussians[etaBin] = new TCanvas(Form("purityGaussians_%i",etaBin),"Electron Purity Gaussians",50,50,1050,1050);
    purityGaussians[etaBin]->Divide(4,3);
    parameterGaussians[etaBin] = new TCanvas(Form("parameterGaussians_%i",etaBin),"Electron parameter Gaussians",50,50,1050,1050);
    parameterGaussians[etaBin]->Divide(4,3);

    for(Int_t q = 0; q < numCanvas; q++)
    {
      nSigE[etaBin][q] = new TCanvas(Form("nSigE_%i_%i",etaBin,q),"nSigma Electron Projections",50,50,1050,1050);
      nSigE[etaBin][q] -> Divide(3,3);
    }
  }

  // Make Projections (first get 2d/3d hists, then project)
  THnSparse* nSigmaESparse = (THnSparse*)f->Get(Form("nSigmaE_%s_%i",Cuts,trig));
  if(wHFT)
  {
    nSigmaESparse = (THnSparse*)f->Get(Form("nSigmaE_%s_HFT_%i",Cuts,trig));
  }

  nSigmaESparse->CalculateErrors();
  TH3D* nSigmaEPtEta;
  nSigmaESparse->GetAxis(3)->SetRange(0, 15); // Only select 0-80% Centrality (centrality16)
  nSigmaEPtEta = nSigmaESparse->Projection(0,1,2);
  if(DEBUG)cout << "got hists trig "<< trig << endl;

  double parStorage[12] = {0.};
  double parHold[12] = {0.};
  double parErrStorage[12] = {0.};
  for(Int_t etaBin=0; etaBin<numEtaBins; etaBin++)
  {
    for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
    {
      projnSigmaE[etaBin][ptbin]  = (TH1D*)nSigmaEPtEta->ProjectionY(Form("projnSigmaE_%i_%i",etaBin,ptbin),nSigmaEPtEta->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEPtEta->GetXaxis()->FindBin(highpt[ptbin])-1,nSigmaEPtEta->GetZaxis()->FindBin(loweta[etaBin]),nSigmaEPtEta->GetZaxis()->FindBin(higheta[etaBin])-1);
    }
  }

  // Analyze the projections
  for(Int_t etaBin= 0; etaBin < numEtaBins; etaBin++)
  {
    parStorage[0] = 0; // reset fitting algorithm for each eta bin
    if(DEBUG) cout << endl << "EtaBin: " << etaBin << endl;

    for(Int_t ptbin = 0; ptbin < numPtBins; ptbin++){
      // Clear the variables of interest
      purity[etaBin][ptbin] = dNdpT[etaBin][ptbin] = dx[etaBin][ptbin] = dy[etaBin][ptbin] = 0;
      for(int ii=0; ii<12; ii++)
      {
        parPlot[ii][etaBin][ptbin] = errPlot[ii][etaBin][ptbin] = 0;
      }

      // Don't analyze below the trigger turn on curve
      if(lowpt[ptbin] < anaConst::trigThreshold[trig]) continue; 

      // Init necessary plotting tools
      lbl[etaBin][ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
      sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etaBin],higheta[etaBin]);
      lbl[etaBin][ptbin]->AddText(textLabel);
      sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
      lbl[etaBin][ptbin]->AddText(textLabel);
      lbl[etaBin][ptbin]->SetFillColor(kWhite);

      int activeCanvas = (int)ptbin/9;
      int activeBin = ptbin - activeCanvas*9;
      if(DEBUG)cout << ptbin << ": " << activeCanvas << " " << activeBin << endl;

      nSigE[etaBin][activeCanvas]->cd(activeBin+1);
      gPad->SetLogy(isLogY);
      projnSigmaE[etaBin][ptbin]->SetMarkerColor(kBlack);
      projnSigmaE[etaBin][ptbin]->SetMarkerStyle(20);
      projnSigmaE[etaBin][ptbin]->SetLineColor(kBlack);
      projnSigmaE[etaBin][ptbin]->SetMarkerSize(0.5);
      projnSigmaE[etaBin][ptbin]->GetXaxis()->SetRangeUser(-11.,16.);
      projnSigmaE[etaBin][ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
      projnSigmaE[etaBin][ptbin]->SetTitle("");
      drawnSigmaE[etaBin][ptbin] = (TH1D*)projnSigmaE[etaBin][ptbin]->Clone();
      drawnSigmaE[etaBin][ptbin]->SetName(Form("drawnSigmaE_%i_%i",etaBin,ptbin));
      drawnSigmaE[etaBin][ptbin]->Write();
      projnSigmaE[etaBin][ptbin]->Draw();
      lbl[etaBin][ptbin]->Draw("same");
      if(DEBUG) cout << "After Clone projnSigmaE" << endl;

      // Do fits on nSig E
      //      nSigE[etaBin][activeCanvas]->cd(activeBin+1);
      gStyle->SetOptFit(1111);
      Double_t par[12];
      Double_t parErr[12];
      Double_t parS[3];
      Double_t kpLow,kpHigh,piLow,piHigh,eLow,eHigh, mpiLow, mpiHigh;
      kpLow = -10.; kpHigh =-6.; 
      piLow = -4. ; piHigh = -2.;
      eLow = -1.0  ; eHigh = 2.0;
      mpiLow = 3.5; mpiHigh = 6;
      fitKP[etaBin][ptbin] = new TF1(Form("fitKP_%i_%i",etaBin,ptbin),"gaus",kpLow,kpHigh);
      fitPi[etaBin][ptbin] = new TF1(Form("fitPi_%i_%i",etaBin,ptbin),"gaus",piLow,piHigh);
      fitE[etaBin][ptbin]  = new TF1(Form("fitE_%i_%i",etaBin,ptbin),"gaus",eLow,eHigh);
      fitmPi[etaBin][ptbin]= new TF1(Form("fitmPi_%i_%i",etaBin,ptbin),"gaus",mpiLow,mpiHigh);
      if(withMergedPion)
        fitCom[etaBin][ptbin]= new TF1(Form("fitCom_%i_%i",etaBin,ptbin),"gaus(0)+gaus(3)+gaus(6)+gaus(9)",-10,10);
      else
        fitCom[etaBin][ptbin]= new TF1(Form("fitCom_%i_%i",etaBin,ptbin),"gaus(0)+gaus(3)+gaus(6)",-10,10);

      fitCom[etaBin][ptbin]->SetParName(0,"#pi C");
      fitCom[etaBin][ptbin]->SetParName(1,"#pi #mu");
      fitCom[etaBin][ptbin]->SetParName(2,"#pi #sigma");
      fitCom[etaBin][ptbin]->SetParName(3,"Kp C");
      fitCom[etaBin][ptbin]->SetParName(4,"Kp #mu");
      fitCom[etaBin][ptbin]->SetParName(5,"Kp #sigma");
      fitCom[etaBin][ptbin]->SetParName(6,"e C");
      fitCom[etaBin][ptbin]->SetParName(7,"e #mu");
      fitCom[etaBin][ptbin]->SetParName(8,"e #sigma");
      fitCom[etaBin][ptbin]->SetParName(9, "mer#pi C");
      fitCom[etaBin][ptbin]->SetParName(10,"mer#pi #mu");
      fitCom[etaBin][ptbin]->SetParName(11,"mer#pi #sigma");

      fitPi[etaBin][ptbin]->SetLineColor(kRed);
      fitKP[etaBin][ptbin]->SetLineColor(kCyan);
      fitE[etaBin][ptbin]->SetLineColor(kBlue);
      fitmPi[etaBin][ptbin]->SetLineColor(kGreen);
      fitCom[etaBin][ptbin]->SetLineColor(kMagenta);
      fitPi[etaBin][ptbin]->SetLineWidth(1);
      fitmPi[etaBin][ptbin]->SetLineWidth(1);
      fitKP[etaBin][ptbin]->SetLineWidth(1);
      fitE[etaBin][ptbin]->SetLineWidth(1);
      fitCom[etaBin][ptbin]->SetLineWidth(2);

      // make versions to draw fits in full
      fitPiD[etaBin][ptbin] = new TF1(Form("drawPi_%i_%i",etaBin,ptbin),"gaus",-10,10);
      fitmPiD[etaBin][ptbin]= new TF1(Form("drawmPi_%i_%i",etaBin,ptbin),"gaus",-10,10);
      fitKPD[etaBin][ptbin] = new TF1(Form("drawKP_%i_%i",etaBin,ptbin),"gaus",-10,10);
      fitED[etaBin][ptbin]  = new TF1(Form("drawE_%i_%i",etaBin,ptbin),"gaus",-10,10);

      // If this is the first bin of the trigger, fit the roughly
      // expected areas for each Gaussian to get starting pars for total fit
      if(parStorage[0] == 0.)
      {
        projnSigmaE[etaBin][ptbin]->Fit(fitPi[etaBin][ptbin],"RQ");
        projnSigmaE[etaBin][ptbin]->Fit(fitKP[etaBin][ptbin],"RQ");
        projnSigmaE[etaBin][ptbin]->Fit(fitE[etaBin][ptbin],"RQ");
        if(withMergedPion)
          projnSigmaE[etaBin][ptbin]->Fit(fitmPi[etaBin][ptbin],"RQ");

        fitPi[etaBin][ptbin]->GetParameters(&par[0]);
        fitKP[etaBin][ptbin]->GetParameters(&par[3]);
        fitE[etaBin][ptbin]->GetParameters(&par[6]);
        if(withMergedPion)
          fitmPi[etaBin][ptbin]->GetParameters(&par[9]);
        fitCom[etaBin][ptbin]->SetParameters(par);
      }
      else // Otherwise, use results from last fit as input
        fitCom[etaBin][ptbin]->SetParameters(parStorage);
      fitCom[etaBin][ptbin]->GetParameters(&parHold[0]);
      // Apply parameter constraints
      //  fitCom[etaBin][ptbin]->SetParLimits(0,0,1000000);
      fitCom[etaBin][ptbin]->SetParLimits(1,-6,-4);
      fitCom[etaBin][ptbin]->SetParLimits(2,0.8,1.2);
      //  fitCom[etaBin][ptbin]->SetParLimits(3,0,1000000);
      fitCom[etaBin][ptbin]->SetParLimits(4,-5,-2);
      fitCom[etaBin][ptbin]->SetParLimits(5,0.8,1.2);
      //  fitCom[etaBin][ptbin]->SetParLimits(6,0,1000000);
      fitCom[etaBin][ptbin]->SetParLimits(7,-0.5,0.0);
      fitCom[etaBin][ptbin]->SetParLimits(8,0.8,1.2);
      //  fitCom[etaBin][ptbin]->SetParLimits(9,0,1000000);
      //  fitCom[etaBin][ptbin]->SetParLimits(10,2,4);
      //  fitCom[etaBin][ptbin]->SetParLimits(11,0.8,1.5);

      // Do fits 1000 times on slightly randomized data to get mean and error of purity
      TH1F* purGaus = new TH1F("purGaus","Purity Results for All Fits",1000,0,2);
      TH1F* dNdpTGaus = new TH1F("dNdpTGaus","dNdpT Results for All Fits",20000,0,20000);
      TH1F* parGaus[12];
      TH1F* chiGaus = new TH1F("chiGaus","Chi2/NDF Results for All Fits",500,-10,10);
      for(int ii=0; ii<12; ii++)
      {
        if(ii == 0 || ii == 3 || ii == 6 || ii == 9)
          parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1600,-8000,8000);
        else if(ii == 1 || ii == 4 || ii == 7 || ii == 10)
          parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1000,-10,10);
        else if(ii == 2 || ii == 5 || ii == 8 || ii == 11)
          parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1000,0,2);
      }  
      TRandom3 *gRnd= new TRandom3(0);
      TH1D* HH;
      for(int it=0; it<1000; it++)
      {
        HH = (TH1D*)projnSigmaE[etaBin][ptbin]->Clone();
        int bins = HH->GetNbinsX();
        for(int b=0; b<bins; b++)
        {
          if(projnSigmaE[etaBin][ptbin]->GetBinContent(b) && projnSigmaE[etaBin][ptbin]->GetBinError(b))
          {
            HH->SetBinContent(b,gRnd->Gaus(HH->GetBinContent(b),HH->GetBinError(b)));
          }
        }
        HH->Fit(fitCom[etaBin][ptbin],"RQ");
        //HH->UseCurrentStyle();
        fitCom[etaBin][ptbin]->GetParameters(&par[0]);
        fitCom[etaBin][ptbin]->GetParameters(&parStorage[0]);
        float chi2 = fitCom[etaBin][ptbin]->GetChisquare();
        int ndf = fitCom[etaBin][ptbin]->GetNDF();
        float chiNDF = chi2/(float)ndf;

        for(int ii=0; ii<12; ii++)
        {
          parGaus[ii] -> Fill(par[ii]);
        }
        if(DEBUG) cout << "finish fit num: " << it << endl;

        Double_t parPi[3],parKP[3],parE[3],parmPi[3];
        for(Int_t i=0; i<12; i++)
        {
          if(i < 3)
            parPi[i] = par[i];
          else if(i < 6)
            parKP[i-3] = par[i];
          else if(i < 9)
            parE[i-6] = par[i];
          else if(i<12)
            parmPi[i-9] = par[i];
        }

        fitPiD[etaBin][ptbin]->SetParameters(parPi);
        fitmPiD[etaBin][ptbin]->SetParameters(parmPi);
        fitKPD[etaBin][ptbin]->SetParameters(parKP);
        fitED[etaBin][ptbin]->SetParameters(parE);
        if(DEBUG) cout << "Parameters Set" << endl;

        // Integrate the fits 
        float eInt = fitED[etaBin][ptbin]->Integral(-1,3);
        float piInt = fitPiD[etaBin][ptbin]->Integral(-1,3);
        float mpiInt = fitmPiD[etaBin][ptbin]->Integral(-1,3);
        float kpInt = fitKPD[etaBin][ptbin]->Integral(-1,3);
        float purNum = eInt;
        if(!withMergedPion)
          mpiInt = 0;
        float purDen = eInt+piInt+kpInt+mpiInt;
        float numPerPt = purDen/(highpt[ptbin]-lowpt[ptbin]);
        float purTemp = purNum/purDen;
        purGaus->Fill(purTemp);
        dNdpTGaus->Fill(numPerPt);
        chiGaus->Fill(chiNDF);
      }

      //Retreive central values and errors
      double tempPar[3];
      TF1* purityFind = new TF1("purityFind","gaus");
      TF1* dndptFind = new TF1("dndptFind","gaus");
      TF1* chi2Find = new TF1("chi2Find","gaus");
      TF1* parFind = new TF1("parFind","gaus");
      if(DEBUG) cout << "Before Par Fits" << endl;

      int aPad = ptbin-9;
      purityGaussians[etaBin]->cd(aPad);
      gPad->SetLogy(isLogY);
      purGaus->Draw();
      purityFind->SetLineColor(kRed);
      purGaus->Fit(purityFind,"Q");
      lbl[etaBin][ptbin]->Draw("same");
      purityFind->GetParameters(&tempPar[0]);
      purity[etaBin][ptbin] = tempPar[1]; // mean
      dy[etaBin][ptbin] = tempPar[2];    // sigma
      
      junk->cd();
      dNdpTGaus->Fit(dndptFind,"Q0");
      dndptFind->GetParameters(&tempPar[0]);
      dNdpT[etaBin][ptbin] = tempPar[1]; //mean
      dy2[etaBin][ptbin] = tempPar[2];    // sigma
      chiGaus->Fit(chi2Find,"Q0");
      chi2Find->GetParameters(&tempPar[0]);
      XnDOF[etaBin][ptbin] = tempPar[1]; //mean
      aYield[etaBin][ptbin] = purity[etaBin][ptbin]* dNdpT[etaBin][ptbin];
      aPur[etaBin][ptbin] = 2.*purity[etaBin][ptbin] - 1.;

      if(DEBUG) cout << "Mid Par fits" << endl;
      parameterGaussians[etaBin]->Divide(4,3);
      for(int ii=0; ii<12; ii++)
      {
       // if(ptbin != 3) continue;
        parameterGaussians[etaBin]->cd(ii+1);
        parGaus[ii] -> Fit(parFind,"Q0");
        parGaus[ii]->Draw();
        parFind->GetParameters(&tempPar[0]);
        parPlot[ii][etaBin][ptbin] = tempPar[1]; //store the mean
        errPlot[ii][etaBin][ptbin] = tempPar[2];
        if(DEBUG)  
          cout << "par " << ii << ": " << tempPar[1] << endl;
        par[ii] = tempPar[1];
        parStorage[ii] = tempPar[1];
      }
      if(DEBUG) cout << "After Par fits" << endl;

      // bypass the parameter selction from multi fit version
      fitCom[etaBin][ptbin]->SetParameters(parHold);
      projnSigmaE[etaBin][ptbin]->Fit(fitCom[etaBin][ptbin],"Q0");
      fitCom[etaBin][ptbin]->GetParameters(&par[0]);
      float chi2 = fitCom[etaBin][ptbin]->GetChisquare();
      int ndf = fitCom[etaBin][ptbin]->GetNDF();
      float chiNDF = chi2/(float)ndf;
      XnDOF[etaBin][ptbin] = chiNDF;

      Double_t parPi[3],parKP[3],parE[3],parmPi[3];
      for(Int_t i=0; i<12; i++)
      {
        if(i < 3)
          parPi[i] = par[i];
        else if(i < 6)
          parKP[i-3] = par[i];
        else if(i < 9)
          parE[i-6] = par[i];
        else if(i<12)
          parmPi[i-9] = par[i];
      }
      if(DEBUG) cout << "Assign par vals" << endl;

      nSigE[etaBin][activeCanvas] -> cd(activeBin+1);
      fitCom[etaBin][ptbin]  -> SetParameters(par);
      fitPiD[etaBin][ptbin]  -> SetParameters(parPi);
      fitmPiD[etaBin][ptbin] -> SetParameters(parmPi);
      fitKPD[etaBin][ptbin]  -> SetParameters(parKP);
      fitED[etaBin][ptbin]   -> SetParameters(parE);
      fitPiD[etaBin][ptbin]  -> SetLineStyle(1);
      fitPiD[etaBin][ptbin]  -> SetLineWidth(1);
      fitPiD[etaBin][ptbin]  -> SetLineColor(kRed);
      fitmPiD[etaBin][ptbin] -> SetLineStyle(1);
      fitmPiD[etaBin][ptbin] -> SetLineWidth(1);
      fitmPiD[etaBin][ptbin] -> SetLineColor(kGreen+3);
      fitKPD[etaBin][ptbin]  -> SetLineStyle(1);
      fitKPD[etaBin][ptbin]  -> SetLineWidth(1);
      fitKPD[etaBin][ptbin]  -> SetLineColor(kCyan);
      fitED[etaBin][ptbin]   -> SetLineStyle(1);
      fitED[etaBin][ptbin]   -> SetLineWidth(1);
      fitED[etaBin][ptbin]   -> SetLineColor(kBlue);
      projnSigmaE[etaBin][ptbin]->Draw();
      //HH->Draw();
      fitCom[etaBin][ptbin]  -> Draw("same");
      fitPiD[etaBin][ptbin]  -> Draw("same");
      if(withMergedPion)     fitmPiD[etaBin][ptbin]->Draw("same");
      fitKPD[etaBin][ptbin]  -> Draw("same");
      fitED[etaBin][ptbin]   -> Draw("same");
      lbl[etaBin][ptbin]     -> Draw("same");
      if(DEBUG) cout << "Finish Draw" << endl; 

      // Prepare for TGraphErrors
      // Integrate the fits 
      /* eInte[etaBin][ptbin] = fitED[etaBin][ptbin]->Integral(-1,3);
         piInte[etaBin][ptbin] = fitPiD[etaBin][ptbin]->Integral(-1,3);
         mpiInte[etaBin][ptbin] = fitmPiD[etaBin][ptbin]->Integral(-1,3);
         kpInte[etaBin][ptbin] = fitKPD[etaBin][ptbin]->Integral(-1,3);
         if(DEBUG) cout << "After Integrate" << endl; 

         if(!withMergedPion)
         mpiInte[etaBin][ptbin] = 0;
         sum[etaBin][ptbin] = eInte[etaBin][ptbin]+piInte[etaBin][ptbin]+mpiInte[etaBin][ptbin]+kpInte[etaBin][ptbin];
         dNdpT[etaBin][ptbin] = sum[etaBin][ptbin]/(highpt[ptbin]-lowpt[ptbin]);
         purity[etaBin][ptbin] = eInte[etaBin][ptbin]/sum[etaBin][ptbin];*/
      pT[etaBin][ptbin] = (lowpt[ptbin]+highpt[ptbin])/2.;
      dx[etaBin][ptbin] = (highpt[ptbin]-lowpt[ptbin])/2.;
      // XnDOF[etaBin][ptbin] = chiNDF;

      // Make Stats box legible
      nSigE[etaBin][activeCanvas]->cd(activeBin+1);
      TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
      s->SetX1NDC(0.65);
      s->SetX2NDC(0.95);
      s->SetY1NDC(0.45);
      s->SetY2NDC(0.95);
      nSigE[etaBin][activeCanvas]->Modified();

      nSigE[etaBin][activeCanvas]->Update();
      if(DEBUG) cout << "Stats Modified" << endl;
    }
  }

  // make graphs that need all pt bins
  TGraphErrors* purGr[numEtaBins];
  TGraphErrors* dNdpTGr[numEtaBins];
  TGraphErrors* adjYield[numEtaBins];
  TGraphErrors* adjPurity[numEtaBins];
  TGraphErrors* drawPur[numEtaBins];
  TGraphErrors* drawdNdpT[numEtaBins];
  TGraphErrors* chi2dof[numEtaBins];
  TGraphErrors* drawX2[numEtaBins];
  TGraphErrors* fitPars[12][numEtaBins];
  TString parNames[12] = {"#pi C","#pi #mu","#pi #sigma","Kp C","Kp #mu","Kp #sigma",
    "e C","e #mu","e #sigma","mer. #pi C","mer. #pi #mu","mer. #pi #sigma"};



  for(Int_t etaBin=0; etaBin < numEtaBins; etaBin++)
  {
    if(DEBUG) cout << "TGraphErrors etaBin: " << etaBin << endl;

    purGr[etaBin] = new TGraphErrors(numPtBins,pT[etaBin],purity[etaBin],dx[etaBin],dy[etaBin]);
    adjYield[etaBin] = new TGraphErrors(numPtBins,pT[etaBin],aYield[etaBin],dx[etaBin],dy2[etaBin]);
    adjPurity[etaBin] = new TGraphErrors(numPtBins,pT[etaBin],aPur[etaBin],dx[etaBin],dy[etaBin]);
    dNdpTGr[etaBin] = new TGraphErrors(numPtBins,pT[etaBin],dNdpT[etaBin],dx[etaBin],dy2[etaBin]);
    chi2dof[etaBin] = new TGraphErrors(numPtBins,pT[etaBin],XnDOF[etaBin],dx[etaBin],dy[etaBin]);
    if(DEBUG) cout << "TGraphErrors Assigned" << endl;

    int numParams = 12; if(!withMergedPion) numParams = 9;
    for(int ii=0; ii<numParams; ii++)
    {
      parameterFitCanvas[etaBin] -> cd(ii+1);
      fitPars[ii][etaBin] =  new TGraphErrors(numPtBins,pT[etaBin],parPlot[ii][etaBin],dx[etaBin],errPlot[ii][etaBin]);
      fitPars[ii][etaBin] -> SetMarkerStyle(20);
      fitPars[ii][etaBin] -> SetMarkerSize(0.7);
      fitPars[ii][etaBin] -> SetTitle(parNames[ii]);
      fitPars[ii][etaBin] -> GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fitPars[ii][etaBin] -> GetYaxis()->SetTitle("Par. Value");
      fitPars[ii][etaBin] -> SetMarkerColor(kRed);
      fitPars[ii][etaBin] -> SetLineColor(kRed);
      fitPars[ii][etaBin] -> Draw("AP");
      fitPars[ii][etaBin] -> GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20.);
      fitPars[ii][etaBin] -> Draw("AP");
      fitPars[ii][etaBin] -> SetName(Form("fitpar_%i_%i",etaBin,ii));//Set object name for write to .root
      fitPars[ii][etaBin] -> Write(); // write to .root
    }
    if(DEBUG) cout << "FitPars Hists Created" << endl;

    purC[etaBin]->cd(1);
    purGr[etaBin]->SetMarkerStyle(20);
    purGr[etaBin]->SetMarkerSize(0.7);
    purGr[etaBin]->SetTitle("Electron Sample Purity");
    purGr[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    purGr[etaBin]->GetYaxis()->SetTitle("Purity");
    purGr[etaBin]->GetYaxis()->SetRangeUser(0.0,1.2);
    purGr[etaBin]->SetMarkerColor(kRed);
    purGr[etaBin]->SetLineColor(kRed);
    purGr[etaBin]->Draw("AP");
    purGr[etaBin]->GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20.);
    purGr[etaBin]->Draw("AP");
    purGr[etaBin]->SetName("purity");//Set object name for write to .root
    purGr[etaBin]->Write(); // write to .root
    drawPur[etaBin] = (TGraphErrors*)purGr[etaBin]->Clone();
    drawPur[etaBin]->SetName(Form("drawPurity_%i",etaBin));
    drawPur[etaBin]->Write();
    if(DEBUG) cout << "Purity Graph Done" << endl;

    // Fit Purity
    //purGr->Fit("pol2","Q");
    //TF1* purFit = new TF1("PurityFit","pol2",lowpt[0],8.);//highpt[numPtBins-1]);
    //purGr->Fit(purFit,"RQ");

    purC[etaBin]->cd(2);
 // Handle drawing for option set before actually draw dn/dpt
    adjYield[etaBin]->SetMarkerStyle(20);
    adjYield[etaBin]->SetMarkerSize(0.7);
    adjYield[etaBin]->SetTitle("Electron Sample");
    adjYield[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    adjYield[etaBin]->GetYaxis()->SetTitle("Purity*dN/dpT");
    adjYield[etaBin]->SetMarkerColor(kRed);
    adjYield[etaBin]->SetLineColor(kRed);
    adjYield[etaBin]->Draw("AP");
    adjYield[etaBin]->GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20);
    adjYield[etaBin]->Draw("AP");
    adjYield[etaBin]->SetName(Form("adjYield_%i",etaBin));
    adjYield[etaBin]->Write();

    adjPurity[etaBin]->SetMarkerStyle(20);
    adjPurity[etaBin]->SetMarkerSize(0.7);
    adjPurity[etaBin]->SetTitle("Electron Sample");
    adjPurity[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    adjPurity[etaBin]->GetYaxis()->SetTitle("2*Purity - 1");
    adjPurity[etaBin]->SetMarkerColor(kRed);
    adjPurity[etaBin]->SetLineColor(kRed);
    adjPurity[etaBin]->Draw("AP");
    adjPurity[etaBin]->GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20);
    adjPurity[etaBin]->Draw("AP");
    adjPurity[etaBin]->SetName(Form("adjPurity_%i",etaBin));
    adjPurity[etaBin]->Write();

    gPad->SetLogy(isLogY);
    dNdpTGr[etaBin]->SetMarkerStyle(20);
    dNdpTGr[etaBin]->SetMarkerSize(0.7);
    dNdpTGr[etaBin]->SetTitle("dN/dpT In Electron Range");
    dNdpTGr[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dNdpTGr[etaBin]->GetYaxis()->SetTitle("dN/dpT");
    dNdpTGr[etaBin]->SetMarkerColor(kRed);
    dNdpTGr[etaBin]->SetLineColor(kRed);
    dNdpTGr[etaBin]->Draw("AP");
    dNdpTGr[etaBin]->GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20.);
    dNdpTGr[etaBin]->Draw("AP");
    dNdpTGr[etaBin]->SetName("dNdpT");//Set object name for write to .root
    dNdpTGr[etaBin]->Write(); // write to .root
    drawdNdpT[etaBin] = (TGraphErrors*)dNdpTGr[etaBin]->Clone();
    drawdNdpT[etaBin]->SetName(Form("drawdNdpT_%i",etaBin));
    drawdNdpT[etaBin]->Write();
    if(DEBUG) cout << "dN/dpT Done" << endl;

    xSquare[etaBin]->cd();
    chi2dof[etaBin]->SetMarkerStyle(20);
    chi2dof[etaBin]->SetMarkerSize(0.7);
    chi2dof[etaBin]->SetTitle("Fit Chi2/DoF");
    chi2dof[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    chi2dof[etaBin]->GetYaxis()->SetTitle("#Chi^{2}/DOF");
    chi2dof[etaBin]->SetMarkerColor(kRed);
    chi2dof[etaBin]->SetLineColor(kRed);
    chi2dof[etaBin]->Draw("AP");
    chi2dof[etaBin]->GetXaxis()->SetRangeUser(anaConst::trigThreshold[trig],20.);
    chi2dof[etaBin]->GetYaxis()->SetRangeUser(0,10);
    chi2dof[etaBin]->Draw("AP");
    chi2dof[etaBin]->SetName("purity");//Set object name for write to .root
    chi2dof[etaBin]->Write(); // write to .root
    drawX2[etaBin] = (TGraphErrors*)chi2dof[etaBin]->Clone();
    drawX2[etaBin]->SetName(Form("drawX2_%i",etaBin));
    drawX2[etaBin]->Write();

  }

  // Make PDF with output canvases
  if(makePDF)
  {
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

    TString titlename = FileName;
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
    sprintf(tlName,"eID:                      0.3 < p/E < 1.5; p    _{T} > 1.5 GeV/c; DCA < 1.5 cm");
    tl.DrawLatex(0.1, 0.75,tlName);
    sprintf(tlName,"Track Quality:    nHitsFit > 20; nHits     #frac{dE}{dx}  > 15;");
    tl.DrawLatex(0.1, 0.7,tlName);
    sprintf(tlName,"                            nHitFit/Max > 0.52;    #left|#eta#right| < 0.7;");
    tl.DrawLatex(0.1, 0.65,tlName);
    if(wHFT)
    {
      sprintf(tlName,"                            isHFTTrack();");
      tl.DrawLatex(0.1, 0.6,tlName);
    }
    sprintf(tlName,"Event:            #left|V_{z}#right| < 6 cm; #left|#DeltaV_{z}#right| < 4 cm");
    tl.DrawLatex(0.1, 0.45,tlName);
    sprintf(tlName,"Triggers:            %s %s",FileLabel,Cuts);
    tl.DrawLatex(0.1, 0.4,tlName);


    // Place canvases in order
    TCanvas* temp = new TCanvas();

    sprintf(name, "%s_Eta_%s_%s.pdf[", FileName,FileLabel,Cuts);
    if(wHFT) sprintf(name, "%s_Eta_%s_%s_wHFT.pdf[", FileName,FileLabel,Cuts);
    temp->Print(name);
    sprintf(name, "%s_Eta_%s_%s.pdf", FileName,FileLabel,Cuts);
    if(wHFT) sprintf(name, "%s_Eta_%s_%s_wHFT.pdf", FileName,FileLabel,Cuts);
    temp = fp; // print front page
    temp->Print(name);
    for(int e=0; e<numEtaBins; e++)
    {
      for(int q=0; q<numCanvas; q++)
      {
        temp = nSigE[e][q]; // print data canvases
        temp->Print(name);
      }
      temp = purC[e];
      temp->Print(name);
      temp = xSquare[e];
      temp->Print(name);
      temp = parameterFitCanvas[e];
      temp->Print(name);
      temp = purityGaussians[e];
      temp->Print(name);
    // Remove when bypassing parameter in multi fits
      //  temp = parameterGaussians[e];
      //  temp->Print(name);
    }
    cout << name << " Made." << endl;
    sprintf(name, "%s_Eta_%s_%s.pdf]", FileName,FileLabel,Cuts);
    if(wHFT) sprintf(name, "%s_Eta_%s_%s_wHFT.pdf]", FileName,FileLabel,Cuts);
    temp->Print(name);
  }

  if(DEBUG) cout << "Make Root: " << makeROOT << endl;
  if(makeROOT)
  {
    file->Write();
    file->Close();
  }
}

void checkBatchMode(){

  // sets batch mode, so don't draw canvas
  Int_t number = 2;
  while(number > 1 || number < 0){
    std::cout << "Batch Mode? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ) {
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
        gROOT->SetBatch(kFALSE);
      if(number == 1)
        gROOT->SetBatch(kTRUE);
    }
    else
    {
      number = 1;
      gROOT->SetBatch(kTRUE);
    }
  }
}

Bool_t checkMakePDF(){

  // Set option for pdf creation
  Int_t number = 2; Bool_t fmakePDF = kTRUE;
  while(number > 1 || number < 0){
    std::cout << "Make PDF? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ){
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
        fmakePDF = kFALSE;
      if(number == 1)
        fmakePDF = kTRUE;
    }
    else
      number = 1; 
  }
  if(!fmakePDF)
    cout << "If not making PDF, need to choose a single trigger sample to fix strange plotting behavior!" << endl;
  return fmakePDF;
}

Bool_t checkMakeRoot(){

  // Set option for .root creation
  Int_t number = 2; Bool_t fmakeROOT = kTRUE;
  while(number > 1 || number < 0){
    std::cout << "Make .root? [default: 1]: ";
    std::string input;
    std::getline( std::cin, input );
    if ( !input.empty() ){
      std::istringstream stream( input );
      stream >> number;
      if(number == 0)
        fmakeROOT = kFALSE;
      if(number == 1){
        fmakeROOT = kTRUE;
      }
    }
    else
      number = 1; 
  }
  return fmakeROOT;
}
