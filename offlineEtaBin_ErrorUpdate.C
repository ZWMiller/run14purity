
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////
//// Offline Plots - Z. Miller Jan 6, 2016
//// --------------------
//// Updated 3/3/16  to allow for all cut sets in a single call - ZWM
//// Updated 3/8/16  to produce eta dependent binning in production - ZWM
//// Updated 3/16/16 to use THnSparse input, project to TH3, then fit - ZWM
//// Updated 6/14/16 reads in parFits for constraint mode, has w and wo
////                 constraint modes - ZYe
//// Updated 6/23/16 produces text files with all the required systematics -ZWM
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
//// For full use follow these directions
//// 1. change wConstraint to kFALSE
//// 2. root -l -b -q 'macros/offlineEtaBin_ErrorUpdate.C++("outputs/purityHists_June9",4,"ALL",kFALSE)'
//// 3. root -l -b -q macros/pl_FitParams.C // this will generate fitpar_output.root
//// 4. cd output
//// 5. mkdir noConstraint
//// 6. mv *.pdf *.root *.txt noConstraint
//// 7. mv noConstraint/purity_June9.root .
//// 8. cd ..
//// 9. change kFALSE to kTRUE
//// 10. repeat step 2-5
//// 12. cd output
//// 11. mkdir wConstraint
//// 12. mv *.pdf *.root *.txt wConstraint
//// 13. mv wConstraint/purity_June9.root .
//// 14. cd ..
//// 15. root -l -b -q macros/pl_EtaComparison.C 
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
#include <fstream>
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
void makeHist(const char* FileName="test", Int_t trig=4,const char* Cut="BEMC");
void offlineEtaBin_ErrorUpdate(const char* FileName="test", Int_t trig=4, const char* Cuts="BEMC",Bool_t ishft=kFALSE);
void checkBatchMode();
Bool_t checkMakePDF();
Bool_t checkMakeRoot();
Int_t getLowCentralityLabel(int);
Int_t getHighCentralityLabel(int);
Double_t singleGaussian(Double_t *x, Double_t *par);
Double_t threeGaussian(Double_t *x, Double_t *par);
Double_t fourGaussian(Double_t *x, Double_t *par);
Bool_t makePDF,makeROOT;
int isLogY = 1;
Bool_t drawAll = kFALSE;
Bool_t withMergedPion = kTRUE;
Bool_t DEBUG = kFALSE;
Bool_t wHFT = kFALSE;
Bool_t uFit = kFALSE;
Bool_t wFitConstraint = kTRUE;

const Int_t numEtaBins = anaConst::nEtaBins;
const Int_t numCentBins = anaConst::nCentBins;
const int numparams = 9;
TF1* fitTotal[3][numCentBins][numEtaBins][numparams];

void offline(const char* FileName="test", Int_t trig=4, const char* Cuts="BEMC",Bool_t ishft=kFALSE)
{
  offlineEtaBin_ErrorUpdate(FileName,trig,Cuts,ishft);
}

void offlineEtaBin_ErrorUpdate(const char* FileName, Int_t trig, const char* Cuts,Bool_t ishft) //0=MB,1=HT1,2=HT2,3=HT3,4=ALL, BEMC, SMD, TOF
{
  TFile *outputfile=new TFile("fitpar_output.root","read");
  for(int ier=0; ier<3; ier++) {
    for(int centbin=0;centbin<numCentBins;centbin++) {
      for(int etabin=0;etabin<numEtaBins;etabin++) {
        for(int parnum=0; parnum<numparams; parnum++) {
          if(parnum%3!=0) fitTotal[ier][centbin][etabin][parnum] = (TF1*)outputfile->Get(Form("Fit_%i_%i_%i_%i",ier,centbin,etabin,parnum));
        }
      }
    }
  }

  wHFT = ishft;
  const char* cutTypes[3] = {"BEMC","SMD","SMD2"};
  if(wHFT) cout << "---Running with HFT Match Requirement!---" << endl;
  if(trig == 4)
  {
    gErrorIgnoreLevel = kInfo;
    //checkBatchMode();
    makePDF = true;//checkMakePDF();
    makeROOT= true;//checkMakeRoot();

    if(!strncmp(Cuts,"ALL",3))
    {
      for(int q=0;q<3;q++)
      {
        for(Int_t i=0;i<trig;i++)
          makeHist(FileName,i,cutTypes[q]);
      }
    }
    else
    {
      for(Int_t i=0;i<trig;i++)
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

void makeHist(const char* FileName, Int_t trig,const char* Cut)
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
  sprintf(name,"%s.root",FileName);
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
    sprintf(fname,"%s_Eta_%s_%s_processed.root",FileName,FileLabel,Cuts);
    if(wHFT)
      sprintf(fname,"%s_Eta_%s_%s_wHFT_processed.root",FileName,FileLabel,Cuts);
    file = new TFile(fname,"RECREATE");
    if (file->IsOpen()==kFALSE)
    {
      std::cout << "!!! Outfile Not Opened !!!" << std::endl;
      makeROOT = kFALSE;
    }
  }

  const Int_t numPtBins  = anaConst::nPtBins;
  const Int_t numEtaBins = anaConst::nEtaBins;
  const Int_t numCentBins = anaConst::nCentBins;
  const Int_t numCanvas = numPtBins/9 + 1;// Each eta gets its own pt bins
  Float_t lowpt[numPtBins],highpt[numPtBins];
  Float_t loweta[numEtaBins],higheta[numEtaBins];
  Float_t lowcent[numCentBins],highcent[numCentBins];
  for(Int_t c=0; c< numPtBins; c++){
    lowpt[c] = anaConst::lpt[c];
    highpt[c] = anaConst::hpt[c];
  }
  for(Int_t c=0; c< numEtaBins; c++){
    loweta[c] = anaConst::etaLow[c];
    higheta[c] = anaConst::etaHigh[c];
  }
  for(Int_t c=0; c< numCentBins; c++){
    lowcent[c] = anaConst::centLow[c];
    highcent[c] = anaConst::centHigh[c];
  }
  Float_t hptCut=anaConst::hptCut;
  const Int_t numTrigs = 4;
  Float_t hptMax=anaConst::hptMax; // Set max above range to allow overflow
  Float_t lowPhi=anaConst::lowPhi, highPhi=anaConst::highPhi;
  // Flat constraints
  //float lowLimit[12] =  {      0,-5.,0.8,       0,-8.,0.8,     0,-0.35,0.9,     0,3.,0.8};
  //float highLimit[12] = {1000000,-1.,1.4,10000000,-2.,1.4,100000, -0.2,1.05,100000,8.,2.5};

  // "No" constraints (just force electron mean to be right of pi and Kp
  //float lowLimit[12] = {-100000.,-6.,0.,-100000., -5.,0.0,-100000,-2,0.0,-100000,2.,0.0};
  //float highLimit[12] = {100000.,-4.,100.,100000.,-2.,100,100000,  2,100,100000, 4.,100};

  TPaveText* lbl[numCentBins][numEtaBins][numPtBins];
  char textLabel[100];

  TH1D* projnSigmaE[numCentBins][numEtaBins][numPtBins];
  TH1D* drawnSigmaE[numCentBins][numEtaBins][numPtBins];
  TF1 *fitPi[numCentBins][numEtaBins][numPtBins];
  TF1 *fitKP[numCentBins][numEtaBins][numPtBins];
  TF1 *fitmPi[numCentBins][numEtaBins][numPtBins];
  TF1 *fitE[numCentBins][numEtaBins][numPtBins];
  TF1 *fitCom[numCentBins][numEtaBins][numPtBins];
  TF1 *fitPiD[numCentBins][numEtaBins][numPtBins];
  TF1 *fitmPiD[numCentBins][numEtaBins][numPtBins];
  TF1 *fitKPD[numCentBins][numEtaBins][numPtBins];
  TF1 *fitED[numCentBins][numEtaBins][numPtBins];
  TF1 *fitComD[numCentBins][numEtaBins][numPtBins];
  double eInte[numCentBins][numEtaBins][numPtBins];
  double piInte[numCentBins][numEtaBins][numPtBins];
  double mpiInte[numCentBins][numEtaBins][numPtBins];
  double kpInte[numCentBins][numEtaBins][numPtBins];
  double pT[numCentBins][numEtaBins][numPtBins],dNdpT[numCentBins][numEtaBins][numPtBins],purityFit[numCentBins][numEtaBins][numPtBins], purityCount[numCentBins][numEtaBins][numPtBins],XnDOF[numCentBins][numEtaBins][numPtBins];
  double pTErr[numCentBins][numEtaBins][numPtBins],dNdpTErr[numCentBins][numEtaBins][numPtBins],purityFitErr[numCentBins][numEtaBins][numPtBins],purityCountErr[numCentBins][numEtaBins][numPtBins];
  double parPlot[12][numCentBins][numEtaBins][numPtBins],
         errPlot[12][numCentBins][numEtaBins][numPtBins];

  // Make Canvas
  TCanvas* nSigE[numCentBins][numEtaBins][numCanvas]; 
  TCanvas* purityFitGaussians[numCentBins][numEtaBins][numCanvas];
  TCanvas* purityCountGaussians[numCentBins][numEtaBins][numCanvas];
  TCanvas* purC[numCentBins][numEtaBins];
  TCanvas* xSquare[numCentBins][numEtaBins];
  TCanvas* parameterGaussians[numCentBins][numEtaBins];
  TCanvas* parameterFitCanvas[numCentBins][numEtaBins];
  TCanvas* junk = new TCanvas("junk","junk",50,50,1050,1050);

  for(int centbin = 0; centbin<numCentBins; centbin++){
    for(Int_t etabin=0; etabin < numEtaBins; etabin++){
      for(Int_t q = 0; q < numCanvas; q++)
      {
        nSigE[centbin][etabin][q] = new TCanvas(Form("nSigE_%i_%i_%i",centbin,etabin,q),"nSigma Electron Projections",50,50,1050,1050);
        nSigE[centbin][etabin][q] -> Divide(3,3);
        purityFitGaussians[centbin][etabin][q] = new TCanvas(Form("purityFitGaussians_%i_%i_%i",centbin,etabin,q),"Electron Fit Purity Gaussians",50,50,1050,1050);
        purityFitGaussians[centbin][etabin][q]->Divide(3,3);
        purityCountGaussians[centbin][etabin][q] = new TCanvas(Form("purityCountGaussians_%i_%i_%i",centbin,etabin,q),"Electron Count Purity Gaussians",50,50,1050,1050);
        purityCountGaussians[centbin][etabin][q]->Divide(3,3);
      }
      purC[centbin][etabin] = new TCanvas(Form("purC_%i_%i",centbin,etabin),"Electron Purity",50,50,1050,1050);
      purC[centbin][etabin]->Divide(1,2);
      xSquare[centbin][etabin] = new TCanvas(Form("xSquare_%i_%i",centbin,etabin),"Chi Square Check",50,50,1050,1050);
      parameterGaussians[centbin][etabin] = new TCanvas(Form("parameterGaussians_%i_%i",centbin,etabin),"Electron parameter Gaussians",50,50,1050,1050);
      parameterGaussians[centbin][etabin]->Divide(4,3);
      parameterFitCanvas[centbin][etabin] = new TCanvas(Form("parFitCanvas_%i_%i",centbin,etabin),"Fit Params vs pT",50,50,1050,1050);
      if(withMergedPion)
        parameterFitCanvas[centbin][etabin]->Divide(3,4);
      else
        parameterFitCanvas[centbin][etabin]->Divide(3,3);
    }
  }

  // Make Projections (first get 2d/3d hists, then project)
  THnSparse* nSigmaESparse = (THnSparse*)f->Get(Form("nSigmaE_%s_%i",Cuts,trig));
  if(wHFT)
  {
    nSigmaESparse = (THnSparse*)f->Get(Form("nSigmaE_%s_HFT_%i",Cuts,trig));
  }

  //nSigmaESparse->CalculateErrors();
  for(int centbin = 0; centbin<numCentBins; centbin++)
  {
    TH3D* nSigmaEPtEta;
    nSigmaESparse->GetAxis(3)->SetRange(lowcent[centbin], highcent[centbin]); // Only select 0-80% Centrality (centrality16)
    nSigmaEPtEta = nSigmaESparse->Projection(0,1,2);
    if(DEBUG)cout << "got hists trig "<< trig << endl;

    double parStorage[12] = {0.};
    double parHold[12] = {0.};
    double parErrStorage[12] = {0.};
    for(Int_t etabin=0; etabin<numEtaBins; etabin++)
    {
      for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
      {
        projnSigmaE[centbin][etabin][ptbin]  = (TH1D*)nSigmaEPtEta->ProjectionY(Form("projnSigmaE_%i_%i_%i",centbin,etabin,ptbin),nSigmaEPtEta->GetXaxis()->FindBin(lowpt[ptbin]+0.001),nSigmaEPtEta->GetXaxis()->FindBin(highpt[ptbin]-0.001),nSigmaEPtEta->GetZaxis()->FindBin(loweta[etabin]+0.001),nSigmaEPtEta->GetZaxis()->FindBin(higheta[etabin]-0.001));
      }
    }

    // Analyze the projections
    for(Int_t etabin = numEtaBins-1; etabin >= 0; etabin--)
    {
      parStorage[0] = -1; // reset fitting algorithm for each eta bin
      if(DEBUG) cout << endl << "EtaBin: " << etabin << endl;

      for(Int_t ptbin = 0; ptbin < numPtBins; ptbin++){
        // Clear the variables of interest
        pT[centbin][etabin][ptbin] = dNdpT[centbin][etabin][ptbin] = purityFit[centbin][etabin][ptbin] = purityCount[centbin][etabin][ptbin] = 0;
        pTErr[centbin][etabin][ptbin] = dNdpTErr[centbin][etabin][ptbin] = purityFitErr[centbin][etabin][ptbin] = purityCountErr[centbin][etabin][ptbin] = 0;
        for(int ii=0; ii<12; ii++)
        {
          parPlot[ii][centbin][etabin][ptbin] = errPlot[ii][centbin][etabin][ptbin] = 0;
        }

        // Don't analyze below the trigger turn on curve
        if(lowpt[ptbin] < anaConst::trigThreshold[trig]) continue; 

        // Init necessary plotting tools
        lbl[centbin][etabin][ptbin] = new TPaveText(.5,.2,.85,.3,Form("NB NDC%i",ptbin));
        sprintf(textLabel,"%i - %i%% Centrality",getLowCentralityLabel(lowcent[centbin]),getHighCentralityLabel(highcent[centbin]));
        lbl[centbin][etabin][ptbin]->AddText(textLabel);
        sprintf(textLabel,"%.2f < #eta < %.2f",loweta[etabin],higheta[etabin]);
        lbl[centbin][etabin][ptbin]->AddText(textLabel);
        sprintf(textLabel,"%.2f < p_{T} < %.2f",lowpt[ptbin],highpt[ptbin]);
        lbl[centbin][etabin][ptbin]->AddText(textLabel);
        lbl[centbin][etabin][ptbin]->SetFillColor(kWhite);
        lbl[centbin][etabin][ptbin]->SetTextSize(0.03);

        int activeCanvas = (int)ptbin/9;
        int activeBin = ptbin - activeCanvas*9;
        if(DEBUG)cout << ptbin << ": " << activeCanvas << " " << activeBin << endl;

        nSigE[centbin][etabin][activeCanvas]->cd(activeBin+1);
        gPad->SetLogy(isLogY);
        projnSigmaE[centbin][etabin][ptbin]->SetMarkerColor(kBlack);
        projnSigmaE[centbin][etabin][ptbin]->SetMarkerStyle(20);
        projnSigmaE[centbin][etabin][ptbin]->SetLineColor(kBlack);
        projnSigmaE[centbin][etabin][ptbin]->SetMarkerSize(0.5);
        projnSigmaE[centbin][etabin][ptbin]->GetXaxis()->SetRangeUser(-11.,16.);
        projnSigmaE[centbin][etabin][ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
        projnSigmaE[centbin][etabin][ptbin]->SetTitle("");
        drawnSigmaE[centbin][etabin][ptbin] = (TH1D*)projnSigmaE[centbin][etabin][ptbin]->Clone();
        drawnSigmaE[centbin][etabin][ptbin]->SetName(Form("drawnSigmaE_%i_%i_%i",centbin,etabin,ptbin));
        drawnSigmaE[centbin][etabin][ptbin]->Write();
        projnSigmaE[centbin][etabin][ptbin]->Draw();
        lbl[centbin][etabin][ptbin]->Draw("same");
        if(DEBUG) cout << "After Clone projnSigmaE" << endl;

        // Do fits on nSig E
        //      nSigE[etabin][activeCanvas]->cd(activeBin+1);
        gStyle->SetOptFit(1111);
        Double_t par[12];
        Double_t parErr[12];
        Double_t parS[3];
        Double_t kpLow,kpHigh,piLow,piHigh,eLow,eHigh, mpiLow, mpiHigh;
        kpLow = -10.; kpHigh =-7.;
        piLow = -5. ; piHigh = -2.;
        eLow = -1.0  ; eHigh = 1.0;
        mpiLow = 3.5; mpiHigh = 6;
        fitKP[centbin][etabin][ptbin] = new TF1(Form("fitKP_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,kpLow,kpHigh,3);
        fitPi[centbin][etabin][ptbin] = new TF1(Form("fitPi_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,piLow,piHigh,3);
        fitE[centbin][etabin][ptbin]  = new TF1(Form("fitE_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,eLow,eHigh,3);
        fitmPi[centbin][etabin][ptbin]= new TF1(Form("fitmPi_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,mpiLow,mpiHigh,3);
        if(withMergedPion)
          fitCom[centbin][etabin][ptbin]= new TF1(Form("fitCom_%i_%i_%i",centbin,etabin,ptbin),fourGaussian,-10,10,12);//"gaus(0)+gaus(3)+gaus(6)+gaus(9)",-10,10);
        else
          fitCom[centbin][etabin][ptbin]= new TF1(Form("fitCom_%i_%i_%i",centbin,etabin,ptbin),threeGaussian,-10,10,9);//"gaus(0)+gaus(3)+gaus(6)",-10,10);

        fitCom[centbin][etabin][ptbin]->SetParName(0,"#pi C");
        fitCom[centbin][etabin][ptbin]->SetParName(1,"#pi #mu");
        fitCom[centbin][etabin][ptbin]->SetParName(2,"#pi #sigma");
        fitCom[centbin][etabin][ptbin]->SetParName(3,"Kp C");
        fitCom[centbin][etabin][ptbin]->SetParName(4,"Kp #mu");
        fitCom[centbin][etabin][ptbin]->SetParName(5,"Kp #sigma");
        fitCom[centbin][etabin][ptbin]->SetParName(6,"e C");
        fitCom[centbin][etabin][ptbin]->SetParName(7,"e #mu");
        fitCom[centbin][etabin][ptbin]->SetParName(8,"e #sigma");
        fitCom[centbin][etabin][ptbin]->SetParName(9, "mer#pi C");
        fitCom[centbin][etabin][ptbin]->SetParName(10,"mer#pi #mu");
        fitCom[centbin][etabin][ptbin]->SetParName(11,"mer#pi #sigma");

        fitPi[centbin][etabin][ptbin]->SetLineColor(kRed);
        fitKP[centbin][etabin][ptbin]->SetLineColor(kCyan);
        fitE[centbin][etabin][ptbin]->SetLineColor(kBlue);
        fitmPi[centbin][etabin][ptbin]->SetLineColor(kGreen);
        fitCom[centbin][etabin][ptbin]->SetLineColor(kMagenta);
        fitPi[centbin][etabin][ptbin]->SetLineWidth(1);
        fitmPi[centbin][etabin][ptbin]->SetLineWidth(1);
        fitKP[centbin][etabin][ptbin]->SetLineWidth(1);
        fitE[centbin][etabin][ptbin]->SetLineWidth(1);
        fitCom[centbin][etabin][ptbin]->SetLineWidth(2);

        // make versions to draw fits in full
        fitPiD[centbin][etabin][ptbin] = new TF1(Form("drawPi_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,-10,10,3);
        fitmPiD[centbin][etabin][ptbin]= new TF1(Form("drawmPi_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,-10,10,3);
        fitKPD[centbin][etabin][ptbin] = new TF1(Form("drawKP_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,-10,10,3);
        fitED[centbin][etabin][ptbin]  = new TF1(Form("drawE_%i_%i_%i",centbin,etabin,ptbin),singleGaussian,-10,10,3);

        // If this is the first bin of the trigger, fit the roughly
        // expected areas for each Gaussian to get starting pars for total fit
        if(fabs(parStorage[0]+1)<1e-6)
        {
          projnSigmaE[centbin][etabin][ptbin]->Fit(fitPi[centbin][etabin][ptbin],"RQ");
          projnSigmaE[centbin][etabin][ptbin]->Fit(fitKP[centbin][etabin][ptbin],"RQ");
          projnSigmaE[centbin][etabin][ptbin]->Fit(fitE[centbin][etabin][ptbin],"RQ");
          if(withMergedPion)
            projnSigmaE[centbin][etabin][ptbin]->Fit(fitmPi[centbin][etabin][ptbin],"RQ");

          fitPi[centbin][etabin][ptbin]->GetParameters(&par[0]);
          fitKP[centbin][etabin][ptbin]->GetParameters(&par[3]);
          fitE[centbin][etabin][ptbin]->GetParameters(&par[6]);
          if(withMergedPion)
            fitmPi[centbin][etabin][ptbin]->GetParameters(&par[9]);
          fitCom[centbin][etabin][ptbin]->SetParameters(par);
          fitCom[centbin][etabin][ptbin]->GetParameters(&parStorage[0]);
        }
        else // Otherwise, use results from last fit as input
          fitCom[centbin][etabin][ptbin]->SetParameters(parStorage);
        fitCom[centbin][etabin][ptbin]->GetParameters(&parHold[0]);

        // Apply parameter constraints
        for(int q=0; q<12; q++)
        {
          fitCom[centbin][etabin][ptbin]->SetParLimits(q,anaConst::lowLimit[q],anaConst::highLimit[q]);
          if(wFitConstraint&&(q==1||q==4||q==7)) {
            cout<<0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])<<" "<<fitTotal[2][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin]))<<" "<<fitTotal[0][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin]))<<endl;
            fitCom[centbin][etabin][ptbin]->SetParLimits(q,fitTotal[2][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])),fitTotal[0][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])));
          }
          if(wFitConstraint&&(q==2||q==8)) {
            cout<<0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])<<" "<<fitTotal[0][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin]))<<" "<<fitTotal[2][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin]))<<endl;
            fitCom[centbin][etabin][ptbin]->SetParLimits(q,fitTotal[0][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])),fitTotal[2][centbin][etabin][q]->Eval(0.5*(anaConst::lpt[ptbin]+anaConst::hpt[ptbin])));
          }
        }
        /*//  fitCom[etabin][ptbin]->SetParLimits(0,0,1000000);
          fitCom[etabin][ptbin]->SetParLimits(1,-6,-4);
          fitCom[etabin][ptbin]->SetParLimits(2,0.8,1.2);
        //  fitCom[etabin][ptbin]->SetParLimits(3,0,1000000);
        fitCom[etabin][ptbin]->SetParLimits(4,-5,-2);
        fitCom[etabin][ptbin]->SetParLimits(5,0.8,1.2);
        //  fitCom[etabin][ptbin]->SetParLimits(6,0,1000000);
        fitCom[etabin][ptbin]->SetParLimits(7,-0.5,0.0);
        fitCom[etabin][ptbin]->SetParLimits(8,0.8,1.2);
        //  fitCom[etabin][ptbin]->SetParLimits(9,0,1000000);
        //  fitCom[etabin][ptbin]->SetParLimits(10,2,4);
        //  fitCom[etabin][ptbin]->SetParLimits(11,0.8,1.5);*/

        // Do fits 1000 times on slightly randomized data to get mean and error of purity
        TH1F* dNdpTGaus = new TH1F("dNdpTGaus","dNdpT Results for All Fits",100,0,2000000);
        TH1F* purFitGaus = new TH1F("purFitGaus","Purity Fit Results for All Fits",1000,0,2);
        TH1F* purCountGaus = new TH1F("purCountGaus","Purity Count Results for All Fits",1000,0,2);
        TH1F* parGaus[12];
        TH1F* chiGaus = new TH1F("chiGaus","Chi2/NDF Results for All Fits",500,-10,10);
        for(int ii=0; ii<12; ii++)
        {
          if(ii == 0 || ii == 3 || ii == 6 || ii == 9)
            parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1600,0,8e6);
          else if(ii == 1 || ii == 4 || ii == 7 || ii == 10)
            parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1000,-10,10);
          else if(ii == 2 || ii == 5 || ii == 8 || ii == 11)
            parGaus[ii] = new TH1F(Form("parGaus_%i",ii),Form("Par %i Results for All Fits",ii),1000,0,4);
        }  
        TRandom3 *gRnd= new TRandom3(0);
        TH1D* HH;
        for(int it=0; it<10; it++)
        {
          fitCom[centbin][etabin][ptbin]->SetParameters(parStorage);
          HH = (TH1D*)projnSigmaE[centbin][etabin][ptbin]->Clone();
          int bins = HH->GetNbinsX();
          for(int b=1; b<=bins; b++)
          {
            if(projnSigmaE[centbin][etabin][ptbin]->GetBinContent(b) && projnSigmaE[centbin][etabin][ptbin]->GetBinError(b))
            {
              HH->SetBinContent(b,gRnd->Gaus(HH->GetBinContent(b),HH->GetBinError(b)));
            }
          }
          HH->Fit(fitCom[centbin][etabin][ptbin],"BRQ");

          //HH->UseCurrentStyle();
          fitCom[centbin][etabin][ptbin]->GetParameters(&par[0]);
          //fitCom[etabin][ptbin]->GetParameters(&parStorage[0]);
          float chi2 = fitCom[centbin][etabin][ptbin]->GetChisquare();
          int ndf = fitCom[centbin][etabin][ptbin]->GetNDF();
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

          fitPiD[centbin][etabin][ptbin]->SetParameters(parPi);
          fitmPiD[centbin][etabin][ptbin]->SetParameters(parmPi);
          fitKPD[centbin][etabin][ptbin]->SetParameters(parKP);
          fitED[centbin][etabin][ptbin]->SetParameters(parE);
          if(DEBUG) cout << "Parameters Set" << endl;

          // Integrate the fits 
          float eInt = fitED[centbin][etabin][ptbin]->Integral(-1,3);
          float piInt = fitPiD[centbin][etabin][ptbin]->Integral(-1,3);
          float mpiInt = fitmPiD[centbin][etabin][ptbin]->Integral(-1,3);
          float kpInt = fitKPD[centbin][etabin][ptbin]->Integral(-1,3);
          if(!withMergedPion) mpiInt = 0;
          float purNum = eInt/HH->GetBinWidth(1);
          float numPerPt = purNum/(highpt[ptbin]-lowpt[ptbin]);
          float purFitDen = (eInt+piInt+kpInt+mpiInt)/HH->GetBinWidth(1);
          float purCountDen = HH->Integral(HH->GetXaxis()->FindBin(-1+0.001),HH->GetXaxis()->FindBin(3-0.001));
          float purFitTemp = purNum/purFitDen;
          float purCountTemp = purNum/purCountDen;
          dNdpTGaus->Fill(numPerPt);
          purFitGaus->Fill(purFitTemp);
          purCountGaus->Fill(purCountTemp);
          chiGaus->Fill(chiNDF);
        }

        //Retreive central values and errors
        double tempPar[3];
        TF1* dNdpTFind = new TF1("dNdpTFind","gaus");
        TF1* purityFitFind = new TF1("purityFitFind","gaus");
        TF1* purityCountFind = new TF1("purityCountFind","gaus");
        TF1* chi2Find = new TF1("chi2Find","gaus");
        TF1* parFind = new TF1("parFind","gaus");
        if(DEBUG) cout << "Before Par Fits" << endl;

        int aCanvas = (int)ptbin/9;
        int aPad = ptbin-aCanvas*9;
        purityFitGaussians[centbin][etabin][aCanvas]->cd(aPad+1);
        gPad->SetLogy(isLogY);
        purFitGaus->Draw();
        purityFitFind->SetLineColor(kRed);
        purFitGaus->Fit(purityFitFind,"Q");
        lbl[centbin][etabin][ptbin]->Draw("same");
        purityFitFind->GetParameters(&tempPar[0]);
        purityFit[centbin][etabin][ptbin] = tempPar[1]; // mean
        purityFitErr[centbin][etabin][ptbin] = tempPar[2];    // sigma
        purityCountGaussians[centbin][etabin][aCanvas]->cd(aPad+1);
        gPad->SetLogy(isLogY);
        purCountGaus->Draw();
        purityCountFind->SetLineColor(kRed);
        purCountGaus->Fit(purityCountFind,"Q");
        lbl[centbin][etabin][ptbin]->Draw("same");
        purityCountFind->GetParameters(&tempPar[0]);
        purityCount[centbin][etabin][ptbin] = tempPar[1]; // mean
        purityCountErr[centbin][etabin][ptbin] = tempPar[2];    // sigma
        junk->cd();
        dNdpTGaus->Fit(dNdpTFind,"Q0");
        dNdpTFind->GetParameters(&tempPar[0]);
        dNdpT[centbin][etabin][ptbin] = tempPar[1]; //mean
        dNdpTErr[centbin][etabin][ptbin] = tempPar[2];    // sigma
        chiGaus->Fit(chi2Find,"Q0");
        chi2Find->GetParameters(&tempPar[0]);
        XnDOF[centbin][etabin][ptbin] = tempPar[1]; //mean
        if(!uFit) {
          purityFit[centbin][etabin][ptbin] = purFitGaus->GetMean(); // mean
          purityFitErr[centbin][etabin][ptbin] = purFitGaus->GetRMS();    // sigma
          purityCount[centbin][etabin][ptbin] = purCountGaus->GetMean(); // mean
          purityCountErr[centbin][etabin][ptbin] = purCountGaus->GetRMS();    // sigma
          dNdpT[centbin][etabin][ptbin] = dNdpTGaus->GetMean(); //mean
          dNdpTErr[centbin][etabin][ptbin] = dNdpTGaus->GetRMS();    // sigma
          XnDOF[centbin][etabin][ptbin] = chiGaus->GetMean(); //mean
        }

        if(DEBUG) cout << "Mid Par fits" << endl;
        parameterGaussians[centbin][etabin]->Divide(4,3);
        for(int ii=0; ii<12; ii++)
        {
          // if(ptbin != 3) continue;
          parameterGaussians[centbin][etabin]->cd(ii+1);
          parGaus[ii] -> Fit(parFind,"Q0");
          parGaus[ii]->Draw();
          parFind->GetParameters(&tempPar[0]);
          parPlot[ii][centbin][etabin][ptbin] = tempPar[1]; //store the mean
          errPlot[ii][centbin][etabin][ptbin] = tempPar[2];
          par[ii] = tempPar[1];
          parStorage[ii] = tempPar[1];
          if(!uFit) {
            parPlot[ii][centbin][etabin][ptbin] = parGaus[ii]->GetMean(); //store the mean
            errPlot[ii][centbin][etabin][ptbin] = parGaus[ii]->GetRMS();
            par[ii] = parGaus[ii]->GetMean();
            parStorage[ii] = parGaus[ii]->GetMean();
          }
          if(DEBUG)
            cout << "par " << ii << ": " << parPlot[ii][centbin][etabin][ptbin] << endl;
        }
        if(DEBUG) cout << "After Par fits" << endl;

        // bypass the parameter selction from multi fit version
        cout << endl << "---PTBIN--- " << lowpt[ptbin] << " -> " << highpt[ptbin] << "GeV/c" << endl;

        fitCom[centbin][etabin][ptbin]->SetParameters(parHold);
        projnSigmaE[centbin][etabin][ptbin]->Fit(fitCom[centbin][etabin][ptbin],"B0R");
        fitCom[centbin][etabin][ptbin]->GetParameters(&par[0]);
        float tempErr;
        for(int ii=0; ii<12; ii++)
        {
          tempErr = fitCom[centbin][etabin][ptbin]->GetParError(ii);
          if(ii%3 == 0) parPlot[ii][centbin][etabin][ptbin] = par[ii]/(highpt[ptbin]-lowpt[ptbin]); //store the fit result
          else parPlot[ii][centbin][etabin][ptbin] = par[ii]; //store the fit result
          errPlot[ii][centbin][etabin][ptbin] = tempErr; //store the fit error
        }
        float chi2 = fitCom[centbin][etabin][ptbin]->GetChisquare();
        int ndf = fitCom[centbin][etabin][ptbin]->GetNDF();
        float chiNDF = chi2/(float)ndf;
        XnDOF[centbin][etabin][ptbin] = chiNDF;

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

        nSigE[centbin][etabin][activeCanvas] -> cd(activeBin+1);
        fitCom[centbin][etabin][ptbin]  -> SetParameters(par);
        fitPiD[centbin][etabin][ptbin]  -> SetParameters(parPi);
        fitmPiD[centbin][etabin][ptbin] -> SetParameters(parmPi);
        fitKPD[centbin][etabin][ptbin]  -> SetParameters(parKP);
        fitED[centbin][etabin][ptbin]   -> SetParameters(parE);
        fitPiD[centbin][etabin][ptbin]  -> SetLineStyle(1);
        fitPiD[centbin][etabin][ptbin]  -> SetLineWidth(1);
        fitPiD[centbin][etabin][ptbin]  -> SetLineColor(kRed);
        fitmPiD[centbin][etabin][ptbin] -> SetLineStyle(1);
        fitmPiD[centbin][etabin][ptbin] -> SetLineWidth(1);
        fitmPiD[centbin][etabin][ptbin] -> SetLineColor(kGreen+3);
        fitKPD[centbin][etabin][ptbin]  -> SetLineStyle(1);
        fitKPD[centbin][etabin][ptbin]  -> SetLineWidth(1);
        fitKPD[centbin][etabin][ptbin]  -> SetLineColor(kCyan);
        fitED[centbin][etabin][ptbin]   -> SetLineStyle(1);
        fitED[centbin][etabin][ptbin]   -> SetLineWidth(1);
        fitED[centbin][etabin][ptbin]   -> SetLineColor(kBlue);
        projnSigmaE[centbin][etabin][ptbin]->Draw();
        //HH->Draw();
        fitCom[centbin][etabin][ptbin]  -> Draw("same");
        fitPiD[centbin][etabin][ptbin]  -> Draw("same");
        if(withMergedPion)     fitmPiD[centbin][etabin][ptbin]->Draw("same");
        fitKPD[centbin][etabin][ptbin]  -> Draw("same");
        fitED[centbin][etabin][ptbin]   -> Draw("same");
        lbl[centbin][etabin][ptbin]     -> Draw("same");
        if(DEBUG) cout << "Finish Draw" << endl; 

        // Prepare for TGraphErrors
        // Integrate the fits 
        /* eInte[etabin][ptbin] = fitED[etabin][ptbin]->Integral(-1,3);
           piInte[etabin][ptbin] = fitPiD[etabin][ptbin]->Integral(-1,3);
           mpiInte[etabin][ptbin] = fitmPiD[etabin][ptbin]->Integral(-1,3);
           kpInte[etabin][ptbin] = fitKPD[etabin][ptbin]->Integral(-1,3);
           if(DEBUG) cout << "After Integrate" << endl; 

           if(!withMergedPion)
           mpiInte[etabin][ptbin] = 0;
           sum[etabin][ptbin] = eInte[etabin][ptbin]+piInte[etabin][ptbin]+mpiInte[etabin][ptbin]+kpInte[etabin][ptbin];
           dNdpT[etabin][ptbin] = sum[etabin][ptbin]/(highpt[ptbin]-lowpt[ptbin]);
           purity[etabin][ptbin] = eInte[etabin][ptbin]/sum[etabin][ptbin];*/
        pT[centbin][etabin][ptbin] = (lowpt[ptbin]+highpt[ptbin])/2.;
        pTErr[centbin][etabin][ptbin] = (highpt[ptbin]-lowpt[ptbin])/2.;
        // XnDOF[etabin][ptbin] = chiNDF;

        // Make Stats box legible
        nSigE[centbin][etabin][activeCanvas]->cd(activeBin+1);
        TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
        s->SetX1NDC(0.65);
        s->SetX2NDC(0.95);
        s->SetY1NDC(0.45);
        s->SetY2NDC(0.95);
        nSigE[centbin][etabin][activeCanvas]->Modified();

        nSigE[centbin][etabin][activeCanvas]->Update();
        if(DEBUG) cout << "Stats Modified" << endl;
      }
    }

    // make graphs that need all pt bins
    TGraphErrors* dNdpTGr[numCentBins][numEtaBins];
    TGraphErrors* purFitGr[numCentBins][numEtaBins];
    TGraphErrors* purCountGr[numCentBins][numEtaBins];
    TGraphErrors* drawdNdpT[numCentBins][numEtaBins];
    TGraphErrors* drawPurFit[numCentBins][numEtaBins];
    TGraphErrors* drawPurCount[numCentBins][numEtaBins];
    TGraphErrors* chi2dof[numCentBins][numEtaBins];
    TGraphErrors* drawX2[numCentBins][numEtaBins];
    TGraphErrors* fitPars[12][numCentBins][numEtaBins];
    TString parNames[12] = {"#pi C","#pi #mu","#pi #sigma","Kp C","Kp #mu","Kp #sigma",
      "e C","e #mu","e #sigma","mer. #pi C","mer. #pi #mu","mer. #pi #sigma"};

    for(Int_t etabin=0; etabin < numEtaBins; etabin++)
    {
      if(DEBUG) cout << "TGraphErrors etabin: " << etabin << endl;

      dNdpTGr[centbin][etabin] = new TGraphErrors(numPtBins,pT[centbin][etabin],dNdpT[centbin][etabin],pTErr[centbin][etabin],dNdpTErr[centbin][etabin]);
      purFitGr[centbin][etabin] = new TGraphErrors(numPtBins,pT[centbin][etabin],purityFit[centbin][etabin],pTErr[centbin][etabin],purityFitErr[centbin][etabin]);
      purCountGr[centbin][etabin] = new TGraphErrors(numPtBins,pT[centbin][etabin],purityCount[centbin][etabin],pTErr[centbin][etabin],purityCountErr[centbin][etabin]);
      chi2dof[centbin][etabin] = new TGraphErrors(numPtBins,pT[centbin][etabin],XnDOF[centbin][etabin],pTErr[centbin][etabin],0);
      if(DEBUG) cout << "TGraphErrors Assigned" << endl;

      // for each eta bin, make a text file with pT, purityCount, purityFit to calculate systematics from
      ofstream systematicsFile;
      systematicsFile.open (Form("outputs/systematicsTextFiles/systematicInformation_%s_%s_CentBin%i_EtaBin%i.txt",FileLabel,Cuts,centbin,etabin));
      for(int ptbin = 0; ptbin<numPtBins; ptbin++)
      {
        systematicsFile << pT[centbin][etabin][ptbin] << " " << pTErr[centbin][etabin][ptbin] << " " << purityFit[centbin][etabin][ptbin] << " " << purityFitErr[centbin][etabin][ptbin] 
          << " " << purityCount[centbin][etabin][ptbin] << " " << purityCountErr[centbin][etabin][ptbin] << endl;
      }
      systematicsFile.close();

      int numParams = 12; if(!withMergedPion) numParams = 9;
      for(int ii=0; ii<numParams; ii++)
      {
        parameterFitCanvas[centbin][etabin] -> cd(ii+1);
        for(int ptbin = 0; ptbin < numPtBins; ptbin++)
        {
          if(DEBUG)cout << "par " << ii << " pt: " << lowpt[ptbin] << " val: " << parPlot[ii][centbin][etabin][ptbin] << endl;
        }
        fitPars[ii][centbin][etabin] =  new TGraphErrors(numPtBins,pT[centbin][etabin],parPlot[ii][centbin][etabin],pTErr[centbin][etabin],errPlot[ii][centbin][etabin]);
        fitPars[ii][centbin][etabin] -> SetMarkerStyle(20);
        fitPars[ii][centbin][etabin] -> SetMarkerSize(0.7);
        fitPars[ii][centbin][etabin] -> SetTitle(parNames[ii]);
        fitPars[ii][centbin][etabin] -> GetXaxis()->SetTitle("p_{T} (GeV/c)");
        if(ii%3==0)fitPars[ii][centbin][etabin] -> GetYaxis()->SetTitle("Par. Value/Bin Width");
        else fitPars[ii][centbin][etabin] -> GetYaxis()->SetTitle("Par. Value");
        fitPars[ii][centbin][etabin] -> SetMarkerColor(kRed);
        fitPars[ii][centbin][etabin] -> SetLineColor(kRed);
        fitPars[ii][centbin][etabin] -> GetXaxis()->SetRangeUser(anaConst::lowPt,anaConst::highPt);
        if(ii==1||ii==4) fitPars[ii][centbin][etabin] -> GetYaxis()->SetRangeUser(-10,0);
        if(ii==7) fitPars[ii][centbin][etabin] -> GetYaxis()->SetRangeUser(-1,1);
        if(ii==10) fitPars[ii][centbin][etabin] -> GetYaxis()->SetRangeUser(2,9);
        if(ii==2||ii==5||ii==8) fitPars[ii][centbin][etabin] -> GetYaxis()->SetRangeUser(0.7,1.5);
        if(ii==11) fitPars[ii][centbin][etabin] -> GetYaxis()->SetRangeUser(0,3.5);
        fitPars[ii][centbin][etabin] -> Draw("AP");
        fitPars[ii][centbin][etabin] -> SetName(Form("fitpar_%i_%i_%i",centbin,etabin,ii));//Set object name for write to .root
        TLine* lineLow = new TLine(anaConst::lowPt,anaConst::lowLimit[ii],anaConst::highPt,anaConst::lowLimit[ii]);
        lineLow->SetLineColor(kBlue);
        TLine* lineHigh = new TLine(anaConst::lowPt,anaConst::highLimit[ii],anaConst::highPt,anaConst::highLimit[ii]);
        lineHigh->SetLineColor(kBlue);
        lineLow->Draw("same");
        lineHigh->Draw("same");
        fitPars[ii][centbin][etabin] -> Write(); // write to .root
      }

      if(DEBUG) cout << "FitPars Hists Created" << endl;

      purC[centbin][etabin]->cd(1);
      purFitGr[centbin][etabin]->SetMarkerStyle(20);
      purFitGr[centbin][etabin]->SetMarkerSize(0.7);
      purFitGr[centbin][etabin]->SetTitle("Electron Sample Purity");
      purFitGr[centbin][etabin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      purFitGr[centbin][etabin]->GetYaxis()->SetTitle("Purity");
      purFitGr[centbin][etabin]->GetYaxis()->SetRangeUser(0.0,1.2);
      purFitGr[centbin][etabin]->SetMarkerColor(kRed);
      purFitGr[centbin][etabin]->SetLineColor(kRed);
      purFitGr[centbin][etabin]->GetXaxis()->SetRangeUser(anaConst::lowPt,anaConst::highPt);
      purFitGr[centbin][etabin]->Draw("AP");
      purFitGr[centbin][etabin]->SetName(Form("purityFit_%i_%i",centbin,etabin));//Set object name for write to .root
      purFitGr[centbin][etabin]->Write(); // write to .root
      drawPurFit[centbin][etabin] = (TGraphErrors*)purFitGr[centbin][etabin]->Clone();
      drawPurFit[centbin][etabin]->SetName(Form("drawPurityFit_%i_%i",centbin,etabin));
      drawPurFit[centbin][etabin]->Write();
      if(DEBUG) cout << "Purity Fit Graph Done" << endl;

      purC[centbin][etabin]->cd(1);
      purCountGr[centbin][etabin]->SetMarkerStyle(20);
      purCountGr[centbin][etabin]->SetMarkerSize(0.7);
      purCountGr[centbin][etabin]->SetTitle("Electron Sample Purity");
      purCountGr[centbin][etabin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      purCountGr[centbin][etabin]->GetYaxis()->SetTitle("Purity");
      purCountGr[centbin][etabin]->GetYaxis()->SetRangeUser(0.0,1.2);
      purCountGr[centbin][etabin]->SetMarkerColor(kBlue);
      purCountGr[centbin][etabin]->SetLineColor(kBlue);
      purCountGr[centbin][etabin]->GetXaxis()->SetRangeUser(anaConst::lowPt,anaConst::highPt);
      purCountGr[centbin][etabin]->Draw("PSAME");
      purCountGr[centbin][etabin]->SetName(Form("purityCount_%i_%i",centbin,etabin));//Set object name for write to .root
      purCountGr[centbin][etabin]->Write(); // write to .root
      drawPurCount[centbin][etabin] = (TGraphErrors*)purCountGr[centbin][etabin]->Clone();
      drawPurCount[centbin][etabin]->SetName(Form("drawPurityCount_%i_%i",centbin,etabin));
      drawPurCount[centbin][etabin]->Write();
      if(DEBUG) cout << "Purity Count Graph Done" << endl;

      // Fit Purity
      //purGr->Fit("pol2","Q");
      //TF1* purFit = new TF1("PurityFit","pol2",lowpt[0],8.);//highpt[numPtBins-1]);
      //purGr->Fit(purFit,"RQ");

      purC[centbin][etabin]->cd(2);
      // Handle drawing for option set before actually draw dn/dpt
      gPad->SetLogy(isLogY);
      dNdpTGr[centbin][etabin]->SetMarkerStyle(20);
      dNdpTGr[centbin][etabin]->SetMarkerSize(0.7);
      dNdpTGr[centbin][etabin]->SetTitle("dN/dpT In Electron Range");
      dNdpTGr[centbin][etabin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      dNdpTGr[centbin][etabin]->GetYaxis()->SetTitle("dN/dpT");
      dNdpTGr[centbin][etabin]->SetMarkerColor(kRed);
      dNdpTGr[centbin][etabin]->SetLineColor(kRed);
      dNdpTGr[centbin][etabin]->GetXaxis()->SetRangeUser(anaConst::lowPt,anaConst::highPt);
      dNdpTGr[centbin][etabin]->Draw("AP");
      dNdpTGr[centbin][etabin]->SetName(Form("dNdpT_%i_%i",centbin,etabin));//Set object name for write to .root
      dNdpTGr[centbin][etabin]->Write(); // write to .root
      drawdNdpT[centbin][etabin] = (TGraphErrors*)dNdpTGr[centbin][etabin]->Clone();
      drawdNdpT[centbin][etabin]->SetName(Form("drawdNdpT_%i_%i",centbin,etabin));
      drawdNdpT[centbin][etabin]->Write();
      if(DEBUG) cout << "dN/dpT Done" << endl;

      xSquare[centbin][etabin]->cd();
      chi2dof[centbin][etabin]->SetMarkerStyle(20);
      chi2dof[centbin][etabin]->SetMarkerSize(0.7);
      chi2dof[centbin][etabin]->SetTitle("Fit Chi2/DoF");
      chi2dof[centbin][etabin]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      chi2dof[centbin][etabin]->GetYaxis()->SetTitle("#Chi^{2}/DOF");
      chi2dof[centbin][etabin]->SetMarkerColor(kRed);
      chi2dof[centbin][etabin]->SetLineColor(kRed);
      chi2dof[centbin][etabin]->GetXaxis()->SetRangeUser(anaConst::lowPt,anaConst::highPt);
      chi2dof[centbin][etabin]->GetYaxis()->SetRangeUser(0,10);
      chi2dof[centbin][etabin]->Draw("AP");
      chi2dof[centbin][etabin]->SetName(Form("chi2dof_%i_%i",centbin,etabin));//Set object name for write to .root
      chi2dof[centbin][etabin]->Write(); // write to .root
      drawX2[centbin][etabin] = (TGraphErrors*)chi2dof[centbin][etabin]->Clone();
      drawX2[centbin][etabin]->SetName(Form("drawX2_%i_%i",centbin,etabin));
      drawX2[centbin][etabin]->Write();
    }
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
    sprintf(tlName,"Triggers:            %s",FileLabel);
    tl.DrawLatex(0.1, 0.8,tlName);
    if(wHFT)
    {
      sprintf(tlName,"                            isHFTTrack();");
      tl.DrawLatex(0.1, 0.75,tlName);
      sprintf(tlName,"Event:       #left|V_{z}#right| < 6 cm; #left|#DeltaV_{z}#right| < 4 cm");
    }
    else {
      sprintf(tlName,"Event:       #left|V_{z}#right| < 30 cm; #left|#DeltaV_{z}#right| < 4 cm");
    }
    tl.DrawLatex(0.1, 0.7,tlName);
    sprintf(tlName,"Kinematic Cut:    pT > 1.5 GeV/c; #left|#eta#right| < 0.7");
    tl.DrawLatex(0.1, 0.65,tlName);
    sprintf(tlName,"Track Quality:     DCA < 1.5 cm; nHitsFit > 20; nHitsdEdx > 15; nHitFit/Max > 0.52");
    tl.DrawLatex(0.1, 0.6,tlName);
    sprintf(tlName,"eID:                      0.3 < p/E < 1.5, %s",Cuts);
    tl.DrawLatex(0.1, 0.55,tlName);


    // Place canvases in order
    TCanvas* temp = new TCanvas();

    sprintf(name, "%s_Eta_%s_%s.pdf[", FileName,FileLabel,Cuts);
    if(wHFT) sprintf(name, "%s_Eta_%s_%s_wHFT.pdf[", FileName,FileLabel,Cuts);
    temp->Print(name);
    sprintf(name, "%s_Eta_%s_%s.pdf", FileName,FileLabel,Cuts);
    if(wHFT) sprintf(name, "%s_Eta_%s_%s_wHFT.pdf", FileName,FileLabel,Cuts);
    temp = fp; // print front page
    temp->Print(name);
    for(int c=0; c<numCentBins; c++)
    {
    for(int e=0; e<numEtaBins; e++)
    {
      for(int q=0; q<numCanvas; q++)
      {
        temp = nSigE[c][e][q]; // print data canvases
        temp->Print(name);
      }
      temp = parameterFitCanvas[c][e];
      temp->Print(name);
      temp = xSquare[c][e];
      temp->Print(name);
      for(int q=0; q<numCanvas; q++)
      {
        temp = purityFitGaussians[c][e][q]; // print data canvases
        temp->Print(name);
      }
      for(int q=0; q<numCanvas; q++)
      {
        temp = purityCountGaussians[c][e][q]; // print data canvases
        temp->Print(name);
      }
      temp = purC[c][e];
      temp->Print(name);
      // Remove when bypassing parameter in multi fits
      //  temp = parameterGaussians[e];
      //  temp->Print(name);
    }
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


Double_t singleGaussian(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f = par[0]/2/3.14159267/par[2]*TMath::Exp(-0.5*TMath::Power((xx-par[1])/par[2],2));
  return f;
}

Double_t threeGaussian(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f = singleGaussian(x,&par[0])+singleGaussian(x,&par[3])+singleGaussian(x,&par[6]);
  return f;
}

Double_t fourGaussian(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f = singleGaussian(x,&par[0])+singleGaussian(x,&par[3])+singleGaussian(x,&par[6])+singleGaussian(x,&par[9]);
  return f;
}

Int_t getLowCentralityLabel(int x)
{
  return 75-x*5;
}
Int_t getHighCentralityLabel(int x)
{
  return 80-x*5;
}
