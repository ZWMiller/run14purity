
#include "anaConst14.h"

void setRangeAndDraw(TGraphErrors*, int,TCanvas*);
void findRejectPoints(int, int, int);
double myFunc(double*, double*);
int numFitPars = 0;
float rejectPointHigh[10], rejectPointLow[10];
int rejectPointLength;
void setNumFitPars();

int pl_FitParams(){
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
  const int numparams = 9;

  // create generic structures for use in all plots
  TPaveText* lbl[numEtaBins][numparams];
  TPaveText* lblE[numEtaBins];
  TPaveText* runInfo[numEtaBins];
  TPaveText* runInfo2[numEtaBins];
  char textLabel[100]; 
  TString parNames[12] = {" #pi C"," #pi #mu"," #pi #sigma"," Kp C"," Kp #mu"," Kp #sigma",
    " e C"," e #mu"," e #sigma"," mer. #pi C"," mer. #pi #mu"," mer. #pi #sigma"};

  // Load all trigger type files
  const char* baseName = "outputs/purityHists_Mar29";
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
    BEMC[q] = new TFile(BEMCName[q],"READ");
    sprintf(SMDLName[q],"%s_Eta_%s_SMD_processed.root",baseName,trigName[q]);
    SMDL[q] = new TFile(SMDLName[q],"READ");
    sprintf(SMDTName[q],"%s_Eta_%s_SMD2_processed.root",baseName,trigName[q]);
    SMDT[q] = new TFile(SMDTName[q],"READ");
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
  TGraphErrors* parBEMC[4][numEtaBins][numparams];
  TGraphErrors* parSMDL[4][numEtaBins][numparams];
  TGraphErrors* parSMDT[4][numEtaBins][numparams];
  TF1* fitBEMC[4][numEtaBins][numparams];
  TF1* fitSMDL[4][numEtaBins][numparams];
  TF1* fitSMDT[4][numEtaBins][numparams];

  float fitRangeLow[4] = {3,2,3.5,4.5};
  float fitRangeHigh[4] = {7,18,18,18};

  setNumFitPars();
  const int nFitPars = numFitPars;
  float parStart[numparams][4] = {
    {0,0,0,0},
    {8.59,-3.276,8.769,12.4},
    {0.058,4.02,1.31,-1.04},
    {0,0,0,0},
    {36.3,-15.77,3.25,38.29},
    {-0.114,3.36,0.94,-1.06},
    {0,0,0,0},
    {0.065,4.74,-0.26,0.435},
    {0.058,4.02,1.31,-1.04},
  };

  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int etabin=0; etabin<numEtaBins; etabin++)
    {
      for(int parnum=0; parnum<numparams; parnum++)
      {
        parBEMC[w][etabin][parnum] = (TGraphErrors*)BEMC[w]->Get(Form("fitpar_%i_%i",etabin,parnum));
        parSMDL[w][etabin][parnum] = (TGraphErrors*)SMDL[w]->Get(Form("fitpar_%i_%i",etabin,parnum));
        parSMDT[w][etabin][parnum] = (TGraphErrors*)SMDT[w]->Get(Form("fitpar_%i_%i",etabin,parnum)); 

        fitBEMC[w][etabin][parnum] = new TF1(Form("BEMCFit_%i_%i_%i",w,etabin,parnum),myFunc,fitRangeLow[w],fitRangeHigh[w],numFitPars);
        fitSMDL[w][etabin][parnum] = new TF1(Form("SMDLFit_%i_%i_%i",w,etabin,parnum),myFunc,fitRangeLow[w],fitRangeHigh[w],numFitPars);
        fitSMDT[w][etabin][parnum] = new TF1(Form("SMDTFit_%i_%i_%i",w,etabin,parnum),myFunc,fitRangeLow[w],fitRangeHigh[w],numFitPars);

        for(int ll=0;ll<4;ll++)
        {
          fitBEMC[w][etabin][parnum]->SetParameter(ll,parStart[parnum][ll]);
          fitSMDL[w][etabin][parnum]->SetParameter(ll,parStart[parnum][ll]);
          fitSMDT[w][etabin][parnum]->SetParameter(ll,parStart[parnum][ll]);
        }
        //fitBEMC[w][etabin][parnum]->SetLineColor(kRed);
        //fitSMDL[w][etabin][parnum]->SetLineColor(kRed);
        //fitSMDT[w][etabin][parnum]->SetLineColor(kRed);
        //fitBEMC[w][etabin][parnum]->DrawClone();
        //fitSMDL[w][etabin][parnum]->DrawClone();
        //fitSMDT[w][etabin][parnum]->DrawClone();
        fitBEMC[w][etabin][parnum]->SetLineColor(kBlack);
        fitSMDL[w][etabin][parnum]->SetLineColor(kBlack);
        fitSMDT[w][etabin][parnum]->SetLineColor(kBlack);
        parBEMC[w][etabin][parnum]->SetMarkerSize(1.5);
        parSMDL[w][etabin][parnum]->SetMarkerSize(1.5);
        parSMDT[w][etabin][parnum]->SetMarkerSize(1.5);


      }
    }
  }

  if(DEBUG)
    cout << "Found Hists." << endl;

  // Make Canvas
  TCanvas* paramCanv[3][4][numEtaBins];
  if(DEBUG)
    cout << "Canvas Made." << endl;

  // Plot Settings
  for(int r=1; r<4; r++)
  {
    if(DEBUG) cout << "plot settings trig: " << r << endl;

    for(int etabin=0;etabin<numEtaBins;etabin++)
    {
      if(DEBUG) cout << "plot settings eta: " << etabin << endl;
      for(int j=0; j<3; j++)
      {
        paramCanv[j][r][etabin]= new TCanvas(Form("paramCanv_%i_%i_%i",j,r,etabin),"Parameter Fits",50,50,1050,1050);
        paramCanv[j][r][etabin]->Divide(3,3);
      }
    }
  }

  if(DEBUG)
    cout << "Settings Assigned." << endl;
// Regular constraints
 // float lowLimit[12] = {0.,-6.,0.8,0.,-5.,0.8,0.,-0.5,0.8,0.,2.,0.8};
 // float highLimit[12] = {100000,-4.,1.2,100000,-2.,1.2,10000,-0.12,1.2,100000,4.,1.5};
  // "No" constraints (just force electron mean to be right of pi and Kp
  float lowLimit[12] = {-100000.,-6.,0.,-100000., -5.,0.0,-100000,-2,0.0,-100000,2.,0.0};
  float highLimit[12] = {100000.,-4.,100.,100000.,-2.,100,100000,  2,100,100000, 4.,100};
  for(int t=1; t<4; t++)
  {
    if(DEBUG) cout << "Trig: " << t << endl;
    for(int etabin=0;etabin<numEtaBins;etabin++)
    {
      for(int parnum=0; parnum<numparams; parnum++)
      {
        TLine* lineLow = new TLine(anaConst::trigThreshold[t],lowLimit[parnum],20,lowLimit[parnum]);
        lineLow->SetLineColor(kBlue);
        TLine* lineHigh = new TLine(anaConst::trigThreshold[t],highLimit[parnum],20,highLimit[parnum]);
        lineHigh->SetLineColor(kBlue);

        TString title = trigName[t] + parNames[parnum];
        if(DEBUG) cout << "BEMC Par Num: " << parnum << endl;
        parBEMC[t][etabin][parnum]->SetTitle("BEMC "+title);
        setRangeAndDraw(parBEMC[t][etabin][parnum],parnum, paramCanv[0][t][etabin]);
        findRejectPoints(t,1,parnum);
        fitBEMC[t][etabin][parnum]->SetLineColor(kMagenta);
        //fitBEMC[t][etabin][parnum]->DrawClone("same");
        fitBEMC[t][etabin][parnum]->SetLineColor(kBlack);
        cout << "Fit: BEMC, Trig: " << t << ", Par: " << parnum << endl; 
        if(!((parnum%3) == 0))
          parBEMC[t][etabin][parnum]->Fit(fitBEMC[t][etabin][parnum],"R");
        lineLow->Draw("same");
        lineHigh->Draw("same");

        if(DEBUG) cout << "SMDL Par Num: " << parnum << endl;
        parSMDL[t][etabin][parnum]->SetTitle("SMD "+title);
        setRangeAndDraw(parSMDL[t][etabin][parnum],parnum,paramCanv[1][t][etabin]);
        findRejectPoints(t,2,parnum);
        fitSMDL[t][etabin][parnum]->SetLineColor(kMagenta);
        //fitSMDL[t][etabin][parnum]->DrawClone("same");
        fitSMDL[t][etabin][parnum]->SetLineColor(kBlack);
        cout << "Fit: SMDL, Trig: " << t << ", Par: " << parnum << endl; 
        if(!((parnum%3) == 0))
          parSMDL[t][etabin][parnum]->Fit(fitSMDL[t][etabin][parnum],"R");
        lineLow->Draw("same");
        lineHigh->Draw("same");

        if(DEBUG) cout << "SMDT Par Num: " << parnum << endl;
        parSMDT[t][etabin][parnum]->SetTitle("SMD2 "+title);
        setRangeAndDraw(parSMDT[t][etabin][parnum],parnum,paramCanv[2][t][etabin]);
        findRejectPoints(t,3,parnum);
        fitSMDT[t][etabin][parnum]->SetLineColor(kMagenta);
        //fitSMDT[t][etabin][parnum]->DrawClone("same");
        fitSMDT[t][etabin][parnum]->SetLineColor(kBlack);
        cout << "Fit: SMDT, Trig: " << t << ", Par: " << parnum << endl; 
        if(!((parnum%3) == 0))
          parSMDT[t][etabin][parnum]->Fit(fitSMDT[t][etabin][parnum],"R");
        lineLow->Draw("same");
        lineHigh->Draw("same");

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
  sprintf(name, "%s_paramFits_EtaBins.pdf[", baseName);
  temp->Print(name);
  sprintf(name, "%s_paramFits_EtaBins.pdf", baseName);
  temp = fp; // print front page
  temp->Print(name);

  for(int j=0; j<3; j++)
  {
    for(int q=1; q<4; q++)
    {
      if(DEBUG) cout << "Place Canvas in PDF trig: " << q << endl;
      for(int e=0; e<numEtaBins; e++)
      {  
        if(DEBUG) cout << "Place Canvas in PDF eta: " << e << endl;
        temp = paramCanv[j][q][e]; // print data canvases
        temp->Print(name);
      }
    }
  }

  sprintf(name, "%s_paramFits_EtaBins.pdf]", baseName);
  temp->Print(name);


  return 1;
}

void setRangeAndDraw(TGraphErrors* gr, int parnum,TCanvas* can)
{
  can->cd(parnum+1);
  gr->Draw();
  switch (parnum)
  {
    case 1:
      gr->GetYaxis()->SetRangeUser(-7,-3);
      break;
    case 2:
      gr->GetYaxis()->SetRangeUser(0.6,1.3);
      break;
    case 4:
      gr->GetYaxis()->SetRangeUser(-6,-1);
      break;
    case 5:
      gr->GetYaxis()->SetRangeUser(0.6,1.3);
      break;
    case 7:
      gr->GetYaxis()->SetRangeUser(-1,0.5);
      break;
    case 8:
      gr->GetYaxis()->SetRangeUser(0.6,1.3);
      break;
  }
  gr->Draw("ape");
}

double myFunc(double* x, double* par)
{
  for(int ii=0; ii<rejectPointLength; ii++)
  {
    if(x[0] > rejectPointLow[ii] && x[0] < rejectPointHigh[ii])
    {
      TF1::RejectPoint();
      return 0;
    }
  }
  double ff = par[0]*TMath::Erf((x[0]-par[1])/par[2])-par[3];
  return ff;
}

void setNumFitPars()
{
  numFitPars = 4;
}

void findRejectPoints(int tr, int cut, int par)
{
  // cut 1=BEMC,2=SMD,3=SMD2
  if(cut < 1 or cut > 3)
  { 
    cout << "Cuts not correct for Reject Points" << endl; 
    return;
  }

  if((par % 3) == 0)
    return;
  for(int ii = 0; ii < rejectPointLength; ii++)
  {
    rejectPointLow[ii] = rejectPointHigh[ii] = 0.;
  }
  rejectPointLength = 0;

  /// Reject point length tells how many regions are being rejected.
  /// Must be set to 1 if any regions are being rejected.
  if(cut == 1 && tr == 1)
  {   
    if(par == 1 || par == 7)
    {
      rejectPointLow[0]  = 10.;
      rejectPointHigh[0]  = 13.;
      rejectPointLength  = 1;
    }
    if(par == 8)
    {
      rejectPointLow[0]   = 0.;
      rejectPointHigh[0]  = 4.;
      rejectPointLength  = 1;
    }  
  }
  if(cut == 1 && tr == 2)
    if(par == 7)
    {
      rejectPointLow[0]   = 12.;
      rejectPointHigh[0]  = 20.;
      rejectPointLength = 1;
    }


  if(cut == 2 && tr == 1)
  {
    if(par == 1 || par == 7)
    {
      rejectPointLow[0]   = 10.;
      rejectPointHigh[0]  = 13.;
      rejectPointLow[1]   = 0.;
      rejectPointHigh[1]  = 2.5;
      rejectPointLength = 2;
    }
  }
  if(cut == 3 && tr == 1)
  {   
    if(par == 1 || par == 7)
    {
      rejectPointLow[0]   = 10.;
      rejectPointHigh[0]  = 13.;
      rejectPointLow[1]   = 0.;
      rejectPointHigh[1]  = 3.5;
      rejectPointLength = 2;
    }
    if(par == 2 || par == 5)
    {
      rejectPointLow[0]   = 10.;
      rejectPointHigh[0]  = 20.;
      rejectPointLow[1]   = 0.;
      rejectPointHigh[1]  = 2.5;
      rejectPointLength = 2;
    }
    if(par == 4)
    {
      rejectPointLow[0]   = 0.;
      rejectPointHigh[0]  = 2.5;
      rejectPointLength = 1;
    }
  }
  if(cut == 3 && tr == 2)
  {
    if(par == 5)
    {
      rejectPointLow[0]   = 6.;
      rejectPointHigh[0]  = 7.;
      rejectPointLength = 1;
    }
  }

  return;
}
