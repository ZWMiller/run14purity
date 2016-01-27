// Offline Plots - Z. Miller Jan 6, 2016
//
// .L offline.C
// offline("FILENAME", triggerType) # Without .root Extension

#include "anaConst14.h"

// Declare functions
void makeHist(const char*,Int_t);
void checkBatchMode();
Bool_t checkMakePDF();
Bool_t checkMakeRoot();
Bool_t makePDF,makeROOT;
Bool_t DEBUG = kFALSE;
int isLogY = 1;

void offline(const char* FileName="test", Int_t trig=4) //0=MB,1=HT1,2=HT2,3=HT3,4=ALL
{
  if(trig == 4)
  {
    checkBatchMode();
    makePDF = checkMakePDF();
    makeROOT= checkMakeRoot();

    for(Int_t i=0;i<trig;i++)
      makeHist(FileName,i);
  }
  else if(trig > 4)
    cout << "Bad Trigger Input, 0-4 only" << endl;
  else
  {
    checkBatchMode();
    makePDF = checkMakePDF();
    makeROOT= checkMakeRoot();
    makeHist(FileName,trig);
  }
}

void makeHist(const char* FileName="test", Int_t trig=4)
{
  TH1F::SetDefaultSumw2();
  // Set Style parameters for this macro
  gStyle->SetOptTitle(1); // Show Title (off by default for cleanliness)
  gErrorIgnoreLevel = kError; // Set Verbosity Level (kPrint shows all)


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
  if(makeROOT){
    sprintf(fname,"/Users/zach/Research/rootFiles/run14NPEpurity/%s_%s_processed.root",FileName,FileLabel);
    file = new TFile(fname,"RECREATE");
    if (file->IsOpen()==kFALSE)
    {
      std::cout << "!!! Outfile Not Opened !!!" << std::endl;
      makeROOT = kFALSE;
    }
  }

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

  TPaveText* lbl[numPtBins];
  TPaveText* pulbl[numPtBins];
  TPaveText* stat[numPtBins];
  char textLabel[100];

  // Make Canvas

  TCanvas* nSigPi  = new TCanvas("nSigPi","nSigma Pion Projections",50,50,1050,1050);
  TCanvas* nSigK   = new TCanvas("nSigP","nSigma Kaon Projections",50,50,1050,1050);
  TCanvas* nSigP   = new TCanvas("nSigK","nSigma Proton Projections",50,50,1050,1050);
  TCanvas* nSigE   = new TCanvas("nSigE","nSigma Electron Projections",50,50,1050,1050);
  TCanvas* nSigEK  = new TCanvas("nSigEK","nSigma Electron K enh Projections",50,50,1050,1050);
  TCanvas* nSigEP  = new TCanvas("nSigEP","nSigma Electron p enh Projections",50,50,1050,1050);
  TCanvas* nSigEPi = new TCanvas("nSigEPi","nSigma Electron pi enh Projections",50,50,1050,1050);
  TCanvas* invBeta = new TCanvas("invBeta","Beta^(-1) Projections",50,50,1050,1050);
  TCanvas* purC = new TCanvas("purC","Electron Purity",50,50,1050,1050);
  nSigPi  -> Divide(3,3);
  nSigK   -> Divide(3,3);
  nSigP   -> Divide(3,3);
  nSigE   -> Divide(3,3);
  invBeta -> Divide(3,3);

  // Make Projections (first get 2d/3d hists, then project)

  TH2F* nSigmaPiPt  = (TH2F*)f->Get(Form("nSigmaPI_Pt_BEMC_%i",trig));
  TH2F* nSigmaKPt   = (TH2F*)f->Get(Form("nSigmaK_Pt_BEMC_%i",trig));
  TH2F* nSigmaPPt   = (TH2F*)f->Get(Form("nSigmaP_Pt_BEMC_%i",trig));
  TH2F* nSigmaEPt   = (TH2F*)f->Get(Form("nSigmaE_Pt_BEMC_%i",trig));
  TH2F* invBetaPt   = (TH2F*)f->Get(Form("invsBeta_Pt_%i",trig));
  TH2F* nSigmaEKPt  = (TH2F*)f->Get(Form("nSigmaE_KEnh_Pt_%i",trig));
  TH2F* nSigmaEPPt  = (TH2F*)f->Get(Form("nSigmaE_PEnh_Pt_%i",trig));
  TH2F* nSigmaEPiPt = (TH2F*)f->Get(Form("nSigmaE_PiEnh_Pt_%i",trig));
  if(DEBUG)cout << "got hists" << endl;


  TH1D* projnSigmaPi[numPtBins];
  TH1D* projnSigmaK[numPtBins];
  TH1D* projnSigmaP[numPtBins];
  TH1D* projnSigmaE[numPtBins];
  TH1D* projnSigmaEK[numPtBins];
  TH1D* projnSigmaEP[numPtBins];
  TH1D* projnSigmaEPi[numPtBins];
  TH1D* projinvBeta[numPtBins];
  TF1 *fitPi[numPtBins];
  TF1 *fitKP[numPtBins];
  TF1 *fitmPi[numPtBins];
  TF1 *fitE[numPtBins];
  TF1 *fitCom[numPtBins];
  TF1 *fitPiD[numPtBins];
  TF1 *fitmPiD[numPtBins];
  TF1 *fitKPD[numPtBins];
  TF1 *fitED[numPtBins];
  TF1 *fitComD[numPtBins];
  double eInte[numPtBins];
  double piInte[numPtBins];
  double mpiInte[numPtBins];
  double kpInte[numPtBins];
  double sum[numPtBins], purity[numPtBins],pT[numPtBins];
  double dx[numPtBins],dy[numPtBins];
  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
  {
    projnSigmaPi[ptbin] = nSigmaPiPt->ProjectionY(Form("projnSigmaPi_%i",ptbin),nSigmaPiPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaPiPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaK[ptbin]  = nSigmaKPt->ProjectionY(Form("projnSigmaK_%i",ptbin),nSigmaKPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaKPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaP[ptbin]  = nSigmaPPt->ProjectionY(Form("projnSigmaP_%i",ptbin),nSigmaPPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaPPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaE[ptbin]  = nSigmaEPt->ProjectionY(Form("projnSigmaE_%i",ptbin),nSigmaEPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projinvBeta[ptbin]  = invBetaPt->ProjectionY(Form("invBeta_%i",ptbin),invBetaPt->GetXaxis()->FindBin(lowpt[ptbin]),invBetaPt->GetXaxis()->FindBin(highpt[ptbin])-1);

    projnSigmaEK[ptbin]  = nSigmaEKPt->ProjectionY(Form("projnSigmaEK_%i",ptbin),nSigmaEKPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEKPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaEP[ptbin]  = nSigmaEPPt->ProjectionY(Form("projnSigmaEP_%i",ptbin),nSigmaEPPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEPPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaEPi[ptbin]  = nSigmaEPiPt->ProjectionY(Form("projnSigmaEPi_%i",ptbin),nSigmaEPiPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEPiPt->GetXaxis()->FindBin(highpt[ptbin])-1);
  }

  // Analyze the projections
  for(Int_t ptbin = 0; ptbin < numPtBins; ptbin++){

    // Init necessary plotting tools
    lbl[ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);

    nSigPi->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaPi[ptbin]->SetMarkerColor(kBlack);
    projnSigmaPi[ptbin]->SetMarkerStyle(20);
    projnSigmaPi[ptbin]->SetMarkerSize(0.5);
    projnSigmaPi[ptbin]->SetLineColor(kBlack);
    projnSigmaPi[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaPi[ptbin]->SetTitle("");
    projnSigmaPi[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigK->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaK[ptbin]->SetMarkerColor(kBlack);
    projnSigmaK[ptbin]->SetMarkerStyle(20);
    projnSigmaK[ptbin]->SetMarkerSize(0.5);
    projnSigmaK[ptbin]->SetLineColor(kBlack);
    projnSigmaK[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaK[ptbin]->SetTitle("");
    projnSigmaK[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigP->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaP[ptbin]->SetMarkerColor(kBlack);
    projnSigmaP[ptbin]->SetMarkerStyle(20);
    projnSigmaP[ptbin]->SetMarkerSize(0.5);
    projnSigmaP[ptbin]->SetLineColor(kBlack);
    projnSigmaP[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaP[ptbin]->SetTitle("");
    projnSigmaP[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigE->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaE[ptbin]->SetMarkerColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerStyle(20);
    projnSigmaE[ptbin]->SetLineColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerSize(0.5);
    projnSigmaE[ptbin]->GetXaxis()->SetRangeUser(-11.,16.);
    projnSigmaE[ptbin]->SetTitle("");
    projnSigmaE[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigEK->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaEK[ptbin]->SetMarkerColor(kBlack);
    projnSigmaEK[ptbin]->SetMarkerStyle(20);
    projnSigmaEK[ptbin]->SetLineColor(kBlack);
    projnSigmaEK[ptbin]->SetMarkerSize(0.5);
    projnSigmaEK[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaEK[ptbin]->SetTitle("");
    projnSigmaEK[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigEP->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaEP[ptbin]->SetMarkerColor(kBlack);
    projnSigmaEP[ptbin]->SetMarkerStyle(20);
    projnSigmaEP[ptbin]->SetLineColor(kBlack);
    projnSigmaEP[ptbin]->SetMarkerSize(0.5);
    projnSigmaEP[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaEP[ptbin]->SetTitle("");
    projnSigmaEP[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigEPi->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projnSigmaEPi[ptbin]->SetMarkerColor(kBlack);
    projnSigmaEPi[ptbin]->SetMarkerStyle(20);
    projnSigmaEPi[ptbin]->SetLineColor(kBlack);
    projnSigmaEPi[ptbin]->SetMarkerSize(0.5);
    projnSigmaEPi[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaEPi[ptbin]->SetTitle("");
    projnSigmaEPi[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    // Do fits on nSig E
    nSigE->cd(ptbin+1);
    gStyle->SetOptFit(1111);
    Double_t par[12];
    Double_t parS[3];
    Double_t kpLow,kpHigh,piLow,piHigh,eLow,eHigh, mpiLow, mpiHigh;
    kpLow = -10.; kpHigh =-6.; 
    piLow = -5. ; piHigh = -2.;
    eLow = -1.  ; eHigh = 3.;
    mpiLow = 3.5; mpiHigh = 5.5;
    fitKP[ptbin] = new TF1(Form("fitKP_%i",ptbin),"gaus",kpLow,kpHigh);
    fitPi[ptbin] = new TF1(Form("fitPi_%i",ptbin),"gaus",piLow,piHigh);
    fitE[ptbin]  = new TF1(Form("fitE_%i",ptbin),"gaus",eLow,eHigh);
    fitmPi[ptbin]= new TF1(Form("fitmPi_%i",ptbin),"gaus",mpiLow,mpiHigh);
    fitCom[ptbin]= new TF1(Form("fitCom_%i",ptbin),"gaus(0)+gaus(3)+gaus(6)+gaus(9)",-10,10);
    fitCom[ptbin]->SetParName(0,"#pi C");
    fitCom[ptbin]->SetParName(1,"#pi #mu");
    fitCom[ptbin]->SetParName(2,"#pi #sigma");
    fitCom[ptbin]->SetParName(3,"Kp C");
    fitCom[ptbin]->SetParName(4,"Kp #mu");
    fitCom[ptbin]->SetParName(5,"Kp #sigma");
    fitCom[ptbin]->SetParName(6,"e C");
    fitCom[ptbin]->SetParName(7,"e #mu");
    fitCom[ptbin]->SetParName(8,"e #sigma");
    fitCom[ptbin]->SetParName(9, "mer#pi C");
    fitCom[ptbin]->SetParName(10,"mer#pi #mu");
    fitCom[ptbin]->SetParName(11,"mer#pi #sigma");
    fitCom[ptbin]->SetParLimits(0,0,1000000);
    fitCom[ptbin]->SetParLimits(2,0,1000000);
    fitCom[ptbin]->SetParLimits(3,0,1000000);
    fitCom[ptbin]->SetParLimits(5,0,1000000);
    fitCom[ptbin]->SetParLimits(6,0,1000000);
    fitCom[ptbin]->SetParLimits(8,0,1000000);
    fitCom[ptbin]->SetParLimits(9,0,1000000);
    fitCom[ptbin]->SetParLimits(11,0,1000000);
    fitPi[ptbin]->SetLineColor(kRed);
    fitKP[ptbin]->SetLineColor(kCyan);
    fitE[ptbin]->SetLineColor(kBlue);
    fitmPi[ptbin]->SetLineColor(kGreen);
    fitCom[ptbin]->SetLineColor(kMagenta);
    fitPi[ptbin]->SetLineWidth(1);
    fitmPi[ptbin]->SetLineWidth(1);
    fitKP[ptbin]->SetLineWidth(1);
    fitE[ptbin]->SetLineWidth(1);
    fitCom[ptbin]->SetLineWidth(2);
    projnSigmaE[ptbin]->Fit(fitPi[ptbin],"R");
    projnSigmaE[ptbin]->Fit(fitKP[ptbin],"R");
    projnSigmaE[ptbin]->Fit(fitE[ptbin],"R");
    projnSigmaE[ptbin]->Fit(fitmPi[ptbin],"R");
    fitPi[ptbin]->GetParameters(&par[0]);
    fitKP[ptbin]->GetParameters(&par[3]);
    fitE[ptbin]->GetParameters(&par[6]);
    fitmPi[ptbin]->GetParameters(&par[9]);
    fitCom[ptbin]->SetParameters(par);
    projnSigmaE[ptbin]->Fit(fitCom[ptbin],"R");
    projnSigmaE[ptbin]->UseCurrentStyle();
    fitCom[ptbin]->GetParameters(&par[0]);

    if(DEBUG) cout << "finish fit" << endl;

    // make versions to draw fits in full
    fitPiD[ptbin] = new TF1(Form("drawPi_%i",ptbin),"gaus",-10,10);
    fitmPiD[ptbin]= new TF1(Form("drawmPi_%i",ptbin),"gaus",-10,10);
    fitKPD[ptbin] = new TF1(Form("drawKP_%i",ptbin),"gaus",-10,10);
    fitED[ptbin] = new TF1(Form("drawE_%i",ptbin),"gaus",-10,10);
    Double_t parPi[3],parKP[3],parE[3],parmPi[3];
    for(Int_t i=0; i<9; i++)
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
    fitPiD[ptbin]->SetParameters(parPi);
    fitmPiD[ptbin]->SetParameters(parmPi);
    fitKPD[ptbin]->SetParameters(parKP);
    fitED[ptbin]->SetParameters(parE);
    if(DEBUG) cout << "Parameters Set" << endl;
    fitPiD[ptbin]->SetLineStyle(2);
    fitPiD[ptbin]->SetLineWidth(1);
    fitPiD[ptbin]->SetLineColor(kRed);
    fitmPiD[ptbin]->SetLineStyle(2);
    fitmPiD[ptbin]->SetLineWidth(1);
    fitmPiD[ptbin]->SetLineColor(kBlack);
    fitKPD[ptbin]->SetLineStyle(2);
    fitKPD[ptbin]->SetLineWidth(1);
    fitKPD[ptbin]->SetLineColor(kCyan);
    fitED[ptbin]->SetLineStyle(2);
    fitED[ptbin]->SetLineWidth(1);
    fitED[ptbin]->SetLineColor(kBlue);
    fitPiD[ptbin]->Draw("same");
    fitmPiD[ptbin]->Draw("same");
    fitKPD[ptbin]->Draw("same");
    fitED[ptbin]->Draw("same");
    lbl[ptbin]->Draw("same");
    if(DEBUG) cout << "Finish Draw" << endl; 

    // Integrate the fits 
    eInte[ptbin] = fitED[ptbin]->Integral(-1,3);
    piInte[ptbin] = fitPiD[ptbin]->Integral(-1,3);
    mpiInte[ptbin] = fitmPiD[ptbin]->Integral(-1,3);
    cout << "mPi Inte: " << mpiInte[ptbin] << endl;
    kpInte[ptbin] = fitKPD[ptbin]->Integral(-1,3);
    sum[ptbin] = eInte[ptbin]+piInte[ptbin]+mpiInte[ptbin]+kpInte[ptbin];
    purity[ptbin] = eInte[ptbin]/sum[ptbin];
    pT[ptbin] = (lowpt[ptbin]+highpt[ptbin])/2.;
    dx[ptbin] = (highpt[ptbin]-lowpt[ptbin])/2.;
    dy[ptbin] = 1/eInte[ptbin]*purity[ptbin];


    // Make Stats box legible
    /*nSigE->cd(ptbin+1);
    TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
    s->SetX1NDC(0.65);
    s->SetX2NDC(0.95);
    s->SetY1NDC(0.45);
    s->SetY2NDC(0.95);
    nSigE->Modified();
    nSigE->Update();
    if(DEBUG) cout << "Stats Modified" << endl;*/

    invBeta->cd(ptbin+1);
    gPad->SetLogy(isLogY);
    projinvBeta[ptbin]->SetMarkerColor(kBlack);
    projinvBeta[ptbin]->SetMarkerStyle(20);
    projinvBeta[ptbin]->SetLineColor(kBlack);
    projinvBeta[ptbin]->SetMarkerSize(0.5);
    projinvBeta[ptbin]->SetTitle("");
    projinvBeta[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

  }
 
  // make graphs that need all pt bins
  TGraphErrors* purGr = new TGraphErrors(numPtBins,pT,purity,dx,dy);
  purC->cd();
  purGr->SetMarkerStyle(20);
  purGr->SetMarkerSize(0.7);
  purGr->SetTitle("Electron Sample Purity");
  purGr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  purGr->GetYaxis()->SetTitle("Purity");
  purGr->GetYaxis()->SetRangeUser(0.5,1.2);
  purGr->SetMarkerColor(kRed);
  purGr->SetLineColor(kRed);
  purGr->Draw("AP");
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
    sprintf(tlName,"                            isHFTTrack();");
    tl.DrawLatex(0.1, 0.6,tlName);

    sprintf(tlName,"Event:            #left|V_{z}#right| < 6 cm; #left|#DeltaV_{z}#right| < 4 cm");
    tl.DrawLatex(0.1, 0.45,tlName);
    sprintf(tlName,"Triggers:            MB");
    tl.DrawLatex(0.1, 0.4,tlName);


    // Place canvases in order
    TCanvas* temp = new TCanvas();
   
    sprintf(name, "%s_%s.pdf[", FileName,FileLabel);
    temp->Print(name);
    sprintf(name, "%s_%s.pdf", FileName,FileLabel);
    temp = fp; // print front page
    temp->Print(name);
    temp = nSigE; // print data canvases
    temp->Print(name);
    temp = purC;
    temp->Print(name);
    temp = nSigPi;
    temp->Print(name);
    temp = nSigK;
    temp->Print(name);
    temp = nSigP;
    temp->Print(name);
    temp = nSigEK;
    temp->Print(name);
    temp = nSigEP;
    temp->Print(name);
    temp = nSigEPi;
    temp->Print(name);
    temp = invBeta;
    temp->Print(name);
    sprintf(name, "%s_%s.pdf]", FileName,FileLabel);
    temp->Print(name);
  }

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
