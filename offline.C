// Offline Plots - Z. Miller July 24, 2015
//
// .L offline.C
// offline("FILENAME") # Without .root Extension

#include "anaConst14.h"

// Declare functions
void checkBatchMode();
Bool_t checkMakePDF();
Bool_t checkMakeRoot();

void offline(const char* FileName="test")
{
  TH1F::SetDefaultSumw2();
  // Set Style parameters for this macro
  gStyle->SetOptTitle(1); // Show Title (off by default for cleanliness)
  gErrorIgnoreLevel = kError; // Set Verbosity Level (kPrint shows all)

  // Set Output options
  Int_t number;
  checkBatchMode();
  Bool_t makePDF = checkMakePDF();
  Bool_t makeROOT= checkMakeRoot();

  // Open ROOT File
  char name[1000];
  sprintf(name,"/Users/zach/Research/rootFiles/run14NPEpurity/%s.root",FileName);
  TFile *f = new TFile(name,"READ");
  if (f->IsOpen()==kFALSE)
  { std::cout << "!!! File Not Found !!!" << std::endl;
    exit(1); }
  // f->ls(); // - DEBUG by printing all objects in ROOT file

  char fname[100];
  TFile* file;
  if(makeROOT){
    sprintf(fname,"/Users/zach/Research/rootFiles/run14NPEpurity/%s_processed.root",FileName);
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
  TCanvas* invBeta = new TCanvas("invBeta","Beta^(-1) Projections",50,50,1050,1050);
  nSigPi  -> Divide(3,3);
  nSigK   -> Divide(3,3);
  nSigP   -> Divide(3,3);
  nSigE   -> Divide(3,3);
  invBeta -> Divide(3,3);

  // Make Projections (first get 2d/3d hists, then project)

  TH2F* nSigmaPiPt = (TH2F*)f->Get("nsigmaPI_Pt");
  TH2F* nSigmaKPt  = (TH2F*)f->Get("nsigmaK_Pt");
  TH2F* nSigmaPPt  = (TH2F*)f->Get("nsigmaP_Pt");
  TH2F* nSigmaEPt  = (TH2F*)f->Get("nsigmaE_Pt");
  TH2F* invBetaPt  = (TH2F*)f->Get("invsBeta_Pt");


  TH1D* projnSigmaPi[numPtBins];
  TH1D* projnSigmaK[numPtBins];
  TH1D* projnSigmaP[numPtBins];
  TH1D* projnSigmaE[numPtBins];
  TH1D* projinvBeta[numPtBins];
  TF1 *fitPi[numPtBins];
  TF1 *fitKP[numPtBins];
  TF1 *fitE[numPtBins];
  TF1 *fitCom[numPtBins];
  TF1 *fitPiD[numPtBins];
  TF1 *fitKPD[numPtBins];
  TF1 *fitED[numPtBins];
  TF1 *fitComD[numPtBins];

  for(Int_t ptbin=0; ptbin<numPtBins; ptbin++)
  {
    projnSigmaPi[ptbin] = nSigmaPiPt->ProjectionY(Form("projnSigmaPi_%i",ptbin),nSigmaPiPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaPiPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaK[ptbin]  = nSigmaKPt->ProjectionY(Form("projnSigmaK_%i",ptbin),nSigmaKPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaKPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaP[ptbin]  = nSigmaPPt->ProjectionY(Form("projnSigmaP_%i",ptbin),nSigmaPPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaPPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projnSigmaE[ptbin]  = nSigmaEPt->ProjectionY(Form("projnSigmaE_%i",ptbin),nSigmaEPt->GetXaxis()->FindBin(lowpt[ptbin]),nSigmaEPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    projinvBeta[ptbin]  = invBetaPt->ProjectionY(Form("invBeta_%i",ptbin),invBetaPt->GetXaxis()->FindBin(lowpt[ptbin]),invBetaPt->GetXaxis()->FindBin(highpt[ptbin])-1);
  }

  // Analyze the projections
  for(Int_t ptbin = 0; ptbin < numPtBins; ptbin++){

    // Init necessary plotting tools
    lbl[ptbin] = new TPaveText(.65,.37,.81,.44,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);

    nSigPi->cd(ptbin+1);
    gPad->SetLogy();
    projnSigmaPi[ptbin]->SetMarkerColor(kBlack);
    projnSigmaPi[ptbin]->SetMarkerStyle(20);
    projnSigmaPi[ptbin]->SetMarkerSize(0.5);
    projnSigmaPi[ptbin]->SetLineColor(kBlack);
    projnSigmaPi[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaPi[ptbin]->SetTitle("");
    projnSigmaPi[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigK->cd(ptbin+1);
    gPad->SetLogy();
    projnSigmaK[ptbin]->SetMarkerColor(kBlack);
    projnSigmaK[ptbin]->SetMarkerStyle(20);
    projnSigmaK[ptbin]->SetMarkerSize(0.5);
    projnSigmaK[ptbin]->SetLineColor(kBlack);
    projnSigmaK[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaK[ptbin]->SetTitle("");
    projnSigmaK[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigP->cd(ptbin+1);
    gPad->SetLogy();
    projnSigmaP[ptbin]->SetMarkerColor(kBlack);
    projnSigmaP[ptbin]->SetMarkerStyle(20);
    projnSigmaP[ptbin]->SetMarkerSize(0.5);
    projnSigmaP[ptbin]->SetLineColor(kBlack);
    projnSigmaP[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaP[ptbin]->SetTitle("");
    projnSigmaP[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigE->cd(ptbin+1);
    gPad->SetLogy();
    projnSigmaE[ptbin]->SetMarkerColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerStyle(20);
    projnSigmaE[ptbin]->SetLineColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerSize(0.5);
    projnSigmaE[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaE[ptbin]->SetTitle("");
    projnSigmaE[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    // Do fits on nSig E
    gStyle->SetOptFit(1111);
    Double_t par[9];
    Double_t parS[3];
    fitKP[ptbin] = new TF1(Form("fitKP_%i",ptbin),"gaus",-10,-5);
    fitPi[ptbin] = new TF1(Form("fitPi_%i",ptbin),"gaus",-5,-1);
    fitE[ptbin]  = new TF1(Form("fitE_%i",ptbin),"gaus",-1,4);
    fitCom[ptbin]= new TF1(Form("fitCom_%i",ptbin),"gaus(0)+gaus(3)+gaus(6)",-10,10);
    fitCom[ptbin]->SetParName(0,"#pi C");
    fitCom[ptbin]->SetParName(1,"#pi #mu");
    fitCom[ptbin]->SetParName(2,"#pi #sigma");
    fitCom[ptbin]->SetParName(3,"Kp C");
    fitCom[ptbin]->SetParName(4,"Kp #mu");
    fitCom[ptbin]->SetParName(5,"Kp #sigma");
    fitCom[ptbin]->SetParName(6,"e C");
    fitCom[ptbin]->SetParName(7,"e #mu");
    fitCom[ptbin]->SetParName(8,"e #sigma");
    fitPi[ptbin]->SetLineColor(kRed);
    fitKP[ptbin]->SetLineColor(kCyan);
    fitE[ptbin]->SetLineColor(kBlue);
    fitCom[ptbin]->SetLineColor(kMagenta);
    fitPi[ptbin]->SetLineWidth(1);
    fitKP[ptbin]->SetLineWidth(1);
    fitE[ptbin]->SetLineWidth(1);
    fitCom[ptbin]->SetLineWidth(2);
    projnSigmaE[ptbin]->Fit(fitPi[ptbin],"R");
    projnSigmaE[ptbin]->Fit(fitKP[ptbin],"R");
    projnSigmaE[ptbin]->Fit(fitE[ptbin],"R");
    fitPi[ptbin]->GetParameters(&par[0]);
    fitKP[ptbin]->GetParameters(&par[3]);
    fitE[ptbin]->GetParameters(&par[6]);
    fitCom[ptbin]->SetParameters(par);
    projnSigmaE[ptbin]->Fit(fitCom[ptbin],"R");
    projnSigmaE[ptbin]->UseCurrentStyle();

    // make versions to draw fits in full
    fitPiD[ptbin] = new TF1(Form("drawPi_%i",ptbin),"gaus",-10,10);
    fitKPD[ptbin] = new TF1(Form("drawKP_%i",ptbin),"gaus",-10,10);
    fitED[ptbin] = new TF1(Form("drawE_%i",ptbin),"gaus",-10,10);
    fitPi[ptbin]->GetParameters(&parS[0]);
    fitPiD[ptbin]->SetParameters(parS);
    fitKP[ptbin]->GetParameters(&parS[0]);
    fitKPD[ptbin]->SetParameters(parS);
    fitE[ptbin]->GetParameters(&parS[0]);
    fitED[ptbin]->SetParameters(parS);
    fitPiD[ptbin]->SetLineStyle(2);
    fitPiD[ptbin]->SetLineWidth(1);
    fitPiD[ptbin]->SetLineColor(kRed);
    fitKPD[ptbin]->SetLineStyle(2);
    fitKPD[ptbin]->SetLineWidth(1);
    fitKPD[ptbin]->SetLineColor(kCyan);
    fitED[ptbin]->SetLineStyle(2);
    fitED[ptbin]->SetLineWidth(1);
    fitED[ptbin]->SetLineColor(kBlue);
    fitPiD[ptbin]->Draw("same");
    fitKPD[ptbin]->Draw("same");
    fitED[ptbin]->Draw("same");

    // Make Stats box legible
    TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
    s->SetX1NDC(0.65);
    s->SetX2NDC(0.95);
    s->SetY1NDC(0.45);
    s->SetY2NDC(0.95);
    nSigE->Modified();
    nSigE->Update();


    invBeta->cd(ptbin+1);
    gPad->SetLogy();
    projinvBeta[ptbin]->SetMarkerColor(kBlack);
    projinvBeta[ptbin]->SetMarkerStyle(20);
    projinvBeta[ptbin]->SetLineColor(kBlack);
    projinvBeta[ptbin]->SetMarkerSize(0.5);
    projinvBeta[ptbin]->SetTitle("");
    projinvBeta[ptbin]->Draw();
    lbl[ptbin]->Draw("same");
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
     sprintf(tlName,"                            isHFTTrack();");
    tl.DrawLatex(0.1, 0.6,tlName);

    sprintf(tlName,"Event:            #left|V_{z}#right| < 6 cm; #left|#DeltaV_{z}#right| < 4 cm");
    tl.DrawLatex(0.1, 0.45,tlName);
    sprintf(tlName,"Triggers:            MB");
    tl.DrawLatex(0.1, 0.4,tlName);


    // Place canvases in order
    TCanvas* temp = new TCanvas();
    sprintf(name, "%s.pdf[", FileName);
    temp->Print(name);
    sprintf(name, "%s.pdf", FileName);
    temp = fp; // print front page
    temp->Print(name);
    temp = nSigE; // print data canvases
    temp->Print(name);
    temp = nSigPi;
    temp->Print(name);
    temp = nSigK;
    temp->Print(name);
    temp = nSigP;
    temp->Print(name);
    temp = invBeta;
    temp->Print(name);
    sprintf(name, "%s.pdf]", FileName);
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
