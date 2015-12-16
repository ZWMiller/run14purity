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
  sprintf(name,"/Users/zach/Research/rootFiles/run14NPEpurity/outputs/%s.root",FileName);
  TFile *f = new TFile(name,"READ");
  if (f->IsOpen()==kFALSE)
  { std::cout << "!!! File Not Found !!!" << std::endl;
    exit(1); }
  // f->ls(); // - DEBUG by printing all objects in ROOT file

  char fname[100];
  TFile* file;
  if(makeROOT){
    sprintf(fname,"/Users/zach/Research/rootFiles/run14NPEpurity/outputs/%s_processed.root",FileName);
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
    lbl[ptbin] = new TPaveText(.61,.75,.81,.83,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.1f < P_{T,e} < %.1f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);

    nSigPi->cd(ptbin+1);
    projnSigmaPi[ptbin]->SetMarkerColor(kBlack);
    projnSigmaPi[ptbin]->SetMarkerStyle(20);
    projnSigmaPi[ptbin]->SetMarkerSize(0.5);
    projnSigmaPi[ptbin]->SetLineColor(kBlack);
    projnSigmaPi[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaPi[ptbin]->SetTitle("");
    projnSigmaPi[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigK->cd(ptbin+1);
    projnSigmaK[ptbin]->SetMarkerColor(kBlack);
    projnSigmaK[ptbin]->SetMarkerStyle(20);
    projnSigmaK[ptbin]->SetMarkerSize(0.5);
    projnSigmaK[ptbin]->SetLineColor(kBlack);
    projnSigmaK[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaK[ptbin]->SetTitle("");
    projnSigmaK[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigP->cd(ptbin+1);
    projnSigmaP[ptbin]->SetMarkerColor(kBlack);
    projnSigmaP[ptbin]->SetMarkerStyle(20);
    projnSigmaP[ptbin]->SetMarkerSize(0.5);
    projnSigmaP[ptbin]->SetLineColor(kBlack);
    projnSigmaP[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaP[ptbin]->SetTitle("");
    projnSigmaP[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    nSigE->cd(ptbin+1);
    projnSigmaE[ptbin]->SetMarkerColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerStyle(20);
    projnSigmaE[ptbin]->SetLineColor(kBlack);
    projnSigmaE[ptbin]->SetMarkerSize(0.5);
    projnSigmaE[ptbin]->GetXaxis()->SetRangeUser(-10.,10.);
    projnSigmaE[ptbin]->SetTitle("");
    projnSigmaE[ptbin]->Draw();
    lbl[ptbin]->Draw("same");
  
    invBeta->cd(ptbin+1);
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
    sprintf(tlName,"-1.5 < n #sigma E < 3; n #sigma P < #left|20#right|; ");
    tl.DrawLatex(0.1, 0.8,tlName);
    sprintf(tlName,"n #sigma K < #left|20#right|; n #sigma Pi < #left|20#right|;");
    tl.DrawLatex(0.1, 0.75,tlName);
    sprintf(tlName,"");
    tl.DrawLatex(0.1, 0.7,tlName);
    sprintf(tlName,"");
    tl.DrawLatex(0.1, 0.6,tlName);
    sprintf(tlName,"Event:  #left|V_{z}#right| < 100 cm;");
    tl.DrawLatex(0.1, 0.5,tlName);
    sprintf(tlName,"Triggers:  HT");
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
