#include "anaConst14.h"

int pl_Purity(){
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

  // create generic structures for use in all plots
  TPaveText* lbl[numPtBins];
  char textLabel[100]; 


  // Load all trigger type files
  const char* baseName = "outputs/fullSample_Feb7";
  char BEMCName[4][100],SMDLName[4][100],SMDTName[4][100],BEMCHFTName[4][100],SMDLHFTName[4][100],SMDTHFTName[4][100];
  char trigName[4][100] = {"MB","BHT1","BHT2","BHT3"};
  TFile* BEMC[4];
  TFile* SMDL[4];
  TFile* SMDT[4];
  for(int q=0; q<4; q++)
  {
    sprintf(BEMCName[q],"%s_%s_BEMC_processed.root",baseName,trigName[q]);
    BEMC[q] = new TFile(BEMCName[q],"READ");
    sprintf(SMDLName[q],"%s_%s_SMD_processed.root",baseName,trigName[q]);
    SMDL[q] = new TFile(SMDLName[q],"READ");
    sprintf(SMDTName[q],"%s_%s_SMD2_processed.root",baseName,trigName[q]);
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
  TH1F* nSigBEMC[4][numPtBins];
  TH1F* nSigSMDL[4][numPtBins];
  TH1F* nSigSMDT[4][numPtBins];
  TH1F* PurityBEMC[4];
  TH1F* PuritySMDL[4];
  TH1F* PuritySMDT[4];
  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      nSigBEMC[w][ptbin] = (TH1F*)BEMC[w]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDL[w][ptbin] = (TH1F*)SMDL[w]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDT[w][ptbin] = (TH1F*)SMDT[w]->Get(Form("drawnSigmaE_%i",ptbin));
    }
    PurityBEMC[w] = (TH1F*)BEMC[w]->Get("drawPurity");
    PuritySMDL[w] = (TH1F*)SMDL[w]->Get("drawPurity");
    PuritySMDT[w] = (TH1F*)SMDT[w]->Get("drawPurity");
  }
  if(DEBUG)
    cout << "Found Hists." << endl;

  // Make Canvas
  TCanvas* purityOL = new TCanvas("purityOL","purityOL",50,50,1050,1050);
  TCanvas* nSigmaOL[numPtBins];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    nSigmaOL[ptbin]= new TCanvas(Form("nSigmaOL_%i",ptbin),Form("%.2f < pT < %.2f",lowpt[ptbin],highpt[ptbin]),50,50,1050,1050);
    nSigmaOL[ptbin]->Divide(2,2);

    sprintf(textLabel,"%.2f < P_{T} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);
  }
  purityOL->Divide(2,2);

  if(DEBUG)
    cout << "Canvas Made." << endl;

  // Plot Settings
  for(int r=0; r<4; r++)
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      nSigBEMC[r][ptbin]->SetLineColor(kRed);
      nSigBEMC[r][ptbin]->SetMarkerColor(kRed);
      nSigBEMC[r][ptbin]->SetMarkerStyle(20);
      nSigBEMC[r][ptbin]->SetTitle(Form("%s nSigmaE, %.2f < pT < %.2f",trigName[r],lowpt[ptbin],highpt[ptbin]));
      nSigBEMC[r][ptbin]->GetXaxis()->SetTitle("nSigmaE");
      nSigBEMC[r][ptbin]->GetYaxis()->SetTitle("Counts");
      if(DEBUG) cout << "after BEMC" << endl;

      nSigSMDL[r][ptbin]->SetLineColor(kBlue);
      nSigSMDL[r][ptbin]->SetMarkerColor(kBlue);
      nSigSMDL[r][ptbin]->SetMarkerStyle(21);
      if(DEBUG) cout << "after SMDL" << endl;

      nSigSMDT[r][ptbin]->SetLineColor(kBlack);
      nSigSMDT[r][ptbin]->SetMarkerColor(kBlack);
      nSigSMDT[r][ptbin]->SetMarkerStyle(22);
    }
    if(DEBUG) cout << "after ptbin loop" << endl;

    PurityBEMC[r]->SetLineColor(kRed);
    PurityBEMC[r]->SetMarkerColor(kRed);
    PurityBEMC[r]->SetMarkerStyle(20);
    PurityBEMC[r]->SetTitle(Form("%s Purity",trigName[r]));

    PuritySMDL[r]->SetLineColor(kBlue);
    PuritySMDL[r]->SetMarkerColor(kBlue);
    PuritySMDL[r]->SetMarkerStyle(21);

    PuritySMDT[r]->SetLineColor(kBlack);
    PuritySMDT[r]->SetMarkerColor(kBlack);
    PuritySMDT[r]->SetMarkerStyle(22);
  }

  if(DEBUG)
    cout << "Settings Assigned." << endl;

  // Actually Draw

  // Use input for only one ptbin in legend. They're all the same, no need to loop.
  TLegend* leg = new TLegend(0.55,0.55,0.85,0.85);
  leg->AddEntry(nSigBEMC[0][0],"BEMC","lpe");
  leg->AddEntry(nSigSMDL[0][0],"BEMC+SMD(Loose)","lpe");
  leg->AddEntry(nSigSMDT[0][0],"BEMC+SMD(Tight)","lpe");
  if(DEBUG)
    cout << "Legend Made." << endl;

  for(int t=0; t<4; t++)
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      nSigmaOL[ptbin]->cd(t+1);
      gPad->SetLogy(isLogY);  
      nSigBEMC[t][ptbin]->Draw();
      nSigSMDL[t][ptbin]->Draw("same");
      nSigSMDT[t][ptbin]->Draw("same");
      leg->Draw("same");
      lbl[ptbin]->Draw("same");
    }

    purityOL->cd(t+1);
    PurityBEMC[t]->Draw("Ape");
    PuritySMDL[t]->Draw("same pe");
    PuritySMDT[t]->Draw("same pe");
    leg->Draw("same");
  }

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
    
    TCanvas* temp = new TCanvas();
    char name[100];
    sprintf(name, "%s_purityOverlays.pdf[", baseName);
    temp->Print(name);
    sprintf(name, "%s_purityOverlays.pdf", baseName);
    temp = fp; // print front page
    temp->Print(name);

    for(int q=0; q<numPtBins; q++)
    {
      temp = nSigmaOL[q]; // print data canvases
      temp->Print(name);
    }
    temp = purityOL;
    temp->Print(name);

    sprintf(name, "%s_purityOverlays.pdf]", baseName);
    temp->Print(name);


  return 1;
}
