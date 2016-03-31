#include "anaConst14.h"

int pl_HFTComparison(){
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
  const char* baseName = "outputs/Mar1_PuritySample";
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
    
    // with HFT
    sprintf(BEMCName[q],"%s_%s_BEMC_wHFT_processed.root",baseName,trigName[q]);
    BEMCh[q] = new TFile(BEMCName[q],"READ");
    sprintf(SMDLName[q],"%s_%s_SMD_wHFT_processed.root",baseName,trigName[q]);
    SMDLh[q] = new TFile(SMDLName[q],"READ");
    sprintf(SMDTName[q],"%s_%s_SMD2_wHFT_processed.root",baseName,trigName[q]);
    SMDTh[q] = new TFile(SMDTName[q],"READ");

  }

  for(int i=0; i<4; i++)
  {
    if(BEMC[i]->IsOpen() && SMDL[i]->IsOpen() && SMDT[i]->IsOpen() &&
       BEMCh[i]->IsOpen() && SMDLh[i]->IsOpen() && SMDTh[i]->IsOpen())
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
  TH1F* dNdpTBEMC[4];
  TH1F* dNdpTSMDL[4];
  TH1F* dNdpTSMDT[4];
  TH1F* nSigBEMCh[4][numPtBins];
  TH1F* nSigSMDLh[4][numPtBins];
  TH1F* nSigSMDTh[4][numPtBins];
  TH1F* PurityBEMCh[4];
  TH1F* PuritySMDLh[4];
  TH1F* PuritySMDTh[4];
  TH1F* dNdpTBEMCh[4];
  TH1F* dNdpTSMDLh[4];
  TH1F* dNdpTSMDTh[4];

  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      nSigBEMC[w][ptbin] = (TH1F*)BEMC[w]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDL[w][ptbin] = (TH1F*)SMDL[w]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDT[w][ptbin] = (TH1F*)SMDT[w]->Get(Form("drawnSigmaE_%i",ptbin)); 
      nSigBEMCh[w][ptbin] = (TH1F*)BEMCh[4]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDLh[w][ptbin] = (TH1F*)SMDLh[4]->Get(Form("drawnSigmaE_%i",ptbin));
      nSigSMDTh[w][ptbin] = (TH1F*)SMDTh[4]->Get(Form("drawnSigmaE_%i",ptbin));

    }
    PurityBEMC[w] = (TH1F*)BEMC[w]->Get("drawPurity");
    PuritySMDL[w] = (TH1F*)SMDL[w]->Get("drawPurity");
    PuritySMDT[w] = (TH1F*)SMDT[w]->Get("drawPurity");
    dNdpTBEMC[w] = (TH1F*)BEMC[w]->Get("drawdNdpT");
    dNdpTSMDL[w] = (TH1F*)SMDL[w]->Get("drawdNdpT");
    dNdpTSMDT[w] = (TH1F*)SMDT[w]->Get("drawdNdpT");
    PurityBEMCh[w] = (TH1F*)BEMCh[w]->Get("drawPurity");
    PuritySMDLh[w] = (TH1F*)SMDLh[w]->Get("drawPurity");
    PuritySMDTh[w] = (TH1F*)SMDTh[w]->Get("drawPurity");
  }
  if(DEBUG)
    cout << "Found Hists." << endl;

  // Make Canvas
  TCanvas* purityOL[4];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    sprintf(textLabel,"%.2f < P_{T} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);
  }

  if(DEBUG)
    cout << "Canvas Made." << endl;

  // Plot Settings
  for(int r=0; r<4; r++)
  {
    purityOL[r] = new TCanvas(Form("purityOL_%i",r),"Purity Overlays",50,50,1050,1050);
    purityOL[r]->Divide(2,2);

    PurityBEMC[r]->SetLineColor(kRed);
    PurityBEMC[r]->SetMarkerColor(kRed);
    PurityBEMC[r]->SetMarkerStyle(20);
    PurityBEMC[r]->SetTitle(Form("%s BEMC Purity",trigName[r]));

    PuritySMDL[r]->SetLineColor(kRed);
    PuritySMDL[r]->SetMarkerColor(kRed);
    PuritySMDL[r]->SetMarkerStyle(21);
    PuritySMDL[r]->SetTitle(Form("%s SMD Purity",trigName[r]));

    PuritySMDT[r]->SetLineColor(kRed);
    PuritySMDT[r]->SetMarkerColor(kRed);
    PuritySMDT[r]->SetMarkerStyle(22);
    PuritySMDT[r]->SetTitle(Form("%s SMD2 Purity",trigName[r]));

    PurityBEMCh[r]->SetLineColor(kBlack);
    PurityBEMCh[r]->SetMarkerColor(kBlack);
    PurityBEMCh[r]->SetMarkerStyle(20);
    PurityBEMCh[r]->SetTitle(Form("%s BEMC Purity",trigName[r]));

    PuritySMDLh[r]->SetLineColor(kBlack);
    PuritySMDLh[r]->SetMarkerColor(kBlack);
    PuritySMDLh[r]->SetMarkerStyle(21);
    PuritySMDLh[r]->SetTitle(Form("%s SMD Purity",trigName[r]));

    PuritySMDTh[r]->SetLineColor(kBlack);
    PuritySMDTh[r]->SetMarkerColor(kBlack);
    PuritySMDTh[r]->SetMarkerStyle(22);
    PuritySMDTh[r]->SetTitle(Form("%s SMD2 Purity",trigName[r]));

  }

  if(DEBUG)
    cout << "Settings Assigned." << endl;

  // Actually Draw

  // Use input for only one ptbin in legend. They're all the same, no need to loop.
  TLegend* leg = new TLegend(0.55,0.55,0.85,0.85);
  leg->AddEntry(PurityBEMC[0],"W/O HFT","lpe");
  leg->AddEntry(PurityBEMCh[0],"W/ HFT","lpe");
  if(DEBUG)
    cout << "Legend Made." << endl;

  for(int t=0; t<4; t++)
  {
    purityOL[t]->cd(1);
    PurityBEMC[t]->Draw("Ape");
    PurityBEMCh[t]->Draw("same pe");
    leg->Draw("same");
    purityOL[t]->cd(2);
    PuritySMDL[t]->Draw("Ape");
    PuritySMDLh[t]->Draw("same pe");
    leg->Draw("same");
    purityOL[t]->cd(3);
    PuritySMDT[t]->Draw("Ape");
    PuritySMDTh[t]->Draw("same pe");
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
     sprintf(tlName,"Study of Various Cut Sets on Inclusive Electron Purity.");
    tl.DrawLatex(0.1, 0.75,tlName);
    sprintf(tlName,"Z.W. Miller - UIC");
    tl.DrawLatex(0.1, 0.67,tlName);

    
    TCanvas* temp = new TCanvas();
    char name[100];
    sprintf(name, "%s_purityOverlays_CompareHFT.pdf[", baseName);
    temp->Print(name);
    sprintf(name, "%s_purityOverlays_CompareHFT.pdf", baseName);
    temp = fp; // print front page
    temp->Print(name);

    for(int q=0; q<4; q++)
    {
      temp = purityOL[q]; // print data canvases
      temp->Print(name);
    }

    sprintf(name, "%s_purityOverlays_CompareHFT.pdf]", baseName);
    temp->Print(name);


  return 1;
}
