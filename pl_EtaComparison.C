#include "anaConst14.h"

int pl_EtaComparison(){
  // Set Flags for Code/Plots
  Bool_t DEBUG = kTRUE;
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
  char textLabel[100]; 


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
  TH1F* nSigBEMC[4][numEtaBins][numPtBins];
  TH1F* nSigSMDL[4][numEtaBins][numPtBins];
  TH1F* nSigSMDT[4][numEtaBins][numPtBins];
  TH1F* PurityBEMC[4][numEtaBins];
  TH1F* PuritySMDL[4][numEtaBins];
  TH1F* PuritySMDT[4][numEtaBins];
  TH1F* dNdpTBEMC[4][numEtaBins];
  TH1F* dNdpTSMDL[4][numEtaBins];
  TH1F* dNdpTSMDT[4][numEtaBins];
  TH1F* Purity2BEMC[4][numEtaBins];
  TH1F* Purity2SMDL[4][numEtaBins];
  TH1F* Purity2SMDT[4][numEtaBins];
  TH1F* dNdpT2BEMC[4][numEtaBins];
  TH1F* dNdpT2SMDL[4][numEtaBins];
  TH1F* dNdpT2SMDT[4][numEtaBins];
  TH1F* adjYieldBEMC[4][numEtaBins];
  TH1F* adjYieldSMDL[4][numEtaBins];
  TH1F* adjYieldSMDT[4][numEtaBins];
  TH1F* nSigBEMCh[4][numEtaBins][numPtBins];
  TH1F* nSigSMDLh[4][numEtaBins][numPtBins];
  TH1F* nSigSMDTh[4][numEtaBins][numPtBins];
  TH1F* PurityBEMCh[4][numEtaBins];
  TH1F* PuritySMDLh[4][numEtaBins];
  TH1F* PuritySMDTh[4][numEtaBins];
  TH1F* dNdpTBEMCh[4][numEtaBins];
  TH1F* dNdpTSMDLh[4][numEtaBins];
  TH1F* dNdpTSMDTh[4][numEtaBins];

  for(int w=0; w<4; w++)
  {
    if(DEBUG) cout << "w: " << w << endl;
    for(int etabin=0; etabin<numEtaBins; etabin++)
    {
      for(int ptbin=0; ptbin<numPtBins; ptbin++)
      {
        nSigBEMC[w][etabin][ptbin] = (TH1F*)BEMC[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDL[w][etabin][ptbin] = (TH1F*)SMDL[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDT[w][etabin][ptbin] = (TH1F*)SMDT[w]->Get(Form("drawnSigmaE_%i%i",etabin,ptbin)); 
       /* nSigBEMCh[w][etabin][ptbin] = (TH1F*)BEMCh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDLh[w][etabin][ptbin] = (TH1F*)SMDLh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));
        nSigSMDTh[w][etabin][ptbin] = (TH1F*)SMDTh[w]->Get(Form("drawnSigmaE_%i_%i",etabin,ptbin));*/
      }

      PurityBEMC[w][etabin] = (TH1F*)BEMC[w]->Get(Form("drawPurity_%i",etabin));
      PuritySMDL[w][etabin] = (TH1F*)SMDL[w]->Get(Form("drawPurity_%i",etabin));
      PuritySMDT[w][etabin] = (TH1F*)SMDT[w]->Get(Form("drawPurity_%i",etabin));
      dNdpTBEMC[w][etabin] = (TH1F*)BEMC[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpTSMDL[w][etabin] = (TH1F*)SMDL[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpTSMDT[w][etabin] = (TH1F*)SMDT[w]->Get(Form("drawdNdpT_%i",etabin));
      Purity2BEMC[w][etabin] = (TH1F*)BEMC[w]->Get(Form("drawPurity_%i",etabin));
      Purity2SMDL[w][etabin] = (TH1F*)SMDL[w]->Get(Form("drawPurity_%i",etabin));
      Purity2SMDT[w][etabin] = (TH1F*)SMDT[w]->Get(Form("drawPurity_%i",etabin));
      dNdpT2BEMC[w][etabin] = (TH1F*)BEMC[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpT2SMDL[w][etabin] = (TH1F*)SMDL[w]->Get(Form("drawdNdpT_%i",etabin));
      dNdpT2SMDT[w][etabin] = (TH1F*)SMDT[w]->Get(Form("drawdNdpT_%i",etabin));
      adjYieldBEMC[w][etabin] = (TH1F*)BEMC[w]->Get(Form("adjYield_%i",etabin));
      adjYieldSMDL[w][etabin] = (TH1F*)SMDL[w]->Get(Form("adjYield_%i",etabin));
      adjYieldSMDT[w][etabin] = (TH1F*)SMDT[w]->Get(Form("adjYield_%i",etabin));

     /* PurityBEMCh[w][etabin] = (TH1F*)BEMCh[w]->Get(Form("drawPurity_%i",etabin));
      PuritySMDLh[w][etabin] = (TH1F*)SMDLh[w]->Get(Form("drawPurity_%i",etabin));
      PuritySMDTh[w][etabin] = (TH1F*)SMDTh[w]->Get(Form("drawPurity_%i",etabin));*/
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

      PurityBEMC[r][etabin]->SetLineColor(1+etabin);
      PurityBEMC[r][etabin]->SetMarkerColor(1+etabin);
      PurityBEMC[r][etabin]->SetMarkerStyle(4+etabin);
      PurityBEMC[r][etabin]->SetTitle(Form("%s BEMC Purity",trigName[r]));

      PuritySMDL[r][etabin]->SetLineColor(1+etabin);
      PuritySMDL[r][etabin]->SetMarkerColor(1+etabin);
      PuritySMDL[r][etabin]->SetMarkerStyle(4+etabin);
      PuritySMDL[r][etabin]->SetTitle(Form("%s SMD Purity",trigName[r]));

      PuritySMDT[r][etabin]->SetLineColor(1+etabin);
      PuritySMDT[r][etabin]->SetMarkerColor(1+etabin);
      PuritySMDT[r][etabin]->SetMarkerStyle(4+etabin);
      PuritySMDT[r][etabin]->SetTitle(Form("%s SMD2 Purity",trigName[r]));

      Purity2BEMC[r][etabin]->SetLineColor(kBlack);
      Purity2SMDL[r][etabin]->SetLineColor(kRed);
      Purity2SMDT[r][etabin]->SetLineColor(kBlue);
      Purity2BEMC[r][etabin]->SetMarkerColor(kBlack);
      Purity2SMDL[r][etabin]->SetMarkerColor(kRed);
      Purity2SMDT[r][etabin]->SetMarkerColor(kBlue);
      Purity2BEMC[r][etabin]->SetMarkerStyle(20);
      Purity2SMDL[r][etabin]->SetMarkerStyle(21);
      Purity2SMDT[r][etabin]->SetMarkerStyle(22);
      Purity2BEMC[r][etabin]->SetTitle(Form("%s Purity Trigger Compare",trigName[r]));
      
      dNdpT2BEMC[r][etabin]->SetLineColor(kBlack);
      dNdpT2SMDL[r][etabin]->SetLineColor(kRed);
      dNdpT2SMDT[r][etabin]->SetLineColor(kBlue);
      dNdpT2BEMC[r][etabin]->SetMarkerColor(kBlack);
      dNdpT2SMDL[r][etabin]->SetMarkerColor(kRed);
      dNdpT2SMDT[r][etabin]->SetMarkerColor(kBlue);
      dNdpT2BEMC[r][etabin]->SetMarkerStyle(20);
      dNdpT2SMDL[r][etabin]->SetMarkerStyle(21);
      dNdpT2SMDT[r][etabin]->SetMarkerStyle(22);
      dNdpT2BEMC[r][etabin]->SetTitle(Form("%s dN/dpT Trigger Compare",trigName[r]));

      adjYieldBEMC[r][etabin]->SetLineColor(kBlack);
      adjYieldSMDL[r][etabin]->SetLineColor(kRed);
      adjYieldSMDT[r][etabin]->SetLineColor(kBlue);
      adjYieldBEMC[r][etabin]->SetMarkerColor(kBlack);
      adjYieldSMDL[r][etabin]->SetMarkerColor(kRed);
      adjYieldSMDT[r][etabin]->SetMarkerColor(kBlue);
      adjYieldBEMC[r][etabin]->SetMarkerStyle(20);
      adjYieldSMDL[r][etabin]->SetMarkerStyle(21);
      adjYieldSMDT[r][etabin]->SetMarkerStyle(22);
      adjYieldBEMC[r][etabin]->SetTitle(Form("%s Adjusted dN/dpT Trigger Compare",trigName[r]));
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
    leg->AddEntry(PurityBEMC[0][etabin],textLabel,"lpe");
  }
  
  TLegend* leg2 = new TLegend(0.15,0.15,0.45,0.45);
  leg2->AddEntry(Purity2BEMC[0][0],"BEMC","lpe");
  leg2->AddEntry(Purity2SMDL[0][0],"SMD","lpe");
  leg2->AddEntry(Purity2SMDT[0][0],"SMD2","lpe");
  TLegend* leg3 = new TLegend(0.55,0.55,0.85,0.85);
  leg3->AddEntry(Purity2BEMC[0][0],"BEMC","lpe");
  leg3->AddEntry(Purity2SMDL[0][0],"SMD","lpe");
  leg3->AddEntry(Purity2SMDT[0][0],"SMD2","lpe");
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
      Purity2SMDL[t][etabin]->Draw("same pe");
      Purity2SMDT[t][etabin]->Draw("same pe");
      leg2->Draw("same");
      lblE[etabin]->Draw("same");
      purityTrigOL[t][etabin]->cd(2);
      gPad->SetLogy(1);
      dNdpT2BEMC[t][etabin]->Draw("Ape");
      dNdpT2SMDL[t][etabin]->Draw("same pe");
      dNdpT2SMDT[t][etabin]->Draw("same pe");
      leg3->Draw("same");
      lblE[etabin]->Draw("same");

      purityTrigOL2[t][etabin]->cd();
      gPad->SetLogy(1);
      adjYieldBEMC[t][etabin]->Draw("Ape");
      adjYieldSMDL[t][etabin]->Draw("same pe");
      adjYieldSMDT[t][etabin]->Draw("same pe");
      leg3->Draw("same");
      lblE[etabin]->Draw("same");
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

  for(int q=1; q<4; q++)
  {
    if(DEBUG) cout << "Place Canvas in PDF trig: " << q << endl;
    temp = purityOL[q]; // print data canvases
    temp->Print(name);
    for(int e=0; e<numEtaBins; e++)
    {  
      if(DEBUG) cout << "Place Canvas in PDF eta: " << e << endl;
      temp = purityTrigOL[q][e]; // print data canvases
      temp->Print(name);
      temp = purityTrigOL2[q][e]; // print data canvases
      temp->Print(name);
    }
  }

  sprintf(name, "%s_purityOverlays_CompareEta.pdf]", baseName);
  temp->Print(name);


  return 1;
}
