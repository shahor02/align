#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <fstream>
#include "AliAlgRes.h"
#include "HistoManager.h"
#include "AliGeomManager.h"
#include "AliITSgeomTGeo.h"
#endif


void ProcRes(const char* inpData, const char* outName=0, const char* outDir=0);
TChain* LoadChain(const char* inpData, const char* chName="res");
void DoOffsets();
void BookHistos();
void BookHistosVTX(HistoManager* hm);
void BookHistosITS(HistoManager* hm);
void BookHistosTRD(HistoManager* hm);
void BookHistosTOF(HistoManager* hm);

void PostProcess(HistoManager* hm, HistoManager* hmpr);
void PostProcessVTX(HistoManager* hm, HistoManager* hmpr);
void PostProcessITS(HistoManager* hm, HistoManager* hmpr);  
void PostProcessTRD(HistoManager* hm, HistoManager* hmpr);  
void PostProcessTOF(HistoManager* hm, HistoManager* hmpr);  


void ProcessRes(Int_t ip);
void FillVTX(Int_t ip);
void FillITS(Int_t ip);
void FillTRD(Int_t ip);
void FillTOF(Int_t ip);

void DrawReport(TObjArray* hmans, const char* psname="algRep");
void DrawReportVTX(TObjArray* hmans, const char* psnm);
void DrawReportITS(TObjArray* hmans, const char* psnm);
void DrawReportTRD(TObjArray* hmans, const char* psnm);
void DrawReportTOF(TObjArray* hmans, const char* psnm);
void DrawHistos(TObjArray* hmans, int id);

Bool_t FitMeanSlope(TObjArray *histos, int minEnt=50);
Bool_t FitProfile(TObjArray *histos, int minEnt=50);

AliAlgRes* res=0;
TChain* ch = 0;
HistoManager* hman=0, *hmanProc=0;
TCanvas* repCanv=0;

//---------------------------------------------------
enum {kHOffsVTX=1000,kHOffsITS=100000,kHOffsTRD=300000,kHOffsTOF=400000};
enum {kDYSnp,kDZSnp,kDYPullSnp,kDZPullSnp};
enum {kNLrITS=6,kNLrTRD=6};
//
int kITSLrOffs[kNLrITS] = {0};
int kTRDLrOffs[kNLrTRD] = {0};
//
const int kNBinPull = 50;
const double kMaxPull = 5;
//
const int kNBinsResVTX = 100;
const int kNBinsAlpVTX = 18;
//
const int kNBinsResITS = 100;
const int kNBinsSnpITS = 10;
//
const int kNBinsResTRD = 100;
const int kNBinsSnpTRD = 10;
//
const int kNBinsResTOF = 100;
const int kNBinsSnpTOF = 10;
//
const double kSnpMaxITS = 0.6;
const double kSnpMaxTRD = 0.6;
const double kSnpMaxTOF = 0.6;
//
//
const double kMaxDYVTX = 0.02;
const double kMaxDZVTX = 0.05;
//
const double kMaxDYITS[kNLrITS] = {0.01,0.01,0.10,0.10,0.01,0.01};
const double kMaxDZITS[kNLrITS] = {0.10,0.10,0.10,0.10,0.20,0.10};
//
const double kMaxDYTRD[kNLrTRD] = {1.00,1.00,1.00,1.00,1.00,1.00};
const double kMaxDZTRD[kNLrTRD] = {10.0,10.0,10.0,10.0,10.0,10.0};
//
const double kMaxDYTOF = 10.;
const double kMaxDZTOF = 10.;
//---------------------------------------------------


void ProcRes(const char* inpData, const char* outName, const char* outDir)
{
  ch = LoadChain(inpData);
  if (!ch) return;
  res = new AliAlgRes();
  ch->SetBranchAddress("t",&res);
  int nent = ch->GetEntries();
  //
  BookHistos();
  printf("Processing %d events\n",nent);
  //
  for (int ient=0;ient<nent;ient++) {
    ch->GetEntry(ient);
    //
    int np = res->GetNPoints();
    for (int ip=np;ip--;) {
      ProcessRes(ip);
    }
  }
  //
  PostProcess(hman,hmanProc);
  //
  if (outName) {
    TString outNS = outName;
    if (outNS.IsNull()) outNS = "algHOut";
    if (!outNS.EndsWith(".root")) outNS += ".root";
    hman->SetFileName(outNS.Data());
    hmanProc->SetFileName(outNS.Data());
    //
    TString outDirR = "hraw";
    TString outDirP = "hpost";
    if (outDir) {
      outDirR += outDir;
      outDirP += outDir;
      hman->AddPrefix(outDir);
      hmanProc->AddPrefix(outDir);
    }
    hman->SetDirName(outDirR.Data());
    hmanProc->SetDirName(outDirP.Data());
    //
    hman->Write();
    hmanProc->Write();
  }
}

//__________________________________________
void ProcessRes(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID);
  if   (lrID==0)                                                       FillVTX(ip);
  else if (lrID>=AliGeomManager::kSPD1 && lrID<=AliGeomManager::kSSD2) FillITS(ip);
  else if (lrID>=AliGeomManager::kTRD1 && lrID<=AliGeomManager::kTRD6) FillTRD(ip);
  else if (lrID==AliGeomManager::kTOF)                                 FillTOF(ip);
  //
}

//__________________________________________
void FillITS(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID);
  int offs = kHOffsITS + (kITSLrOffs[lrID-AliGeomManager::kSPD1]+modID)*10;
  hman->GetHisto2F(offs+kDYSnp)->Fill(res->GetSnp(ip),res->GetDY(ip));
  hman->GetHisto2F(offs+kDZSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip));
  hman->GetHisto2F(offs+kDYPullSnp)->Fill(res->GetSnp(ip),res->GetDY(ip)/res->GetSigmaY(ip));
  hman->GetHisto2F(offs+kDZPullSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip)/res->GetSigmaZ(ip));
  //
}

//__________________________________________
void FillTRD(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID);
  int offs = kHOffsTRD + (kTRDLrOffs[lrID-AliGeomManager::kTRD1]+modID)*10;
  hman->GetHisto2F(offs+kDYSnp)->Fill(res->GetSnp(ip),res->GetDY(ip));
  hman->GetHisto2F(offs+kDZSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip));
  hman->GetHisto2F(offs+kDYPullSnp)->Fill(res->GetSnp(ip),res->GetDY(ip)/res->GetSigmaY(ip));
  hman->GetHisto2F(offs+kDZPullSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip)/res->GetSigmaZ(ip));
  //
}

//__________________________________________
void FillTOF(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  AliGeomManager::VolUIDToLayer(volID,modID);
  int offs = kHOffsTOF + modID*10;
  hman->GetHisto2F(offs+kDYSnp)->Fill(res->GetSnp(ip),res->GetDY(ip));
  hman->GetHisto2F(offs+kDZSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip));
  hman->GetHisto2F(offs+kDYPullSnp)->Fill(res->GetSnp(ip),res->GetDY(ip)/res->GetSigmaY(ip));
  hman->GetHisto2F(offs+kDZPullSnp)->Fill(res->GetSnp(ip),res->GetDZ(ip)/res->GetSigmaZ(ip));
  //
}

//__________________________________________
void FillVTX(int ip)
{
  int offs = kHOffsVTX;
  hman->GetHisto2F(offs+kDYSnp)->Fill(res->GetAlpha(ip),res->GetDY(ip));
  hman->GetHisto2F(offs+kDZSnp)->Fill(res->GetAlpha(ip),res->GetDZ(ip));
  hman->GetHisto2F(offs+kDYPullSnp)->Fill(res->GetAlpha(ip),res->GetDY(ip)/res->GetSigmaY(ip));
  hman->GetHisto2F(offs+kDZPullSnp)->Fill(res->GetAlpha(ip),res->GetDZ(ip)/res->GetSigmaZ(ip));
  //
}


//________________________________________________________
TChain* LoadChain(const char* inpData, const char* chName)
{
  TChain* chain = new TChain(chName);
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    chain->AddFile(inpData);
  }
  else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return 0;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while ( !flName.IsNull() ) {
      flName = flName.Strip(TString::kBoth,' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
      flName = flName.Strip(TString::kBoth,',');
      flName = flName.Strip(TString::kBoth,'"');
      printf("Adding %s\n",flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n<1) {
    printf("Obtained chain is empty\n");
    return 0;
  }
  else printf("Opened %s chain with %d entries\n",chName,n);
  return chain;
}

//___________________________________________
void BookHistos()
{
  DoOffsets();
  hman = new HistoManager();
  hmanProc = new HistoManager();
  BookHistosVTX(hman);
  BookHistosITS(hman);
  BookHistosTRD(hman);
  BookHistosTOF(hman);
}


//___________________________________________
void BookHistosITS(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsITS;
  int cntM = 0;
  for (int ilr=0;ilr<kNLrITS;ilr++) {
    int nmod = AliGeomManager::LayerSize(ilr+AliGeomManager::kSPD1);
    for (int imod=0;imod<nmod;imod++) {
      int moffs = hoffs + cntM*10;
      //
      hnm = Form("ITSDYvsSnp_L%d_M%d",ilr,imod);
      htl = Form("ITS #DeltaY vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpITS,-kSnpMaxITS,kSnpMaxITS,
		    kNBinsResITS,-kMaxDYITS[ilr],kMaxDYITS[ilr]);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaY");
      //
      hm->AddHisto(h2, moffs + kDYSnp);
      //
      hnm = Form("ITSDZvsSnp_L%d_M%d",ilr,imod);
      htl = Form("ITS #DeltaZ vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpITS,-kSnpMaxITS,kSnpMaxITS,
		    kNBinsResITS,-kMaxDZITS[ilr],kMaxDZITS[ilr]);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaZ");
      //
      hm->AddHisto(h2, moffs + kDZSnp);
      //
      //--------------------
      //
      hnm = Form("ITSDYPullvsSnp_L%d_M%d",ilr,imod);
      htl = Form("ITS #DeltaY Pull vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpITS,-kSnpMaxITS,kSnpMaxITS,
		    kNBinPull,-kMaxPull,kMaxPull);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaY Pull");
      //
      hm->AddHisto(h2, moffs + kDYPullSnp);
      //
      hnm = Form("ITSDZPullvsSnp_L%d_M%d",ilr,imod);
      htl = Form("ITS #DeltaZ Pull vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpITS,-kSnpMaxITS,kSnpMaxITS,
		    kNBinPull,-kMaxPull,kMaxPull);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaZ Pull");
      //
      hm->AddHisto(h2, moffs + kDZPullSnp);
      //
      //
      cntM++;
    }
  }
}


//___________________________________________
void BookHistosTRD(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsTRD;
  int cntM = 0;
  for (int ilr=0;ilr<kNLrTRD;ilr++) {
    int nmod = AliGeomManager::LayerSize(ilr+AliGeomManager::kTRD1);
    for (int imod=0;imod<nmod;imod++) {
      int moffs = hoffs + cntM*10;
      //
      hnm = Form("TRDYvsSnp_L%d_M%d",ilr,imod);
      htl = Form("TRD #DeltaY vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpTRD,-kSnpMaxTRD,kSnpMaxTRD,
		    kNBinsResTRD,-kMaxDYTRD[ilr],kMaxDYTRD[ilr]);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaY");
      //
      hm->AddHisto(h2, moffs + kDYSnp);
      //
      hnm = Form("TRDDZvsSnp_L%d_M%d",ilr,imod);
      htl = Form("TRD #DeltaZ vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpTRD,-kSnpMaxTRD,kSnpMaxTRD,
		    kNBinsResTRD,-kMaxDZTRD[ilr],kMaxDZTRD[ilr]);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaZ");
      //
      hm->AddHisto(h2, moffs + kDZSnp);
      //
      //---------------------------------
      hnm = Form("TRDYPullvsSnp_L%d_M%d",ilr,imod);
      htl = Form("TRD #DeltaY Pull vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpTRD,-kSnpMaxTRD,kSnpMaxTRD,
		    kNBinPull,-kMaxPull,kMaxPull);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaY Pull");
      //
      hm->AddHisto(h2, moffs + kDYPullSnp);
      //
      hnm = Form("TRDZPullvsSnp_L%d_M%d",ilr,imod);
      htl = Form("TRD #DeltaZ Pull vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpTRD,-kSnpMaxTRD,kSnpMaxTRD,
		    kNBinPull,-kMaxPull,kMaxPull);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaZ Pull");
      //
      hm->AddHisto(h2, moffs + kDZPullSnp);
      //
      cntM++;
    }
  }
}

//___________________________________________
void BookHistosTOF(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsTOF;
  int cntM = 0;
  int nmod = AliGeomManager::LayerSize(AliGeomManager::kTOF);
  for (int imod=0;imod<nmod;imod++) {
    int moffs = hoffs + cntM*10;
    //
    hnm = Form("TOFYvsSnp_M%d",imod);
    htl = Form("TOF #DeltaY vs Snp Mod%d",imod);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  kNBinsSnpTOF,-kSnpMaxTOF,kSnpMaxTOF,
		  kNBinsResTOF,-kMaxDYTOF,kMaxDYTOF);
    h2->SetXTitle("snp");
    h2->SetYTitle("#DeltaY");
    //
    hm->AddHisto(h2, moffs + kDYSnp);
    //
    hnm = Form("TOFDZvsSnp_M%d",imod);
    htl = Form("TOF #DeltaZ vs Snp Mod%d",imod);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  kNBinsSnpTOF,-kSnpMaxTOF,kSnpMaxTOF,
		  kNBinsResTOF,-kMaxDZTOF,kMaxDZTOF);
    h2->SetXTitle("snp");
    h2->SetYTitle("#DeltaZ");
    //
    hm->AddHisto(h2, moffs + kDZSnp);
    //
    hnm = Form("TOFYPullvsSnp_M%d",imod);
    htl = Form("TOF #DeltaY Pull vs Snp Mod%d",imod);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  kNBinsSnpTOF,-kSnpMaxTOF,kSnpMaxTOF,
		  kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("snp");
    h2->SetYTitle("#DeltaY Pull");
    //
    hm->AddHisto(h2, moffs + kDYPullSnp);
    //
    hnm = Form("TOFDZPullvsSnp_M%d",imod);
    htl = Form("TOF #DeltaZ Pull vs Snp Mod%d",imod);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  kNBinsSnpTOF,-kSnpMaxTOF,kSnpMaxTOF,
		  kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("snp");
    h2->SetYTitle("#DeltaZ Pull");
    //
    hm->AddHisto(h2, moffs + kDZPullSnp);
    //
    cntM++;
  }
}

//___________________________________________
void BookHistosVTX(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsVTX;
  //
  hnm = Form("VTXYvsAlp");
  htl = Form("VTX #DeltaY vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),	       
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinsResVTX,-kMaxDYVTX,kMaxDYVTX);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaY");
  //
  hm->AddHisto(h2, hoffs + kDYSnp);
  //
  hnm = Form("VTXDZvsAlp");
  htl = Form("VTX #DeltaZ vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinsResVTX,-kMaxDZVTX,kMaxDZVTX);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaZ");
  //
  hm->AddHisto(h2, hoffs + kDZSnp);
  //
  hnm = Form("VTXYPullvsAlp");
  htl = Form("VTX #DeltaY Pull vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),	       
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinPull,-kMaxPull,kMaxPull);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaY Pull");
  //
  hm->AddHisto(h2, hoffs + kDYPullSnp);
  //
  hnm = Form("VTXDZPullvsAlp");
  htl = Form("VTX #DeltaZ Pull vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinPull,-kMaxPull,kMaxPull);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaZ Pull");
  //
  hm->AddHisto(h2, hoffs + kDZPullSnp);
  //
}

//___________________________________________
void PostProcess(HistoManager* hm, HistoManager* hmpr)
{
  // process filled histos
  PostProcessVTX(hm,hmpr);
  PostProcessITS(hm,hmpr);  
  PostProcessTRD(hm,hmpr);  
  PostProcessTOF(hm,hmpr);  
  //
}

//___________________________________________
void PostProcessVTX(HistoManager* hm,HistoManager* hmProc)
{
  // process filled histos
  int hoffs = kHOffsVTX;
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h;
  TH1* h1;
  //
  h = hm->GetHisto2F(hoffs + kDYSnp);
  if (h) arrY.Add(h);
  h = hm->GetHisto2F(hoffs + kDZSnp);
  if (h) arrZ.Add(h);
  h = hm->GetHisto2F(hoffs + kDYPullSnp);
  if (h) arrYP.Add(h);
  h = hm->GetHisto2F(hoffs + kDZPullSnp);
  if (h) arrZP.Add(h);
  //
  if (FitProfile(&arrY)) {
    h1 = (TH1*)arrY.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_Y_Mean","DCA_{Y} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDYSnp*10+0);
    }
    h1 = (TH1*)arrY.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_Y_Sigm","DCA_{Y} sigma");
      hmProc->AddHisto(h1, kHOffsVTX + kDYSnp*10+1);
    }
  }
  //
  if (FitProfile(&arrZ)) {
    h1 = (TH1*)arrZ.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_Z_Mean","DCA_{Z} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDZSnp*10+0);
    }
    h1 = (TH1*)arrZ.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_Z_Sigm","DCA_{Z} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDZSnp*10+1);
    }
  }
  //
  if (FitProfile(&arrYP)) {
    h1 = (TH1*)arrYP.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_YPull_Mean","DCA_{YPull} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDYPullSnp*10+0);
    }
    h1 = (TH1*)arrYP.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_YPull_Sigm","DCA_{YPull} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDYPullSnp*10+1);
    }
  }
  //
  if (FitProfile(&arrZP)) {
    h1 = (TH1*)arrZP.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_ZPull_Mean","DCA_{ZPull} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDZPullSnp*10+0);
    }
    h1 = (TH1*)arrZP.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_ZPull_Sigm","DCA_{ZPull} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDZPullSnp*10+1);
    }
  }
  //
}


//___________________________________________
void PostProcessTOF(HistoManager* hm,HistoManager* hmProc)
{
  // process filled histos
  int hoffs = kHOffsTOF;
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  // 
  int cntM=0;
  int nmod = AliGeomManager::LayerSize(AliGeomManager::kTOF);
  int nmodSM = nmod/18;
  //
  for (int isc=0;isc<18;isc++) {
    for (int im=0;im<nmodSM;im++) {
      h = hm->GetHisto2F(hoffs+cntM*10+kDYSnp);      if (h) arrY.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZSnp);      if (h) arrZ.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDYPullSnp);  if (h) arrYP.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZPullSnp);  if (h) arrZP.Add(h);
      cntM++;
    }
    //
    //    arrY.Print();
    if (FitMeanSlope(&arrY)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrY[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TOF_DY_SM%d_%s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TOF DY SM%d %s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+isc*100 + kDYSnp*10 + i*2+j);
	}
      }
    }
    arrY.Clear();
    //
    if (FitMeanSlope(&arrZ)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZ[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TOF_DZ_SM%d_%s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TOF DZ SM%d %s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+isc*100 + kDZSnp*10 + i*2+j);
	}
      }
    }
    arrZ.Clear();
    //
    //
    if (FitMeanSlope(&arrYP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrYP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TOF_DYPull_SM%d_%s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TOF DY Pull SM%d %s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+isc*100 + kDYPullSnp*10 + i*2+j);
	}
      }
    }
    arrYP.Clear();
    //
    if (FitMeanSlope(&arrZP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TOF_DZPull_SM%d_%s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TOF DZ Pull SM%d %s%s",isc,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+isc*100 + kDZPullSnp*10 + i*2+j);
	}
      }
    }
    arrZP.Clear();
    //
  }
  //  
  //
}

//___________________________________________
void PostProcessITS(HistoManager* hm, HistoManager* hmProc) 
{
  // process filled histos
  int hoffs = kHOffsITS;
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  // 
  int cntM=0;
  for (int ilr=0;ilr<6;ilr++) {
    int nmod = AliGeomManager::LayerSize(AliGeomManager::kSPD1+ilr);
    for (int im=0;im<nmod;im++) {
      h = hm->GetHisto2F(hoffs+cntM*10+kDYSnp);      if (h) arrY.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZSnp);      if (h) arrZ.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDYPullSnp);  if (h) arrYP.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZPullSnp);  if (h) arrZP.Add(h);
      cntM++;
    }
    //
    //    arrY.Print();
    if (FitMeanSlope(&arrY)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrY[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("ITS_DY_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("ITS DY Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDYSnp*10 + i*2+j);
	}
      }
    }
    arrY.Clear();
    //
    if (FitMeanSlope(&arrZ)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZ[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("ITS_DZ_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("ITS DZ Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDZSnp*10 + i*2+j);
	}
      }
    }
    arrZ.Clear();
    //
    //
    if (FitMeanSlope(&arrYP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrYP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("ITS_DYPull_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("ITS DY Pull Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDYPullSnp*10 + i*2+j);
	}
      }
    }
    arrYP.Clear();
    //
    if (FitMeanSlope(&arrZP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("ITS_DZPull_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("ITS DZ Pull Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDZPullSnp*10 + i*2+j);
	}
      }
    }
    arrZP.Clear();
    //
  }
  //  
  //
}


//_______________________________________________________
void PostProcessTRD(HistoManager* hm, HistoManager* hmProc)
{
  // postprocess histos
  //
  int hoffs = kHOffsTRD;
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  // 
  int cntM=0;
  for (int ilr=0;ilr<6;ilr++) {
    int nmod = AliGeomManager::LayerSize(AliGeomManager::kTRD1+ilr);
    for (int im=0;im<nmod;im++) {
      h = hm->GetHisto2F(hoffs+cntM*10+kDYSnp);      if (h) arrY.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZSnp);      if (h) arrZ.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDYPullSnp);  if (h) arrYP.Add(h);
      h = hm->GetHisto2F(hoffs+cntM*10+kDZPullSnp);  if (h) arrZP.Add(h);
      cntM++;
    }
    //
    //    arrY.Print();
    if (FitMeanSlope(&arrY)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrY[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TRD_DY_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TRD DY Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDYSnp*10 + i*2+j);
	}
      }
    }
    arrY.Clear();
    //
    if (FitMeanSlope(&arrZ)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZ[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TRD_DZ_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TRD DZ Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDZSnp*10 + i*2+j);
	}
      }
    }
    arrZ.Clear();
    //
    //
    if (FitMeanSlope(&arrYP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrYP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TRD_DYPull_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TRD DY Pull Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDYPullSnp*10 + i*2+j);
	}
      }
    }
    arrYP.Clear();
    //
    if (FitMeanSlope(&arrZP)) {
      for (int i=0;i<2;i++) {
	for (int j=0;j<2;j++) {
	  TH1* h1 = (TH1*)arrZP[i*2+j]; if (!h1) continue;
	  h1->SetName(Form("TRD_DZPull_Lr%d_%s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  h1->SetTitle(Form("TRD DZ Pull Lr%d %s%s",ilr,i ? "Sigma":"Mean", j ? "P1":"P0"));
	  hmProc->AddHisto(h1, hoffs+ilr*100 + kDZPullSnp*10 + i*2+j);
	}
      }
    }
    arrZP.Clear();
    //
  }
  //  
}


//___________________________________________
void DoOffsets()
{
  for (int i=0,n=0;i<kNLrITS;i++) {kITSLrOffs[i] = n; n += AliGeomManager::LayerSize(i+AliGeomManager::kSPD1);}
  for (int i=0,n=0;i<kNLrTRD;i++) {kTRDLrOffs[i] = n; n += AliGeomManager::LayerSize(i+AliGeomManager::kTRD1);}
}

//___________________________________________
Bool_t FitMeanSlope(TObjArray *histos, int minEnt)
{
  // fit mean and slope
  int nhist = histos->GetEntriesFast();
  TH2* h = (TH2*)histos->At(0);
  double range = h->GetYaxis()->GetXmax();
  TF1* gs = new TF1("gs","gaus",-range,range);
  TObjArray hfarr;
  TH1F* hOutMeanP0 = new TH1F("_dumMeanP0","",nhist,0,nhist);
  TH1F* hOutMeanP1 = new TH1F("_dumMeanP1","",nhist,0,nhist);
  TH1F* hOutSigP0  = new TH1F("_dumSigP0","",nhist,0,nhist);
  TH1F* hOutSigP1  = new TH1F("_dumSigP1","",nhist,0,nhist);
  for (int ih=0;ih<nhist;ih++) {
    h = (TH2*)histos->At(ih);
    TString nm = h->GetName();
    int idm = nm.Last('_');
    if (idm>0) {
      hOutMeanP0->GetXaxis()->SetBinLabel(ih+1,nm.Data()+idm+1);
      hOutMeanP1->GetXaxis()->SetBinLabel(ih+1,nm.Data()+idm+1);
      hOutSigP0->GetXaxis()->SetBinLabel(ih+1,nm.Data()+idm+1);
      hOutSigP1->GetXaxis()->SetBinLabel(ih+1,nm.Data()+idm+1);
    }
    if (!h || h->GetEntries()<minEnt) continue;
    h->FitSlicesY(gs,-1,-1,0,"qmrl",&hfarr);
    TH1* hmean = (TH1*)hfarr[1];
    TH1* hsig  = (TH1*)hfarr[2];
    if (!hmean || !hsig) {hfarr.Delete(); continue;}
    hmean->Fit("pol1","qmr");
    TF1* plm = (TF1*)hmean->GetListOfFunctions()->FindObject("pol1");
    if (plm) {
      hOutMeanP0->SetBinContent(ih+1,plm->GetParameter(0));
      hOutMeanP0->SetBinError(ih+1,plm->GetParError(0));
      hOutMeanP1->SetBinContent(ih+1,plm->GetParameter(1));
      hOutMeanP1->SetBinError(ih+1,plm->GetParError(1));
    }
    //
    hsig->Fit("pol1","qmr");
    TF1* pls = (TF1*)hsig->GetListOfFunctions()->FindObject("pol1");
    if (pls) {
      hOutSigP0->SetBinContent(ih+1,pls->GetParameter(0));
      hOutSigP0->SetBinError(ih+1,pls->GetParError(0));
      hOutSigP1->SetBinContent(ih+1,pls->GetParameter(1));
      hOutSigP1->SetBinError(ih+1,pls->GetParError(1));
    }
    //
    hfarr.Delete();
  }
  //
  delete gs;
  histos->Clear();
  histos->AddAtAndExpand(hOutMeanP0,0);
  histos->AddAtAndExpand(hOutMeanP1,1);
  histos->AddAtAndExpand(hOutSigP0,2);
  histos->AddAtAndExpand(hOutSigP1,3);
  //
  return kTRUE;
}

//___________________________________________
Bool_t FitProfile(TObjArray *histos, int minEnt)
{
  // fit profile
  TH2* h = (TH2*)histos->At(0);
  if (h->GetEntries()<minEnt) return kFALSE;
  double range = h->GetYaxis()->GetXmax();
  TF1* gs = new TF1("gs","gaus",-range,range);
  TObjArray hfarr;
  h->FitSlicesY(gs,-1,-1,0,"qmrl",&hfarr);
  TH1* hmean = (TH1*)hfarr[1];
  TH1* hsig  = (TH1*)hfarr[2];
  if (!hmean || !hsig) {hfarr.Delete(); return kFALSE;}
  delete gs;
  histos->Clear();
  histos->AddAtAndExpand(hfarr.RemoveAt(1),0);
  histos->AddAtAndExpand(hfarr.RemoveAt(2),1);
  //
  return kTRUE;
}

//============================================
void DrawReport(TObjArray* hmans, const char* psname)
{
  gStyle->SetOptStat(0);
  //
  repCanv = new TCanvas("algCanv","algCanv",700,900);
  TString psnm1 = psname;
  if (psnm1.IsNull()) psnm1 = Form("algRep");
  if (!psnm1.EndsWith(".ps")) psnm1 += ".ps";
  TString psnm0 = psnm1 + "["; 
  TString psnm2 = psnm1 + "]";
  repCanv->Print(psnm0.Data());
  //
  DrawReportVTX(hmans,psnm1.Data());
  DrawReportITS(hmans,psnm1.Data());  
  DrawReportTRD(hmans,psnm1.Data());    
  DrawReportTOF(hmans,psnm1.Data());  
  //
  repCanv->cd();
  repCanv->Print(psnm2.Data());
}

//____________________________________________
void DrawReportVTX(TObjArray* hmans, const char* psnm)
{
  repCanv->Clear();
  repCanv->Divide(2,2);
  //
  repCanv->cd(1);
  DrawHistos(hmans, kHOffsVTX + kDYSnp*10+0);
  repCanv->cd(2);
  DrawHistos(hmans, kHOffsVTX + kDYSnp*10+1);
  //
  repCanv->cd(3);
  DrawHistos(hmans, kHOffsVTX + kDYPullSnp*10+0);
  repCanv->cd(4);
  DrawHistos(hmans, kHOffsVTX + kDYPullSnp*10+1);
  //
  repCanv->cd();
  repCanv->Print(psnm); // Y
  //
  repCanv->Clear();
  repCanv->Divide(2,2);
  //  
  repCanv->cd(1);
  DrawHistos(hmans, kHOffsVTX + kDZSnp*10+0);
  repCanv->cd(2);
  DrawHistos(hmans, kHOffsVTX + kDZSnp*10+1);
  //
  repCanv->cd(3);
  DrawHistos(hmans, kHOffsVTX + kDZPullSnp*10+0);
  repCanv->cd(4);
  DrawHistos(hmans, kHOffsVTX + kDZPullSnp*10+1);
  //
  repCanv->cd();
  repCanv->Print(psnm); // Z
  //
}

//____________________________________________
void DrawReportITS(TObjArray* hmans, const char* psnm)
{
  int icn;
  for (int ilr=0;ilr<6;ilr++) {
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsITS + ilr*100 + kDYSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsITS + ilr*100 + kDYPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsITS + ilr*100 + kDZSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsITS + ilr*100 + kDZPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
  }
  //
}

//____________________________________________
void DrawReportTRD(TObjArray* hmans, const char* psnm)
{
  int icn;
  for (int ilr=0;ilr<6;ilr++) {
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTRD + ilr*100 + kDYSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTRD + ilr*100 + kDYPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
  }
  //
}

//____________________________________________
void DrawReportTOF(TObjArray* hmans, const char* psnm)
{
  int icn;
  for (int isc=0;isc<18;isc++) {
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTOF + isc*100 + kDYSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTOF + isc*100 + kDYPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTOF + isc*100 + kDZSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
    //
    icn = 1;
    repCanv->Clear();
    repCanv->Divide(1,4);
    for (int i=0;i<2;i++) 
      for (int j=0;j<2;j++) {
	repCanv->cd(icn++);
	DrawHistos(hmans, kHOffsTOF + isc*100 + kDZPullSnp*10 + i*2+j);
      }
    repCanv->Print(psnm);
    //
  }
  //
}

//_________________________________________________
void DrawHistos(TObjArray* hmans, int id)
{
  // draw histo ID from set of histo managers
  int nhm = hmans->GetEntriesFast();
  double mn=1e9,mx=-1e9;
  int nhAcc = 0;
  for (int ih=0;ih<nhm;ih++) {
    HistoManager* hm = (HistoManager*)hmans->At(ih);
    if (!hm) continue;
    TH1* h = hm->GetHisto1F(id);
    if (!h || h->GetEntries()==0) continue;
    if (mn>h->GetMinimum()) mn = h->GetMinimum();
    if (mx<h->GetMaximum()) mx = h->GetMaximum();
    //
    nhAcc++;
  }
  if (nhAcc<1 || mn>=mx) return;
  nhAcc = 0;
  gStyle->SetTitleW(0.9);
  for (int ih=0;ih<nhm;ih++) {
    HistoManager* hm = (HistoManager*)hmans->At(ih);
    if (!hm) continue;
    TH1* h = hm->GetHisto1F(id);
    if (!h || h->GetEntries()==0) continue;
    h->SetMinimum( mn - 0.2*(mx-mn) );
    h->SetMaximum( mx + 0.2*(mx-mn) );
    h->Draw(nhAcc==0 ? "":"same");
    nhAcc++;
  }
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();
}
