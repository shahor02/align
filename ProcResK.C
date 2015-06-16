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


void ProcResK(const char* inpData, const char* outName="halgResK", const char* outDir="r");
TChain* LoadChain(const char* inpData, const char* chName="res");
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
void DrawHistos(TObjArray* hmans, int id, float range=-1, float rangeM=-9999);

Bool_t FitMeanSlope(TObjArray *histos, int minEnt=50);
Bool_t FitProfile(TObjArray *histos, int minEnt=50);

AliAlgRes* res=0;
TChain* ch = 0;
HistoManager* hman=0, *hmanProc=0;
TCanvas* repCanv=0;

//---------------------------------------------------
enum {kHOffsVTX=1000,kHOffsITS=100000,kHOffsTRD=300000,kHOffsTOF=400000};
enum {kDY,kDZ,kDYPull,kDZPull};
enum {kNLrITS=6,kNLrTRD=6};
//
const char* kFitOptGS = "qmr"; // "qmrl"
const char* kFitOptP1 = "qmr"; //
//
const int kNBinPull = 50;
const double kMaxPull = 5;
//
const int kNBinsResVTX = 100;
const int kNBinsAlpVTX = 36;
//
const int kNBinsResITS = 100;
//
const int kNBinsResTRD = 100;
//
const int kNBinsResTOF = 100;
//
const double kMaxDYVTX = 0.05;
const double kMaxDZVTX = 0.1;
//
const double kMaxDYITS[kNLrITS] = {0.40,0.40,0.40,0.40,0.50,0.50};
const double kMaxDZITS[kNLrITS] = {0.20,0.20,0.40,0.40,0.50,0.50};
//
const double kMaxDYTRD[kNLrTRD] = {4.00,2.00,2.00,2.00,2.00,2.00};
const double kMaxDZTRD[kNLrTRD] = {10.0,10.0,10.0,10.0,10.0,10.0};
//
const double kMaxDYTOF = 10.;
const double kMaxDZTOF = 10.;
//---------------------------------------------------


void ProcResK(const char* inpData, const char* outName, const char* outDir)
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
    if (!res->GetKalmanDone()) continue;
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
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID)-AliGeomManager::kSPD1;
  int offs = kHOffsITS + lrID*10;
  hman->GetHisto2F(offs+kDY)->Fill(modID,res->GetDYK(ip));
  hman->GetHisto2F(offs+kDZ)->Fill(modID,res->GetDZK(ip));
  hman->GetHisto2F(offs+kDYPull)->Fill(modID,res->GetDYK(ip)/
				       TMath::Sqrt(res->GetSigY2(ip)+res->GetSigY2K(ip)));
  hman->GetHisto2F(offs+kDZPull)->Fill(modID,res->GetDZK(ip)/
				       TMath::Sqrt(res->GetSigZ2(ip)+res->GetSigZ2K(ip)));
  //
}

//__________________________________________
void FillTRD(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID)-AliGeomManager::kTRD1;
  int offs = kHOffsTRD + lrID*10;
  hman->GetHisto2F(offs+kDY)->Fill(modID,res->GetDYK(ip));
  hman->GetHisto2F(offs+kDZ)->Fill(modID,res->GetDZK(ip));
  hman->GetHisto2F(offs+kDYPull)->Fill(modID,res->GetDYK(ip)/
				       TMath::Sqrt(res->GetSigY2(ip)+
						   res->GetSigY2K(ip)));
  hman->GetHisto2F(offs+kDZPull)->Fill(modID,res->GetDZK(ip)/
				       TMath::Sqrt(res->GetSigZ2(ip)+
						   res->GetSigZ2K(ip)));
  //
}

//__________________________________________
void FillTOF(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  AliGeomManager::VolUIDToLayer(volID,modID);
  int offs = kHOffsTOF + (modID/91)*100;
  modID = modID%91;
  hman->GetHisto2F(offs+kDY)->Fill(modID,res->GetDYK(ip));
  hman->GetHisto2F(offs+kDZ)->Fill(modID,res->GetDZK(ip));
  hman->GetHisto2F(offs+kDYPull)->Fill(modID,res->GetDYK(ip)/
				       TMath::Sqrt(res->GetSigY2(ip)+res->GetSigY2K(ip)));
  hman->GetHisto2F(offs+kDZPull)->Fill(modID,res->GetDZK(ip)/
				       TMath::Sqrt(res->GetSigZ2(ip)+res->GetSigZ2K(ip)));
  //
}

//__________________________________________
void FillVTX(int ip)
{
  int offs = kHOffsVTX;
  hman->GetHisto2F(offs+kDY)->Fill(res->GetAlpha(ip),res->GetDYK(ip));
  hman->GetHisto2F(offs+kDZ)->Fill(res->GetAlpha(ip),res->GetDZK(ip));
  hman->GetHisto2F(offs+kDYPull)->Fill(res->GetAlpha(ip),res->GetDYK(ip)/
					  TMath::Sqrt(res->GetSigY2(ip)+res->GetSigY2K(ip)));
  hman->GetHisto2F(offs+kDZPull)->Fill(res->GetAlpha(ip),res->GetDZ(ip)/
					  TMath::Sqrt(res->GetSigZ2(ip)+res->GetSigZ2K(ip)));
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
  for (int ilr=0;ilr<kNLrITS;ilr++) {
    int moffs = hoffs + ilr*10;
    int nmod = AliGeomManager::LayerSize(ilr+AliGeomManager::kSPD1);
    hnm = Form("ITSDY_L%d",ilr);
    htl = Form("ITS #DeltaY Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,kNBinsResITS,-kMaxDYITS[ilr],kMaxDYITS[ilr]);
    h2->SetXTitle("module");
    h2->SetYTitle("#DeltaY");
    hm->AddHisto(h2, moffs + kDY);
    //
    hnm = Form("ITSDZ_L%d",ilr);
    htl = Form("ITS #DeltaZ Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,kNBinsResITS,-kMaxDZITS[ilr],kMaxDZITS[ilr]);
    h2->SetXTitle("module");
    h2->SetYTitle("#DeltaZ");
    hm->AddHisto(h2, moffs + kDZ);
    //
    //--------------------
    //
    hnm = Form("ITSDYPull_L%d",ilr);
    htl = Form("ITS #DeltaY Pull Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("module");
    h2->SetYTitle("#DeltaY Pull");
    hm->AddHisto(h2, moffs + kDYPull);
    //
    hnm = Form("ITSDZPull_L%d",ilr);
    htl = Form("ITS #DeltaZ Pull Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("module");
    h2->SetYTitle("#DeltaZ Pull");
    hm->AddHisto(h2, moffs + kDZPull);
  }
}


//___________________________________________
void BookHistosTRD(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsTRD;
  for (int ilr=0;ilr<kNLrTRD;ilr++) {
    int nmod = AliGeomManager::LayerSize(ilr+AliGeomManager::kTRD1);
    int moffs = hoffs + ilr*10;
    //
    hnm = Form("TRDY_L%d",ilr);
    htl = Form("TRD #DeltaY Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTRD,-kMaxDYTRD[ilr],kMaxDYTRD[ilr]);
    h2->SetXTitle("ch");
    h2->SetYTitle("#DeltaY");
    hm->AddHisto(h2, moffs + kDY);
    //
    hnm = Form("TRDZ_L%d",ilr);
    htl = Form("TRD #DeltaZ Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTRD,-kMaxDZTRD[ilr],kMaxDZTRD[ilr]);
    h2->SetXTitle("ch");
    h2->SetYTitle("#DeltaZ");
    hm->AddHisto(h2, moffs + kDZ);
    //
    //---------------------------------
    //
    hnm = Form("TRDYPull_L%d",ilr);
    htl = Form("TRD #DeltaY Pull Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("ch");
    h2->SetYTitle("#DeltaY Pull");
    hm->AddHisto(h2, moffs + kDYPull);
    //
    hnm = Form("TRDZPull_L%d",ilr);
    htl = Form("TRD #DeltaZ Pull Lr%d",ilr);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinPull,-kMaxPull,kMaxPull);
    h2->SetXTitle("ch");
    h2->SetYTitle("#DeltaZ Pull");
    hm->AddHisto(h2, moffs + kDZPull);
  }
}

//___________________________________________
void BookHistosTOF(HistoManager* hm)
{
  //
  TH2* h2=0;
  TString hnm,htl;
  int hoffs = kHOffsTOF;
  int nmod = AliGeomManager::LayerSize(AliGeomManager::kTOF)/18;
  for (int isc=0;isc<18;isc++) {
    int moffs = hoffs + isc*100;
    //
    hnm = Form("TOFY_SM%d",isc);
    htl = Form("TOF #DeltaY SM%d",isc);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTOF,-kMaxDYTOF,kMaxDYTOF);
    h2->SetXTitle("#strip");
    h2->SetYTitle("#DeltaY");
    hm->AddHisto(h2, moffs + kDY);
    //
    hnm = Form("TOFZ_SM%d",isc);
    htl = Form("TOF #DeltaZ SM%d",isc);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTOF,-kMaxDZTOF,kMaxDZTOF);
    h2->SetXTitle("#strip");
    h2->SetYTitle("#DeltaZ");
    hm->AddHisto(h2, moffs + kDZ);
    //
    //-----------------
    //
    hnm = Form("TOFYPull_SM%d",isc);
    htl = Form("TOF #DeltaYPull SM%d",isc);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTOF,-kMaxPull,kMaxPull);
    h2->SetXTitle("#strip");
    h2->SetYTitle("#DeltaYPull");
    hm->AddHisto(h2, moffs + kDYPull);
    //
    hnm = Form("TOFZPull_SM%d",isc);
    htl = Form("TOF #DeltaZPull SM%d",isc);
    h2 = new TH2F(hnm.Data(),htl.Data(),
		  nmod,0,nmod,
		  kNBinsResTOF,-kMaxPull,kMaxPull);
    h2->SetXTitle("#strip");
    h2->SetYTitle("#DeltaZPull");
    hm->AddHisto(h2, moffs + kDZPull);
    //
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
  hm->AddHisto(h2, hoffs + kDY);
  //
  hnm = Form("VTXDZvsAlp");
  htl = Form("VTX #DeltaZ vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinsResVTX,-kMaxDZVTX,kMaxDZVTX);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaZ");
  //
  hm->AddHisto(h2, hoffs + kDZ);
  //
  hnm = Form("VTXYPullvsAlp");
  htl = Form("VTX #DeltaY Pull vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),	       
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinPull,-kMaxPull,kMaxPull);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaY Pull");
  //
  hm->AddHisto(h2, hoffs + kDYPull);
  //
  hnm = Form("VTXDZPullvsAlp");
  htl = Form("VTX #DeltaZ Pull vs #alpha");
  h2 = new TH2F(hnm.Data(),htl.Data(),
		kNBinsAlpVTX,-TMath::Pi(),TMath::Pi(),
		kNBinPull,-kMaxPull,kMaxPull);
  h2->SetXTitle("#alpha");
  h2->SetYTitle("#DeltaZ Pull");
  //
  hm->AddHisto(h2, hoffs + kDZPull);
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
  h = hm->GetHisto2F(hoffs + kDY);
  if (h) arrY.Add(h);
  h = hm->GetHisto2F(hoffs + kDZ);
  if (h) arrZ.Add(h);
  h = hm->GetHisto2F(hoffs + kDYPull);
  if (h) arrYP.Add(h);
  h = hm->GetHisto2F(hoffs + kDZPull);
  if (h) arrZP.Add(h);
  //
  if (FitProfile(&arrY)) {
    h1 = (TH1*)arrY.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_Y_Mean","DCA_{Y} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDY*10+0);
    }
    h1 = (TH1*)arrY.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_Y_Sigm","DCA_{Y} sigma");
      hmProc->AddHisto(h1, kHOffsVTX + kDY*10+1);
    }
  }
  //
  if (FitProfile(&arrZ)) {
    h1 = (TH1*)arrZ.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_Z_Mean","DCA_{Z} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDZ*10+0);
    }
    h1 = (TH1*)arrZ.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_Z_Sigm","DCA_{Z} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDZ*10+1);
    }
  }
  //
  if (FitProfile(&arrYP)) {
    h1 = (TH1*)arrYP.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_YPull_Mean","DCA_{YPull} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDYPull*10+0);
    }
    h1 = (TH1*)arrYP.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_YPull_Sigm","DCA_{YPull} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDYPull*10+1);
    }
  }
  //
  if (FitProfile(&arrZP)) {
    h1 = (TH1*)arrZP.RemoveAt(0);
    if (h1) {
      h1->SetNameTitle("DCA_ZPull_Mean","DCA_{ZPull} mean");
      hmProc->AddHisto(h1, kHOffsVTX + kDZPull*10+0);
    }
    h1 = (TH1*)arrZP.RemoveAt(1);
    if (h1) {
      h1->SetNameTitle("DCA_ZPull_Sigm","DCA_{ZPull} sigma");      
      hmProc->AddHisto(h1, kHOffsVTX + kDZPull*10+1);
    }
  }
  //
}


//___________________________________________
void PostProcessTOF(HistoManager* hm,HistoManager* hmProc)
{
  // process filled histos
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  TH1* h1=0;
  // 
  for (int isc=0;isc<18;isc++) {
    //
    int hoffs = kHOffsTOF + isc*100;
    h = hm->GetHisto2F(hoffs + kDY);
    if (h) arrY.Add(h);
    h = hm->GetHisto2F(hoffs + kDZ);
    if (h) arrZ.Add(h);
    h = hm->GetHisto2F(hoffs + kDYPull);
    if (h) arrYP.Add(h);
    h = hm->GetHisto2F(hoffs + kDZPull);
    if (h) arrZP.Add(h);
    //
    if (FitProfile(&arrY)) {
      h1 = (TH1*)arrY.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDY*10+0);
      }
      h1 = (TH1*)arrY.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDY*10+1);
      }
    }
    //
    if (FitProfile(&arrZ)) {
      h1 = (TH1*)arrZ.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDZ*10+0);
      }
      h1 = (TH1*)arrZ.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDZ*10+1);
      }
    }
    //
    if (FitProfile(&arrYP)) {
      h1 = (TH1*)arrYP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDYPull*10+0);
      }
      h1 = (TH1*)arrYP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDYPull*10+1);
      }
    }
    //
    if (FitProfile(&arrZP)) {
      h1 = (TH1*)arrZP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDZPull*10+0);
      }
      h1 = (TH1*)arrZP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffs + kDZPull*10+1);
      }
    }
    //
    arrY.Clear();
    arrZ.Clear();
    arrYP.Clear();
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
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  TH1* h1=0;
  // 
  for (int ilr=0;ilr<6;ilr++) {
    //
    int hoffs = kHOffsITS + ilr*10;
    int hoffsP= kHOffsITS + ilr*100;
    h = hm->GetHisto2F(hoffs + kDY);
    if (h) arrY.Add(h);
    h = hm->GetHisto2F(hoffs + kDZ);
    if (h) arrZ.Add(h);
    h = hm->GetHisto2F(hoffs + kDYPull);
    if (h) arrYP.Add(h);
    h = hm->GetHisto2F(hoffs + kDZPull);
    if (h) arrZP.Add(h);
    //
    if (FitProfile(&arrY)) {
      h1 = (TH1*)arrY.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDY*10+0);
      }
      h1 = (TH1*)arrY.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDY*10+1);
      }
    }
    //
    if (FitProfile(&arrZ)) {
      h1 = (TH1*)arrZ.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZ*10+0);
      }
      h1 = (TH1*)arrZ.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZ*10+1);
      }
    }
    //
    if (FitProfile(&arrYP)) {
      h1 = (TH1*)arrYP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDYPull*10+0);
      }
      h1 = (TH1*)arrYP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDYPull*10+1);
      }
    }
    //
    if (FitProfile(&arrZP)) {
      h1 = (TH1*)arrZP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZPull*10+0);
      }
      h1 = (TH1*)arrZP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZPull*10+1);
      }
    }
    //
    arrY.Clear();
    arrZ.Clear();
    arrYP.Clear();
    arrZP.Clear();
    //
  }
  //
}

//_______________________________________________________
void PostProcessTRD(HistoManager* hm, HistoManager* hmProc)
{
  // postprocess histos
  //
  TObjArray arrY,arrZ,arrYP,arrZP;
  TH2* h=0;
  TH1* h1=0;
  // 
  for (int ilr=0;ilr<6;ilr++) {
    //
    int hoffs = kHOffsTRD + ilr*10;
    int hoffsP = kHOffsTRD + ilr*100;
    h = hm->GetHisto2F(hoffs + kDY);
    if (h) arrY.Add(h);
    h = hm->GetHisto2F(hoffs + kDZ);
    if (h) arrZ.Add(h);
    h = hm->GetHisto2F(hoffs + kDYPull);
    if (h) arrYP.Add(h);
    h = hm->GetHisto2F(hoffs + kDZPull);
    if (h) arrZP.Add(h);
    //
    if (FitProfile(&arrY)) {
      h1 = (TH1*)arrY.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDY*10+0);
      }
      h1 = (TH1*)arrY.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDY*10+1);
      }
    }
    //
    if (FitProfile(&arrZ)) {
      h1 = (TH1*)arrZ.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZ*10+0);
      }
      h1 = (TH1*)arrZ.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZ*10+1);
      }
    }
    //
    if (FitProfile(&arrYP)) {
      h1 = (TH1*)arrYP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDYPull*10+0);
      }
      h1 = (TH1*)arrYP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDYPull*10+1);
      }
    }
    //
    if (FitProfile(&arrZP)) {
      h1 = (TH1*)arrZP.RemoveAt(0);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZPull*10+0);
      }
      h1 = (TH1*)arrZP.RemoveAt(1);
      if (h1) {
	hmProc->AddHisto(h1, hoffsP + kDZPull*10+1);
      }
    }
    //
    arrY.Clear();
    arrZ.Clear();
    arrYP.Clear();
    arrZP.Clear();
    //
  }
  //  
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
    h->FitSlicesY(gs,0,-1,0,kFitOptGS,&hfarr);
    TH1* hmean = (TH1*)hfarr[1];
    TH1* hsig  = (TH1*)hfarr[2];
    if (!hmean || !hsig) {hfarr.Delete(); continue;}
    hmean->Fit("pol1",kFitOptP1);
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
  h->FitSlicesY(gs,0,-1,0,kFitOptGS,&hfarr);
  hfarr.SetOwner(kFALSE);
  TH1* hmean = (TH1*)hfarr[1];
  TH1* hsig  = (TH1*)hfarr[2];
  if (!hmean || !hsig) {hfarr.Delete(); return kFALSE;}
  delete gs;
  histos->Clear();
  hmean->SetName(Form("%s_%s",h->GetName(),"Mean"));
  hmean->SetTitle(Form("%s %s",h->GetTitle(),"Mean"));
  hsig->SetName(Form("%s_%s",h->GetName(),"Sigma"));
  hsig->SetTitle(Form("%s %s",h->GetTitle(),"Sigma"));
  histos->AddAtAndExpand(hmean,0);
  histos->AddAtAndExpand(hsig,1);
  hfarr.Clear();
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
  if (psnm1.IsNull()) psnm1 = Form("algRepK");
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
  double rangeY = 50e-4;
  double rangeZ = 50e-4;  
  double rangeYSig = 150e-4;
  double rangeZSig = 250e-4;  
  repCanv->Clear();
  repCanv->Divide(2,2);
  //
  repCanv->cd(1);
  DrawHistos(hmans, kHOffsVTX + kDY*10+0,rangeY);
  repCanv->cd(2);
  DrawHistos(hmans, kHOffsVTX + kDY*10+1,rangeYSig,0);
  //
  repCanv->cd(3);
  DrawHistos(hmans, kHOffsVTX + kDYPull*10+0);
  repCanv->cd(4);
  DrawHistos(hmans, kHOffsVTX + kDYPull*10+1);
  //
  repCanv->cd();
  repCanv->Print(psnm); // Y
  //
  repCanv->Clear();
  repCanv->Divide(2,2);
  //  
  repCanv->cd(1);
  DrawHistos(hmans, kHOffsVTX + kDZ*10+0,rangeZ);
  repCanv->cd(2);
  DrawHistos(hmans, kHOffsVTX + kDZ*10+1,rangeZSig);
  //
  repCanv->cd(3);
  DrawHistos(hmans, kHOffsVTX + kDZPull*10+0);
  repCanv->cd(4);
  DrawHistos(hmans, kHOffsVTX + kDZPull*10+1);
  //
  repCanv->cd();
  repCanv->Print(psnm); // Z
  //
}

//____________________________________________
void DrawReportITS(TObjArray* hmans, const char* psnm)
{
  int icn;
  float drRangeY[2][6] = {
    { 50e-4 ,300e-4, 300e-4,300e-4,400e-4, 400e-4},
    { 0.03, 0.01,   0.1, 0.1 ,0.01, 0.01}
  };
  float drRangeZ[2][6] = {
    {100e-4 ,300e-4, 500e-4, 500e-4, 1000e-4, 1000e-4},
    {500e-4, 500e-4, 0.2 , 0.2,  0.2, 0.2}
  };
  
  for (int ilr=0;ilr<6;ilr++) {
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDY*10 + 0, drRangeY[0][ilr]); // meanY
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDY*10 + 1, drRangeY[1][ilr],0); // sigma
      //
      repCanv->Print(psnm);
    }
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDYPull*10 + 0, 1); // pullY
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDYPull*10 + 1, 2,0); // pullY sigma
      //
      repCanv->Print(psnm);
    }

    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDZ*10 + 0, drRangeZ[0][ilr]); // meanZ
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDZ*10 + 1, drRangeZ[1][ilr],0); // sigma 
      //
      repCanv->Print(psnm);
    }
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDZPull*10 + 0, 1); // pullZ
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsITS + ilr*100 + kDZPull*10 + 1, 2,0); // sigma 
      //
      repCanv->Print(psnm);
    }
  }
  //
}

//____________________________________________
void DrawReportTRD(TObjArray* hmans, const char* psnm)
{
  float drRangeY[2][6] = {
    { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
    { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4}
  };
  float drRangeZ[2][6] = {
    { 2. , 2., 2., 2., 2., 2.},
    { 5. , 5. , 5., 5., 5., 5.}
  };
  //
  int icn;
  for (int ilr=0;ilr<6;ilr++) {
    { // Y
      icn=0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDY*10 + 0,drRangeY[0][ilr]);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDY*10 + 1,drRangeY[1][ilr],0);   
      //
      repCanv->Print(psnm);
    }
    
    { // Y pull
      icn=0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDYPull*10 + 0, 1);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDYPull*10 + 1, 2,0);   
      //
      repCanv->Print(psnm);
    }
    
    
    { // Z
      icn=0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZ*10 + 0,drRangeZ[0][ilr]);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZ*10 + 1,drRangeZ[1][ilr],0);   
      //
      repCanv->Print(psnm);
    }
    
    { // Z pull
      icn=0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZPull*10 + 0, 1);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTRD + ilr*100 + kDZPull*10 + 1, 2,0);   
      //
      repCanv->Print(psnm);
    }
  }
  //
}

//____________________________________________
void DrawReportTOF(TObjArray* hmans, const char* psnm)
{
  float drRangeY[2] = {5., 10.};
  float drRangeZ[2] = {5., 10.};
  // 
  int icn=0;
  for (int isc=0;isc<18;isc++) {
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDY*10 + 0,drRangeY[0]);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDY*10 + 1,drRangeY[1],0);
      repCanv->Print(psnm);
    }
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDYPull*10 + 0, 1);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDYPull*10 + 1, 2,0);
      repCanv->Print(psnm);
    }
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDZ*10 + 0,drRangeZ[0]);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDZ*10 + 1,drRangeZ[1],0);
      repCanv->Print(psnm);
    }
    //
    {
      icn = 0;
      repCanv->Clear();
      repCanv->Divide(1,2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDZPull*10 + 0, 2);
      repCanv->cd(++icn);
      DrawHistos(hmans, kHOffsTOF + isc*100 + kDZPull*10 + 1, 2,0);
      repCanv->Print(psnm);
    }
    //
  }
  //
}

//_________________________________________________
void DrawHistos(TObjArray* hmans, int id, float range, float rangeM)
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
  float mnd = mn - 0.2*(mx-mn);
  float mxd = mx + 0.2*(mx-mn);
  if (range>0) {
    mxd = range; 
    mnd = rangeM>-999 ? rangeM : -range;
  }
  nhAcc = 0;
  gStyle->SetTitleW(0.9);
  for (int ih=0;ih<nhm;ih++) {
    HistoManager* hm = (HistoManager*)hmans->At(ih);
    if (!hm) continue;
    TH1* h = hm->GetHisto1F(id);
    if (!h || h->GetEntries()==0) continue;
    h->SetMinimum( mnd );
    h->SetMaximum( mxd );
    h->Draw(nhAcc==0 ? "":"same");
    nhAcc++;
  }
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();
}
