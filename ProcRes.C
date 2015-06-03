#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TString.h>
#include <fstream>
#include "AliAlgRes.h"
#include "HistoManager.h"
#include "AliGeomManager.h"
#include "AliITSgeomTGeo.h"
#endif

TChain* LoadChain(const char* inpData, const char* chName="res");
void DoOffsets();
void BookHistos();
void BookHistosITS(HistoManager* hm);
void BookHistosTRD(HistoManager* hm);
void BookHistosTOF(HistoManager* hm);

void ProcessRes(Int_t ip);
void FillITS(Int_t ip);
void FillTRD(Int_t ip);
void FillTOF(Int_t ip);

AliAlgRes* res=0;
TChain* ch = 0;
HistoManager* hman=0;

//---------------------------------------------------
enum {kHOffsITS=100000,kHOffsTRD=300000,kHOffsTOF=400000};
enum {kDYSnp,kDZSnp};
enum {kNLrITS=6,kNLrTRD=6};
//
int kITSLrOffs[kNLrITS] = {0};
int kTRDLrOffs[kNLrTRD] = {0};
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
const double kMaxDYITS[kNLrITS] = {0.01,0.01,0.10,0.10,0.01,0.01};
const double kMaxDZITS[kNLrITS] = {0.10,0.10,0.10,0.10,0.20,0.10};
//
const double kMaxDYTRD[kNLrTRD] = {1.00,1.00,1.00,1.00,1.00,1.00};
const double kMaxDZTRD[kNLrTRD] = {10.0,10.0,10.0,10.0,10.0,10.0};
//
const double kMaxDYTOF = 10.;
const double kMaxDZTOF = 10.;
//---------------------------------------------------


void ProcRes(const char* inpData)
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
  
}

//__________________________________________
void ProcessRes(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID);
  if (lrID>=AliGeomManager::kSPD1 && lrID<=AliGeomManager::kSSD2) FillITS(ip);
  if (lrID>=AliGeomManager::kTRD1 && lrID<=AliGeomManager::kTRD6) FillTRD(ip);
  if (lrID==AliGeomManager::kTOF)                                 FillTOF(ip);
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
    cntM++;
  }
}

//___________________________________________
void DoOffsets()
{
  for (int i=0,n=0;i<kNLrITS;i++) {kITSLrOffs[i] = n; n += AliGeomManager::LayerSize(i+AliGeomManager::kSPD1);}
  for (int i=0,n=0;i<kNLrTRD;i++) {kTRDLrOffs[i] = n; n += AliGeomManager::LayerSize(i+AliGeomManager::kTRD1);}
}
