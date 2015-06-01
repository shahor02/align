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
void BookHistos();
void BookHistosITS(HistoManager* hm);

void ProcessRes(Int_t ip);
void FillITS(Int_t ip);

AliAlgRes* res=0;
TChain* ch = 0;
HistoManager* hman=0;

//---------------------------------------------------
enum {kHOffsITS=100000};
enum {kDYSnp,kDZSnp};
const int kNLrITS = 6;
const int kNBinsResITS = 100;
const int kNBinsSnpITS = 10;
const double kSnpMaxITS = 0.6;
const double kMaxDYITS[kNLrITS] = {0.01,0.01,0.10,0.10,0.01,0.01};
const double kMaxDZITS[kNLrITS] = {0.10,0.10,0.10,0.10,0.20,0.10};
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
  if (lrID>=AliGeomManager::kSPD1 && lrID>=AliGeomManager::kSSD2) FillITS(ip);
  //
}

//__________________________________________
void FillITS(int ip)
{
  int volID = res->GetVolID(ip);
  int modID = 0;
  int lrID = AliGeomManager::VolUIDToLayer(volID,modID);
  int offs = kHOffsITS + modID*10;
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
  hman = new HistoManager();
  BookHistosITS(hman);
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
      hm->AddHisto(h2, hoffs+cntM*10+kDYSnp);
      //
      hnm = Form("ITSDZvsSnp_L%d_M%d",ilr,imod);
      htl = Form("ITS #DeltaZ vs Snp Lr%d Mod%d",ilr,imod);
      h2 = new TH2F(hnm.Data(),htl.Data(),
		    kNBinsSnpITS,-kSnpMaxITS,kSnpMaxITS,
		    kNBinsResITS,-kMaxDZITS[ilr],kMaxDZITS[ilr]);
      h2->SetXTitle("snp");
      h2->SetYTitle("#DeltaZ");
      //
      hm->AddHisto(h2, hoffs+cntM*10+kDZSnp);
      //
      cntM++;
    }
  }
}
