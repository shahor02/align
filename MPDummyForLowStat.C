#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TAxis.h>
#include "Mille.h"
#endif

// Generate fake MP records for DOFs below statistics threshold

const char* mpDummy = "mpDummy.mille";
const Float_t kDummyDer = 1e-6;
const Float_t kDummyRes = 0;
const Float_t kDummyErr = 999.;

void MPDummyForLowStat(const char* stfile, int thr=30, int nGen=40)
{
  // show degrees of freedom with low stat
  TFile* fl = TFile::Open(stfile);
  if (!fl) {printf("Failed to open %s\n",stfile); return;}
  TList* lst = (TList*)fl->Get("clist");
  if (!lst) {printf("No clist in %s\n",stfile); return;}
  TH1* hstdof = (TH1*)lst->FindObject("DOFstat");
  if (!hstdof) {printf("No DOFstat histo in %s\n",stfile); return;}
  //
  int ndof = hstdof->GetNbinsX();
  TAxis* xax = hstdof->GetXaxis();
  printf("%4s\t%-50s\t%s","cnt"," DOF ID_name","entries");
  Mille ml(mpDummy);
  int   labDum[1] = {0}, cnt=0;
  float locDum[1] = {0}, gloDum[1] = {kDummyDer};
  //
  for (int i=1;i<=ndof;i++) {
    if (hstdof->GetBinContent(i)>thr) continue;
    TString labS = xax->GetBinLabel(i);
    printf("%4d\t%-50s\t%7d\n",cnt++,labS.Data(),(int)hstdof->GetBinContent(i));
    int indL = labS.Index("-");
    if (indL>0) labS.Resize(indL);
    else {
      printf("Failed to extract label from %s\n",labS.Data());
      exit(1);
    }
    int lab = labS.Atoi();
    for (int j=nGen;j--;) {
      labDum[i] = lab;
      ml.mille(0, locDum, 1, gloDum, labDum, kDummyRes, kDummyErr);
    }
  }
  //
  lst->SetOwner();
  delete lst;
  fl->Close();
  delete fl;
}

