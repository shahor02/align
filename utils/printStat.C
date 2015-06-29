#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TAxis.h"
#endif

void printStat(const char* fl)
{
  printf("Stat for %s\n",fl);
  TFile* ff = TFile::Open(fl);
  TList* lst = (TList*)ff->Get("clist");
  if (!lst) {printf("no clist\n");return;}
  TH1* hstat = (TH1*)lst->FindObject("stat");
  if (!hstat) {printf("no hstat\n");return;}
  //
  TAxis* ax = hstat->GetXaxis();
  for (int ib=1;ib<ax->GetNbins();ib++) {
    double val = hstat->GetBinContent(ib);
    if (val) printf("%-20s\t%9d\n",ax->GetBinLabel(ib),int(val));
  }
  ff->Close();
  delete ff;
  return;
}
