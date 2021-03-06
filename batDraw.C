{
  hv0bm = new HistoManager("hpostv0BM","halgV0_newTRDnewRC.root",1,"v0");
  hv0bp = new HistoManager("hpostv0BP","halgV0_newTRDnewRC.root",1,"v0");
  hv0b0 = new HistoManager("hpostv0B0","halgV0_newTRDnewRC.root",1,"v0");
  hv0bm->SetColor(kBlue);
  hv0bp->SetColor(kBlue);
  hv0b0->SetColor(kBlue);
  hv0bm->SetMarkerStyle(24);
  hv0bp->SetMarkerStyle(24);
  hv0b0->SetMarkerStyle(24);
  hpfb0 = new HistoManager("hpostpfB0","halgPF3_newTRDnewRC.root",1);
  hpfbp = new HistoManager("hpostpfBP","halgPF3_newTRDnewRC.root",1);
  hpfbm = new HistoManager("hpostpfBM","halgPF3_newTRDnewRC.root",1);
  hpfb0->SetColor(kRed);
  hpfbp->SetColor(kRed);
  hpfbm->SetColor(kRed);
  hpfbm->SetMarkerStyle(20);
  hpfbp->SetMarkerStyle(20);
  hpfb0->SetMarkerStyle(20);
  //
  gROOT->ProcessLine(".L ProcRes.C++g");

  TObjArray harr;
  harr.Add(hv0b0);
  harr.Add(hpfb0);
  for (int i=harr.GetEntriesFast();i--;) {
    HistoManager* hm = (HistoManager*)harr[i];
    hm->SetMarkerSize(0.6);
  }
  DrawReport(&harr,"algRepB0_newTRDnewRCPF3");
  //
  
  harr.Clear();
  harr.Add(hv0bp);
  harr.Add(hpfbp);
  for (int i=harr.GetEntriesFast();i--;) {
    HistoManager* hm = (HistoManager*)harr[i];
    hm->SetMarkerSize(0.6);
  }
  DrawReport(&harr,"algRepBP_newTRDnewRCPF3");
  //
  harr.Clear();
  harr.Add(hv0bm);
  harr.Add(hpfbm);
  for (int i=harr.GetEntriesFast();i--;) {
    HistoManager* hm = (HistoManager*)harr[i];
    hm->SetMarkerSize(0.6);
  }
  DrawReport(&harr,"algRepBM_newTRDnewRCPF3");
  //
  //
}
