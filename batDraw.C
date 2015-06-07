{
  hv0bm = new HistoManager("hpostv0BM","halgV0_trdNew.root",1);
  hv0bp = new HistoManager("hpostv0BP","halgV0_trdNew.root",1);
  hv0b0 = new HistoManager("hpostv0B0","halgV0_trdNew.root",1);
  hv0bm->SetColor(kRed);
  hv0bp->SetColor(kBlue);
  hv0b0->SetColor(kGreen+2);
  hv0bm->SetMarkerStyle(24);
  hv0bp->SetMarkerStyle(25);
  hv0b0->SetMarkerStyle(26);
  hsnb0 = new HistoManager("hpostsnB0","halgSN.root",1);
  hsnbp = new HistoManager("hpostsnBP","halgSN.root",1);
  hsnbm = new HistoManager("hpostsnBM","halgSN.root",1);
  hsnb0->SetColor(kGreen+2);
  hsnbp->SetColor(kBlue);
  hsnbm->SetColor(kRed);
  hsnbm->SetMarkerStyle(20);
  hsnbp->SetMarkerStyle(21);
  hsnb0->SetMarkerStyle(22);
  //
  TObjArray harr;
  harr.Add(hv0b0);
  harr.Add(hv0bp);
  harr.Add(hv0bm);
  harr.Add(hsnb0);
  harr.Add(hsnbp);
  harr.Add(hsnbm);

  //
  for (int i=harr.GetEntriesFast();i--;) {
    HistoManager* hm = (HistoManager*)harr[i];
    hm->SetMarkerSize(0.4);
  }
  gROOT->ProcessLine(".L ProcRes.C++g");
  DrawReport(&harr,"algRep");
  //
}
