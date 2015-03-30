
{
  printf("root logon: Loading custom settings\n");
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1001111);
  gStyle->SetPalette(1);
  gStyle->SetTitleX(0.20);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.1);
  gStyle->SetTitleAlign(13);
  gStyle->SetTickLength(-.03,"x");
  gStyle->SetTickLength(-.02,"y");
  gStyle->SetTickLength(-.02,"z");
  gStyle->SetLabelOffset(0.03,"x");
  gStyle->SetLabelOffset(0.02,"y");
  gStyle->SetLabelOffset(0.02,"z");
  //
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  //
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetTitleOffset(1.35);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameFillColor(kWhite);
  //
  gROOT->ForceStyle();
  //gStyle->SetTitleColor(0);
  //gStyle->SetStatColor(0);
  //  TString cmd = gSystem->GetMakeSharedLib();      
  //  cmd.ReplaceAll("$Opt","-g $Opt");
  //  gSystem->SetMakeSharedLib(cmd.Data());
  //
  // LOAD NECESSARY ADDITIONAL LIBRARIES
  gSystem->AddIncludePath(" -I./ -I$ALICE_ROOT/include");
  //
  if (gClassTable->GetID("HistoManager")<0) {
    gROOT->ProcessLine(".L HistoManager.cxx+");
    if (gClassTable->GetID("HistoManager")<0) {
      printf("Failed to load HistoManager class\n");
    }
  }
  if (gClassTable->GetID("RStopwatch")<0) {
    gROOT->ProcessLine(".L RStopwatch.cxx+");
    if (gClassTable->GetID("RStopwatch")<0) {
      printf("Failed to load RStopwatch class\n");
    }
  }
  //
  //
  gROOT->ProcessLine(".L SaveCanvas.C+");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  //
  //  gROOT->ProcessLine(".x LoadLibs.C");
  if (gClassTable->GetID("AliAlgPoint")<0) {
    gROOT->ProcessLine(".x LoadAlgLibs.C");
  }
  //
}
