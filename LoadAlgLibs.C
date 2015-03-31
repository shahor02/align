Bool_t LoadAlgLibs()
{
  printf("LoadAlgLibs >>\n");
  //
  gROOT->ProcessLine(".L AliAlgAux.cxx+g");
  gROOT->ProcessLine(".L AliAlgPoint.cxx+g");
  if (gClassTable->GetID("AliAlgPoint")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgTrack.cxx+g");
  if (gClassTable->GetID("AliAlgTrack")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgVol.cxx+g");  
  if (gClassTable->GetID("AliAlgVol")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSens.cxx+g");  
  if (gClassTable->GetID("AliAlgSens")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDet.cxx+g");
  if (gClassTable->GetID("AliAlgDet")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSensITS.cxx+g");
  if (gClassTable->GetID("AliAlgSensITS")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDetITS.cxx+g");
  if (gClassTable->GetID("AliAlgDetITS")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSensTPC.cxx+g");
  if (gClassTable->GetID("AliAlgSensTPC")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDetTPC.cxx+g");
  if (gClassTable->GetID("AliAlgDetTPC")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSensTRD.cxx+g");
  if (gClassTable->GetID("AliAlgSensTRD")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDetTRD.cxx+g");
  if (gClassTable->GetID("AliAlgDetTRD")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSensTOF.cxx+g");
  if (gClassTable->GetID("AliAlgSensTOF")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDetTOF.cxx+g");
  if (gClassTable->GetID("AliAlgDetTOF")<0) return kFALSE;
   //
  gROOT->ProcessLine(".L AliAlgSteer.cxx+g");
  if (gClassTable->GetID("AliAlgSteer")<0) return kFALSE;
  //
  printf("LoadAlgLibs <<\n");
  return kTRUE;
}
