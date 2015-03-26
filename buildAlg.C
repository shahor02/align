Bool_t LoadAlgLibs();

void buildAlg()
{
  if (!LoadAlgLibs()) return;
  AliGeomManager::LoadGeometry("geometry.root");
  AliAlgDet* its = new AliAlgDetITS("its");
  its->Init();
}





Bool_t LoadAlgLibs()
{
  gROOT->ProcessLine(".L AliAlgAux.cxx+");
  gROOT->ProcessLine(".L AliAlgPoint.cxx+");
  if (gClassTable->GetID("AliAlgPoint")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgTrack.cxx+");
  if (gClassTable->GetID("AliAlgTrack")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgVol.cxx+");  
  if (gClassTable->GetID("AliAlgVol")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgSens.cxx+");  
  if (gClassTable->GetID("AliAlgSens")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDet.cxx+");
  if (gClassTable->GetID("AliAlgDet")<0) return kFALSE;
  //
  gROOT->ProcessLine(".L AliAlgDetITS.cxx+");
  if (gClassTable->GetID("AliAlgDetITS")<0) return kFALSE;
  //
  return kTRUE;
}
