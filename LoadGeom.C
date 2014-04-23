void LoadGeom(const char* cdb="local:///data/CDBMirror/alice/data/2011/OCDB",
	      Int_t run=156000,
	      const char* misalignMod=0)
{
  //
  man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdb);
  man->SetRun(run);
  AliGeomManager::LoadGeometry();
  TString algMod = misalignMod;
  //
  if (!algMod.IsNull()) AliGeomManager::ApplyAlignObjsFromCDB(algMod.Data());
  //
}
