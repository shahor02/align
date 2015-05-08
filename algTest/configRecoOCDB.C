// this is a macro to setup the OCDB whose objects will be used
// to undo the alignment/calibration applied during reconstruction

void configRecoOCDB(int run = 188503) 
{
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->ClearCache();
  if (man->IsDefaultStorageSet()) man->UnsetDefaultStorage();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetSnapshotMode("OCDB.root");
  }
  //
}
