// this is a macro to setup the OCDB whose objects will be used
// as a reference for the alignment/calibration, i.e. coorections
// will be evaluated wrt these objects

void configRefOCDB(int run = 188503) 
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
