// this is a macro to setup the OCDB whose objects will be used
// to undo the alignment/calibration applied during reconstruction

void configRecoOCDB(int run = 188503) 
{
  //
  AliAlgAux::CleanOCDB();
  AliCDBManager* man = AliCDBManager::Instance();
  //
  if (gSystem->AccessPathName("data/OCDB.root", kFileExists)==0) {
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetSnapshotMode("data/OCDB.root");
  }
  else {
    man->SetRaw(1);
  }
  //
  man->SetRun(run);

}
