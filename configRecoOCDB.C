// this is a macro to setup the OCDB whose objects will be used
// to undo the alignment/calibration applied during reconstruction

void configRecoOCDB(int run = 117220) 
{
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->UnsetDefaultStorage();
  man->UnsetSnapshotMode();
  //
  if (gSystem->AccessPathName("data/OCDB.root", kFileExists)==0) {
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(run);
    man->SetSnapshotMode("data/OCDB.root");
  }
  else {
    AliFatalGeneralF("","Failed to setup Reco-Time OCDB for run %d",run);
  }
  //
}
