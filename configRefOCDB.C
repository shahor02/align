// this is a macro to setup the OCDB whose objects will be used
// as a reference for the alignment/calibration, i.e. coorections
// will be evaluated wrt these objects

const char* specOCDB = "/alice/cern.ch/user/s/shahoian/tstOCDB/outOCDB_LHC12tst0";

void configRefOCDB(int run = 188503) 
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
  if (specOCDB) {
    man->SetSpecificStorage("ITS/Align/Data",Form("alien://folder=%s",specOCDB));
    man->SetSpecificStorage("TRD/Align/Data",Form("alien://folder=%s",specOCDB));
    man->SetSpecificStorage("TOF/Align/Data",Form("alien://folder=%s",specOCDB));
  }
  //
  man->SetRun(run);
  //
}
