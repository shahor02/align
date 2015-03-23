void ExtDetVols(const char* det="ITS",
		const char* ocdb="local:///home/shahoian/alicesw/aliroot/master/src/OCDB")
{
  AliGeomManager::LoadGeometry("geometry.root");
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  AliCDBEntry* entry = man->Get(Form("%s/Align/Data",det),0);
  TObjArray* arr = (TObjArray*)entry->GetObject();
  int nobj = arr->GetEntriesFast();
  for (int i=0;i<nobj;i++) {
    AliAlignObjParams* par = (AliAlignObjParams*)arr->At(i);
    int vID = par->GetVolUID();
    TString snam = par->GetSymName();
    //
    if (!vID) { // passive volume
      vid = AliITSAlignMilleModule::GetVolumeIDFromSymname(snam.Data());
      printf("vol = new AliAlgVol(\"%s\",%d);\n",snam.Data(),vID);
      printf("AddVolume(vol);\n\\\\\n");
    }
    else {
      printf("sens = new AliAlgSens(\"%s\",%d);\n",snam.Data(),vID);
      printf("AddSensor(sens);\n");
    }
  }
}
