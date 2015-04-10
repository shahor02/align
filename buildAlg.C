#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
//
#include "AliAlgSteer.h"
#include "AliAlgDet.h"
#endif


TFile* flIn=0;
AliESDEvent *esdEv=0;
AliESDfriend *esdFr=0;
TTree *esdTree = 0;

Int_t LoadESD(const char* path="data/AliESDs.root",Bool_t friends=kTRUE);
void PrintTrack(Int_t i);
void PrintTracks();
void ConnectFriends();


AliAlgSteer * algSteer = 0;

void buildAlg(int evID=4062, int trID=0, Bool_t cosm=kTRUE) // for cosm: data -> LHC15c_000218623_cosmics_15000218623020_10
//  void buildAlg(int evID=4, int trID=2, Bool_t cosm=kFALSE) // for beam: data -> LHC10b_000117220_vpass1_pass4_10000117220022_30
{
  LoadESD();
  LoadEvent(evID);
  esdEv->InitMagneticField();
  int run = esdEv->GetRunNumber();
  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("data/OCDB.root", kFileExists)==0) {        
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");    
    man->SetRun(run);
    man->SetSnapshotMode("data/OCDB.root");
  }
  else {
    man->SetRaw(1);
    man->SetRun(run);    
  }
  AliGeomManager::LoadGeometry("data/geometry.root");
  //
  algSteer = new AliAlgSteer();
  //
  algSteer->AddDetector(AliAlgSteer::kITS);
  algSteer->AddDetector(AliAlgSteer::kTPC);
  algSteer->AddDetector(AliAlgSteer::kTRD);
  algSteer->AddDetector(AliAlgSteer::kTOF);
  //
  TString detalg = "GRP";
  for (int id=0;id<AliAlgSteer::kNDetectors;id++) {
    if (algSteer->GetDetectorByDetID(id)) {
      detalg += " ";
      detalg += algSteer->GetDetNameByDetID(id);
    }
  }
  // apply alignment
  AliGeomManager::ApplyAlignObjsFromCDB(detalg.Data());
  //

  //  AliAlgDet* its = new AliAlgDetITS("its");
  //  its->Init();
  algSteer->Init();
  //
  //
  PrintTracks();
  //
  algSteer->SetESDEvent(esdEv);
  if (!cosm) {
    AliESDtrack* tr = esdEv->GetTrack(trID);
    algSteer->ProcessTrack(tr);
  }
  else {
    AliESDCosmicTrack* trC = esdEv->GetCosmicTrack(trID);
    algSteer->ProcessTrack(trC);
  }
  //
}

//-----------------------------------------------------------------
Int_t LoadESD(const char* path,Bool_t friends)
{
  flIn = TFile::Open(path);
  if (!flIn) {
    printf("Failed to open %s\n",path);
    return -1;
  }
  //
  esdTree = (TTree*) flIn->Get("esdTree");
  if (!esdTree) {
    printf("No ESDtree found in %s\n",path);
    return -1;
  }
  //
  printf("Loaded esdTree with %d entries\n",(int)esdTree->GetEntries());
  esdEv = new AliESDEvent();
  esdEv->ReadFromTree(esdTree);
  if (friends) ConnectFriends();
  //
  return 0;
}

Int_t LoadEvent(Int_t iev)
{
  if (!esdEv) return -1; 
  esdEv->Reset();
  esdTree->GetEntry(iev);
  if (esdFr) esdEv->SetESDfriend(esdFr);
  esdEv->ConnectTracks();
  return 0;
}


void PrintTracks()
{
  for (int i=0;i<esdEv->GetNumberOfTracks();i++) PrintTrack(i);
  printf("NCosm: %d\n",esdEv->GetNumberOfCosmicTracks());
}

void PrintTrack(Int_t i)
{
  AliESDtrack* tr = esdEv->GetTrack(i);
  if (!tr) return;
  AliESDfriendTrack* trf = tr->GetFriendTrack();
  if (!trf || !trf->GetTrackPointArray()) trf=0;
  printf("%3d: its:%d tpc:%d trd:%d tof:%d | P:%6.2f Fr:%s\n",
	 i,
	 tr->IsOn(AliESDtrack::kITSrefit) ? tr->GetNcls(0) : -tr->GetNcls(0),
	 tr->IsOn(AliESDtrack::kTPCrefit) ? tr->GetNcls(1) : -tr->GetNcls(1),
	 tr->IsOn(AliESDtrack::kTRDout)   ? tr->GetNcls(2) : -tr->GetNcls(2),
	 tr->IsOn(AliESDtrack::kTOFout),
	 tr->GetP(),trf ? "ON":"OFF");
}


void ConnectFriends()
{
  // Connect the friends tree as soon as available.
  //
  // Handle the friends first
  //
  if (!esdTree->FindBranch("ESDfriend.")) {
    // Try to add ESDfriend. branch as friend
      TString esdFriendTreeFName;
      esdFriendTreeFName = (esdTree->GetCurrentFile())->GetName();    
      TString basename = gSystem->BaseName(esdFriendTreeFName);
      Int_t index = basename.Index("#")+1;
      basename.Remove(index);
      basename += "AliESDfriends.root";
      TString dirname = gSystem->DirName(esdFriendTreeFName);
      dirname += "/";
      esdFriendTreeFName = dirname + basename;
      //
      TTree* cTree = esdTree->GetTree();
      if (!cTree) cTree = esdTree;      
      cTree->AddFriend("esdFriendTree", esdFriendTreeFName.Data());
      cTree->SetBranchStatus("ESDfriend.", 1);
      esdFr = (AliESDfriend*)(esdEv->FindListObject("AliESDfriend"));
      if (esdFr) cTree->SetBranchAddress("ESDfriend.", &esdFr);
  }
}
