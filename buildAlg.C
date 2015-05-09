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
//
void ConfigAlign(AliAlgSteer* algSteer);
void ConfigITS(AliAlgSteer* algSteer);
void ConfigTPC(AliAlgSteer* algSteer);
void ConfigTRD(AliAlgSteer* algSteer);
void ConfigTOF(AliAlgSteer* algSteer);
//
AliAlgSteer * algSTEER = 0;

//void buildAlg(int evID=4062, int trID=0) // for cosm: data -> LHC15c_000218623_cosmics_15000218623020_10
//void buildAlg(int evID=6594, int trID=0) // for cosm: data -> LHC15c_000218623_cosmics_15000218623020_10
void buildAlg(int evID=4, int trID=2) // for beam: data -> LHC10b_000117220_vpass1_pass4_10000117220022_30
{
  LoadESD();
  LoadEvent(evID>=0 ? evID : 0);
  //
  int nEv = esdTree->GetEntries();
  //
  algSTEER = new AliAlgSteer();
  //
  ConfigAlign(algSTEER);
  //
  algSTEER->SetESDTree(esdTree);
  int evFirst=0,evLast = esdTree->GetEntries()-1;
  if (evID>0) {
    evFirst = evID;
    evLast  = evID;    
  }

  for (int iev=evFirst;iev<=evLast;iev++) {
    LoadEvent(iev);
    algSTEER->ProcessEvent(esdEv);
  }
  //
  algSTEER->Terminate();

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

//======================================================================
//======================================================================
//
// Configure alignment steering
//
//======================================================================
//======================================================================
//======================================================================
void ConfigAlign(AliAlgSteer* algSteer)
{
  //
  algSteer->AddDetector(AliAlgSteer::kITS);
  algSteer->AddDetector(AliAlgSteer::kTPC);
  algSteer->AddDetector(AliAlgSteer::kTRD);
  algSteer->AddDetector(AliAlgSteer::kTOF);
  algSteer->InitDetectors();
  //
  algSteer->GetDetectorByDetID(AliAlgSteer::kTPC)->Disable();
  //  algSteer->GetDetectorByDetID(AliAlgSteer::kTOF)->Disable();
  //  algSteer->GetDetectorByDetID(AliAlgSteer::kTRD)->Disable();

  ConfigITS(algSteer);
  ConfigTPC(algSteer);
  ConfigTRD(algSteer);
  ConfigTOF(algSteer);
  //
  algSteer->SetCosmicSelStrict();
  algSteer->SetMinDetAcc(1);
  algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kMPRec | AliAlgSteer::kContR);
  //
  algSteer->InitDOFs();   
  //  
  //  algSteer->SetMilleTXT(1);

}

//======================================================================
void ConfigITS(AliAlgSteer* algSteer)
{
  //
  AliAlgDetITS* det = (AliAlgDetITS*)algSteer->GetDetectorByDetID(AliAlgSteer::kITS);
  if (!det||det->IsDisabled()) return;
  det->SetObligatory(kTRUE);
  det->SetUseErrorParam(kTRUE);
  det->SetTrackFlagSel(AliESDtrack::kITSin);
  det->SetNPointsSel(2);
  //
  det->SetAddErrorLr(0,30e-4,200e-4);
  det->SetAddErrorLr(1,30e-4,200e-4);
  det->SetAddErrorLr(2,500e-4,80e-4);
  det->SetAddErrorLr(3,500e-4,80e-4);
  det->SetAddErrorLr(4,50e-4,500e-4);
  det->SetAddErrorLr(5,50e-4,500e-4);   
}

//======================================================================
void ConfigTPC(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTPC* det = (AliAlgDetTPC*)algSteer->GetDetectorByDetID(AliAlgSteer::kTPC);
  if (!det||det->IsDisabled()) return;
  det->SetObligatory(kTRUE);
  det->SetTrackFlagSel(AliESDtrack::kTPCin);
  det->SetNPointsSel(50);
  det->SetAddError(3,10.); // HUGE errors
}

//======================================================================
void ConfigTRD(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTRD* det = (AliAlgDetTRD*)algSteer->GetDetectorByDetID(AliAlgSteer::kTRD);
  if (!det||det->IsDisabled()) return;
  det->SetObligatory(kFALSE);
  det->SetTrackFlagSel(AliESDtrack::kTRDout);
  det->SetNPointsSel(2);
}

//======================================================================
void ConfigTOF(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTOF* det = (AliAlgDetTOF*)algSteer->GetDetectorByDetID(AliAlgSteer::kTOF);
  if (!det||det->IsDisabled()) return;
  det->SetObligatory(kTRUE);
  det->SetTrackFlagSel(AliESDtrack::kTOFout);
  det->SetNPointsSel(1);
  //
}
