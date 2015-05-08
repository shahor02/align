#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
//
#include "AliAlgSteer.h"
#include "AliAlgDet.h"
#include "AliAlgDetITS.h"
#include "AliAlgDetTPC.h"
#include "AliAlgDetTRD.h"
#include "AliAlgDetTOF.h"
#endif

void alignConf(AliAlgSteer* algSteer);
void ConfigITS(AliAlgSteer* algSteer);
void ConfigTPC(AliAlgSteer* algSteer);
void ConfigTRD(AliAlgSteer* algSteer);
void ConfigTOF(AliAlgSteer* algSteer);
//

void alignConf(AliAlgSteer* algSteer)
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
  algSteer->SetMinDetAcc(2);
  algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kMPRec | AliAlgSteer::kContR);
  //
  algSteer->InitDOFs();   
  //  
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
