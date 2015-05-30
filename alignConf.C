#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
//
#include "AliAlgSteer.h"
#include "AliAlgDet.h"
#include "AliAlgSens.h"
#include "AliAlgDetITS.h"
#include "AliAlgDetTPC.h"
#include "AliAlgDetTRD.h"
#include "AliAlgDetTOF.h"
#include "AliAlgVtx.h"
#endif

void alignConf(AliAlgSteer* algSteer);
void ConfigITS(AliAlgSteer* algSteer);
void ConfigTPC(AliAlgSteer* algSteer);
void ConfigTRD(AliAlgSteer* algSteer);
void ConfigTOF(AliAlgSteer* algSteer);
void ConfigVTX(AliAlgSteer* algSteer);


void alignConf(AliAlgSteer* algSteer)
{
  //
  algSteer->SetRefOCDBConfigMacro("configRefOCDB.C");
  //  algSteer->SetRecoOCDBConfigMacro("configRecoOCDB.C");  
  algSteer->SetRecoOCDBConfigMacro(""); // Use ESD info
  //
  algSteer->AddDetector(AliAlgSteer::kITS);
  algSteer->AddDetector(AliAlgSteer::kTPC);
  algSteer->AddDetector(AliAlgSteer::kTRD);
  algSteer->AddDetector(AliAlgSteer::kTOF);
  algSteer->InitDetectors();
  //
  //algSteer->GetDetectorByDetID(AliAlgSteer::kITS)->SetDisabled();
  algSteer->GetDetectorByDetID(AliAlgSteer::kTPC)->SetDisabled();
  //algSteer->GetDetectorByDetID(AliAlgSteer::kTOF)->SetDisabled();
  //algSteer->GetDetectorByDetID(AliAlgSteer::kTRD)->SetDisabled();
  //
  /*
  for (int i=algSteer->GetNDetectors();i--;) { //!!!!RS
    AliAlgDet* det = algSteer->GetDetector(i);
    for (int iv=det->GetNVolumes();iv--;) 
      //  det->GetVolume(iv)->SetFreeDOFPattern( (0x1<<AliAlgVol::kDOFTY)|(0x1<<AliAlgVol::kDOFTZ)|(0x1<<AliAlgVol::kDOFTX) );
      //      det->GetVolume(iv)->SetFreeDOFPattern( (0x1<<AliAlgVol::kDOFTX) );
      det->FixNonSensors();
  }
  */

  ConfigITS(algSteer);
  ConfigTPC(algSteer);
  ConfigTRD(algSteer);
  ConfigTOF(algSteer);
  //
  algSteer->SetVtxMinCont(5);   // accept events with min number of vertexTracks contributors
  algSteer->SetVtxMinContVC(10); // use for vertex constraint only those with min number of contributors
  algSteer->SetMaxDCAforVC(0.1,0.6); // dcaR/Z primary selection to allow vertex constraint
  algSteer->SetMaxChi2forVC(10);     // track-vertex chi2 primary selection to allow vertex constraint

  algSteer->SetCosmicSelStrict(kTRUE); // apply track selection to each leg separately
  //
  algSteer->SetMinDetAccColl(2);       // min number of detectors in track
  algSteer->SetMinDetAccCosm(2);
  //
  algSteer->SetMinPointsColl(6,6);     // min number of points per track Boff/Bon
  algSteer->SetMinPointsCosm(4,4);
  //
  algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kMPRec | AliAlgSteer::kContR);
  //algSteer->SetMPOutType(AliAlgSteer::kMille | AliAlgSteer::kContR);
  //  
  //  algSteer->SetMilleTXT(1);

  AliAlgSens* vtxS = (AliAlgSens*)algSteer->GetVertexSensor();
  printf("VTS: %p\n",vtxS);
  vtxS->SetFreeDOFPattern(0); // fix vertex movements

  algSteer->InitDOFs();
  //  
}

//======================================================================
void ConfigVTX(AliAlgSteer* algSteer)
{
  //
  AliAlgVtx *vtx = algSteer->GetVertexSensor();
  vtx->SetFreeDOFPattern(0);
}

//======================================================================
void ConfigITS(AliAlgSteer* algSteer)
{
  //
  const double kCondSig[AliAlgVol::kNDOFGeom] = {0.1,0.1,0.2,1,1,1}; // precondition sigmas
  //
  AliAlgDetITS* det = (AliAlgDetITS*)algSteer->GetDetectorByDetID(AliAlgSteer::kITS);
  if (!det||det->IsDisabled()) return;
  det->SetUseErrorParam(kTRUE);
  //
  det->SetObligatoryColl(kTRUE);
  det->SetObligatoryCosm(kTRUE);
  //
  det->SetTrackFlagSelCosm(AliESDtrack::kITSin);
  det->SetTrackFlagSelColl(AliESDtrack::kITSrefit | AliESDtrack::kTPCrefit);
  //
  det->SetNPointsSelCosm(2);
  det->SetNPointsSelColl(3);
  //
  det->SetITSSelPatternColl(AliAlgDetITS::kSPDAny);
  det->SetITSSelPatternCosm(AliAlgDetITS::kSPDNoSel);
  //
  // precondition
  for (int iv=det->GetNVolumes();iv--;) {
    AliAlgVol* vol = det->GetVolume(iv);
    for (int idf=AliAlgVol::kNDOFGeom;idf--;) {
      if (TMath::Abs(vol->GetParErr(idf))>1e-6) continue; // there is already condition
      vol->SetParErr(idf, kCondSig[idf] );    // set condition
    }
  }  
  //
  //  det->SetAddError(2,2);
  /*
  det->SetAddErrorLr(0,20e-4,100e-4);
  det->SetAddErrorLr(1,20e-4,100e-4);
  det->SetAddErrorLr(2,2000e-4,20e-4);
  det->SetAddErrorLr(3,2000e-4,20e-4);
  det->SetAddErrorLr(4,20e-4,500e-4);
  det->SetAddErrorLr(5,20e-4,500e-4);   
  */
}

//======================================================================
void ConfigTPC(AliAlgSteer* algSteer)
{
  //
  AliAlgDetTPC* det = (AliAlgDetTPC*)algSteer->GetDetectorByDetID(AliAlgSteer::kTPC);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kFALSE);
  det->SetObligatoryCosm(kFALSE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);
  det->SetTrackFlagSelCosm(AliESDtrack::kTPCin);
  //
  det->SetNPointsSelColl(70);
  det->SetNPointsSelCosm(50);
  //
  //  det->SetAddError(3,10.); // HUGE errors
  det->SetAddError(0.1,0.1); // HUGE errors
}

//======================================================================
void ConfigTRD(AliAlgSteer* algSteer)
{
  //
  const double kCondSigSMD[AliAlgVol::kNDOFGeom] = {2,2,2,1,1,1}; // precondition sigmas
  const double kCondSigCHA[AliAlgVol::kNDOFGeom] = {0.5,0.5,5.0,1,1,1}; // precondition sigmas

  AliAlgDetTRD* det = (AliAlgDetTRD*)algSteer->GetDetectorByDetID(AliAlgSteer::kTRD);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kFALSE);
  det->SetObligatoryCosm(kFALSE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTRDout);
  det->SetTrackFlagSelCosm(AliESDtrack::kTRDout);
  //
  det->SetNPointsSelColl(2);
  det->SetNPointsSelCosm(2);
  //
  // precondition
  // precondition
  for (int idf=0;idf<AliAlgVol::kNDOFGeom;idf++) {
    det->SetDOFCondition(idf,kCondSigSMD[idf],0); // level 0 - supermodules
    det->SetDOFCondition(idf,kCondSigCHA[idf],1); // level 1 - chambers
  }
  //  det->SetAddError(0.1,0.1);
  //
}

//======================================================================
void ConfigTOF(AliAlgSteer* algSteer)
{
  //
  const double kCondSigSMD[AliAlgVol::kNDOFGeom] = {1,5.,5.,1,1,1}; // precondition sigmas
  const double kCondSigSTR[AliAlgVol::kNDOFGeom] = {0.1,0.1,0.1, -1,-1,-1}; // precondition sigmas
  //
  AliAlgDetTOF* det = (AliAlgDetTOF*)algSteer->GetDetectorByDetID(AliAlgSteer::kTOF);
  if (!det||det->IsDisabled()) return;
  //
  det->SetObligatoryColl(kTRUE);
  det->SetObligatoryCosm(kTRUE);
  //
  det->SetTrackFlagSelColl(AliESDtrack::kTOFout);
  det->SetTrackFlagSelCosm(AliESDtrack::kTOFout);
  //
  det->SetNPointsSelColl(1);
  det->SetNPointsSelCosm(1);
  //
  // precondition
  for (int idf=0;idf<AliAlgVol::kNDOFGeom;idf++) {
    det->SetDOFCondition(idf,kCondSigSMD[idf],0); // level 0 - supermodules
    det->SetDOFCondition(idf,kCondSigSTR[idf],1); // level 1 - strips
  }
  //
  for (int is=det->GetNSensors();is--;) det->GetSensor(is)->SetVarFrame(AliAlgVol::kLOC);
  //
  //  det->SetAddError(12,12);

  //
}
