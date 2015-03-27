/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAlgSteer.h"
#include "AliLog.h"
#include "AliAlgDet.h"
#include "AliAlgDetITS.h"

const char* AliAlgSteer::fgkDetectorName[AliAlgSteer::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID" };
const int AliAlgSteer::fgkSkipLayers[AliAlgSteer::kNLrSkip] = {AliGeomManager::kPHOS1,AliGeomManager::kPHOS2,
							       AliGeomManager::kMUON,AliGeomManager::kEMCAL};

ClassImp(AliAlgSteer)

//________________________________________________________________
AliAlgSteer::AliAlgSteer()
:  fNDet(0)
  ,fRunNumber(-1)
  ,fAlgTrack(0)
{
  // def c-tor
  for (int i=kNDetectors;i--;) {
    fDetectors[i] = 0;
    fDetPos[i] = -1;
  }
}

//________________________________________________________________
AliAlgSteer::~AliAlgSteer()
{
  // d-tor
  delete fAlgTrack;
  for (int i=0;i<fNDet;i++) delete fDetectors[i];
}

//________________________________________________________________
void AliAlgSteer::Init()
{
  // init all detectors
  //
  static Bool_t done = kFALSE;
  if (done) return;
  done = kTRUE;
  //
  fAlgTrack = new AliAlgTrack();
  //
  for (int i=0;i<fNDet;i++) fDetectors[i]->Init();
}

//________________________________________________________________
  void AliAlgSteer::AddDetector(UInt_t id, AliAlgDet* det)
{
  // add detector participating in the alignment, optionally constructed externally
  if (id>=kNDetectors)  AliFatalF("Detector type %d exceeds allowed range %d:%d",
				  id,0,kNDetectors-1);
  //
  if (fDetPos[id]!=-1) AliFatalF("Detector %d was already added",id);
  if (!det) {
    switch(id) {
    case kITS: det = new AliAlgDetITS("ITS"); break;
      //  case kTPC: det = new AliAlgDetTPC(); break;
      //  case kTRD: det = new AliAlgDetTRD(); break;
      //  case kTOF: det = new AliAlgDetTOF(); break;
    default: AliErrorF("%d not implemented yet",id); break;
    };
  }
  //
  fDetectors[fNDet] = det;
  fDetPos[id] = fNDet;
  fNDet++;
  //
}

//_________________________________________________________
void AliAlgSteer::AddDetector(AliAlgDet* det)
{
  // add detector constructed externally to alignment framework
  AddDetector(det->GetType(), det);
}


//_________________________________________________________
Bool_t AliAlgSteer::AcceptTrack(const AliESDtrack* /*esdTr*/) const
{
  // decide if the track should be processed
  return kTRUE;
  //
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessTrack(const AliESDtrack* esdTr)
{
  // process single track
  //
  int nPnt = 0;
  const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
  if (!trF) return kFALSE;
  const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
  if (!trPoints || (nPnt=trPoints->GetNPoints())<1) return kFALSE;
  //
  if (!AcceptTrack(esdTr)) return kFALSE;
  //
  fAlgTrack->Clear();
  //
  // process the track points for each detector, fill needed points (tracking frame) in the fAlgTrack
  AliAlgDet* det = 0;
  for (int idet=0;idet<kNDetectors;idet++) {
    if (!(det=GetDetectorByType(idet))) continue;
    if (!det->PresentInTrack(esdTr) ) continue;
    //
    GetDetector(idet)->ProcessPoints(esdTr, fAlgTrack);
  }
  //
  return kTRUE;
}

//_________________________________________________________
void AliAlgSteer::SetRunNumber(Int_t run)
{
  if (run==fRunNumber) return;  // nothing to do
  //
  AcknowledgeNewRun(run); 
}

//_________________________________________________________
void AliAlgSteer::AcknowledgeNewRun(Int_t run)
{
  // load needed info for new run
  fRunNumber = run;
  for (int idet=0;idet<fNDet;idet++) GetDetector(idet)->AcknowledgeNewRun(run);
  //
}

//_________________________________________________________
AliAlgDet* AliAlgSteer::GetDetectorByVolID(Int_t vid) const
{
  // get detector by sensor volid
  for (int i=fNDet;i--;) if (fDetectors[i]->SensorOfDetector(vid)) return fDetectors[i];
  return 0;
}

//____________________________________________
void AliAlgSteer::Print(const Option_t *opt) const
{
  // print info
  for (int idt=0;idt<kNDetectors;idt++) {
    AliAlgDet* det = GetDetectorByType(idt);
    if (det) det->Print(opt);
    else printf("Detector:%5s is not defined\n",fgkDetectorName[idt]);
  }
  
}
