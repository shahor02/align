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
#include "AliAlgAux.h"
#include "AliAlgDet.h"
#include "AliAlgDetITS.h"
#include "AliAlgDetTPC.h"
#include "AliAlgDetTRD.h"
#include "AliAlgDetTOF.h"
#include "AliTrackerBase.h"
#include "AliESDCosmicTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include <TMath.h>

using namespace TMath;
using namespace AliAlgAux;

ClassImp(AliAlgSteer)


const Char_t* AliAlgSteer::fgkDetectorName[AliAlgSteer::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID" };
const Int_t   AliAlgSteer::fgkSkipLayers[AliAlgSteer::kNLrSkip] = {AliGeomManager::kPHOS1,AliGeomManager::kPHOS2,
								   AliGeomManager::kMUON,AliGeomManager::kEMCAL};

//________________________________________________________________
AliAlgSteer::AliAlgSteer()
  :fNDet(0)
  ,fRunNumber(-1)
  ,fAlgTrack(0)
  ,fRefPoint()
  ,fESDEvent()
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
  if (id>=kNDetectors)  AliFatalF("Detector typeID %d exceeds allowed range %d:%d",
				  id,0,kNDetectors-1);
  //
  if (fDetPos[id]!=-1) AliFatalF("Detector %d was already added",id);
  if (!det) {
    switch(id) {
    case kITS: det = new AliAlgDetITS(GetDetNameByDetID(kITS)); break;
    case kTPC: det = new AliAlgDetTPC(GetDetNameByDetID(kTPC)); break;
    case kTRD: det = new AliAlgDetTRD(GetDetNameByDetID(kTRD)); break;
    case kTOF: det = new AliAlgDetTOF(GetDetNameByDetID(kTOF)); break;
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
  AddDetector(det->GetDetID(), det);
}


//_________________________________________________________
Bool_t AliAlgSteer::AcceptTrack(const AliESDtrack* /*esdTr*/) const
{
  // decide if the track should be processed
  return kTRUE;
  //
}

//_________________________________________________________
Bool_t AliAlgSteer::AcceptTrack(const AliESDtrack* esdPairCosm[kNCosmLegs]) const
{
  // decide if the pair of tracks making cosmic track should be processed
  for (int i=kNCosmLegs;i--;) if (!AcceptTrack(esdPairCosm[i])) return kFALSE; // TODO
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
  ResetDetectors();
  fAlgTrack->Clear();
  //
  // process the track points for each detector, 
  // fill needed points (tracking frame) in the fAlgTrack
  fRefPoint.SetContainsMeasurement(kFALSE);
  fRefPoint.SetContainsMaterial(kFALSE);  
  fAlgTrack->AddPoint(&fRefPoint); // reference point which the track will refer to
  //
  AliAlgDet* det = 0;
  for (int idet=0;idet<kNDetectors;idet++) {
    if (!(det=GetDetectorByDetID(idet))) continue;
    if (!det->PresentInTrack(esdTr) ) continue;
    //
    det->ProcessPoints(esdTr, fAlgTrack);
  }
  //
  fAlgTrack->Set(esdTr->GetX(),esdTr->GetAlpha(),esdTr->GetParameter(),esdTr->GetCovariance());
  fAlgTrack->CopyFrom(esdTr);
  fAlgTrack->SetFieldON( Abs(AliTrackerBase::GetBz())>kAlmost0Field );
  fAlgTrack->SortPoints();
  //
  // at this stage the points are sorted from maxX to minX, the latter corresponding to
  // reference point (e.g. vertex) with X~0. The fAlgTrack->GetInnerPointID() points on it,
  // hence fAlgTrack->GetInnerPointID() is the 1st really measured point. We will set the 
  // alpha of the reference point to alpha of the barrel sector corresponding to this 
  // 1st measured point
  int pntMeas = fAlgTrack->GetInnerPointID()-1;
  if (pntMeas<0) { // this should not happen
    fAlgTrack->Print("p meas");
    AliFatal("AliAlgTrack->GetInnerPointID() cannot be 0");
  }
  fRefPoint.SetAlphaSens(Sector2Alpha(fAlgTrack->GetPoint(pntMeas)->GetAliceSector()));
  //
  if (!fAlgTrack->IniFit()) return kFALSE;
  if (!fAlgTrack->ProcessMaterials()) return kFALSE;
  fAlgTrack->DefineDOFs();
  //
  if (!fAlgTrack->CalcResiduals()) return kFALSE;
  if (!fAlgTrack->CalcResidDeriv()) return kFALSE;
  //
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessTrack(const AliESDCosmicTrack* cosmTr)
{
  // process single cosmic track
  //
  int nPnt = 0;
  const AliESDtrack*        trCosmE[kNCosmLegs]  = {0};
  //  const AliESDfriendTrack*  trCosmF[kNCosmLegs] = {0};
  //  const AliTrackPointArray* trCosmP[kNCosmLegs] = {0};
  //
  for (int leg=kNCosmLegs;leg--;) {
    const AliESDtrack* esdTr = 
      fESDEvent->GetTrack(leg==kCosmLow ? 
			  cosmTr->GetESDLowerTrackIndex():cosmTr->GetESDUpperTrackIndex());
    const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
    if (!trF) return kFALSE;
    const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
    if (!trPoints || (nPnt+=trPoints->GetNPoints())<1) return kFALSE;
    //
    trCosmE[leg] = esdTr;
    //    trCosmF[leg] = trF;
    //    trCosmP[leg] = trPoints;
  }
  //
  if (!AcceptTrack(trCosmE)) return kFALSE;
  //
  ResetDetectors();
  fAlgTrack->Clear();
  fAlgTrack->SetCosmic(kTRUE);
  //
  // process the track points for each detector, 
  // fill needed points (tracking frame) in the fAlgTrack
  fRefPoint.SetContainsMeasurement(kFALSE);
  fRefPoint.SetContainsMaterial(kFALSE);  
  fAlgTrack->AddPoint(&fRefPoint); // reference point which the track will refer to
  //
  AliAlgDet* det = 0;
  for (int leg=kNCosmLegs;leg--;) {
    for (int idet=0;idet<kNDetectors;idet++) {
      if (!(det=GetDetectorByDetID(idet))) continue;
      if (!det->PresentInTrack(trCosmE[leg]) ) continue;
      // upper leg points marked as the track goes in inverse direction
      det->ProcessPoints(trCosmE[leg],fAlgTrack, leg==kCosmUp);
    }
  }
  fAlgTrack->CopyFrom(cosmTr);
  fAlgTrack->SetFieldON( Abs(AliTrackerBase::GetBz())>kAlmost0Field );
  fAlgTrack->SortPoints();
   //
  // at this stage the points are sorted from maxX to minX, the latter corresponding to
  // reference point (e.g. vertex) with X~0. The fAlgTrack->GetInnerPointID() points on it,
  // hence fAlgTrack->GetInnerPointID() is the 1st really measured point. We will set the 
  // alpha of the reference point to alpha of the barrel sector corresponding to this 
  // 1st measured point
  int pntMeas = fAlgTrack->GetInnerPointID()-1;
  if (pntMeas<0) { // this should not happen
    fAlgTrack->Print("p meas");
    AliFatal("AliAlgTrack->GetInnerPointID() cannot be 0");
  }
  fRefPoint.SetAlphaSens(Sector2Alpha(fAlgTrack->GetPoint(pntMeas)->GetAliceSector()));
  // 
  if (!fAlgTrack->IniFit()) return kFALSE;
  //
  if (!fAlgTrack->ProcessMaterials()) return kFALSE;
  fAlgTrack->DefineDOFs();
  if (!fAlgTrack->CalcResiduals()) return kFALSE;
  if (!fAlgTrack->CalcResidDeriv()) return kFALSE;
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
    AliAlgDet* det = GetDetectorByDetID(idt);
    if (det) det->Print(opt);
    else printf("Detector:%5s is not defined\n",GetDetNameByDetID(idt));
  }
  //  
}

void AliAlgSteer::ResetDetectors()
{
  // reset detectors for next track
  for (int idet=fNDet;idet--;) {
    AliAlgDet* det = GetDetector(idet);
    det->ResetPool();   // reset used alignment points
  }
}

//____________________________________________
AliSymMatrix* AliAlgSteer::BuildMatrix(TVectorD &vec)
{
  int npnt = fAlgTrack->GetNPoints();
  int nlocpar = fAlgTrack->GetNLocPar();
  int nlocparETP = fAlgTrack->GetNLocExtPar(); // parameters of external track param
  //
  vec.ResizeTo(nlocpar);
  memset(vec.GetMatrixArray(),0,nlocpar*sizeof(double));
  AliSymMatrix* matp = new AliSymMatrix(nlocpar);
  AliSymMatrix& mat = *matp;
  //
  for (int ip=npnt;ip--;) {
    AliAlgPoint* pnt = fAlgTrack->GetPoint(ip);
    //
    if (pnt->ContainsMeasurement()) {
      //      pnt->Print("meas");
      for (int idim=2;idim--;) { // each point has 2 position residuals
	double  sigma2 = pnt->GetErrDiag(idim);              // residual error
	double  resid  = fAlgTrack->GetResidual(idim,ip);    // residual
	double* deriv  = fAlgTrack->GetDerivative(idim,ip);  // array of Dresidual/Dparams
	//
	double sg2inv = 1./sigma2;
	for (int parI=nlocpar;parI--;) { 
	  vec[parI] -= deriv[parI]*resid*sg2inv;
	  //  printf("%d %d %d %+e %+e %+e -> %+e\n",ip,idim,parI,sg2inv,deriv[parI],resid,vec[parI]);
	  //	  for (int parJ=nlocpar;parJ--;) {
	  for (int parJ=parI+1;parJ--;) {
	    mat(parI,parJ) += deriv[parI]*deriv[parJ]*sg2inv;	  
	  }
	}
      } // loop over 2 orthogonal measurements at the point 
    } // derivarives at measured points
    //
    // if the point contains material, consider its expected kinks, eloss 
    // as measurements
    if (pnt->ContainsMaterial()) {
      // at least 4 parameters: 2 spatial + 2 angular kinks with 0 expectaction
      int npm = pnt->GetNMatPar();
      const float* expMatCorr = pnt->GetMatCorrExp(); // expected correction (diagonalized)
      const float* expMatCov  = pnt->GetMatCorrCov(); // its error
      int offs  = pnt->GetMaxLocVarID() - npm;
      for (int ipar=0;ipar<npm;ipar++) {
	int parI = offs + ipar;
	vec[parI] += expMatCorr[ipar]/expMatCov[ipar]; // consider expectation as measurement
	mat(parI,parI) += 1./expMatCov[ipar]; // this measurement is orthogonal to all others
	printf("Pnt:%3d MatVar:%d DOF %3d | ExpVal: %+e Cov: %+e\n",ip,ipar,parI, expMatCorr[ipar], expMatCov[ipar]);
      }
    } // material effect descripotion params
    //
  } // loop over track points
  //
  return matp;
}

