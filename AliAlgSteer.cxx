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
#include "AliAlgMPRecord.h"
#include "AliTrackerBase.h"
#include "AliESDCosmicTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliRecoParam.h"
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>

using namespace TMath;
using namespace AliAlgAux;

ClassImp(AliAlgSteer)


const Char_t* AliAlgSteer::fgkDetectorName[AliAlgSteer::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID" };
const Int_t   AliAlgSteer::fgkSkipLayers[AliAlgSteer::kNLrSkip] = {AliGeomManager::kPHOS1,AliGeomManager::kPHOS2,
								   AliGeomManager::kMUON,AliGeomManager::kEMCAL};

const Char_t* AliAlgSteer::fgkStatClName[AliAlgSteer::kNStatCl] = {"Inp: ","Acc: "};
const Char_t* AliAlgSteer::fgkStatName[AliAlgSteer::kMaxStat] =   
  {"runs","Ev.Coll", "Ev.Cosm", "Trc.Coll", "Trc.Cosm"};

//________________________________________________________________
AliAlgSteer::AliAlgSteer()
  :fNDet(0)
  ,fNDOFs(0)
  ,fRunNumber(-1)
  ,fFieldOn(kFALSE)
  ,fCosmicEvent(kFALSE)
  ,fAlgTrack(0)
  ,fSelEventSpecii(AliRecoParam::kCosmic|AliRecoParam::kLowMult|AliRecoParam::kHighMult)
  ,fObligatoryDetPattern(0)
  ,fCosmicSelStrict(kFALSE)
  ,fMinDetAcc(0)
  ,fPtMin(0.3)
  ,fEtaMax(1.5)
  ,fVtxMinCont(-1)
  ,fVtxMaxCont(-1)
  ,fVtxMinContCS(10)
  ,fGloParVal(0)
  ,fGloParErr(0)
  ,fRefPoint()
  ,fESDEvent(0)
  ,fVertex(0)
  ,fMPRecord(0)
  ,fMPRecTree(0)
  ,fMPRecFile(0)
  ,fMPFileName("mpRecord.root")
{
  // def c-tor
  for (int i=kNDetectors;i--;) {
    fDetectors[i] = 0;
    fDetPos[i] = -1;
  }
  for (int i=kNCosmLegs;i--;) fESDTrack[i] = 0;
  memset(fStat,0,kNStatCl*kMaxStat*sizeof(float));
}

//________________________________________________________________
AliAlgSteer::~AliAlgSteer()
{
  // d-tor
  if (fMPRecFile) CloseMPOutput();
  delete fMPRecord;
  //
  delete fAlgTrack;
  delete[] fGloParVal;
  delete[] fGloParErr;
  for (int i=0;i<fNDet;i++) delete fDetectors[i];
  //
}

//________________________________________________________________
void AliAlgSteer::InitDetectors()
{
  // init all detectors geometry
  //
  static Bool_t done = kFALSE;
  if (done) return;
  done = kTRUE;
  //
  fAlgTrack = new AliAlgTrack();
  //
  int dofCnt = 0;
  for (int i=0;i<fNDet;i++) dofCnt += fDetectors[i]->InitGeom();
  if (!dofCnt) AliFatal("No DOFs found");
  //
  fGloParVal = new Float_t[dofCnt];
  fGloParErr = new Float_t[dofCnt];
  memset(fGloParVal,0,dofCnt*sizeof(Float_t));
  memset(fGloParErr,0,dofCnt*sizeof(Float_t));
  AliInfoF("Booked %d global parameters",dofCnt);
  //
  fNDOFs = 0;
  for (int idt=0;idt<kNDetectors;idt++) {
    AliAlgDet* det = GetDetectorByDetID(idt);
    if (det) fNDOFs += det->AssignDOFs();
  }
  
}

//________________________________________________________________
void AliAlgSteer::InitDOFs()
{
  // scan all free global parameters, link detectors to array of params
  //
  static Bool_t done = kFALSE;
  if (done) return;
  done = kTRUE;
  //
  for (int i=0;i<fNDet;i++) fDetectors[i]->InitDOFs();
  //
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
  det->SetAlgSteer(this);
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
Bool_t AliAlgSteer::CheckDetectorPattern(UInt_t patt) const 
{
  //validate detector pattern
  return (patt&fObligatoryDetPattern)==fObligatoryDetPattern && NumberOfBitsSet(patt)>=fMinDetAcc;
}

//_________________________________________________________
UInt_t AliAlgSteer::AcceptTrack(const AliESDtrack* esdTr, Bool_t strict) const
{
  // decide if the track should be processed
  AliAlgDet* det = 0;
  UInt_t detAcc = 0;
  if (fFieldOn && esdTr->Pt()<fPtMin) return 0;
  if (Abs(esdTr->Eta())>fEtaMax)      return 0;
  //
  for (int idet=0;idet<kNDetectors;idet++) {
    if (!(det=GetDetectorByDetID(idet))) continue;
    if (!det->AcceptTrack(esdTr)) {
      if (strict && det->IsObligatory()) return 0;
      else continue;
    }
    //
    detAcc |= 0x1<<idet;
  }
  if (NumberOfBitsSet(detAcc)<fMinDetAcc) return 0;
  return detAcc;
  //
}

//_________________________________________________________
UInt_t AliAlgSteer::AcceptTrackCosmic(const AliESDtrack* esdPairCosm[kNCosmLegs]) const
{
  // decide if the pair of tracks making cosmic track should be processed
  UInt_t detAcc=0,detAccLeg;
  for (int i=kNCosmLegs;i--;) {
    detAccLeg = AcceptTrack(esdPairCosm[i],fCosmicSelStrict); // missing obligatory detectors in one leg might be allowed
    if (!detAccLeg) return 0;
    detAcc |= detAccLeg;
  }
  if (fCosmicSelStrict) return detAcc;
  //
  // for non-stric selection check convolution of detector presence
  if (!CheckDetectorPattern(detAcc)) return 0;
  return detAcc;
  //
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessEvent(const AliESDEvent* esdEv)
{
  // process event
  AliInfoF("Processing event %d of ev.specie %d",esdEv->GetEventNumberInFile(),esdEv->GetEventSpecie());
  if (!(esdEv->GetEventSpecie()&fSelEventSpecii)) {
    AliInfoF("Reject: specie does not match, allowed 0x%0x",fSelEventSpecii);
    return kFALSE;
  }
  //
  SetESDEvent(esdEv);
  SetCosmicEvent(esdEv->GetEventSpecie()==AliRecoParam::kCosmic);
  SetFieldOn(Abs(esdEv->GetMagneticField())>kAlmostZeroF);
  if (!IsCosmicEvent() && !CheckSetVertex(esdEv->GetPrimaryVertexTracks())) return kFALSE;
  //
  int accTr = 0;
  if (IsCosmicEvent()) {
    fStat[kInpStat][kEventCosm]++;
    int ntr = esdEv->GetNumberOfCosmicTracks();
    for (int itr=0;itr<ntr;itr++) {
      accTr += ProcessTrack(esdEv->GetCosmicTrack(itr));      
    }
    if (accTr) fStat[kAccStat][kEventCosm]++;
  }
  else {
    fStat[kInpStat][kEventColl]++;
    int ntr = esdEv->GetNumberOfTracks();
    for (int itr=0;itr<ntr;itr++) {
      accTr += ProcessTrack(esdEv->GetTrack(itr));      
    }
    if (accTr) fStat[kAccStat][kEventColl]++;
  }    
  //
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessTrack(const AliESDtrack* esdTr)
{
  // process single track
  //
  fStat[kInpStat][kTrackColl]++;
  fESDTrack[0] = esdTr;
  fESDTrack[1] = 0;
  //
  int nPnt = 0;
  const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
  if (!trF) return kFALSE;
  const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
  if (!trPoints || (nPnt=trPoints->GetNPoints())<1) return kFALSE;
  //
  UInt_t detAcc = AcceptTrack(esdTr);
  if (!detAcc) return kFALSE;
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
    if (!(detAcc&(0x1<<idet))) continue;
    det=GetDetectorByDetID(idet);
    if (!det->ProcessPoints(esdTr, fAlgTrack)) {
      detAcc &= ~(0x1<<idet); // did not survive, suppress detector in the track
      if (det->IsObligatory()) return kFALSE;
    }
    if (NumberOfBitsSet(detAcc)<fMinDetAcc) return kFALSE; // abandon track
  }
  //
  fAlgTrack->Set(esdTr->GetX(),esdTr->GetAlpha(),esdTr->GetParameter(),esdTr->GetCovariance());
  fAlgTrack->CopyFrom(esdTr);
  fAlgTrack->SetFieldON( GetFieldOn() );
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
  if (!StoreProcessedTrack()) return kFALSE;
  //
  fStat[kAccStat][kTrackColl]++;
  //
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::CheckSetVertex(const AliESDVertex *vtx)
{ 
  // vertex selection/constraint check
  if (!vtx) {
    fVertex = 0;
    return kTRUE;
  }
  int ncont = vtx->GetNContributors();
  if (fVtxMinCont>0 && fVtxMinCont>ncont) {
#if DEBUG>2    
    AliInfoF("Rejecting event with %d vertex contributors (min %d asked)",ncont,fVtxMinCont);
#endif
    return kFALSE;
  }
  if (fVtxMaxCont>0 && ncont>fVtxMaxCont) {
#if DEBUG>2    
    AliInfoF("Rejecting event with %d vertex contributors (max %d asked)",ncont,fVtxMaxCont);
#endif
    return kFALSE;
  }
  fVertex = (ncont>=fVtxMinContCS) ? vtx : 0; // use vertex as a constraint
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessTrack(const AliESDCosmicTrack* cosmTr)
{
  // process single cosmic track
  //
  fStat[kInpStat][kTrackCosm]++;
  int nPnt = 0;
  fESDTrack[0] = 0;
  fESDTrack[1] = 0;
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
    fESDTrack[leg] = esdTr;
  }
  //
  UInt_t detAcc = AcceptTrackCosmic(fESDTrack);
  if (!detAcc) return kFALSE;
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
  UInt_t detAccLeg[2] = {detAcc,detAcc};
  for (int leg=kNCosmLegs;leg--;) {
    for (int idet=0;idet<kNDetectors;idet++) {
      if (!(detAccLeg[leg]&(0x1<<idet))) continue;
      det = GetDetectorByDetID(idet);
      //
      // upper leg points marked as the track goes in inverse direction
      if (!det->ProcessPoints(fESDTrack[leg],fAlgTrack, leg==kCosmUp)) {
	if (fCosmicSelStrict && det->IsObligatory()) return kFALSE;
	detAccLeg[leg] &= ~(0x1<<idet); // did not survive, suppress detector
      }
      if (NumberOfBitsSet(detAccLeg[leg])<fMinDetAcc) return kFALSE; // abandon track
    }
  }
  // last check on legs-combined patter
  if (!CheckDetectorPattern(detAccLeg[0]&detAccLeg[1])) return kFALSE;
  //
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
  if (!StoreProcessedTrack()) return kFALSE;
  //
  fStat[kAccStat][kTrackCosm]++;
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::StoreProcessedTrack()
{
  // write alignment track
  if (!fMPRecord) InitMPOutput();
  //
  fMPRecord->Clear();
  if (!fMPRecord->FillTrack(fAlgTrack)) return kFALSE;
  fMPRecord->SetRun(fRunNumber);
  fMPRecord->SetTimeStamp(fESDEvent->GetTimeStamp());
  UInt_t tID = 0xffff & UInt_t(fESDTrack[0]->GetID());
  if (IsCosmicEvent()) tID |= (0xffff & UInt_t(fESDTrack[1]->GetID()))<<16; 
  fMPRecord->SetTrackID(tID);
  fMPRecTree->Fill();
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
  if (run==fRunNumber) return;  // nothing to do
  fRunNumber = run;
  for (int idet=0;idet<fNDet;idet++) GetDetector(idet)->AcknowledgeNewRun(run);
  fStat[kInpStat][kRun]++;
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
  TString opts = opt; 
  opts.ToLower();
  printf("%5d DOFs in %d detectors\n",fNDOFs,fNDet);
  for (int idt=0;idt<kNDetectors;idt++) {
    AliAlgDet* det = GetDetectorByDetID(idt);
    if (det) det->Print(opt);
    else printf("Detector:%5s is not defined\n",GetDetNameByDetID(idt));
  }
  //
  if (opts.Contains("stat")) PrintStatistics();
}

//________________________________________________________
void AliAlgSteer::PrintStatistics() const
{
  // print processing stat
  printf("\nProcessing Statistics\n");
  printf("Type: ");
  for (int i=0;i<kMaxStat;i++) printf("%s ",fgkStatName[i]); printf("\n");
  for (int icl=0;icl<kNStatCl;icl++) {
    printf("%s ",fgkStatClName[icl]);
    for (int i=0;i<kMaxStat;i++) printf(Form("%%%dd ",(int)strlen(fgkStatName[i])),(int)fStat[icl][i]); printf("\n");
  }
}

//________________________________________________________
void AliAlgSteer::ResetDetectors()
{
  // reset detectors for next track
  for (int idet=fNDet;idet--;) {
    AliAlgDet* det = GetDetector(idet);
    det->ResetPool();   // reset used alignment points
  }
}

//____________________________________________
Bool_t AliAlgSteer::TestLocalSolution()
{
  // test track local solution
  TVectorD rhs;
  AliSymMatrix* mat = BuildMatrix(rhs);
  if (!mat) return kFALSE;
  TVectorD vsl(rhs);
  if (!mat->SolveChol(rhs,vsl,kTRUE)) {
   delete mat;
   return kFALSE;
  }
  //
  // print solution vector
  int nlocpar = fAlgTrack->GetNLocPar();
  int nlocparETP = fAlgTrack->GetNLocExtPar(); // parameters of external track param
  printf("ETP Update: ");
  for (int i=0;i<nlocparETP;i++) printf("%+.2e(%+.2e) ",vsl[i],Sqrt((*mat)(i,i))); printf("\n");
  //
  if (nlocpar>nlocparETP) printf("Mat.Corr. update:\n");
  for (int ip=fAlgTrack->GetNPoints();ip--;) {
    AliAlgPoint* pnt = fAlgTrack->GetPoint(ip);
    int npm = pnt->GetNMatPar();
    const float* expMatCov  = pnt->GetMatCorrCov(); // its error
    int offs  = pnt->GetMaxLocVarID() - npm;
    for (int ipar=0;ipar<npm;ipar++) {
      int parI = offs + ipar;
      double err = Sqrt(expMatCov[ipar]);
      printf("Pnt:%3d MatVar:%d DOF %3d | %+.3e(%+.3e) -> sig:%+.3e -> pull: %+.2e\n",
	     ip,ipar,parI,vsl[parI],Sqrt((*mat)(parI,parI)), err,vsl[parI]/err);
    }
  }
  //
  fAlgTrack->CalcResiduals(vsl.GetMatrixArray());
  fAlgTrack->SetLocPars(vsl.GetMatrixArray());
  delete mat;
  //
  return kTRUE;
}

//____________________________________________
AliSymMatrix* AliAlgSteer::BuildMatrix(TVectorD &vec)
{
  // build matrix/vector for local track
  int npnt = fAlgTrack->GetNPoints();
  int nlocpar = fAlgTrack->GetNLocPar();
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
	double* deriv  = fAlgTrack->GetDResDLoc(idim,ip);  // array of Dresidual/Dparams
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
	vec[parI] -= expMatCorr[ipar]/expMatCov[ipar]; // consider expectation as measurement
	mat(parI,parI) += 1./expMatCov[ipar]; // this measurement is orthogonal to all others
	//printf("Pnt:%3d MatVar:%d DOF %3d | ExpVal: %+e Cov: %+e\n",ip,ipar,parI, expMatCorr[ipar], expMatCov[ipar]);
      }
    } // material effect descripotion params
    //
  } // loop over track points
  //
  return matp;
}

//____________________________________________
void AliAlgSteer::InitMPOutput()
{
  // prepare MP record output
  if (!fMPRecord) fMPRecord = new AliAlgMPRecord();
  //
  fMPRecFile = TFile::Open(fMPFileName.Data(),"recreate");
  if (!fMPRecFile) AliFatalF("Failed to create output file %s",fMPFileName.Data());
  //
  fMPRecTree = new TTree("mpTree","MPrecord Tree");
  fMPRecTree->Branch("mprec","AliAlgMPRecord",&fMPRecord);
  //
}

//____________________________________________
void AliAlgSteer::CloseMPOutput()
{
  // close output
  if (!fMPRecFile) return;
  fMPRecFile->cd();
  fMPRecTree->Write();
  delete fMPRecTree;
  fMPRecTree = 0;
  fMPRecFile->Close();
  delete fMPRecFile;
  fMPRecFile = 0;
  delete fMPRecord;
  fMPRecord = 0;
}

//____________________________________________
void AliAlgSteer::SetMPRecFileName(const char* name) 
{
  // set output file name
  fMPFileName = name; 
  if (fMPFileName.IsNull()) fMPFileName = "mpRecord.root"; 
  //
}

//____________________________________________
void AliAlgSteer::SetObligatoryDetector(Int_t detID, Bool_t v)
{
  // mark detector presence obligatory in the track
  AliAlgDet* det = GetDetectorByDetID(detID);
  if (!det) {
    AliErrorF("Detector %d is not defined",detID);
  }
  if (v) fObligatoryDetPattern |=  0x1<<v;
  else   fObligatoryDetPattern &=~(0x1<<v);
  if (det->IsObligatory()!=v) det->SetObligatory(v);
  //
}
