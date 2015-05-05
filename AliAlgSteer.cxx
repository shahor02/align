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
#include "AliAlgPoint.h"
#include "AliAlgDet.h"
#include "AliAlgDetITS.h"
#include "AliAlgDetTPC.h"
#include "AliAlgDetTRD.h"
#include "AliAlgDetTOF.h"
#include "AliAlgVtx.h"
#include "AliAlgMPRecord.h"
#include "AliTrackerBase.h"
#include "AliESDCosmicTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliRecoParam.h"
#include "AliCDBRunRange.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "Mille.h"
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <stdio.h>

using namespace TMath;
using namespace AliAlgAux;

ClassImp(AliAlgSteer)

const Char_t* AliAlgSteer::fgkMPDataExt[AliAlgSteer::kMilleMPRec] = {".mille",".root"};
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
  ,fVtxSens(0)
  ,fSelEventSpecii(AliRecoParam::kCosmic|AliRecoParam::kLowMult|AliRecoParam::kHighMult)
  ,fObligatoryDetPattern(0)
  ,fCosmicSelStrict(kFALSE)
  ,fMinDetAcc(0)
  ,fPtMin(0.3)
  ,fEtaMax(1.5)
  ,fVtxMinCont(-1)
  ,fVtxMaxCont(-1)
  ,fVtxMinContVC(10)
  ,fMinITSClforVC(3)
  ,fITSPattforVC(kSPDAny)
  ,fMaxChi2forVC(10)
   //
  ,fGloParVal(0)
  ,fGloParErr(0)
  ,fRefPoint(0)
  ,fESDTree(0)
  ,fESDEvent(0)
  ,fVertex(0)
  ,fMPOutType(kMille)
  ,fMille(0)
  ,fMPRecord(0)
  ,fMPRecTree(0)
  ,fMPRecFile(0)
  ,fMilleDBuffer()
  ,fMilleIBuffer()
  ,fMPDatFileName("mpData")
  ,fMPParFileName("mpParam.txt")
  //
  ,fOutCDBPath("local://outOCDB")
  ,fOutCDBComment("AliAlgSteer")
  ,fOutCDBResponsible("")
   //
  ,fRefOCDBConf("configRefOCDB.C")
  ,fRefOCDBLoaded(0)
  ,fUseRecoOCDB(kTRUE)
{
  // def c-tor
  for (int i=kNDetectors;i--;) {
    fDetectors[i] = 0;
    fDetPos[i] = -1;
  }
  for (int i=kNCosmLegs;i--;) fESDTrack[i] = 0;
  memset(fStat,0,kNStatCl*kMaxStat*sizeof(float));
  fMaxDCAforVC[0] = 0.1;
  fMaxDCAforVC[1] = 0.6;
  SetOutCDBRunRange();
}

//________________________________________________________________
AliAlgSteer::~AliAlgSteer()
{
  // d-tor
  if (fMPRecFile) CloseMPRecOutput();
  if (fMille)     CloseMilleOutput();
  //
  delete fAlgTrack;
  delete[] fGloParVal;
  delete[] fGloParErr;
  for (int i=0;i<fNDet;i++) delete fDetectors[i];
  delete fVtxSens;
  delete fRefPoint;
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
  // special fake sensor for vertex constraint point
  fVtxSens = new AliAlgVtx();
  fVtxSens->PrepareMatrixL2G();
  fVtxSens->PrepareMatrixL2GIdeal();
  //
  fRefPoint = new AliAlgPoint();
  fRefPoint->SetSensor(fVtxSens);
  //
  for (int idt=0;idt<kNDetectors;idt++) {
    AliAlgDet* det = GetDetectorByDetID(idt);
    if (det) det->CacheReferenceOCDB();
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
  //
  if (!fRefOCDBLoaded) LoadRefOCDB();
  //
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
  SetESDEvent(esdEv);
  AliInfoF("Processing event %d of ev.specie %d",esdEv->GetEventNumberInFile(),esdEv->GetEventSpecie());
  //
  if (esdEv->GetRunNumber() != GetRunNumber()) SetRunNumber(esdEv->GetRunNumber());
  //
  if (!(esdEv->GetEventSpecie()&fSelEventSpecii)) {
    AliInfoF("Reject: specie does not match, allowed 0x%0x",fSelEventSpecii);
    return kFALSE;
  }
  //
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
  // fill needed points (tracking frame) in the fAlgTrack
  fRefPoint->SetContainsMeasurement(kFALSE);
  fRefPoint->SetContainsMaterial(kFALSE);
  fAlgTrack->AddPoint(fRefPoint); // reference point which the track will refer to
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
  // do we want to add the vertex as a measured point ?
  if (!AddVertexConstraint()) { // no constrain, just reference point w/o measurement
    fRefPoint->SetXYZTracking(0,0,0);
    fRefPoint->SetAlphaSens(Sector2Alpha(fAlgTrack->GetPoint(pntMeas)->GetAliceSector()));
  }
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
  fVertex = (ncont>=fVtxMinContVC) ? vtx : 0; // use vertex as a constraint
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
  fRefPoint->SetContainsMeasurement(kFALSE);
  fRefPoint->SetContainsMaterial(kFALSE);  
  fAlgTrack->AddPoint(fRefPoint); // reference point which the track will refer to
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
  fRefPoint->SetAlphaSens(Sector2Alpha(fAlgTrack->GetPoint(pntMeas)->GetAliceSector()));
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
  Bool_t res = kTRUE;
  if ((fMPOutType==kMille)||(fMPOutType==kMilleMPRec)) res &= FillMilleData();
  if ((fMPOutType==kMPRec)||(fMPOutType==kMilleMPRec)) res &= FillMPRecData();
  //
  return res;
}

//_________________________________________________________
Bool_t AliAlgSteer::FillMilleData()
{
  // store MP2 in Mille format
  if (!fMille) {
    TString mo = Form("%s%s",fMPDatFileName.Data(),fgkMPDataExt[kMille]);
    fMille = new Mille(mo.Data());
    if (!fMille) AliFatalF("Failed to create output file %s",mo.Data());
  }
  //
  if (!fAlgTrack->GetDerivDone()) {
    AliError("Track derivatives are not yet evaluated");
    return kFALSE;
  }
  int np(fAlgTrack->GetNPoints()),nDGloTot(0); // total number global derivatives stored
  int nParETP(fAlgTrack->GetNLocExtPar()); // numnber of local parameters for reference track param
  int nVarLoc(fAlgTrack->GetNLocPar());    // number of local degrees of freedom in the track
  float* buffDL(0),*buffDG(0);               // faster acces arrays
  int *buffI(0);
  //
  const int* gloParID(fAlgTrack->GetGloParID()); // IDs of global DOFs this track depends on
  for (int ip=0;ip<np;ip++) {
    AliAlgPoint* pnt = fAlgTrack->GetPoint(ip);
    if (pnt->ContainsMeasurement()) {
      int gloOffs = pnt->GetDGloOffs(); // 1st entry of global derivatives for this point
      int nDGlo   = pnt->GetNGloDOFs(); // number of global derivatives (number of DOFs it depends on)
      // check buffer sizes
      {
	if (fMilleDBuffer.GetSize()<nVarLoc+nDGlo) fMilleDBuffer.Set(100+nVarLoc+nDGlo);
	if (fMilleIBuffer.GetSize()<nDGlo)         fMilleIBuffer.Set(100+nDGlo);
	buffDL = fMilleDBuffer.GetArray(); // faster acces
	buffDG = buffDL + nVarLoc;         // faster acces	
	buffI  = fMilleIBuffer.GetArray(); // faster acces
      }
      // local der. array cannot be 0-suppressed by Mille construction, need to reset all to 0
      //
      for (int idim=0;idim<2;idim++) { // 2 dimensional orthogonal measurement
	memset(buffDL,0,nVarLoc*sizeof(float));
	const double* deriv  = fAlgTrack->GetDResDLoc(idim,ip);  // array of Dresidual/Dparams_loc
	// derivatives over reference track parameters
	for (int j=0;j<nParETP;j++) buffDL[j] = (IsZeroAbs(deriv[j])) ? 0:deriv[j];
	//
	// point may depend on material variables within this limits
	int lp0 = pnt->GetMinLocVarID(), lp1 = pnt->GetMaxLocVarID();
	for (int j=lp0;j<lp1;j++)   buffDL[j] = (IsZeroAbs(deriv[j])) ? 0:deriv[j];
	//
	// derivatives over global params: this array can be 0-suppressed, no need to reset
	int nGlo(0);
	deriv = fAlgTrack->GetDResDGlo(idim, gloOffs);
	const int* gloIDP(gloParID + gloOffs);
	for (int j=0;j<nDGlo;j++) {
	  if (IsZeroAbs(deriv[j])) {
	    buffDG[nGlo]  = deriv[j];        // value of derivative
	    buffI[nGlo++] = gloIDP[j];       // global DOF ID
	  }
	}
	fMille->mille(nVarLoc,buffDL, nGlo,buffDG, buffI, 
		      fAlgTrack->GetResidual(idim,ip),Sqrt(pnt->GetErrDiag(idim)));
	nDGloTot += nGlo;
	//
      }
    }
    if (pnt->ContainsMaterial()) {     // material point can add 4 or 5 otrhogonal pseudo-measurements
      memset(buffDL,0,nVarLoc*sizeof(float));      
      int nmatpar = pnt->GetNMatPar();  // residuals (correction expectation value)
      const float* expMatCorr = pnt->GetMatCorrExp(); // expected corrections (diagonalized)
      const float* expMatCov  = pnt->GetMatCorrCov(); // their diagonalized error matrix
      int offs  = pnt->GetMaxLocVarID() - nmatpar;    // start of material variables
      // here all derivatives are 1 = dx/dx
      for (int j=0,j1=j+offs;j<nmatpar;j++) { // mat. "measurements" don't depend on global params
	buffDL[j1] = 1.0;                     // only 1 non-0 derivative	
	fMille->mille(nVarLoc,buffDL,0,buffDG,buffI,expMatCorr[j],Sqrt(expMatCov[j]));
	buffDL[j1] = 0.0;                     // reset buffer
      }
    } // material "measurement"
  } // loop over points
  //
  if (!nDGloTot) {
    AliInfo("Track does not depend on free global parameters, discard");
    fMille->kill();
    return kFALSE;
  }
  fMille->end(); // store the record
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::FillMPRecData()
{
  // store MP2 in MPRecord format
  if (!fMPRecord) InitMPRecOutput();
  //
  fMPRecord->Clear();
  if (!fMPRecord->FillTrack(fAlgTrack)) return kFALSE;
  fMPRecord->SetRun(fRunNumber);
  fMPRecord->SetTimeStamp(fESDEvent->GetTimeStamp());
  UInt_t tID = 0xffff & UInt_t(fESDTrack[0]->GetID());
  if (IsCosmicEvent()) tID |= (0xffff & UInt_t(fESDTrack[1]->GetID()))<<16; 
  fMPRecord->SetTrackID(tID);
  fMPRecTree->Fill();
  return kTRUE;
}

//______________________________________________
void  AliAlgSteer::SetESDTree(const TTree* tr)
{
  // attach tree with input events
  fESDTree = tr;
  // go the the real ESD tree
  while (fESDTree->GetTree() && fESDTree->GetTree()!=fESDTree) fESDTree = fESDTree->GetTree();
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
  if (fRunNumber>0) fStat[kAccStat][kRun]++;
  fRunNumber = run;
  AliInfoF("Processing new run %d",fRunNumber);
  //
  if (!fUseRecoOCDB) {
    AliWarning("Reco-time OCDB will NOT be preloaded");
    return;
  }
  LoadRecoTimeOCDB();
  //
  for (int idet=0;idet<fNDet;idet++) GetDetector(idet)->AcknowledgeNewRun(run);
  //
  AliCDBManager::Destroy();
  //
  //  LoadRefOCDB(); //??? we need to get back reference OCDB ???
  //
  fStat[kInpStat][kRun]++;
  //
}

//_________________________________________________________
Bool_t AliAlgSteer::LoadRecoTimeOCDB()
{
  // Load OCDB paths used for the reconstruction of data being processed
  // In order to avoid unnecessary uploads, the objects are not actually 
  // loaded/cached but just added as specific paths with version
  AliInfoF("Preloading Reco-Time OCDB for run %d from ESD UserInfo list",fRunNumber);
  if (!fESDTree) {
    AliFatal("Cannot preload Reco-Time OCDB since the ESD tree is not set");
  }
  const TList* userInfo = const_cast<TTree*>(fESDTree)->GetUserInfo();
  TMap* cdbMap = (TMap*)userInfo->FindObject("cdbMap");
  TList* cdbList = (TList*)userInfo->FindObject("cdbList");
  //
  if (!cdbMap || !cdbList) {
    userInfo->Print();
    AliFatal("Failed to extract cdbMap and cdbList from UserInfo list");
  }
  //
  return PreloadOCDB(fRunNumber,cdbMap,cdbList);
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
  printf("MPData output :\t");
  if ((fMPOutType==kMille)||(fMPOutType==kMilleMPRec)) printf("%s%s ",fMPDatFileName.Data(),fgkMPDataExt[kMille]);
  if ((fMPOutType==kMPRec)||(fMPOutType==kMilleMPRec)) printf("%s%s ",fMPDatFileName.Data(),fgkMPDataExt[kMPRec]);
  printf("\n");
  printf("MP Params     :\t%s\n",fMPParFileName.Data());
  printf("\n");
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
void AliAlgSteer::InitMPRecOutput()
{
  // prepare MP record output
  if (!fMPRecord) fMPRecord = new AliAlgMPRecord();
  //
  TString mo = Form("%s%s",fMPDatFileName.Data(),fgkMPDataExt[kMPRec]);
  fMPRecFile = TFile::Open(mo.Data(),"recreate");
  if (!fMPRecFile) AliFatalF("Failed to create output file %s",mo.Data());
  //
  fMPRecTree = new TTree("mpTree","MPrecord Tree");
  fMPRecTree->Branch("mprec","AliAlgMPRecord",&fMPRecord);
  //
}

//____________________________________________
void AliAlgSteer::CloseMPRecOutput()
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
void AliAlgSteer::CloseMilleOutput()
{
  // close output
  delete fMille;
  fMille = 0;
}

//____________________________________________
void AliAlgSteer::SetMPDatFileName(const char* name) 
{
  // set output file name
  fMPDatFileName = name;
  // strip root or mille extensions, they will be added automatically later
  if      (fMPDatFileName.EndsWith(fgkMPDataExt[kMille])) 
    fMPDatFileName.Remove(fMPDatFileName.Length()-strlen(fgkMPDataExt[kMille]));
  else if (fMPDatFileName.EndsWith(fgkMPDataExt[kMPRec])) 
      fMPDatFileName.Remove(fMPDatFileName.Length()-strlen(fgkMPDataExt[kMPRec]));
  //
  if (fMPDatFileName.IsNull()) fMPDatFileName = "mpData"; 
  //
}

//____________________________________________
void AliAlgSteer::SetMPParFileName(const char* name) 
{
  // set output file name
  fMPParFileName = name; 
  if (fMPParFileName.IsNull()) fMPParFileName = "mpParam.txt"; 
  //
}

//____________________________________________
void AliAlgSteer::SetOutCDBPath(const char* name)
{
  // set output storage name
  fOutCDBPath = name; 
  if (fOutCDBPath.IsNull()) fOutCDBPath = "local://outOCDB"; 
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

//____________________________________________
Bool_t AliAlgSteer::AddVertexConstraint()
{
  // if vertex is set and if particle is primary, add vertex as a meared point
  //
  const AliESDtrack* esdTr = fESDTrack[0];
  if (!fVertex || !esdTr) return kFALSE;
  //
  if (esdTr->GetNcls(0)<fMinITSClforVC) return kFALSE; // not enough its clusters
  switch (fITSPattforVC) {
  case kSPDBoth: 
    if (!esdTr->HasPointOnITSLayer(0) || !esdTr->HasPointOnITSLayer(1)) return kFALSE;
    break;
  case kSPDAny:
    if (!esdTr->HasPointOnITSLayer(0) && !esdTr->HasPointOnITSLayer(1)) return kFALSE;
    break;
  case kSPD0:
    if (!esdTr->HasPointOnITSLayer(0)) return kFALSE;
    break;
  case kSPD1:    
    if (!esdTr->HasPointOnITSLayer(1)) return kFALSE;
    break;
  default: break;
  }
  //
  AliExternalTrackParam trc = *esdTr;
  Double_t dz[2],dzCov[3];
  if (!trc.PropagateToDCA(fVertex,AliTrackerBase::GetBz(),2*fMaxDCAforVC[0],dz,dzCov)) return kFALSE;
  //
  // check if primary candidate
  if (Abs(dz[0])>fMaxDCAforVC[0] || Abs(dz[1])>fMaxDCAforVC[1]) return kFALSE;
  Double_t covar[6]; fVertex->GetCovMatrix(covar);
  Double_t p[2] = {trc.GetParameter()[0]-dz[0],trc.GetParameter()[1]-dz[1]};
  Double_t c[3] = {0.5*(covar[0]+covar[2]),0.,covar[5]};
  Double_t chi2 = trc.GetPredictedChi2(p,c);
  if (chi2>fMaxChi2forVC) return kFALSE;
  //
  // assing measured vertex rotated to VtxSens frame as reference point
  double xyz[3],xyzT[3]; 
  fVertex->GetXYZ(xyz);
  fVtxSens->SetAlpha(trc.GetAlpha());
  // usually translation from GLO to TRA frame should go via matrix T2G
  // but for the VertexSensor Local and Global are the same frames
  fVtxSens->GetMatrixT2L().MasterToLocal(xyz,xyzT);
  fRefPoint->SetAlphaSens(fVtxSens->GetAlpTracking());
  fRefPoint->SetXYZTracking(xyzT);
  fRefPoint->SetYZErrTracking(c);
  fRefPoint->SetContainsMeasurement(kTRUE);
  fRefPoint->Init();
  //
  return kTRUE;
}

//______________________________________________________
void AliAlgSteer::WriteCalibrationResults() const
{
  // writes output calibration
  AliCDBManager::Destroy();
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(fOutCDBPath.Data());
  //
  AliAlgDet* det;
  for (int idet=0;idet<kNDetectors;idet++) {
    if (!(det=GetDetectorByDetID(idet))) continue;
    det->WriteCalibrationResults();
  }
  //
}

//______________________________________________________
void AliAlgSteer::SetOutCDBRunRange(int rmin,int rmax)
{
  // set output run range
  fOutCDBRunRange[0] = rmin >=0 ? rmin : 0;
  fOutCDBRunRange[1] = rmax>fOutCDBRunRange[0] ? rmax : AliCDBRunRange::Infinity();  
}

//______________________________________________________
Bool_t AliAlgSteer::LoadRefOCDB()
{
  // setup OCDB whose objects will be used as a reference with respect to which the
  // alignment/calibration will prodice its corrections.
  // Detectors which need some reference calibration data must use this one
  //
  //
  AliInfo("Loading reference OCDB");
  AliCDBManager::Destroy();
  AliCDBManager* man = AliCDBManager::Instance();
  if (!fRefOCDBConf.IsNull() && !gSystem->AccessPathName(fRefOCDBConf.Data(), kFileExists)) {
    AliInfoF("Executing reference OCDB setup macro %s",fRefOCDBConf.Data());
    gROOT->ProcessLine(Form(".x %s",fRefOCDBConf.Data()));
  }
  else {
    AliWarningF("No reference OCDB config macro %s is found, assume raw:// with run %d",
		fRefOCDBConf.Data(),AliCDBRunRange::Infinity());
    man->SetRaw(kTRUE);
    man->SetRun(AliCDBRunRange::Infinity());
  }
  //
  if (AliGeomManager::GetGeometry()) {
    AliInfo("Destroying current geometry before loading reference one");
    AliGeomManager::Destroy();
  }
  AliGeomManager::LoadGeometry("geometry.root");
  if (!AliGeomManager::GetGeometry()) AliFatal("Failed to load geometry, cannot run");
  //
  TString detList = "";
  for (int i=0;i<kNDetectors;i++) {detList += GetDetNameByDetID(i); detList += " ";}
  AliGeomManager::ApplyAlignObjsFromCDB(detList.Data());
  //
  fRefOCDBLoaded++;
  //
  return kTRUE;
}


//********************* interaction with PEDE **********************

//______________________________________________________
void AliAlgSteer::GenPedeParamFile(const Option_t *opt) const
{
  // produce steering file template for PEDE
  //
  TString opts = opt;
  opts.ToLower();
  AliInfoF("Generating MP2 parameters template file %s",fMPParFileName.Data());
  //
  FILE* flOut = fopen (fMPParFileName.Data(),"w+");
  //
  for (int idt=0;idt<kNDetectors;idt++) {
    AliAlgDet* det = GetDetectorByDetID(idt);
    if (!det) continue;
    det->WritePedeParamFile(flOut,opt);
    //
  }
  //
  fclose(flOut);
}

