#include "AliAlgTrack.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliAlgAux.h"
#include "TMatrixD.h"
#include "TVectorD.h"


using namespace AliAlgAux;
using namespace TMath;

// RS: this is not good: we define constants outside the class, but it is to
// bypass the CINT limitations on static arrays initializations 
const Int_t kRichardsonOrd = 1;              // Order of Richardson extrapolation for derivative (min=1)
const Int_t kRichardsonN = kRichardsonOrd+1; // N of 2-point symmetric derivatives needed for requested order
const Int_t kNRDClones = kRichardsonN*2     ;// number of variations for derivative of requested order

//____________________________________________________________________________
AliAlgTrack::AliAlgTrack() :
  fNLocPar(0)
  ,fNLocExtPar(0)
  ,fInnerPointID(0)
  //  ,fMinX2X0Pt2Account(0.5e-5/1.0)
  ,fMinX2X0Pt2Account(0.5e-3/1.0)
  ,fMass(0.14)
  ,fChi2(0)
  ,fPoints(0)
{
  // def c-tor
  for (int i=0;i<2;i++) fResidA[i] = fDerivA[i] = 0;
  fNeedInv[0] = fNeedInv[1] = kFALSE;
  //
}

//____________________________________________________________________________
AliAlgTrack::~AliAlgTrack()
{
  // d-tor
}

//____________________________________________________________________________
void AliAlgTrack::Clear(Option_t *)
{
  // reset the track
  ResetBit(0xffffffff);
  fPoints.Clear();
  fChi2 = 0;
  fInnerPointID = -1;
  fNeedInv[0] = fNeedInv[1] = kFALSE;
  //
}

//____________________________________________________________________________
void AliAlgTrack::DefineDOFs()
{
  // define varied DOF's (local parameters) for the track: 
  // 1) kinematic params (5 or 4 depending on Bfield)
  // 2) mult. scattering angles (2)
  // 3) if requested by point: energy loss
  //
  fNLocPar = fNLocExtPar = GetFieldON() ? kNKinParBON : kNKinParBOFF;
  int np = GetNPoints();
  //
  // the points are sorted in order opposite to track direction -> outer points come 1st,
  // but for the 2-leg cosmic track the innermost points are in the middle (1st lower leg, then upper one)
  //
  // start along track direction, i.e. last point in the ordered array
  int minPar = fNLocPar;
  for (int ip=GetInnerPointID()+1;ip--;) { // collision track or cosmic lower leg
    AliAlgPoint* pnt = GetPoint(ip);
    pnt->SetMinLocVarID(minPar);
    if (pnt->ContainsMaterial()) fNLocPar += (pnt->GetELossVaried()&&GetFieldON()) ? kNMSPar+kNELosPar:kNMSPar;
    pnt->SetMaxLocVarID(fNLocPar); // flag up to which parameted ID this points depends on
  }
  //
  if (IsCosmic()) {
    minPar = fNLocPar;
    for (int ip=GetInnerPointID()+1;ip<np;ip++) { // collision track or cosmic lower leg
      AliAlgPoint* pnt = GetPoint(ip);
      pnt->SetMinLocVarID(minPar);
      if (pnt->ContainsMaterial()) fNLocPar += (pnt->GetELossVaried()&&GetFieldON()) ? kNMSPar+kNELosPar:kNMSPar;
      pnt->SetMaxLocVarID(fNLocPar); // flag up to which parameted ID this points depends on
    }
  }
  //
  if (fLocPar.GetSize()<fNLocPar) fLocPar.Set(fNLocPar);
  fLocPar.Reset();
  fLocParA = fLocPar.GetArray();
  //  
  if (fResid[0].GetSize()<np) {
    fResid[0].Set(np);
    fResid[1].Set(np);
  }
  if (fDeriv[0].GetSize()<fNLocPar*np) {
    fDeriv[0].Set(fNLocPar*np);
    fDeriv[1].Set(fNLocPar*np);
  }
  for (int i=2;i--;) {
    fResid[i].Reset();
    fDeriv[i].Reset();
    fResidA[i] = fResid[i].GetArray();
    fDerivA[i] = fDeriv[i].GetArray();
  }
  //
  memcpy(fLocParA,GetParameter(),fNLocExtPar*sizeof(Double_t));
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResidDeriv(double *params)
{
  // Propagate for given local params and calculate residuals and their derivatives derivatives.
  // The 1st 4 or 5 elements of params vector should be the reference externalTrackParam
  // Then parameters of material corrections for each point
  // marked as having materials should come (4 or 5 dependending if ELoss is varied or fixed)
  //
  // If params are not provided, use internal params array
  //
  if (!params) params = fLocParA;
  //
  if (!GetResidDone()) CalcResiduals(params);
  //
  int np = GetNPoints();
  //
  // collision track or cosmic lower leg
  if (!CalcResidDeriv(params,fNeedInv[0],GetInnerPointID(),0)) return kFALSE;
  //
  if (IsCosmic()) { // cosmic upper leg
    if (!CalcResidDeriv(params,fNeedInv[1],GetInnerPointID()+1,np-1)) return kFALSE;
  } 
  //
  SetDerivDone();
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResidDeriv(double *params,Bool_t invert,int pFrom,int pTo)
{
  // Calculate derivatives of residuals vs params for points pFrom to pT. For cosmic upper leg 
  // track parameter may require inversion
  //
  // The derivatives are calculated using Richardson extrapolation 
  // (like http://root.cern.ch/root/html/ROOT__Math__RichardsonDerivator.html)
  //
  AliExternalTrackParam probD[kNRDClones];     // use this to vary supplied param for derivative calculation
  double varDelta[kRichardsonN];
  //
  const double kDelta[kNKinParBON]  = {0.02,0.02, 0.001,0.001, 0.01}; // variations for ExtTrackParam 
  const double kDeltaMat[kNMatDOFs] = {0.02,0.02, 0.001,0.001, 0.01}; // variations for ETP delta at points with material
  //
  int np,pinc;
  if (pTo>pFrom) { // fit in points decreasing order: cosmics upper leg
    pTo++;
    np = pTo - pFrom;
    pinc = 1;
  }
  else {           // fit in points increasing order: collision track or cosmics lower leg
    pTo--;
    np = pFrom - pTo;
    pinc = -1;
  }
  // 1) derivative wrt AliExternalTrackParam parameters
  for (int ipar=fNLocExtPar;ipar--;) {
    SetParams(probD,kNRDClones, GetX(),GetAlpha(),params);
    if (invert) for (int ic=kNRDClones;ic--;) probD[ic].Invert();
    double del = kDelta[ipar];
    if (ipar==kParq2Pt) { // for 1/pt variation use fractional increment
      del *= Abs(Get1P());
      if (del<1e-4) del = 1e-4;
    }
    for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
      varDelta[icl] = del;
      ModParam(probD[(icl<<1)+0], ipar, del);
      ModParam(probD[(icl<<1)+1], ipar,-del);
      del *= 0.5;
    }
    // propagate varied tracks to each point
    for (int ip=pFrom;ip!=pTo;ip+=pinc) { // points are ordered against track direction
      AliAlgPoint* pnt = GetPoint(ip);
      if ( !PropagateParamToPoint(probD, kNRDClones, pnt) ) return kFALSE;
      // account for materials
      if (pnt->ContainsMaterial()) { // apply material corrections
	Bool_t eLossFree = pnt->GetELossVaried() && GetFieldON();
	int nParFree = eLossFree ? kNMSPar + kNELosPar : kNMSPar;
	if (!ApplyMatCorr(probD, kNRDClones, &params[pnt->GetMaxLocVarID()-nParFree], eLossFree)) return kFALSE;
	if (!eLossFree && !ApplyELoss(probD, kNRDClones,pnt)) return kFALSE; // apply precalculated eloss
      }    
      //
      if (pnt->ContainsMeasurement()) {  
	int offsDer = ip*fNLocPar + ipar;
	RichardsonDeriv(probD, varDelta, pnt, fDerivA[0][offsDer], fDerivA[1][offsDer]); // calculate derivatives
      }
    } // loop over points
  } // loop over ExtTrackParam parameters
  //
  // 2) now vary material effect related parameters: MS and eventually ELoss
  //
  for (int ip=pFrom;ip!=pTo;ip+=pinc) { // points are ordered against track direction
    AliAlgPoint* pnt = GetPoint(ip);
    if (!pnt->ContainsMaterial()) continue;
    //
    Bool_t eLossFree = pnt->GetELossVaried() && GetFieldON();
    int nParFree = eLossFree ? kNMatDOFs : kNMatDOFs-kNELosPar;
    int offs = pnt->GetMaxLocVarID() - nParFree; // the parameters for this point start with this offset
    double *currPar = (double*)&params[offs];
    //
    for (int ipar=0;ipar<nParFree;ipar++) { // loop over DOFs related to MS and ELoss are point ip
      double del = kDeltaMat[ipar];
      if (ipar==kParq2Pt) {
	if (eLossFree) { // for eloss variation variation use fractional increment
	  del *= Abs(Get1P());
	  if (del<1e-4) del = 1e-4;
	}
	else continue;
      }
      //
      // we will vary the tracks starting from the original parameters propagated to given point and stored there
      SetParams(probD,kNRDClones, pnt->GetXPoint(),pnt->GetAlphaSens(),pnt->GetTrParamWS());
      // no need for eventual track inversion here: if needed, this is already done in ParamWS
      //
      for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
	varDelta[icl] = del;
	double parOrig = currPar[ipar];
	currPar[ipar] += del;
	//
	// apply varied material effects : incremented by delta
	if (!ApplyMatCorr(probD[(icl<<1)+0], currPar, eLossFree)) return kFALSE;
	//
	// apply varied material effects : decremented by delta
	currPar[ipar] = parOrig - del;
	if (!ApplyMatCorr(probD[(icl<<1)+1], currPar, eLossFree)) return kFALSE;
	//
	currPar[ipar] = parOrig;
	del *= 0.5;
      }     
      if (pnt->ContainsMeasurement()) {   // calculate derivatives at the scattering point itself
	int offsDer = ip*fNLocPar + offs + ipar;
	RichardsonDeriv(probD, varDelta, pnt, fDerivA[0][offsDer], fDerivA[1][offsDer]); // calculate derivatives for ip
      }
      //
      // loop over points whose residuals can be affected by the material effects on point ip
      for (int jp=ip+pinc;jp!=pTo;jp+=pinc) {
	AliAlgPoint* pntJ = GetPoint(jp);
	if ( !PropagateParamToPoint(probD, kNRDClones, pntJ) ) return kFALSE;
	//
	if (pntJ->ContainsMaterial()) { // apply MS and ELoss
	  Bool_t eLossFreeJ = pntJ->GetELossVaried()&&GetFieldON();
	  int nParFreeJ = eLossFreeJ ? kNMatDOFs : kNMatDOFs-kNELosPar;
	  int offsJ = pntJ->GetMaxLocVarID() - nParFreeJ;
	  double *currParJ = &params[offsJ];	  
	  if (!ApplyMatCorr(probD,kNRDClones,&params[offsJ],eLossFreeJ)) return kFALSE;
	  if (!eLossFreeJ && !ApplyELoss(probD, kNRDClones,pntJ)) return kFALSE; // apply precalculated eloss
	}
	if (pntJ->ContainsMeasurement()) {  
	  int offsDer = jp*fNLocPar + offs + ipar;
	  RichardsonDeriv(probD, varDelta, pntJ, fDerivA[0][offsDer], fDerivA[1][offsDer]); // calculate derivatives
	}
	//
      } // << loop over points whose residuals can be affected by the material effects on point ip
    } // << loop over DOFs related to MS and ELoss are point ip
  }  // << loop over all points of the track
  //  
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResiduals(const double *params)
{
  // Propagate for given local params and calculate residuals
  // The 1st 4 or 5 elements of params vector should be the reference externalTrackParam
  // Then parameters of material corrections for each point
  // marked as having materials should come (4 or 5 dependending if ELoss is varied or fixed)
  //
  // If params are not provided, use internal params array
  //
  if (!params) params = fLocParA;
  int np = GetNPoints();
  fChi2 = 0;
  //
  // collision track or cosmic lower leg
  if (!CalcResiduals(params,fNeedInv[0],GetInnerPointID(),0)) return kFALSE;
  //
  if (IsCosmic()) { // cosmic upper leg
    if (!CalcResiduals(params,fNeedInv[1],GetInnerPointID()+1,np-1)) return kFALSE;
  }
  //
  SetResidDone();
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResiduals(const double *params,Bool_t invert,int pFrom,int pTo)
{
  // calculate residuals for the single leg from points pFrom to pT
  //
  AliExternalTrackParam probe;
  SetParams(probe,GetX(),GetAlpha(),params);
  if (invert) probe.Invert();
  int np,pinc;
  if (pTo>pFrom) { // fit in points decreasing order: cosmics upper leg
    pTo++;
    np = pTo - pFrom;
    pinc = 1;
  }
  else {           // fit in points increasing order: collision track or cosmics lower leg
    pTo--;
    np = pFrom - pTo;
    pinc = -1;
  }
  //
  for (int ip=pFrom;ip!=pTo;ip+=pinc) { // points are ordered against track direction
    AliAlgPoint* pnt = GetPoint(ip);
    if (!PropagateParamToPoint(probe, pnt)) return kFALSE;
    //
    // account for materials
    if (pnt->ContainsMaterial()) { // apply material corrections
      Bool_t eLossFree = pnt->GetELossVaried() && GetFieldON();
      int nParFree = eLossFree ? kNMSPar + kNELosPar : kNMSPar;
      if (!ApplyMatCorr(probe, &params[pnt->GetMaxLocVarID()-nParFree], eLossFree)) return kFALSE;
      if (!eLossFree && !ApplyELoss(probe,pnt)) return kFALSE; // apply precalculated eloss
    }
    //
    // store the current track kinematics at the point
    pnt->SetTrParamWS(probe.GetParameter());
    if (pnt->ContainsMeasurement()) { // need to calculate residuals in the frame where errors are orthogonal
      pnt->GetResidualsDiag(probe.GetParameter(),fResidA[0][ip],fResidA[1][ip]);
      fChi2 += fResidA[0][ip]*fResidA[0][ip]/pnt->GetErrDiag(0);
      fChi2 += fResidA[1][ip]*fResidA[1][ip]/pnt->GetErrDiag(1);
    }
    //
  }
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::PropagateParamToPoint(AliExternalTrackParam* tr, int nTr, const AliAlgPoint* pnt)
{
  // Propagate set of tracks to the point  (only parameters, no error matrix)
  // VECTORIZE this
  //
  for (int itr=nTr;itr--;) {
    if (!PropagateParamToPoint(tr[itr],pnt)) {
#if DEBUG>3
      AliErrorF("Failed on clone %d",itr);
#endif
      return kFALSE;
    }
  }
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::PropagateParamToPoint(AliExternalTrackParam &tr, const AliAlgPoint* pnt)
{
  // propagate tracks to the point (only parameters, no error matrix)
  double xyz[3],bxyz[3];
  //
  if (!tr.RotateParamOnly(pnt->GetAlphaSens())) {
#if DEBUG>3
    AliErrorF("Failed to rotate to alpha=%f",pnt->GetAlphaSens());
    tr.Print();
    pnt->Print();
#endif
    return kFALSE;
  }
  tr.GetXYZ(xyz);
  //
  if (GetFieldON()) {
    if (pnt->GetUseBzOnly()) {
      if (!tr.PropagateParamOnlyTo(pnt->GetXPoint(),AliTrackerBase::GetBz(xyz))) {
#if DEBUG>3
	AliErrorF("Failed to propagate(BZ) to X=%f",pnt->GetXPoint());
	tr.Print();
	pnt->Print();
#endif
	return kFALSE;
      }
    }
    else {
      AliTrackerBase::GetBxByBz(xyz,bxyz);
      if (!tr.PropagateParamOnlyBxByBzTo(pnt->GetXPoint(),bxyz)) {
#if DEBUG>3
	AliErrorF("Failed to propagate(BXYZ) to X=%f",pnt->GetXPoint());
	tr.Print();
	pnt->Print();
#endif
	return kFALSE;
      }
    }
  }    
  else { // straigth line propagation
    if ( !tr.PropagateParamOnlyTo(pnt->GetXPoint(),0) ) {
#if DEBUG>3
      AliErrorF("Failed to propagate(B=0) to X=%f",pnt->GetXPoint());
      tr.Print();
      pnt->Print();
#endif
      return kFALSE;
    }
  }
  //
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::PropagateToPoint(AliExternalTrackParam &tr, const AliAlgPoint* pnt, 
				     int minNSteps, double maxStep, Bool_t matCor, double *matPar)
{
  // propagate tracks to the point. If matCor is true, then material corrections will be applied.
  // if matPar pointer is provided, it will be filled by total x2x0 and signed xrho
  if (!tr.Rotate(pnt->GetAlphaSens())) {
    AliDebug(5,Form("Failed to rotate to alpha=%f",pnt->GetAlphaSens()));
    return kFALSE;
  }
  //
  double xyz0[3],xyz1[3],bxyz[3],matarr[7];
  double xPoint=pnt->GetXPoint(),dx=xPoint-tr.GetX(),dxa=Abs(dx),step=dxa/minNSteps;
  if (matPar) matPar[0]=matPar[1]=0;
  if (dxa<kTinyDist) return kTRUE;
  if (step>maxStep) step = maxStep;
  int nstep = dxa/step;
  step = dxa/nstep;
  if (dx<0) step = -step;
  //
  //  printf("-->will go from X:%e to X:%e in %d steps of %f\n",tr.GetX(),xPoint,nstep,step);

  // do we go along or against track direction
  Bool_t alongTrackDir = (dx>0&&!pnt->IsInvDir()) || (dx<0&&pnt->IsInvDir());
  Bool_t queryXYZ = matCor||GetFieldON();
  if (queryXYZ) tr.GetXYZ(xyz0);
  //
  double x2X0Tot=0,xrhoTot=0;
  for (int ist=nstep;ist--;) { // single propagation step >>
    double xToGo = xPoint - step*ist;
    //
    if (GetFieldON()) {
      if (pnt->GetUseBzOnly()) {
	if (!tr.PropagateTo(xToGo,AliTrackerBase::GetBz(xyz0))) {
	  AliDebug(5,Form("Failed to propagate(BZ) to X=%f",xToGo));
	  return kFALSE;
	}
      }
      else {
	AliTrackerBase::GetBxByBz(xyz0,bxyz);
	if (!tr.PropagateToBxByBz(xToGo,bxyz)) {
	  AliDebug(5,Form("Failed to propagate(BXYZ) to X=%f",xToGo));
	  return kFALSE;
	}
      }
    }    
    else { // straigth line propagation
      if ( !tr.PropagateTo(xToGo,0) ) {
	AliDebug(5,Form("Failed to propagate(B=0) to X=%f",xToGo));
	return kFALSE;
      }
    }
    //
    if (queryXYZ) {
      //      printf("-> At %f (%f)\n", tr.GetX(), xToGo);      
      tr.GetXYZ(xyz1);
      if (matCor) {
	AliTrackerBase::MeanMaterialBudget(xyz0,xyz1,matarr);
	Double_t xrho=matarr[0]*matarr[4], xx0=matarr[1];
        if (alongTrackDir) xrho = -xrho; // if we go along track direction, energy correction is negative
	x2X0Tot += xx0;
	xrhoTot += xrho;
	//if (xx0>1e-5) {
	//  printf("    MB: {%.2f %.2f %.2f}->{%.2f %.2f %.2f} [%e %e] -> [%e %e]\n",
	//	 xyz0[0],xyz0[1],xyz0[2], xyz1[0],xyz1[1],xyz1[2],  xx0,xrho, x2X0Tot,xrhoTot);
	//}
	if (!tr.CorrectForMeanMaterial(xx0,xrho,fMass)) return kFALSE;
      }
      for (int l=3;l--;) xyz0[l] = xyz1[l];
    }
  } // single propagation step <<
  //
  if (matPar) {
    matPar[0] = x2X0Tot;
    matPar[1] = xrhoTot;
  }
  return kTRUE;
}

/*
//______________________________________________________
Bool_t AliAlgTrack::ApplyMS(AliExternalTrackParam& trPar, double tms,double pms)
{
  //------------------------------------------------------------------------------
  // Modify track par (e.g. AliExternalTrackParam) in the tracking frame 
  // (dip angle lam, az. angle phi) 
  // by multiple scattering defined by polar and azumuthal scattering angles in 
  // the track collinear frame (tms and pms resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  | Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] |   | Cos[tms]         |
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  | Sin[lam]	       -Cos[lam]	0     |   | Sin[pms]*Sin[tms]|
  //
  //------------------------------------------------------------------------------
  //
  double *par = (double*) trPar.GetParameter();
  //
  if (Abs(tms)<1e-7) return kTRUE;
  //
  double snTms = Sin(tms), csTms = Cos(tms);
  double snPms = Sin(pms), csPms = Cos(pms);  
  double snPhi = par[2],  csPhi = Sqrt((1.-snPhi)*(1.+snPhi));
  double csLam = 1./Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
  //
  double  r00 = csLam*csPhi, r01 = snLam*csPhi, &r02 = snPhi;
  double  r10 = csLam*snPhi, r11 = snLam*snPhi, &r12 = csPhi;
  double &r20 = snLam      ,&r21 = csLam;
  //
  double &v0 = csTms, v1 = snTms*csPms, v2 = snTms*snPms;
  //
  double px = r00*v0 + r01*v1 - r02*v2;
  double py = r10*v0 + r11*v1 + r12*v2;
  double pz = r20*v0 - r21*v1;
  //
  double pt = Sqrt(px*px + py*py);
  par[2] = py/pt;
  par[3] = pz/pt;
  par[4]*= csLam/pt;
  //
  return kTRUE;
}
*/

//______________________________________________________
Bool_t AliAlgTrack::ApplyMatCorr(AliExternalTrackParam& trPar, const Double_t *corrPar, Bool_t eloss)
{
  // Modify track param (e.g. AliExternalTrackParam) in the tracking frame 
  // by delta accounting for material effects
  //
  const double kMaxSnp = 0.95;
  double* par = (double*)trPar.GetParameter();
  double snpNew = par[kParSnp]+corrPar[kParSnp];
  if (Abs(snpNew)>kMaxSnp) {    
#if DEBUG>3
      AliErrorF("Snp is too large: %f",snpNew);
      printf("DeltaPar: "); 
      for (int i=0;i<kNMSPar+(eloss?kNELosPar:0);i++) printf("%+.3e ",corrPar[i]); printf("\n");
      trPar.Print();
#endif
    return kFALSE;
  }
  par[kParY]   += corrPar[kParY];
  par[kParZ]   += corrPar[kParZ];
  par[kParSnp]  = snpNew;
  par[kParTgl] += corrPar[kParTgl];
  if (eloss&&GetFieldON()) par[kParq2Pt] += corrPar[kParq2Pt];
  return kTRUE;
}


//______________________________________________________
Bool_t AliAlgTrack::ApplyMS(AliExternalTrackParam& trPar, double ms1,double ms2)
{
  //------------------------------------------------------------------------------
  // Modify track par (e.g. AliExternalTrackParam) in the tracking frame 
  // (dip angle lam, az. angle phi) 
  // by multiple scattering defined by scattering angles in 2 orthogonal directions in
  // the track collinear frame (ms1 and ms2 resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  | Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] |   | Cos[tms]         |
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  | Sin[lam]	       -Cos[lam]	0     |   | Sin[pms]*Sin[tms]|
  //
  //------------------------------------------------------------------------------
  //
  // convert 2 orthogonal scattering angles to point in cyl. cooordinates
  double snt2 = ms1*ms1 + ms2*ms2;
  if (snt2>kAlmostOneD) return kFALSE;
  if (IsZeroPos(snt2)) return kTRUE; // no scattering
  double snTms = Sqrt(snt2);
  double csTms = Sqrt(1.-snt2);
  double phiMS = ATan2(ms2,ms1);
  double snPms = Sin(phiMS), csPms = Cos(phiMS);  
  //
  double *par = (double*) trPar.GetParameter();
  //
  double snPhi = par[2],  csPhi = Sqrt((1.-snPhi)*(1.+snPhi));
  double csLam = 1./Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
  //
  double  r00 = csLam*csPhi, r01 = snLam*csPhi, &r02 = snPhi;
  double  r10 = csLam*snPhi, r11 = snLam*snPhi, &r12 = csPhi;
  double &r20 = snLam      ,&r21 = csLam;
  //
  double &v0 = csTms, v1 = snTms*csPms, v2 = snTms*snPms;
  //
  double px = r00*v0 + r01*v1 - r02*v2;
  double py = r10*v0 + r11*v1 + r12*v2;
  double pz = r20*v0 - r21*v1;
  //
  double pt = Sqrt(px*px + py*py);
  par[2] = py/pt;
  par[3] = pz/pt;
  par[4]*= csLam/pt;
  //
  return kTRUE;
}


//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam& trPar, const AliAlgPoint* pnt)
{
  // apply eloss effect to q/pt term
  if (!GetFieldON()) return kTRUE;
  double &p4 = ((double*)trPar.GetParameter())[kParq2Pt];
  p4 += pnt->GetMatCorrPar()[kParq2Pt];
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam& trPar, double dE)
{
  // apply eloss according to externally supplied dE
  double &p4 = ((double*)trPar.GetParameter())[4];
  double p = trPar.GetP();
  double p2 = p*p;
  Double_t e = Sqrt(p2 + fMass*fMass);
  if ( Abs(dE) > 0.3*e ) return kFALSE; //30% energy loss is too much!
  if ( (1.+ dE/p2*(dE + 2*e)) < 0. ) return kFALSE;
  double cP4 = 1./Sqrt(1.+ dE/p2*(dE + 2*e));  //A precise formula by Ruben !
  if (Abs(p4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c
  p4 *= cP4;
  return kTRUE;
}

/*
//______________________________________________________
Bool_t AliAlgTrack::ApplyMS(AliExternalTrackParam* trSet, int ntr, double tms,double pms)
{
  //------------------------------------------------------------------------------
  // Modify params for SET of tracks (e.g. AliExternalTrackParam) in the tracking frame 
  // (dip angle lam, az. angle phi) 
  // by multiple scattering defined by polar and azumuthal scattering angles in 
  // the track collinear frame (tms and pms resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  | Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] |   | Cos[tms]         |
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  | Sin[lam]	       -Cos[lam]	0     |   | Sin[pms]*Sin[tms]|
  //
  //------------------------------------------------------------------------------
  //
  // VECTORIZE THIS
  //
  if (Abs(tms)<1e-7) return kTRUE;
  //
  double snTms = Sin(tms), csTms = Cos(tms);
  double snPms = Sin(pms), csPms = Cos(pms);  
  //
  for (int itr=ntr;itr--;) { // at the moment just loop
    double *par = (double*) trSet[itr].GetParameter();
    //
    double snPhi = par[2],  csPhi = Sqrt((1.-snPhi)*(1.+snPhi));
    double csLam = 1./Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
    //
    double  r00 = csLam*csPhi, r01 = snLam*csPhi, &r02 = snPhi;
    double  r10 = csLam*snPhi, r11 = snLam*snPhi, &r12 = csPhi;
    double &r20 = snLam      ,&r21 = csLam;
    //
    double &v0 = csTms, v1 = snTms*csPms, v2 = snTms*snPms;
    //
    double px = r00*v0 + r01*v1 - r02*v2;
    double py = r10*v0 + r11*v1 + r12*v2;
    double pz = r20*v0 - r21*v1;
    //
    double pt = Sqrt(px*px + py*py);
    par[2] = py/pt;
    par[3] = pz/pt;
    par[4]*= csLam/pt;
  }
  //
  return kTRUE;
}
*/

//______________________________________________________
Bool_t AliAlgTrack::ApplyMatCorr(AliExternalTrackParam* trSet, int ntr, const Double_t *corrPar, Bool_t eloss)
{
  // Modify set of track params (e.g. AliExternalTrackParam) in the tracking frame 
  // by delta accounting for material effects
  for (int itr=ntr;itr--;) {
    if (!ApplyMatCorr(trSet[itr],corrPar,eloss)) {
#if DEBUG>3
      AliErrorF("Failed on clone %d",itr);
#endif      
      return kFALSE;
    }
  }
  return kTRUE;
}


//______________________________________________________
Bool_t AliAlgTrack::ApplyMS(AliExternalTrackParam* trSet, int ntr, double ms1,double ms2)
{
  //------------------------------------------------------------------------------
  // Modify params for SET of tracks (e.g. AliExternalTrackParam) in the tracking frame 
  // (dip angle lam, az. angle phi) 
  // by multiple scattering defined by polar and azumuthal scattering angles in 
  // the track collinear frame (tms and pms resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  | Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] |   | Cos[tms]         |
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  | Sin[lam]	       -Cos[lam]	0     |   | Sin[pms]*Sin[tms]|
  //
  //------------------------------------------------------------------------------
  //
  // VECTORIZE THIS
  //
  double snt2 = ms1*ms1 + ms2*ms2;
  if (snt2>kAlmostOneD) return kFALSE;
  if (IsZeroPos(snt2)) return kTRUE; // no scattering
  double snTms = Sqrt(snt2);
  double csTms = Sqrt(1.-snt2);
  double phiMS = ATan2(ms2,ms1);
  double snPms = Sin(phiMS), csPms = Cos(phiMS);  
  //
  for (int itr=ntr;itr--;) { // at the moment just loop
    double *par = (double*) trSet[itr].GetParameter();
    //
    double snPhi = par[2],  csPhi = Sqrt((1.-snPhi)*(1.+snPhi));
    double csLam = 1./Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
    //
    double  r00 = csLam*csPhi, r01 = snLam*csPhi, &r02 = snPhi;
    double  r10 = csLam*snPhi, r11 = snLam*snPhi, &r12 = csPhi;
    double &r20 = snLam      ,&r21 = csLam;
    //
    double &v0 = csTms, v1 = snTms*csPms, v2 = snTms*snPms;
    //
    double px = r00*v0 + r01*v1 - r02*v2;
    double py = r10*v0 + r11*v1 + r12*v2;
    double pz = r20*v0 - r21*v1;
    //
    double pt = Sqrt(px*px + py*py);
    par[2] = py/pt;
    par[3] = pz/pt;
    par[4]*= csLam/pt;
  }
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam* trSet, int ntr, const AliAlgPoint* pnt)
{
  // apply eloss according to x*rho of the point to SET of tracks
  //
  // VECTORIZE THIS
  //
  // at the moment just loop
  for (int itr=ntr;itr--;) {
    if (!ApplyELoss(trSet[itr],pnt)) {
#if DEBUG>3
      AliErrorF("Failed on clone %d",itr);
#endif
      return kFALSE;
    }
  }
  return kTRUE;
  //
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam* trSet, int ntr, double dE)
{
  // apply eloss according to externally supplied dE to SET of tracks
  //
  // VECTORIZE THIS
  //
  // at the moment just loop
  for (int itr=ntr;itr--;) if (!ApplyELoss(trSet[itr],dE)) return kFALSE;
  return kTRUE;
  //
}


//______________________________________________
Double_t AliAlgTrack::RichardsonExtrap(double *val, int ord)
{
  // Calculate Richardson extrapolation of order ord (starting from 1)
  // The array val should contain estimates ord+1 of derivatives with variations
  // d, d/2 ... d/2^ord.
  // The array val is overwritten
  //
  if (ord==1) return (4.*val[1] - val[0])*(1./3);
  do {for (int i=0;i<ord;i++) val[i] = (4.*val[i+1] - val[i])*(1./3);} while(--ord);
  return val[0];
}

//______________________________________________
Double_t AliAlgTrack::RichardsonExtrap(const double *val, int ord)
{
  // Calculate Richardson extrapolation of order ord (starting from 1)
  // The array val should contain estimates ord+1 of derivatives with variations
  // d, d/2 ... d/2^ord.
  // The array val is not overwritten
  //
  if (ord==1) return (4.*val[1] - val[0])*(1./3);
  double* buff = new double[ord+1];
  memcpy(buff,val,(ord+1)*sizeof(double));
  do {for (int i=0;i<ord;i++) buff[i] = (4.*buff[i+1] - buff[i])*(1./3);} while(--ord);
  return buff[0];
}

//______________________________________________
void AliAlgTrack::RichardsonDeriv(const AliExternalTrackParam* trSet, const double *delta, const AliAlgPoint* pnt, double& derY, double& derZ)
{
  // Calculate Richardson derivatives for diagonalized Y and Z from a set of kRichardsonN pairs 
  // of tracks with same parameter of i-th pair varied by +-delta[i]
  static double derRichY[kRichardsonN],derRichZ[kRichardsonN];
  //
  for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
    double resYVP=0,resYVN=0,resZVP=0,resZVN=0;
    pnt->GetResidualsDiag(trSet[(icl<<1)+0].GetParameter(), resYVP, resZVP); // variation with +delta
    pnt->GetResidualsDiag(trSet[(icl<<1)+1].GetParameter(), resYVN, resZVN); // variation with -delta
    derRichY[icl] = 0.5*(resYVP-resYVN)/delta[icl];   // 2-point symmetric derivatives
    derRichZ[icl] = 0.5*(resZVP-resZVN)/delta[icl];
  }
  derY = RichardsonExtrap(derRichY,kRichardsonOrd);   // dY/dPar
  derZ = RichardsonExtrap(derRichZ,kRichardsonOrd);  // dZ/dPar
  //
}

//______________________________________________
void AliAlgTrack::Print(Option_t *opt) const
{
  // print track data
  printf("%s ",IsCosmic() ? "  Cosmic  ":"Collision ");
  AliExternalTrackParam::Print();
  printf("N Free Params: %d, for kinematics: %d | Npoints: %d (Inner:%d) | M : %.3f | Chi2: %.3f\n",fNLocPar,fNLocExtPar,
	 GetNPoints(),GetInnerPointID(),fMass,fChi2);
  //
  TString optS = opt;
  optS.ToLower();
  Bool_t res = optS.Contains("r") && GetResidDone();
  Bool_t der = optS.Contains("d") && GetDerivDone();
  if (optS.Contains("p") || res || der) { 
    for (int ip=0;ip<GetNPoints();ip++) {
      printf("#%3d ",ip);
      GetPoint(ip)->Print(opt);  
      if (res) printf("Residuals  : %+.3e %+.3e\n",GetResidual(0,ip),GetResidual(1,ip));
      if (der) {
	for (int ipar=0;ipar<fNLocPar;ipar++) {
	  printf("Dres/dp%03d : %+.3e %+.3e\n",ipar,GetDerivative(0,ip)[ipar], GetDerivative(1,ip)[ipar]);
	}
      }
    }
  } // print points
}

//______________________________________________
Bool_t AliAlgTrack::IniFit() 
{
  // perform initial fit of the track
  //
  const int    kMinNStep = 3;
  const double kMaxDefStep = 3.0; 
  //
  AliExternalTrackParam trc = *this;
  //
  fChi2 = 0;
  //
  // the points are ranged from outer to inner for collision tracks, 
  // and from outer point of lower leg to outer point of upper leg for the cosmic track 
  //
  // the fit will always start from the outgoing track in inward direction
  if (!FitLeg(trc,0,GetInnerPointID(),fNeedInv[0])) return kFALSE; // collision track or cosmic lower leg
  //
  printf("Lower leg: %d %d\n",0,GetInnerPointID()); trc.Print();
  //
  if (IsCosmic()) {
    AliExternalTrackParam trcU = trc;
    if (!FitLeg(trcU,GetNPoints()-1,GetInnerPointID()+1,fNeedInv[1])) return kFALSE; //fit upper leg of cosmic track
    //
    // propagate to reference point, which is the inner point of lower leg
    const AliAlgPoint* refP = GetPoint(GetInnerPointID());
    if (!PropagateToPoint(trcU,refP,kMinNStep,kMaxDefStep,kTRUE)) return kFALSE;
    //
    printf("Upper leg: %d %d\n",GetInnerPointID()+1,GetNPoints()-1); trcU.Print();
    //
    if (!CombineTracks(trc,trcU)) return kFALSE;
    printf("Combined\n"); 
    trc.Print();    
  }
  CopyFrom(&trc);
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::CombineTracks(AliExternalTrackParam& trcL, const AliExternalTrackParam& trcU)
{
  // Assign to trcL the combined tracks (Kalman update of trcL by trcU)
  // The trcL and trcU MUST be defined at same X,Alpha
  //
  // Update equations: tracks described by vectors vL and vU and coviriances CL and CU resp.
  // then the gain matrix K = CL*(CL+CU)^-1
  // Updated vector and its covariance:
  // CL' = CL - K*CL
  // vL' = vL + K(vU-vL)
  //
  if (Abs(trcL.GetX()-trcU.GetX())>kTinyDist || Abs(trcL.GetAlpha()-trcU.GetAlpha())>kTinyDist) {
    AliError("Tracks must be defined at same reference X and Alpha");
    trcL.Print();
    trcU.Print();
    return kFALSE;
  }
  //
  const double* covU=trcU.GetCovariance(),*parU=trcU.GetParameter();
  double* covL=(double*)trcL.GetCovariance(),*parL=(double*)trcL.GetParameter();
  //
  int mtSize = GetFieldON() ? kNKinParBON : kNKinParBOFF;
  TMatrixD matCL(mtSize,mtSize),matCLplCU(mtSize,mtSize);
  TVectorD vl(mtSize),vUmnvL(mtSize);
  //
  //
  for (int i=mtSize;i--;) {
    vUmnvL[i] = parU[i] - parL[i];     // y = residual of 2 tracks
    vl[i]  = parL[i];
    for (int j=i+1;j--;) {
      int indIJ = ((i*(i+1))>>1)+j; // position of IJ cov element in the AliExternalTrackParam covariance array
      matCL(i,j) = matCL(j,i) = covL[indIJ];
      matCLplCU(i,j) = matCLplCU(j,i) = covL[indIJ] + covU[indIJ];
    }
  }
  matCLplCU.Invert();                      // S^-1 = (Cl + Cu)^-1
  if (!matCLplCU.IsValid()) { 
#if DEBUG>3
    AliError("Failed to invert summed cov.matrix of cosmic track");
    matCLplCU.Print();
#endif    
    return kFALSE; // inversion failed
  }
  TMatrixD matK(matCL,TMatrixD::kMult,matCLplCU); // gain K = Cl*(Cl+Cu)^-1
  TMatrixD matKdotCL(matK,TMatrixD::kMult,matCL); // K*Cl
  TVectorD vlUp = matK*vUmnvL;                   // K*(vl - vu)
  for (int i=mtSize;i--;) {
    parL[i] += vlUp[i];     // updated param: vL' = vL + K(vU-vL)
    for (int j=i+1;j--;) covL[((i*(i+1))>>1)+j] -= matKdotCL(i,j); // updated covariance: Cl' = Cl - K*Cl
  } 
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::FitLeg(AliExternalTrackParam& trc, int pFrom,int pTo, Bool_t &inv) 
{
  // perform initial fit of the track
  // the fit will always start from the outgoing track in inward direction (i.e. if cosmics - bottom leg)
  const int    kMinNStep = 3;
  const double kMaxDefStep = 3.0; 
  const double kErrSpace=50.;
  const double kErrAng = 0.5;
  const double kErrRelPtI = 0.7;
  const double kIniErr[15] = { // initial error
    kErrSpace*kErrSpace,
    0                  , kErrSpace*kErrSpace,
    0                  ,                   0, kErrAng*kErrAng,
    0                  ,                   0,               0, kErrAng*kErrAng,
    0                  ,                   0,               0,               0, kErrRelPtI*kErrRelPtI
  };
  const Double_t kOverShootX = 5;//kMaxDefStep*0.7;
  //
  // prepare seed at outer point
  AliAlgPoint* p0 = GetPoint(pFrom);
  double phi = trc.Phi(),alp=p0->GetAlphaSens();
  BringTo02Pi(phi);
  BringTo02Pi(alp);
  double dphi = DeltaPhiSmall(phi,alp); // abs delta angle
  if (dphi>Pi()/2.) { // need to invert the track to new frame
    inv = kTRUE;
    //    printf("Fit in %d %d Delta: %.3f -> Inverting for\n",pFrom,pTo,dphi); 
    //    p0->Print("meas");
    //    printf("BeforeInv "); trc.Print();
    trc.Invert();
    //    printf("After Inv "); trc.Print();
  }
  if (!trc.RotateParamOnly(p0->GetAlphaSens())) return kFALSE;
  if (!trc.PropagateParamOnlyTo(p0->GetXPoint()+kOverShootX,AliTrackerBase::GetBz())) return kFALSE;
  double* cov = (double*)trc.GetCovariance();
  memcpy(cov,kIniErr,15*sizeof(double));
  cov[14] *= trc.GetSigned1Pt()*trc.GetSigned1Pt();
  //
  int np,pinc;
  if (pTo>pFrom) { // fit in points increasing order: collision track or cosmics lower leg
    pTo++;
    np = pTo - pFrom;
    pinc = 1;
  }
  else {          // fit in points decreasing order: cosmics upper leg
    pTo--;
    np = pFrom - pTo;
    pinc = -1;
  }
  Bool_t res = 0;
  //
  for (int ip=pFrom;ip!=pTo;ip+=pinc) { // inward fit from outer point
    AliAlgPoint* pnt = GetPoint(ip);
    //
    if (!PropagateToPoint(trc,pnt,kMinNStep, kMaxDefStep, kTRUE)) return kFALSE;
    if (pnt->ContainsMeasurement()) {
      const double* yz    = pnt->GetYZTracking();
      const double* errYZ = pnt->GetYZErrTracking();
      double chi = trc.GetPredictedChi2(yz,errYZ);
      if (!trc.Update(yz,errYZ)) return kFALSE;
      fChi2 += chi;
    }
  }
  //
  if (inv) {
    //    printf("Before inverting back "); trc.Print();
    trc.Invert();
  }
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::ProcessMaterials() 
{
  // attach material effect info to alignment points
  AliExternalTrackParam trc = *this;

  // collision track of cosmic lower leg: move along track direction from last (middle for cosmic lower leg) 
  // point (inner) to 1st one (outer)
  if (fNeedInv[0]) trc.Invert(); // track needs to be inverted ? (should be for upper leg)
  if (!ProcessMaterials(trc, GetInnerPointID(),0)) {
#if DEBUG>3
    AliError("Failed to process materials for leg along the track");
#endif    
    return kFALSE;
  }
  if (IsCosmic()) {
    // cosmic upper leg: move againg the track direction from middle point (inner) to last one (outer)
    trc = *this;
    if (fNeedInv[1]) trc.Invert(); // track needs to be inverted ? 
    if (!ProcessMaterials(trc, GetInnerPointID()+1,GetNPoints()-1)) {
#if DEBUG>3
      AliError("Failed to process materials for leg against the track");
#endif    
      return kFALSE;
    }
  }
  return kTRUE;
}

/*
//______________________________________________
Bool_t AliAlgTrack::ProcessMaterials() 
{
  // attach material effect info to alignment points
  const int    kMinNStep = 3;
  const double kMaxDefStep = 3.0; 
  const double kErrSpcT = 1e-6;
  const double kErrAngT = 1e-6;
  const double kErrPtIT = 1e-12;
  const double kErrSpcH = 10.0;
  const double kErrAngH = 0.5;
  const double kErrPtIH = 0.5;
  const double kErrTiny[15] = { // initial tiny error
    kErrSpcT*kErrSpcT,
    0                  , kErrSpcT*kErrSpcT,
    0                  ,                   0, kErrAngT*kErrAngT,
    0                  ,                   0,               0, kErrAngT*kErrAngT,
    0                  ,                   0,               0,               0, kErrPtIT*kErrPtIT
  };
  const double kErrHuge[15] = { // initial tiny error
    kErrSpcH*kErrSpcH,
    0                  , kErrSpcH*kErrSpcH,
    0                  ,                   0, kErrAngH*kErrAngH,
    0                  ,                   0,               0, kErrAngH*kErrAngH,
    0                  ,                   0,               0,               0, kErrPtIH*kErrPtIH
  };
  //
  // 2 copies of the track, one will be propagated accounting for materials, other - w/o
  AliExternalTrackParam tr1 = *this, tr0;
  double x2X0xRho[2] = {0,0};
  //
  // here we move in track direction
  for (int ip=GetInnerPointID()+1;ip--;) { // point are order against track direction
    AliAlgPoint* pnt = GetPoint(ip);
    memcpy((double*)tr1.GetCovariance(),kErrTiny,15*sizeof(double)); // assign tiny errors to both tracks
    tr0 = tr1;
    //
    //printf(">>MatP%d ",ip); pnt->Print();
    if (!PropagateToPoint(tr1,pnt,kMinNStep, kMaxDefStep, kTRUE ,x2X0xRho)) return kFALSE; // with material corrections
    //
    // is there enough material to consider the point as a scatterer?
    if (x2X0xRho[0]*Abs(tr1.GetSigned1Pt()) < GetMinX2X0Pt2Account()) { // ignore materials
      pnt->SetContainsMaterial(kFALSE);
      continue;
    }
    //
    if (!PropagateToPoint(tr0,pnt,kMinNStep, kMaxDefStep, kFALSE,0)    ) return kFALSE; // no material corrections
    double cov1[15];
    memcpy(cov1,tr1.GetCovariance(),15*sizeof(double));                           // save errors with mat.effect
    if (pnt->ContainsMeasurement()) {
      memcpy((double*)tr1.GetCovariance(),kErrHuge,15*sizeof(double));  // assign large errors
      const double* yz    = pnt->GetYZTracking();
      const double* errYZ = pnt->GetYZErrTracking();
      if (!tr1.Update(yz,errYZ)) return kFALSE;                         // adjust to measurement      
    }
    //
    double *cov0=(double*)tr0.GetCovariance(),*par0=(double*)tr0.GetParameter(),*par1=(double*)tr1.GetParameter();
    double *covP=pnt->GetMatCorrCov(),*parP=pnt->GetMatCorrPar();
    for (int l=15;l--;) covP[l] = cov1[l] - cov0[l];
    for (int l=5; l--;) parP[l] = par1[l] - par0[l];
    pnt->SetContainsMaterial(kTRUE);
    pnt->SetX2X0(x2X0xRho[0]);
    pnt->SetXTimesRho(x2X0xRho[1]);
    //printf("Add mat%d %e %e\n",ip, x2X0xRho[0],x2X0xRho[1]);
    //
  }
  //
  if (IsCosmic()) {
    AliFatal(" TODO cosmic");
  }
  //
  return kTRUE;
}
*/

//______________________________________________
Bool_t AliAlgTrack::ProcessMaterials(AliExternalTrackParam& trc, int pFrom,int pTo) 
{
  // attach material effect info to alignment points
  const int    kMinNStep = 3;
  const double kMaxDefStep = 3.0; 
  const double kErrSpcT = 1e-6;
  const double kErrAngT = 1e-6;
  const double kErrPtIT = 1e-12;
  const double kErrSpcH = 10.0;
  const double kErrAngH = 0.5;
  const double kErrPtIH = 0.5;
  const double kErrTiny[15] = { // initial tiny error
    kErrSpcT*kErrSpcT,
    0                  , kErrSpcT*kErrSpcT,
    0                  ,                   0, kErrAngT*kErrAngT,
    0                  ,                   0,               0, kErrAngT*kErrAngT,
    0                  ,                   0,               0,               0, kErrPtIT*kErrPtIT
  };
  const double kErrHuge[15] = { // initial tiny error
    kErrSpcH*kErrSpcH,
    0                  , kErrSpcH*kErrSpcH,
    0                  ,                   0, kErrAngH*kErrAngH,
    0                  ,                   0,               0, kErrAngH*kErrAngH,
    0                  ,                   0,               0,               0, kErrPtIH*kErrPtIH
  };
  //
  // 2 copies of the track, one will be propagated accounting for materials, other - w/o
  AliExternalTrackParam tr0;
  double x2X0xRho[2] = {0,0};
  //
  int np,pinc;
  if (pTo>pFrom) { // fit in points decreasing order: cosmics upper leg
    pTo++;
    np = pTo - pFrom;
    pinc = 1;
  }
  else {           // fit in points increasing order: collision track or cosmics lower leg
    pTo--;
    np = pFrom - pTo;
    pinc = -1;
  }
  //
  for (int ip=pFrom;ip!=pTo;ip+=pinc) { // points are ordered against track direction
    AliAlgPoint* pnt = GetPoint(ip);
    memcpy((double*)trc.GetCovariance(),kErrTiny,15*sizeof(double)); // assign tiny errors to both tracks
    tr0 = trc;
    //
    if (!PropagateToPoint(trc,pnt,kMinNStep, kMaxDefStep, kTRUE ,x2X0xRho)) {  // with material corrections
#if DEBUG>3
      AliErrorF("Failed to take track to point %d (dir: %d -> %d) with mat.corr.",ip,pFrom,pTo);
      trc.Print();
      pnt->Print("meas");
#endif      
      return kFALSE;
    }
    //
    // is there enough material to consider the point as a scatterer?
    if (x2X0xRho[0]*Abs(trc.GetSigned1Pt()) < GetMinX2X0Pt2Account()) { // ignore materials
      pnt->SetContainsMaterial(kFALSE);
      continue;
    }
    //
    if (!PropagateToPoint(tr0,pnt,kMinNStep, kMaxDefStep, kFALSE,0)) { // no material corrections
#if DEBUG>3
      AliErrorF("Failed to take track to point %d (dir: %d -> %d) w/o mat.corr.",ip,pFrom,pTo);
      tr0.Print();
      pnt->Print("meas");
#endif      
      return kFALSE; 
    }
    double cov1[15];
    memcpy(cov1,trc.GetCovariance(),15*sizeof(double));                           // save errors with mat.effect
    if (pnt->ContainsMeasurement()) {
      memcpy((double*)trc.GetCovariance(),kErrHuge,15*sizeof(double));  // assign large errors
      const double* yz    = pnt->GetYZTracking();
      const double* errYZ = pnt->GetYZErrTracking();
      if (!trc.Update(yz,errYZ)) return kFALSE;                         // adjust to measurement      
    }
    //
    double *cov0=(double*)tr0.GetCovariance(),*par0=(double*)tr0.GetParameter(),*par1=(double*)trc.GetParameter();
    double *covP=pnt->GetMatCorrCov(),*parP=pnt->GetMatCorrPar();
    for (int l=15;l--;) covP[l] = cov1[l] - cov0[l];
    for (int l=5; l--;) parP[l] = par1[l] - par0[l];
    pnt->SetContainsMaterial(kTRUE);
    pnt->SetX2X0(x2X0xRho[0]);
    pnt->SetXTimesRho(x2X0xRho[1]);
    //printf("Add mat%d %e %e\n",ip, x2X0xRho[0],x2X0xRho[1]);
    //
  }
  //
  return kTRUE;
}

//______________________________________________
void AliAlgTrack::SortPoints()
{
  // sort points in order against track direction: innermost point is last
  // for collision tracks. 
  // For 2-leg cosmic tracks: 1st points of outgoing (lower) leg are added from large to
  // small radii, then the points of incomint (upper) leg are added in increasing R direction
  //
  // The fInnerPointID will mark the id of the innermost point, i.e. the last one for collision-like
  // tracks and in case of cosmics - the point of lower leg with smallest R
  //
  fPoints.Sort();
  int np = GetNPoints();
  fInnerPointID = np-1;
  if (IsCosmic()) {
    for (int ip=np;ip--;) {
      AliAlgPoint* pnt = GetPoint(ip);
      if (pnt->IsInvDir()) continue;   // this is a point of upper leg
      fInnerPointID = ip;
      break;
    }
  }
  //
}

