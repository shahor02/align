#include "AliAlgTrack.h"

//____________________________________________________________________________
AliAlgTrack::AliAlgTrack() :
  fPoints(0)
{
  // def c-tor
  
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
  for (int ip=0;ip<np;ip++) {
    AliAlgPoint* pnt = GetPoint(ip);
    pnt->SetMaxLocVarID(fNLocPar); // flag up to which parameted ID this points depends on
    if (pnt->ContainsMaterial()) {
      fNLocPar += kNMSPar;
      if (pnt->GetVaryELoss()) fNLocPar += kNELosPar;
    }
  }
  //
  int nOld = fResid[0].GetSize();
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
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResidDeriv(const double *params)
{
  // Propagate for given local params and calculate residuals and their derivatives derivatives.
  // The 1st 4 or 5 elements of params vector should be the externalTrackParam 
  //
  // The derivatives are calculated using Richardson extrapolation 
  // (like http://root.cern.ch/root/html/ROOT__Math__RichardsonDerivator.html)
  //
  static AliExternalTrackParam probD[kNRDClones];     // use this to vary supplied param for derivative calculation
  static double varDelta[kRichardsonN];
  //
  const double kDelta[kNKinParBON] = {0.02,0.02, 0.001,0.001, 0.01}; // variations for ExtTrackParam 
  const double kDeltaMat[kNMatDOFs] = {0.001, 0.1, 0.01}; // variations for MS theta, phi and ELoass
  //
  if (!GetResidDone()) CalcResiduals(params);
  //
  int np = GetNPoints();
  //
  // 1) derivative wrt AliExternalTrackParam parameters
  for (int ipar=fNLocExtPar;ipar--;) {
    SetParams(probD,kNRDClones, GetX(),GetAlpha(),params);
    double del = kDelta[ipar];
    if (ipar==kNKinParBON-1) { // for 1/pt variation use fractional increment
      del *= TMath::Abs(Get1P());
      if (del<1e-4) del = 1e-4;
    }
    for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
      varDelta[icl] = del;
      ModParam(probD[(icl<<1)+0], ipar, del);
      ModParam(probD[(icl<<1)+1], ipar,-del);
      del *= 0.5;
    }
    // propagate varied tracks to each point
    for (int ip=0;ip<np;ip++) {
      AliAlgPoint* pnt = GetPoint(ip);
      if ( !PropagateToPoint(probD, kNRDClones, pnt) ) return kFALSE;
      //
      if (pnt->ContainsMeasurement()) {  
	int offsDer = ip*fNLocPar + ipar;
	RichardsonDeriv(probD, varDelta, pnt, derRichY[offsDer], derRichZ[offsDer]); // calculate derivatives
      }
      //
      if (pnt->ContainsMaterial()) { // apply MS and ELoss
	int offs = pnt->GetMaxLocVarID();
	if (!ApplyMS(probD,params[offs+kMSTheta],params[offs+kMSPhi])) return kFALSE;
	if (GetFieldON()) {	
	  if (pnt->GetELossVaried()) {if (!ApplyELoss(probD,kNRDClones,params[offs+kELoss])) return kFALSE;}
	  else                       {if (!ApplyELoss(probD,kNRDClones,pnt)) return kFALSE;}
	}
      }      
      //
    } // loop over points
    //    
  } // loop over ExtTrackParam parameters
  //
  // now vary material effect related parameters: MS and eventually ELoss
  //
  for (int ip=0;ip<np;ip++) { // loop over all points of the track
    AliAlgPoint* pnt = GetPoint(ip);
    if (!pnt->ContainsMaterial()) continue;
    int offs = pnt->GetMaxLocVarID(); // the parameters for this point start with this offset
    //
    double *vpars = (double*)&params[offs];
    for (int ipar=0;ipar<kNMatDOFs;ipar++) { // loop over DOFs related to MS and ELoss are point ip
      double del = kDeltaMat[ipar];
      if (ipar==kELoss) {
	if(!pnt->GetELossVaried()) continue;      
	// for eloss variation variation use fractional increment
	del *= TMath::Abs(Get1P());
	if (del<1e-4) del = 1e-4;
      }
      //
      for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
	varDelta[icl] = del;
	vpar[ipar] += del;
	// apply varied material effects
	if (!ApplyMS(probD[(icl<<1)+0], vpar[kMSTheta], vpar[kMSPhi])) return kFALSE;
	if (GetFieldON()) {	
	  if (pnt->GetELossVaried()) {if (!ApplyELoss(probD[(icl<<1)+0],kNRDClones,vpar[kELoss])) return kFALSE;}
	  else                       {if (!ApplyELoss(probD[(icl<<1)+0],kNRDClones,pnt)) return kFALSE;}
	}
	vpar[ipar] -= del+del;
	if (!ApplyMS(probD[(icl<<1)+1], vpar[kMSTheta], vpar[kMSPhi])) return kFALSE;
	if (GetFieldON()) {	
	  if (pnt->GetELossVaried()) {if (!ApplyELoss(probD[(icl<<1)+1],kNRDClones,vpar[kELoss])) return kFALSE;}
	  else                       {if (!ApplyELoss(probD[(icl<<1)+1],kNRDClones,pnt)) return kFALSE;}
	}
	//
	vpar[ipar] += del;
	del *= 0.5;
      }
      //
      for (int jp=ip+1;jp<np;jp++) { // loop over points whose residuals can be affected by the material effects on point ip
	AliAlgPoint* pntj = GetPoint(jp);
	if ( !PropagateToPoint(probD, kNRDClones, pntj) ) return kFALSE;
	//
	if (pntj->ContainsMeasurement()) {  
	  int offsDer = ip*fNLocPar + offs + ipar;
	  RichardsonDeriv(probD, varDelta, pntj, derRichY[offsDer], derRichZ[offsDer]); // calculate derivatives
	}
	//
	if (pntj->ContainsMaterial()) { // apply MS and ELoss
	  int offsj = pntj->GetMaxLocVarID();
	  if (!ApplyMS(probD,params[offsj+kMSTheta],params[offsj+kMSPhi])) return kFALSE;
	  if (GetFieldON()) {	
	    if (pntj->GetELossVaried()) {if (!ApplyELoss(probD,kNRDClones,params[offsj+kELoss])) return kFALSE;}
	    else                        {if (!ApplyELoss(probD,kNRDClones,pntj)) return kFALSE;}
	  }
	} 
	//
      } // << loop over points whose residuals can be affected by the material effects on point ip
      //
    } // << loop over DOFs related to MS and ELoss are point ip
  }  // << loop over all points of the track
  //  
  return kTRUE;
}
    // init variation clones by the track propagated to this point during residuals calculation
    SetParams(probD,kNRDClones, pnt->GetXTracking(),pnt->GetAlpha(),pnt->GetTrParamWS());
    
kDeltaMat

    if ( !PropagateToPoint(probe, pnt, params) ) return kFALSE;
    pnt->SetTrParamWS(probe.GetParameter());    // store the current track kinematics at the point
    //
    if (pnt->ContainsMeasurement()) { // need to calculate residuals in the frame where errors are orthogonal
      pnt->GetResidualsDiag(probe.GetParameter(),fResidA[0][ip],fResidA[1][ip]);
    }
  }
  //

  SetDerivDone();
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::CalcResiduals(const double *params)
{
  // Propagate for given local params and calculate residuals
  // The 1st 4 or 5 elements of params vector should be the externalTrackParam 
  static AliExternalTrackParam probe;
  SetParam(probe,GetX(),GetAlpha(),params);
  //
  int np = GetNPoints();
  for (int ip=0;ip<np;ip++) {
    AliAlgPoint* pnt = GetPoint(ip);
    if ( !PropagateToPoint(probe, pnt) ) return kFALSE;
    pnt->SetTrParamWS(probe.GetParameter());    // store the current track kinematics at the point (RS: do we need it here?)
    //
    if (pnt->ContainsMeasurement()) { // need to calculate residuals in the frame where errors are orthogonal
      pnt->GetResidualsDiag(probe.GetParameter(),fResidA[0][ip],fResidA[1][ip]);
    }
    //
    // account for materials
    if (pnt->ContainsMaterial()) { // apply MS and ELoss
      int offs = pnt->GetMaxLocVarID();
      if (!ApplyMS(probe,params[offs+kMSTheta],params[offs+kMSPhi])) return kFALSE;
      if (GetFieldON()) {	
	if (pnt->GetELossVaried()) {if (!ApplyELoss(probe,params[offs+kELoss])) return kFALSE;}
	else                       {if (!ApplyELoss(probe,pnt)) return kFALSE;}
      }
    }
  }
  //
  SetResidDone();
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::PropagateToPoint(AliExternalTrackParam** tr, int nTr, const AliAlgPoint* pnt)
{
  // Propagate set of tracks to the point
  // VECTORIZE this
  //
  double xyz[3],bxyz[3],bz;
  //
  for (int itr=nTr;itr--;) {
    if (!PropagateToPointMulti(*tr[itr],pnt)) return kFALSE;
  }
  return kTRUE;
}

//______________________________________________________
Bool_t AliAlgTrack::PropagateToPoint(AliExternalTrackParam &tr, const AliAlgPoint* pnt)
{
  // propagate tracks to the point
  double xyz[3],bxyz[3],bz;
  //
  if (!tr.RotateParamOnly(pnt->GetAlpha())) return kFALSE;
  tr.GetXYZ(xyz);
  //
  if (GetFieldON()) {
    if (pnt->GetUseBzOnly()) {
      if (!tr.PropagateParamOnlyTo(pnt->GetXTracking(),AliTrackerBase::GetBz(xyz))) return kFALSE;
    }
    else {
      AliTrackerBase::GetBxByBz(xyz,bxyz);
      if (!tr->PropagateParamOnlyBxByBzTo(pnt->GetXTracking(),bxyz)) return kFALSE;
    }
  }    
  else { // straigth line propagation
    if ( !tr.PropagateParamOnlyTo(pnt->GetXTracking(),0) ) return kFALSE;
  }
  //
  return kTRUE;
}

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
  //  / Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] \   / Cos[tms]         \
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  \ Sin[lam]	       -Cos[lam]	0     /   \ Sin[pms]*Sin[tms]/
  //
  //------------------------------------------------------------------------------
  //
  double *par = (double*) trPar.GetParameter();
  //
  if (TMath::Abs(tms)<1e-7) return kTRUE;
  //
  double snTms = TMath::Sin(tms), csTms = TMath::Cos(tms);
  double snPms = TMath::Sin(pms), csPms = TMath::Cos(pms);  
  double snPhi = par[2],  csPhi = TMath::Sqrt((1.-snPhi)*(1.+snPhi));
  double csLam = 1./TMath::Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
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
  double pt = TMath::Sqrt(px*px + py*py);
  par[2] = py/pt;
  par[3] = pz/pt;
  par[4]*= csLam/pt;
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam& trPar, const AliAlgPoint* pnt)
 {
   // apply eloss according to x*rho of the point
   double &p4 = ((double*)trPar.GetParameter())[4];
   double p = trPar.GetP();
   double p2 = p*p;
   Double_t e = TMath::Sqrt(p2 + fMass*fMass);
   Double_t bg = p/fMass;
   double dE = Bethe(bg)*pnt->SetXTimesRho();
   if ( TMath::Abs(dE) > 0.3*e ) return kFALSE; //30% energy loss is too much!
   if ( (1.+ dE/p2*(dE + 2*e)) < 0. ) return kFALSE;
   cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*e));  //A precise formula by Ruben !
   if (TMath::Abs(p4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c
   p4 *= cP4;
   return kTRUE;
 }

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam& trPar, double dE)
 {
   // apply eloss according to externally supplied dE
   double &p4 = ((double*)trPar.GetParameter())[4];
   double p = trPar.GetP();
   double p2 = p*p;
   Double_t e = TMath::Sqrt(p2 + fMass*fMass);
   if ( TMath::Abs(dE) > 0.3*e ) return kFALSE; //30% energy loss is too much!
   if ( (1.+ dE/p2*(dE + 2*e)) < 0. ) return kFALSE;
   cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*e));  //A precise formula by Ruben !
   if (TMath::Abs(p4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c
   p4 *= cP4;
   return kTRUE;
 }


//______________________________________________________
Bool_t AliAlgTrack::ApplyMS(AliExternalTrackParam** trSet, int ntr, double tms,double pms)
{
  //------------------------------------------------------------------------------
  // Modify params for SET of tracks (e.g. AliExternalTrackParam) in the tracking frame 
  // (dip angle lam, az. angle phi) 
  // by multiple scattering defined by polar and azumuthal scattering angles in 
  // the track collinear frame (tms and pms resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  / Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] \   / Cos[tms]         \
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  \ Sin[lam]	       -Cos[lam]	0     /   \ Sin[pms]*Sin[tms]/
  //
  //------------------------------------------------------------------------------
  //
  // VECTORIZE THIS
  //
  if (TMath::Abs(tms)<1e-7) return kTRUE;
  //
  double snTms = TMath::Sin(tms), csTms = TMath::Cos(tms);
  double snPms = TMath::Sin(pms), csPms = TMath::Cos(pms);  
  //
  for (int itr=ntr;itr--;) { // at the moment just loop
    double *par = (double*) trSet[itr]->GetParameter();
    //
    double snPhi = par[2],  csPhi = TMath::Sqrt((1.-snPhi)*(1.+snPhi));
    double csLam = 1./TMath::Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
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
    double pt = TMath::Sqrt(px*px + py*py);
    par[2] = py/pt;
    par[3] = pz/pt;
    par[4]*= csLam/pt;
  }
  //
  return kTRUE;
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam** trSet, int ntr, const AliAlgPoint* pnt)
{
  // apply eloss according to x*rho of the point to SET of tracks
  //
  // VECTORIZE THIS
  //
  // at the moment just loop
  for (int itr=ntr;itr--;) if (!ApplyELoss(*trSet[itr],pnt)) return kFASLE;
  return kTRUE;
  //
}

//______________________________________________
Bool_t AliAlgTrack::ApplyELoss(AliExternalTrackParam** trSet, int ntr, double dE)
{
  // apply eloss according to externally supplied dE to SET of tracks
  //
  // VECTORIZE THIS
  //
  // at the moment just loop
  for (int itr=ntr;itr--;) if (!ApplyELoss(*trSet[itr],dE)) return kFASLE;
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
  if (ord==1) return (4.*fval[1] - val[0])*(1./3);
  do {for (int i=0;i<nord;i++) val[i] = (4.*val[i+1] - val[i])*(1./3);} while(--ord);
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
  if (ord==1) return (4.*fval[1] - val[0])*(1./3);
  double* buff = new double[ord+1];
  memcpy(buff,val,(ord+1)*sizeof(double));
  do {for (int i=0;i<nord;i++) buff[i] = (4.*buff[i+1] - buff[i])*(1./3);} while(--ord);
  return buff[0];
}

//______________________________________________
void AliAlgTrack::RichardsonDeriv(const AliExternalTrackParam** trSet, const double *delta, const AliAlgPoint* pnt, double& derY, double& derZ)
{
  // Calculate Richardson derivatives for diagonalized Y and Z from a set of kRichardsonN pairs 
  // of tracks with same parameter of i-th pair varied by +-delta[i]
  static double varDelta[kRichardsonN],derRichY[kRichardsonN],derRichZ[kRichardsonN];
  //
  for (int icl=0;icl<kRichardsonN;icl++) { // calculate kRichardsonN variations with del, del/2, del/4...
    double resYVP=0,resYVP=0,resZVP=0,resZVN=0;
    pnt->GetResidualsDiag(probD[(icl<<1)+0].GetParameter(), resYVP, resZVP); // variation with +delta
    pnt->GetResidualsDiag(probD[(icl<<1)+1].GetParameter(), resYVN, resZVN); // variation with -delta
    derRichY[icl] = 0.5*(resYVP-resYVN)/varDelta[icl];   // 2-point symmetric derivatives
    derRichZ[icl] = 0.5*(resZVP-resZVN)/varDelta[icl];
  }
  derY = RichardsonExtrap(derRichY,kRichardsonOrd);   // dY/dPar
  detZ = RichardsonExtrap(derRichZ,kRichardsonOrd);  // dZ/dPar
  //
}
