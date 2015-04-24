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

#include "AliAlgSens.h"
#include "AliAlgAux.h"
#include "AliLog.h"
#include "AliGeomManager.h"
ClassImp(AliAlgSens)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSens::AliAlgSens(const char* name,Int_t vid, Int_t iid) 
  : AliAlgVol(name)
  ,fIntID(iid)
{
  // def c-tor
  SetVolID(vid);
  fAddError[0] = fAddError[1] = 0;
}

//_________________________________________________________
AliAlgSens::~AliAlgSens()
{
  // d-tor
}
/*
//__________________________________________________________________???
void AliAlgSens::GetDeltaMatrixTra(TGeoHMatrix& deltaM, const Double_t *delta) const
{
  // prepare delta matrix for the sensitive volume from delta vector
  // of tracking frame: shift by dx,dy,dz THEN rotation by phi,theta,psi
  // Note that this is opposite to TGeo convention (1st rotation then translation)
  //
  // RS: do we need to shift by fX ? 
  const double *tr=&delta[0],*rt=&delta[3]; // translation(cm) and rotation(degree) 
  TGeoHMatrix trM;
  TGeoRotation rtM;
  deltaM.Clear();
  deltaM.SetTranslation(tr);
  //  deltaM.SetDx(delta[0]+fX);
  //  deltaM.SetDy(delta[1]);
  //  deltaM.SetDy(delta[2]);
  rtM.SetAngles(rt[0],rt[1],rt[2]);
  deltaM.MultiplyLeft(&rtM);
  //  
}

//_________________________________________________________???
void AliAlgSens::DeltaTra2DeltaLoc(const TGeoHMatrix& deltaTra, TGeoHMatrix& deltaLoc) const
{
  // convert delta matrix for tracking frame (obtained by GetDeltaMatrixTra)
  // to delta matrix in local frame (like the one from Delta2Matrix)
  deltaLoc = deltaTra;
  deltaLoc.MultiplyLeft(&fMatT2L);
  deltaLoc.Multiply(&fMatT2L.Inverse());
  //  
}
*/

//_________________________________________________________
void AliAlgSens::DPosTraDParLOC(const double *tra, double* deriv) const
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor LOCAL frame 
  // parameters in TGeoHMatrix convention.
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],loc[3],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  for (int ip=kNDOFGeom;ip--;) delta[ip] = 0;
  for (int ip=kNDOFGeom;ip--;) {
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    // variation matrix in tracking frame for variation in sensor LOCAL frame
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(loc,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(loc,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(loc,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(loc,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParLOC(const AliAlgVol* parent, const double *tra, double* deriv) const
{
  // Jacobian of position in sensor tracking frame (tra) vs parent volume LOCAL frame parameters.
  // NO check of parentship is done!
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],loc[3],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  TGeoHMatrix matRel = parent->GetMatrixL2G().Inverse(); // 
  // this is the matrix for transition from sensor to parent volume local frames: LOC=matRel*loc
  matRel *= GetMatrixL2G();
  //
  for (int ip=kNDOFGeom;ip--;) delta[ip] = 0;
  for (int ip=kNDOFGeom;ip--;) {
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(loc,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(loc,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(loc,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(loc,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParTra(const double *tra, double* deriv) const
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor TRACKING 
  // frame parameters in TGeoHMatrix convention, i.e. the modified parameter is
  // tra' = tau*tra
  //
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],loc[3],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  for (int ip=kNDOFGeom;ip--;) delta[ip] = 0;
  for (int ip=kNDOFGeom;ip--;) {
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(loc,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(loc,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(loc,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(loc,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//__________________________________________________________________
void AliAlgSens::GetModifiedMatrixT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the sensitive module tracking2local matrix from its current T2L matrix 
  // by applying local delta of modification of LOCAL frame:
  // loc' = delta*loc = T2L'*tra = T2L'*T2L^-1*loc   ->  T2L' = delta*T2L
  Delta2Matrix(matMod, delta);
  matMod.Multiply(&GetMatrixT2L());
}

//__________________________________________________________________
void AliAlgSens::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in sensor TRACKING frame by applying 
  // local delta of modification of LOCAL frame:
  // tra' = tau*tra = tau*T2L^-1*loc = T2L^-1*loc' = T2L^-1*delta*loc 
  // tau = T2L^-1*delta*T2L
  Delta2Matrix(matMod, delta);
  matMod.Multiply(&GetMatrixT2L());
  matMod.MultiplyLeft(&GetMatrixT2L().Inverse());
}

//__________________________________________________________________
void AliAlgSens::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const
{
  // prepare the variation matrix tau in sensor TRACKING frame by applying 
  // local delta of modification of LOCAL frame of its PARENT; 
  // The relMat is matrix for transformation from child to parent frame: LOC = relMat*loc
  //
  // tra' = tau*tra = tau*T2L^-1*loc = T2L^-1*loc' = T2L^-1*relMat^-1*Delta*relMat*loc
  // tau = (relMat*T2L)^-1*Delta*(relMat*T2L)
  Delta2Matrix(matMod, delta);
  TGeoHMatrix tmp = relMat;
  tmp *= GetMatrixT2L();
  matMod.Multiply(&tmp);
  matMod.MultiplyLeft(&tmp.Inverse());
}


//__________________________________________________________________
void AliAlgSens::GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in sensor TRACKING frame by applying 
  // local delta of modification of the same TRACKING frame:
  // tra' = tau*tra
  Delta2Matrix(matMod, delta);
}

//__________________________________________________________________
void AliAlgSens::GetModifiedMatrixT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the sensitive module tracking2local matrix from its current T2L matrix 
  // by applying local delta of modification of TRACKING frame:
  // loc' = T2L'*tra = T2L*delta*tra    ->  T2L' = T2L*delta
  Delta2Matrix(matMod, delta);
  matMod.MultiplyLeft(&GetMatrixT2L());
}

//__________________________________________________________________
void AliAlgSens::AddChild(AliAlgVol*)
{
  AliFatalF("Sensor volume cannot have childs: id=%d %s",GetVolID(),GetName());
}

//__________________________________________________________________
Int_t AliAlgSens::Compare(const TObject* b) const
{
  // compare VolIDs
  return GetUniqueID()<b->GetUniqueID() ? -1 : 1;
}

//__________________________________________________________________
void AliAlgSens::SetTrackingFrame()
{
  // define tracking frame of the sensor
  AliWarningF("Generic method called for %s",GetSymName());
  double tra[3]={0},loc[3],glo[3];
  const TGeoHMatrix &t2l = GetMatrixT2L();
  const double* t = t2l.GetTranslation();
  double r = TMath::Sqrt(t[0]*t[0]+t[1]*t[1]);
  // ITS defines tracking frame with origin in sensor, others at 0
  if (r>1) tra[0] = r;
  //
  t2l.LocalToMaster(tra,loc);
  GetMatrixL2GOrig().LocalToMaster(loc,glo);
  fX = Sqrt(glo[0]*glo[0]+glo[1]*glo[1]);
  fAlp = ATan2(glo[1],glo[0]);
  AliAlgAux::BringToPiPM(fAlp);
}

//____________________________________________
void AliAlgSens::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("Lev:%2d %s VId:%6d (IntID:%4d) X:%8.4f Alp:%+.4f | Err: %.4e %.4e\n",
	 CountParents(), GetSymName(), GetVolID(), GetInternalID(),fX, fAlp, 
	 fAddError[0],fAddError[1]);
  //
  if (opts.Contains("mat")) { // print matrices
    printf("L2G original: "); 
    GetMatrixL2GOrig().Print();
    printf("L2G misalign: "); 
    GetMatrixL2G().Print();
    printf("T2L         : "); 
    GetMatrixT2L().Print();
  }
  //
}

//____________________________________________
void AliAlgSens::PrepareMatrixT2L()
{
  // extract from geometry T2L matrix
  const TGeoHMatrix* t2l = AliGeomManager::GetTracking2LocalMatrix(GetVolID());
  if (!t2l) AliFatalF("Failed to find T2L matrix for VID:%d %s",GetVolID(),GetSymName());
  SetMatrixT2L(*t2l);  
  //
}
