#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoMatrix.h>
#include <TMath.h>
#include "AliGeomManager.h"
#include "AliAlignObjParams.h"
#endif

void DPosTraDParLoc(const double *tra, double* deriv);
void GetDeltaMatrixLoc(TGeoHMatrix& deltaM, const Double_t *delta);
void GetModifiedMatrixG2L(TGeoHMatrix& matMod, const Double_t *delta);
void GetModifiedMatrixT2L(TGeoHMatrix& matMod, const Double_t *delta);


double fXTracking,alpTra,dtdp[6*3];
TGeoHMatrix fMatT2L,fMatG2T,fMatG2L;


void tstal(double xg,double yg,double zg, int vid)
{
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry("geometry.root");
  //
  fMatT2L = *AliGeomManager::GetTracking2LocalMatrix(vid);
  fMatG2L = *AliGeomManager::GetMatrix(vid);
  AliGeomManager::GetTrackingMatrix(vid,fMatG2T);
  //
  double glo[3] = {xg,yg,zg};
  //
  double tmp[3],tra[3],loc[3]={0,0,0};

  double* trans = fMatG2T.GetTranslation();
  fXTracking = TMath::Sqrt(trans[0]*trans[0]+trans[1]*trans[1]);
  loc[0] = 100; loc[1] = 0;
  fMatG2T.LocalToMaster(loc,tmp);
  alpTra = TMath::ATan2(tmp[1],tmp[0]);
  printf("Xtracking:%e alpTracking: %e\n",fXTracking, alpTra);
  //
  printf("Global:  %+e %+e %+e\n",glo[0],glo[1],glo[2]);
  fMatG2L.MasterToLocal(glo,loc);
  printf("Local :  %+e %+e %+e\n",loc[0],loc[1],loc[2]);
  //
  fMatT2L.MasterToLocal(loc,tra);
  printf("TrackF:  %+e %+e %+e\n",tra[0]+fXTracking,tra[1],tra[2]);
  //
  fMatG2T.MasterToLocal(glo,tmp);
  printf("TrackFF: %+e %+e %+e\n",tmp[0]+fXTracking,tmp[1],tmp[2]);
  //
  // derivative of tracking coordinates vs local parameters
  DPosTraDParLoc(tra, dtdp);
  for (int i=0;i<6;i++) {
    double *curd = dtdp + i*3;
    printf("dXYZtr/dP_%d:\t%+e\t%+e\t%+e\n",i,curd[0],curd[1],curd[2]);
  }
}

void DPosTraDParLoc(const double *tra, double* deriv)
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor local parameters.
  // Result is stored in array deriv as linearized matrix 6x3 
  const int kMaxParGeom = 6;
  const double kDelta[kMaxParGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kMaxParGeom],loc[3],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  fMatT2L.MasterToLocal(tra,loc);
  for (int ip=kMaxParGeom;ip--;) delta[ip] = 0;
  for (int ip=kMaxParGeom;ip--;) {
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetModifiedMatrixT2L(matMod, delta);
    matMod.LocalToMaster(loc,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetModifiedMatrixT2L(matMod, delta);
    matMod.LocalToMaster(loc,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetModifiedMatrixT2L(matMod, delta);
    matMod.LocalToMaster(loc,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetModifiedMatrixT2L(matMod, delta);
    matMod.LocalToMaster(loc,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}


//__________________________________________________________________
void GetDeltaMatrixLoc(TGeoHMatrix& deltaM, const Double_t *delta)
{
  // prepare delta matrix for the sensitive volume from 
  // local delta vector (TGeoMatrix convension): dx,dy,dz,phi,theta,psi
  const double *tr=&delta[0],*rt=&delta[3]; // translation(cm) and rotation(degree) 
  AliAlignObjParams tempAlignObj;
  tempAlignObj.SetRotation(rt[0],rt[1],rt[2]);
  tempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  tempAlignObj.GetMatrix(deltaM);
  //  
}

//__________________________________________________________________
void GetDeltaMatrixTra(TGeoHMatrix& deltaM, const Double_t *delta)
{
  // prepare delta matrix for the sensitive volume from delta vector
  // of tracking frame: shift by dx,dy,dz THEN rotation by phi,theta,psi
  // Note that this is opposite to TGeo convention (1st rotation then translation)
  //
  // RS: do we need to shift by fXTracking ? 
  const double *tr=&delta[0],*rt=&delta[3]; // translation(cm) and rotation(degree) 
  TGeoHMatrix trM,rtM;
  deltaM.Clear();
  deltaM.SetTranlation(tr);
  //  deltaM.SetDx(delta[0]+fXTracking);
  //  deltaM.SetDy(delta[1]);
  //  deltaM.SetDy(delta[2]);
  rtM.SetRotation(rt);
  deltaM.MultiplyLeft(rtM);
  //  
}

//__________________________________________________________________
void DeltaTra2DeltaLoc(const TGeoHMatrix& deltaTra, TGeoHMatrix& deltaLoc)
{
  // convert delta matrix for tracking frame (obtained by GetDeltaMatrixTra)
  // to delta matrix in local frame (like the one from GetDeltaMatrixLoc)
  deltaLoc = deltaTra;
  deltaLoc.MultiplyLoc(fMatT2L);
  //  
}



//__________________________________________________________________
void GetModifiedMatrixG2L(TGeoHMatrix& matMod, const Double_t *delta)
{
  // prepare for the sensitive module global2local matrix from its current matrix 
  // by applying local delta
  GetDeltaMatrixLoc(matMod, delta);
  matMod.MultiplyLeft(&fMatG2L);
}

//__________________________________________________________________
void GetModifiedMatrixT2L(TGeoHMatrix& matMod, const Double_t *delta)
{
  // prepare the sensitive module tracking2local matrix from its current T2L matrix 
  // by applying local delta
  GetDeltaMatrixLoc(matMod, delta);
  matMod.MultiplyLeft(&fMatT2L);
}
