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

#include "AliAlgVol.h"
#include "AliAlgSens.h"
#include "AliAlignObjParams.h"
#include "AliGeomManager.h"
#include "AliAlgAux.h"
#include "AliLog.h"
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

ClassImp(AliAlgVol)

using namespace TMath;

//_________________________________________________________
AliAlgVol::AliAlgVol(const char* symname) :
  TNamed(symname,"")
  ,fX(0)
  ,fAlp(0)
  ,fFirstParOffs(-1)
  ,fParOffs(0)
  ,fDOF(0)
  ,fNDOF(0)
  ,fNDOFGeomFree(0)
  ,fNDOFFree(0)
  //
  ,fParent(0)
  ,fChildren(0)

  ,fNProcPoints(0)
  ,fParVals(0)
  ,fParErrs(0)
  ,fParCstr(0)
  //
  ,fMatL2G()
  ,fMatL2GOrig()
  ,fMatT2L()
{
  // def c-tor
}

//_________________________________________________________
AliAlgVol::~AliAlgVol()
{
  // d-tor
  delete[] fParOffs;
  delete[] fParVals;
  delete[] fParErrs;
  delete[] fParCstr;
  delete fChildren;
}

//_________________________________________________________
void AliAlgVol::GetDeltaMatrixLoc(TGeoHMatrix& deltaM, const Double_t *delta) const
{
  // prepare delta matrix for the volume from its
  // local delta vector (TGeoMatrix convension): dx,dy,dz,phi,theta,psi
  const double *tr=&delta[0],*rt=&delta[3]; // translation(cm) and rotation(degree) 
  AliAlignObjParams tempAlignObj;
  tempAlignObj.SetRotation(rt[0],rt[1],rt[2]);
  tempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  tempAlignObj.GetMatrix(deltaM);
  //  
}

//_________________________________________________________
void AliAlgVol::GetDeltaMatrixLoc(const AliAlgVol* parent, TGeoHMatrix& deltaM, 
				  const Double_t *delta, const TGeoHMatrix* relMat) const
{
  // prepare delta matrix for the child volume from 
  // local delta vector of the parent (TGeoMatrix convension): dx,dy,dz,phi,theta,psi
  // since it requires calculation of child->parent transition matrix, it can be provided
  // as an optional parameter relMat
  // The calculation is done as deltaM = L2G * DeltaM * L2G^-1 l2g
  // where DeltaM is the loca variation matrix of parent volume, 
  // L2G parent local-global matrix and l2g is sensor local-global matrix
  // 
  parent->GetDeltaMatrixLoc(deltaM,delta);
  deltaM.MultiplyLeft(&parent->GetMatrixL2G());
  if (relMat) deltaM.Multiply(relMat);
  else {
    deltaM.Multiply(&parent->GetMatrixL2G().Inverse());
    deltaM.Multiply(&GetMatrixL2G());
  }
  //  
}

//_________________________________________________________
Int_t AliAlgVol::CountParents() const
{
  // count parents in the chain
  int cnt = 0;
  const AliAlgVol* p = this;
  while( (p=p->GetParent()) ) cnt++;
  return cnt;
}

//____________________________________________
void AliAlgVol::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("Lev:%2d %s | %2d nodes | Effective X:%8.4f Alp:%+.4f\n",
	 CountParents(),GetSymName(),GetNChildren(),fX,fAlp);
  if (opts.Contains("mat")) { // print matrices
    printf("L2G original: "); 
    GetMatrixL2GOrig().Print();
    printf("L2G misalign: "); 
    GetMatrixL2G().Print();
    printf("T2L (fake)  : "); 
    GetMatrixT2L().Print();
 }
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixL2G()
{
  // extract from geometry L2G matrix
  const char* path = GetSymName();
  if (gGeoManager->GetAlignableEntry(path)) {
    const TGeoHMatrix* l2g = AliGeomManager::GetMatrix(path);
    if (!l2g) AliFatalF("Failed to find L2G matrix for alignable %s",path);
    SetMatrixL2G(*l2g);
  }
  else { // extract from path
    if (!gGeoManager->CheckPath(path)) AliFatalF("Volume path %s not valid!",path);
    TGeoPhysicalNode* node = (TGeoPhysicalNode*)gGeoManager->GetListOfPhysicalNodes()->FindObject(path);
    TGeoHMatrix l2g;
    if (!node) {
      AliWarningF("Attention: volume %s was not misaligned, extracting original matrix",path);
      if (!AliGeomManager::GetOrigGlobalMatrix(path,l2g)) {
	AliFatalF("Failed to find ideal L2G matrix for %s",path);
      }      
    } 
    else {
      l2g = *node->GetMatrix();
    }
    SetMatrixL2G(l2g);
  }
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixL2GOrig()
{
  // extract from geometry ideal L2G matrix
  TGeoHMatrix mtmp;
  if (!AliGeomManager::GetOrigGlobalMatrix(GetSymName(),mtmp)) 
    AliFatalF("Failed to find ideal L2G matrix for %s",GetSymName());
  SetMatrixL2GOrig(mtmp);
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixT2L()
{
  // for non-sensors we define the fake tracking frame with the alpha angle being
  // the average angle of centers of its children
  //
  double tot[3]={0,0,0},loc[3]={0,0,0},glo[3];
  int nch = GetNChildren();
  for (int ich=nch;ich--;) {
    AliAlgVol* vol = GetChild(ich);
    vol->GetMatrixL2GOrig().LocalToMaster(loc,glo);
    for (int j=3;j--;) tot[j] += glo[j];
  }
  if (nch) for (int j=3;j--;) tot[j] /= nch;
  //
  fAlp = TMath::ATan2(tot[1],tot[0]);
  AliAlgAux::BringToPiPM(fAlp);
  //
  fX = TMath::Sqrt(tot[0]*tot[0]+tot[1]*tot[1]);
  //
  // 1st create Tracking to Global matrix
  fMatT2L.Clear();
  fMatT2L.RotateZ(fAlp*RadToDeg());
  fMatT2L.SetDx(fX);
  // then convert it to Tracking to Local  matrix
  fMatT2L.MultiplyLeft(&GetMatrixL2GOrig().Inverse());
  //
}

//____________________________________________
void AliAlgVol::SetMatrixT2L(const TGeoHMatrix& m)
{
  // set the T2L matrix and define tracking frame
  // Note that this method is used for the externally set matrices
  // (in case of sensors). For other volumes the tracking frame and matrix
  // is defined in the PrepareMatrixT2L method
  fMatT2L = m;
  SetTrackingFrame();
}

//__________________________________________________________________
void AliAlgVol::SetTrackingFrame()
{
  // Define tracking frame of the sensor
  // This method should be implemented for sensors, which receive the T2L
  // matrix from the geometry
  AliErrorF("Volume %s was supposed to implement its own method",GetName());
}
