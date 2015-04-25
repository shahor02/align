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

const char* AliAlgVol::fgkFrameName[AliAlgVol::kNVarFrames] = {"LOC","TRA"};
UInt_t      AliAlgVol::fgDefGeomFree = 
  BIT(AliAlgVol::kDOFTX)|BIT(AliAlgVol::kDOFTY)|BIT(AliAlgVol::kDOFTZ)|
  BIT(AliAlgVol::kDOFPH)|BIT(AliAlgVol::kDOFTH)|BIT(AliAlgVol::kDOFPS);

//_________________________________________________________
AliAlgVol::AliAlgVol(const char* symname) :
  TNamed(symname,"")
  ,fVarFrame(kTRA)
  ,fX(0)
  ,fAlp(0)
  ,fNDOFTot(0)
  ,fDOF(0)
  ,fNDOFGeomFree(0)
  ,fNDOFFree(0)
  //
  ,fParent(0)
  ,fChildren(0)
  //
  ,fNProcPoints(0)
  ,fFirstParOffs(-1)
  ,fParOffs(0)
  ,fParVals(0)
  ,fParErrs(0)
  ,fParCstr(0)
  //
  ,fMatL2G()
  ,fMatL2GOrig()
  ,fMatT2L()
{
  // def c-tor
  if (symname) { // real volumes have at least geometric degrees of freedom
    SetNDOFTot(kNDOFGeom);
  }
  SetFreeDOFPattern(fgDefGeomFree);
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
void AliAlgVol::Delta2Matrix(TGeoHMatrix& deltaM, const Double_t *delta) const
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

//__________________________________________________________________
void AliAlgVol::GetModifiedMatrixL2G(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare for modified L2G matrix from its current matrix 
  // by applying local delta, i.e. glo = L2G*loc = L2G*delta*loc -> L2G = L2G*delta
  Delta2Matrix(matMod, delta);
  matMod.MultiplyLeft(&GetMatrixL2G());
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of LOCAL frame:
  // tra' = tau*tra = tau*T2L^-1*loc = T2L^-1*loc' = T2L^-1*delta*loc 
  // tau = T2L^-1*delta*T2L
  Delta2Matrix(matMod, delta);
  matMod.Multiply(&GetMatrixT2L());
  matMod.MultiplyLeft(&GetMatrixT2L().Inverse());
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
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
void AliAlgVol::GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of the same TRACKING frame:
  // tra' = tau*tra
  Delta2Matrix(matMod, delta);
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of TRACKING frame of its PARENT; 
  // The relMat is matrix for transformation from child to parent frame: TRA = relMat*tra
  // (see DPosTraDParTRA)
  //
  // tra' = tau*tra = tau*relMat^-1*TRA = relMat^-1*TAU*TRA = relMat^-1*TAU*relMat*tra
  // tau = relMat^-1*TAU*relMat
  Delta2Matrix(matMod, delta); // TAU
  matMod.Multiply(&relMat);
  matMod.MultiplyLeft(&relMat.Inverse());
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
  printf("     DOFs: Tot: %d Free: %d (offs: %5d) Geom: %d {",fNDOFTot,fNDOFFree,fFirstParOffs,fNDOFGeomFree);
  for (int i=0;i<kNDOFGeom;i++) printf("%d",IsFreeDOFGeom(DOFGeom_t(i)) ? 1:0); 
  printf("} in %s frame\n",fgkFrameName[fVarFrame]);
  //
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

//__________________________________________________________________
Int_t AliAlgVol::InitDOFs()
{
  // initialize degrees of freedom
  //
  // index allowed DOFs
  fNDOFFree = fNDOFGeomFree = 0;
  for (int i=0;i<kNDOFGeom;i++) { // start with standard DOF geoms
    fParOffs[i] = -1;
    if (!IsFreeDOFGeom(DOFGeom_t(i))) continue;
    fParOffs[i] = fNDOFGeomFree++;
  }
  fNDOFFree += fNDOFGeomFree;
  // TODO loop over other DOFs
  return fNDOFFree;
}

//__________________________________________________________________
void AliAlgVol::SetNDOFTot(Int_t n)
{
  // book global degrees of freedom
  if (n<kNDOFGeom) n = kNDOFGeom;
  fNDOFTot = n;
  fParOffs = new Char_t[fNDOFTot];
  fParVals = new Float_t[fNDOFTot];
  fParErrs = new Float_t[fNDOFTot];
  fParCstr = new Float_t[fNDOFTot];
  for (int i=fNDOFTot;i--;) {
    fParOffs[i] = -1;
    fParVals[i] = 0;
    fParErrs[i] = 0;
    fParCstr[i] = 0;
  }
}
