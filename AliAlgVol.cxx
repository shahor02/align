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
#include "AliAlignObjParams.h"

ClassImp(AliAlgVol)

//_________________________________________________________
AliAlgVol::AliAlgVol(const char* symname) :
  TNamed(symname,"")
  ,fFirstParOffs(-1)
  ,fParOffs(0)
  ,fDOF(0)
  ,fNDOF(0)
  ,fNDOFGeomFree(0)
  ,fNDOFFree(0)

  ,fParent(0)
  ,fChildren(0)

  ,fNProcPoints(0)
  ,fParVals(0)
  ,fParErrs(0)
  ,fParCstr(0)
  //
  ,fMatG2L()
  ,fMatG2LIdeal()
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
  // The calculation is done as deltaM = G2L * DeltaM * G2L^-1 l2g
  // where DeltaM is the loca variation matrix of parent volume, 
  // G2L parent local-global matrix and l2g is sensor local-global matrix
  // 
  parent->GetDeltaMatrixLoc(deltaM,delta);
  deltaM.MultiplyLeft(&parent->GetMatrixG2L());
  if (relMat) deltaM.Multiply(relMat);
  else {
    deltaM.Multiply(&parent->GetMatrixG2L().Inverse());
    deltaM.Multiply(&GetMatrixG2L());
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
