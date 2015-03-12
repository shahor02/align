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

ClassImp(AliAlgVol)

//------------------------------------------------------------
AliAlgVol::AliAlgVol(const char* name,const char* title) :
  TNamed(name,title)
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
{
  // def c-tor
}

//------------------------------------------------------------
AliAlgVol::~AliAlgVol()
{
  // d-tor
  delete[] fParOffs;
  delete[] fParVals;
  delete[] fParErrs;
  delete[] fParCstr;
  delete fChildren;
}

//-------------------------------------------------------------
TGeoHMatrix *AliAlgVol::GetSensitiveVolumeModifiedMatrix(UShort_t voluid, const Double_t *delta,Bool_t local)
{
  // modify the original TGeoHMatrix of the sensitive module 'voluid' according
  // with a delta transform. applied to the supermodule matrix
  // return NULL if error

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager)  return NULL;

  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  tr[0]=delta[0]; // in centimeter
  tr[1]=delta[1]; 
  tr[2]=delta[2];
  ang[0]=delta[3]; // psi   (X)  in deg
  ang[1]=delta[4]; // theta (Y)
  ang[2]=delta[5]; // phi   (Z)
  //
  static AliAlignObjParams tempAlignObj;
  tempAlignObj.SetRotation(ang[0],ang[1],ang[2]);
  tempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  TGeoHMatrix hm;
  tempAlignObj.GetMatrix(hm);
  //printf("\n0: delta matrix\n");hm.Print();

  // 1) start setting fSensVolModif = fSensVol
  if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
  //
  if (local) {
    // 2) set fSensVolModif = SensVolRel
    fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() );
    // 3) multiply left by delta
    fSensVolModifMatrix->MultiplyLeft( &hm );
    // 4) multiply left by fMatrix
    fSensVolModifMatrix->MultiplyLeft( fMatrix );
  }
  else fSensVolModifMatrix->MultiplyLeft( &hm );
  //
  return fSensVolModifMatrix;
}
