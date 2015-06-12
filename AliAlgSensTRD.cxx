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

#include "AliAlgSensTRD.h"
#include "AliTRDgeometry.h"
#include "AliAlgAux.h"
#include "AliLog.h"
ClassImp(AliAlgSensTRD)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSensTRD::AliAlgSensTRD(const char* name,Int_t vid, Int_t iid, Int_t isec) 
  :AliAlgSens(name,vid,iid)
  ,fSector(isec)
{
  // def c-tor
}

//_________________________________________________________
AliAlgSensTRD::~AliAlgSensTRD()
{
  // d-tor
}
/*
//__________________________________________________________________
void AliAlgSensTRD::SetTrackingFrame()
{
  // define tracking frame of the sensor: just rotation by sector angle
  fAlp = Sector2Alpha(fSector);
  fX = 0;
}
*/

//____________________________________________
void AliAlgSensTRD::PrepareMatrixT2L()
{
  // extract from geometry T2L matrix
  double alp = Sector2Alpha(fSector);
  double loc[3]={0,0,0},glo[3];
  GetMatrixL2GIdeal().LocalToMaster(loc,glo);
  double x = Sqrt(glo[0]*glo[0]+glo[1]*glo[1]);
  TGeoHMatrix t2l;
  t2l.SetDx(x);
  t2l.RotateZ(alp*RadToDeg());
  t2l.MultiplyLeft(&GetMatrixL2GIdeal().Inverse());
  /*
  const TGeoHMatrix* t2l = AliGeomManager::GetTracking2LocalMatrix(GetVolID());
  if (!t2l) {
    Print("long");
    AliFatalF("Failed to find T2L matrix for VID:%d %s",GetVolID(),GetSymName());
  }
  */
  SetMatrixT2L(t2l);
  //
}
