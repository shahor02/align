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
#include "AliAlgAux.h"
#include "AliLog.h"
ClassImp(AliAlgSensTRD)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSensTRD::AliAlgSensTRD(const char* name,Int_t vid, Int_t iid) 
  : AliAlgSens(name,vid,iid)
{
  // def c-tor
}

//_________________________________________________________
AliAlgSensTRD::~AliAlgSensTRD()
{
  // d-tor
}

//__________________________________________________________________
void AliAlgSensTRD::SetTrackingFrame()
{
  // define tracking frame of the sensor
  double tra[3]={0},loc[3],glo[3];
  // TRD defines tracking frame with origin in sensor, others at 0
  GetMatrixT2L().LocalToMaster(tra,loc);
  GetMatrixL2GOrig().LocalToMaster(loc,glo);
  fX = Sqrt(glo[0]*glo[0]+glo[1]*glo[1]);
  fAlp = ATan2(glo[1],glo[0]);
  AliAlgAux::BringTo02Pi(fAlp);
}

