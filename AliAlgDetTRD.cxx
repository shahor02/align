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

#include "AliAlgDetTRD.h"
#include "AliAlgVol.h"
#include "AliAlgSensTRD.h"
#include "AliAlgSteer.h"
#include "AliGeomManager.h"
#include "AliESDtrack.h"
#include "AliTRDgeometry.h"
#include <TGeoManager.h>
#include <TMath.h>

using namespace TMath;

ClassImp(AliAlgDetTRD);

//____________________________________________
AliAlgDetTRD::AliAlgDetTRD(const char* title)
  : AliAlgDet()
  ,fNonRCCorrDzDtgl(0)
{
  // default c-tor
  SetNameTitle(AliAlgSteer::GetDetNameByDetID(AliAlgSteer::kTRD),title);
  SetDetID(AliAlgSteer::kTRD);
  fExtraErrRC[0] = fExtraErrRC[1] = 0;
  //
  // ad hoc correction
  SetNonRCCorrDzDtgl();
  SetExtraErrRC();
}


//____________________________________________
AliAlgDetTRD::~AliAlgDetTRD()
{
  // d-tor
}

//____________________________________________
void AliAlgDetTRD::DefineVolumes()
{
  // define TRD volumes
  //
  const int kNSect = 18, kNStacks = 5, kNLayers = 6;
  AliAlgSensTRD *chamb=0;
  //
  int labDet = (GetDetID()+1)*1000000;
  //  AddVolume( volTRD = new AliAlgVol("TRD") ); // no main volume, why?
  AliAlgVol *sect[kNSect] = {0};
  //
  for (int ilr=0;ilr<kNLayers;ilr++) { // layer
    for (int ich=0;ich<kNStacks*kNSect;ich++) { // chamber
      Int_t isector   = ich/AliTRDgeometry::Nstack();
      Int_t istack    = ich%AliTRDgeometry::Nstack();
      //Int_t lid       = AliTRDgeometry::GetDetector(ilr,istack,isector);
      int iid = labDet + (1+ilr)*10000 + (1+isector)*100 + (1+istack);
      const char *symname = Form("TRD/sm%02d/st%d/pl%d",isector,istack,ilr);
      if (!gGeoManager->GetAlignableEntry(symname)) continue;
      UShort_t vid    = AliGeomManager::LayerToVolUID(AliGeomManager::kTRD1+ilr,ich);
      AddVolume( chamb = new AliAlgSensTRD(symname,vid,iid/*lid*/,isector) );
      iid =  labDet + (1+isector)*100;
      if (!sect[isector]) sect[isector] = new AliAlgVol(Form("TRD/sm%02d",isector),iid);
      chamb->SetParent(sect[isector]);
    } // chamber
  } // layer
  //
  for (int isc=0;isc<kNSect;isc++) {
    if (sect[isc]) AddVolume(sect[isc]);
  }
  //
}

//____________________________________________
Bool_t AliAlgDetTRD::AcceptTrack(const AliESDtrack* trc,Int_t trtype) const 
{
  // test if detector had seed this track
  if (!CheckFlags(trc,trtype)) return kFALSE;
  if (trc->GetTRDntracklets()<fNPointsSel[trtype]) return kFALSE;
  return kTRUE;
}

//____________________________________________
AliAlgPoint* AliAlgDetTRD::TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trp, const AliESDtrack* tr)
{
  // custom correction for non-crossing points, to account for the dependence of Z bias
  // on track inclination (prob. to not cross pads)
  //
  const double kTilt = 2.*TMath::DegToRad();
  AliAlgPoint* pnt = AliAlgDet::TrackPoint2AlgPoint(pntId,trp,tr); // process in usual way
  if (!pnt) return 0;
  //
  // is it pad crrossing?
  double* errYZ = (double*) pnt->GetYZErrTracking();
  double sgYZ = errYZ[1];
  if (TMath::Abs(sgYZ)<0.01) {  // crossing
    // increase errors since the error 
    errYZ[0] += fExtraErrRC[0]*fExtraErrRC[0];
    errYZ[2] += fExtraErrRC[1]*fExtraErrRC[1];
  }
  else { // account for probability to not cross the row
    double* pYZ = (double*)pnt->GetYZTracking();
    double corrZ = 1.055*tr->GetTgl();
    pYZ[1] += corrZ; 
    pYZ[0] += corrZ*Sign(kTilt,sgYZ);  // Y and Z are correlated
  }
  return pnt;
  //
}

//__________________________________________
//____________________________________________
void AliAlgDetTRD::Print(const Option_t *opt) const
{
  // print info
  AliAlgDet::Print(opt);
  printf("Correction dZ/dTgl NonRC: %f\n",fNonRCCorrDzDtgl);
  printf("Extra error for RC tracklets: Y:%e Z:%e\n",fExtraErrRC[0],fExtraErrRC[1]);
}
 
