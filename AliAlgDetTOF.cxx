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

#include "AliAlgDetTOF.h"
#include "AliAlgVol.h"
#include "AliAlgSensTOF.h"
#include "AliAlgSteer.h"
#include "AliGeomManager.h"
#include "AliTOFGeometry.h"
#include "AliESDtrack.h"
#include <TGeoManager.h>

ClassImp(AliAlgDetTOF);

//____________________________________________
AliAlgDetTOF::AliAlgDetTOF(const char* title)
{
  // default c-tor
  SetNameTitle(AliAlgSteer::GetDetNameByDetID(AliAlgSteer::kTOF),title);
  SetDetID(AliAlgSteer::kTOF);
}

//____________________________________________
AliAlgDetTOF::~AliAlgDetTOF()
{
  // d-tor
}

//____________________________________________
void AliAlgDetTOF::DefineVolumes()
{
  // define TOF volumes
  //
  const int kNSect = 18, kNStrips = AliTOFGeometry::NStripA()+2*AliTOFGeometry::NStripB()+2*AliTOFGeometry::NStripC();
  int labDet = GetDetLabel();
  AliAlgSensTOF *strip=0;
  //
  //  AddVolume( volTOF = new AliAlgVol("TOF") ); // no main volume, why?
  AliAlgVol *sect[kNSect] = {0};
  //
  for (int isc=0;isc<kNSect;isc++) {
    int iid = labDet + (1+isc)*100;
    AddVolume(sect[isc] = new AliAlgVol(Form("TOF/sm%02d",isc),iid));
  }
  //
  int cnt = 0;
  for (int isc=0;isc<kNSect;isc++) {
    for (int istr=1;istr<=kNStrips;istr++) { // strip
      int iid = labDet + (1+isc)*100 + (1+istr);
      int vid = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF, cnt++);
      const char *symname = Form("TOF/sm%02d/strip%02d",isc,istr);
      if (!gGeoManager->GetAlignableEntry(symname)) continue;
      AddVolume( strip=new AliAlgSensTOF(symname,vid,iid,isc) );
      strip->SetParent(sect[isc]);
    } // strip
  } // layer
  //
}

//____________________________________________
Bool_t AliAlgDetTOF::AcceptTrack(const AliESDtrack* trc,Int_t trtype) const 
{
  // test if detector had seed this track
  return CheckFlags(trc,trtype);
}

//____________________________________________
AliAlgPoint* AliAlgDetTOF::TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trpArr, const AliESDtrack* tr)
{
  // convert the pntId-th point to AliAlgPoint, detectors may override this method
  //
  // TOF stores in the trackpoints X,Y with alignment applied but Z w/o alignment!!!
  // -> need special treatment
  //
  // convert to detector tracking frame
  UShort_t vid = trpArr->GetVolumeID()[pntId];
  Int_t sid = VolID2SID(vid); // sensor index within the detector
  if (!sid<0) return 0;
  AliAlgSens* sens = GetSensor(sid);
  if (sens->GetSkip()) return 0;
  AliAlgPoint* pnt = GetPointFromPool();
  pnt->SetSensor(sens);
  //
  double tra[3],traId[3],loc[3],glo[3] = {trpArr->GetX()[pntId], trpArr->GetY()[pntId], trpArr->GetZ()[pntId]};
  //
  const TGeoHMatrix& matL2Grec = sens->GetMatrixL2GReco(); // local to global matrix used for reconstruction
  //const TGeoHMatrix& matL2G    = sens->GetMatrixL2G();     // local to global orig matrix used as a reference 
  const TGeoHMatrix& matT2L    = sens->GetMatrixT2L();     // matrix for tracking to local frame translation
  //
  // >>>------------- here we fix the z by emulating Misalign action in the tracking frame ------>>>
  TGeoHMatrix mT2G;
  const TGeoHMatrix &mClAlg = sens->GetMatrixClAlg();
  sens->GetMatrixT2G(mT2G);
  mT2G.MasterToLocal(glo,tra);
  mClAlg.MasterToLocal(tra,traId); // here we have almost idead X,Y and wrong Z
  const double *trans = mClAlg.GetTranslation();
  const double *rotmt = mClAlg.GetRotationMatrix();  
  tra[2] = trans[2] + traId[0]*rotmt[6]+traId[1]*rotmt[7]+tra[2]*rotmt[8]; //we got misaligned Z
  mT2G.LocalToMaster(tra,glo);
  // now continue as usual
  // <<<------------- here we fix the z by emulating Misalign action in the tracking frame ------<<<

  // undo reco-time alignment
  matL2Grec.MasterToLocal(glo,loc); // go to local frame using reco-time matrix 
  matT2L.MasterToLocal(loc,traId); // go to tracking frame 
  //
  sens->GetMatrixClAlg().LocalToMaster(traId,tra);   // apply alignment
  //
  if (!fUseErrorParam) {
    // convert error
    TGeoHMatrix hcov;
    Double_t hcovel[9];
    const Float_t *pntcov = trpArr->GetCov()+pntId*6; // 6 elements per error matrix
    hcovel[0] = double(pntcov[0]);
    hcovel[1] = double(pntcov[1]);
    hcovel[2] = double(pntcov[2]);
    hcovel[3] = double(pntcov[1]);
    hcovel[4] = double(pntcov[3]);
    hcovel[5] = double(pntcov[4]);
    hcovel[6] = double(pntcov[2]);
    hcovel[7] = double(pntcov[4]);
    hcovel[8] = double(pntcov[5]);
    hcov.SetRotation(hcovel);
    hcov.Multiply(&matL2Grec);                
    hcov.MultiplyLeft(&matL2Grec.Inverse());    // errors in local frame
    hcov.Multiply(&matT2L);
    hcov.MultiplyLeft(&matT2L.Inverse());       // errors in tracking frame
    //
    Double_t *hcovscl = hcov.GetRotationMatrix();
    const double *sysE = sens->GetAddError(); // additional syst error
    pnt->SetYZErrTracking(hcovscl[4]+sysE[0]*sysE[0],hcovscl[5],hcovscl[8]+sysE[1]*sysE[1]);
  }
  else { // errors will be calculated just before using the point in the fit, using track info
    pnt->SetYZErrTracking(0,0,0);
    pnt->SetNeedUpdateFromTrack();
  }
  pnt->SetXYZTracking(tra[0],tra[1],tra[2]);
  pnt->SetAlphaSens(sens->GetAlpTracking());
  pnt->SetXSens(sens->GetXTracking());
  pnt->SetDetID(GetDetID());
  pnt->SetSID(sid);
  //
  pnt->SetContainsMeasurement();
  //
  pnt->Init();
  //
  return pnt;
  //
}
