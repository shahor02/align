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

#include "AliAlgSensTOF.h"
#include "AliAlgAux.h"
#include "AliAlgDetTOF.h"
#include "AliLog.h"
#include "AliAlgPoint.h"
#include "AliTrackPointArray.h"
#include "AliESDtrack.h"

ClassImp(AliAlgSensTOF)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSensTOF::AliAlgSensTOF(const char* name,Int_t vid, Int_t iid, Int_t isec) 
  :AliAlgSens(name,vid,iid)
  ,fSector(isec)
{
  // def c-tor
}

//_________________________________________________________
AliAlgSensTOF::~AliAlgSensTOF()
{
  // d-tor
}

/*
//__________________________________________________________________
void AliAlgSensTOF::SetTrackingFrame()
{
  // define tracking frame of the sensor: just rotation by sector angle
  fAlp = Sector2Alpha(fSector);
  fX = 0;
}
*/

//____________________________________________
void AliAlgSensTOF::PrepareMatrixT2L()
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

//____________________________________________
AliAlgPoint* AliAlgSensTOF::TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trpArr, const AliESDtrack* tr)
{
  // convert the pntId-th point to AliAlgPoint, detectors may override this method
  //
  // TOF stores in the trackpoints X,Y with alignment applied but Z w/o alignment!!!
  // -> need special treatment
  //
  AliAlgDetTOF* det = (AliAlgDetTOF*)GetDetector();
  AliAlgPoint* pnt = det->GetPointFromPool();
  pnt->SetSensor(this);
  //
  double tra[3],traId[3],loc[3],
    glo[3] = {trpArr->GetX()[pntId], trpArr->GetY()[pntId], trpArr->GetZ()[pntId]};
  const TGeoHMatrix& matL2Grec = GetMatrixL2GReco(); // local to global matrix used for reconstruction
  const TGeoHMatrix& matT2L    = GetMatrixT2L();     // matrix for tracking to local frame translation
  //
  // >>>------------- here we fix the z by emulating Misalign action in the tracking frame ------>>>
  TGeoHMatrix mT2G;
  const TGeoHMatrix &mClAlg = GetMatrixClAlg();
  GetMatrixT2G(mT2G);
  mT2G.MasterToLocal(glo,tra);
  mClAlg.MasterToLocal(tra,traId); // here we have almost idead X,Y and wrong Z
  const double *trans = mClAlg.GetTranslation();
  const double *rotmt = mClAlg.GetRotationMatrix();  
  tra[2] = trans[2] + traId[0]*rotmt[6]+traId[1]*rotmt[7]+tra[2]*rotmt[8]; //we got misaligned Z
  mT2G.LocalToMaster(tra,glo);
  // now continue as usual
  // <<<------------- here we fix the z by emulating Misalign action in the tracking frame ------<<<
  //
  // undo reco-time alignment
  matL2Grec.MasterToLocal(glo,loc); // go to local frame using reco-time matrix, here we recover ideal measurement 
  //
  matT2L.MasterToLocal(loc,traId);  // go to tracking frame 
  //
  GetMatrixClAlg().LocalToMaster(traId,tra);   // apply alignment
  //
  if (!det->GetUseErrorParam()) {
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
    const double *sysE = GetAddError(); // additional syst error
    pnt->SetYZErrTracking(hcovscl[4]+sysE[0]*sysE[0],hcovscl[5],hcovscl[8]+sysE[1]*sysE[1]);
  }
  else { // errors will be calculated just before using the point in the fit, using track info
    pnt->SetYZErrTracking(0,0,0);
    pnt->SetNeedUpdateFromTrack();
  }
  pnt->SetXYZTracking(tra[0],tra[1],tra[2]);
  pnt->SetAlphaSens(GetAlpTracking());
  pnt->SetXSens(GetXTracking());
  pnt->SetDetID(det->GetDetID());
  pnt->SetSID(GetSID());
  //
  pnt->SetContainsMeasurement();
  //
  pnt->Init();
  //
  return pnt;
  //
}
