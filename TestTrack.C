#include "AliAlgPoint.h"
#include "AliAlgTrack.h"
#include "AliExternalTrackParam.h"
#include "AliMagF.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliTrackerBase.h"
#include <TGeoGlobalMagField.h>
#include <TMath.h>


void Load();

AliAlgTrack* algTrack=0;

const int kNITS = 6;
double rITS[kNITS] = {3.9,7.6,15.0,23.9,38.0,43.0};
double pars[200];

Bool_t TestTrack(const AliExternalTrackParam& trSrc)
{
  AliLog::SetClassDebugLevel("AliAlgTrack",10);
  Load();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld->SolenoidField();
  //
  AliExternalTrackParam tr0(trSrc);
  if (!tr0.RotateParamOnly( tr0.Phi() )) return kFALSE;
  //
  algTrack = new AliAlgTrack();
  algTrack->AliExternalTrackParam::operator=(tr0);
  //
  if (TMath::Abs(bz)>0) algTrack->SetFieldON();
  //
  // add points
  double xyz[3];
  for (int i=0;i<kNITS;i++) {
    tr0.GetXYZ(xyz);
    bz = AliTrackerBase::GetBz(xyz);
    if (!tr0.PropagateTo(rITS[i],bz)) return kFALSE;
    //
    //    tr0.Print();
    AliAlgPoint* pnt = new AliAlgPoint();
    //
    pnt->SetAlpha(tr0.GetAlpha());
    pnt->SetXYZTracking(tr0.GetX(),tr0.GetY(),tr0.GetZ());
    pnt->SetYZErrTracking( tr0.GetSigmaY2(), tr0.GetSigmaZY(), tr0.GetSigmaZ2());
    pnt->SetX2X0(1e-2);
    pnt->SetXTimesRho(0.05*5);
    pnt->SetUseBzOnly();
    pnt->Init();
    algTrack->AddPoint(pnt);
  }
  //
  algTrack->DefineDOFs();
  int nDOF   = algTrack->GetNLocPar();
  int nDOFtr = algTrack->GetNLocExtPar();
  memset(pars,0,nDOF*sizeof(double));
  for (int i=nDOFtr;i--;) pars[i] = algTrack->GetParameter()[i];
  //
  algTrack->CalcResidDeriv(pars);

  return kTRUE;
}


//________________________________________________________________________________
void Load()
{
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry("geometry.root");
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliMagF* fld = new AliMagF("bmap","bmap");
    TGeoGlobalMagField::Instance()->SetField( fld );
    TGeoGlobalMagField::Instance()->Lock();
  }
}
