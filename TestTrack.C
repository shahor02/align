#include "AliAlgPoint.h"
#include "AliAlgTrack.h"
#include "AliExternalTrackParam.h"
#include "AliMagF.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliTrackerBase.h"

#include "AliSymMatrix.h"
#include "AliAlgAux.h"

#include <TGeoGlobalMagField.h>
#include <TMath.h>
#include <TRandom.h>

using namespace AliAlgAux;

void Load();
Bool_t BuildEq();

AliAlgTrack* algTrack=0;
AliSymMatrix* matr = 0;
TVectorD* rhs = 0;
//
const int kNITS = 6;
double rITS[kNITS] = {3.9,7.6,15.0,23.9,38.0,43.0};
double pars[200];

const double kSclValTrk = 1;
const double kSclValMS = 1;
TVectorD var(5);

Bool_t TestTrack(const AliExternalTrackParam& trSrc)
{
  AliLog::SetClassDebugLevel("AliAlgTrack",10);
  Load();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld->SolenoidField();
  //
  AliExternalTrackParam tr0(trSrc);
  if (!tr0.RotateParamOnly( tr0.Phi() )) return kFALSE;
  AliExternalTrackParam tr1(tr0);
  double *par = (double*)tr0.GetParameter();
  par[0] += var[0] = kSclValTrk*gRandom->Gaus(0,1)*tr0.GetSigmaY2();
  par[1] += var[1] = kSclValTrk*gRandom->Gaus(0,1)*tr0.GetSigmaZ2();
  par[2] += var[2] = kSclValTrk*gRandom->Gaus(0,1)*tr0.GetSigmaSnp2();  
  par[3] += var[3] = kSclValTrk*gRandom->Gaus(0,1)*tr0.GetSigmaTgl2();  
  par[4] += var[4] = kSclValTrk*gRandom->Gaus(0,1)*tr0.GetSigma1Pt2();  
  algTrack = new AliAlgTrack();
  algTrack->AliExternalTrackParam::operator=(tr1);
  //
  if (TMath::Abs(bz)>0) algTrack->SetFieldON();
  //
  // add points
  double xyz[3];
  int nMSp = 0;
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
    pnt->SetYZErrTracking( TMath::Power(20e-4,2), 0 , TMath::Power(100e-4,2));
    pnt->SetX2X0(1e-2);
    pnt->SetXTimesRho(0);//0.05*5);
    pnt->SetUseBzOnly();
    pnt->Init();
    algTrack->AddPoint(pnt);
    if (pnt->ContainsMaterial()) {
      double x2x0 = pnt->GetX2X0();
      double p = tr0.P();
      double sigMS = 0.014*TMath::Sqrt(x2x0)/p;
      var.ResizeTo(5+(nMSp+1)*2);
      var[5+2*nMSp+0] = kSclValMS*sigMS*gRandom->Gaus(0,1);
      var[5+2*nMSp+1] = kSclValMS*sigMS*gRandom->Gaus(0,1);
      algTrack->ApplyMS(tr0, var[5+2*nMSp+0],var[5+2*nMSp+1]);
      nMSp++;
    }
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

//________________________________________________________________________________
Bool_t BuildEq()
{
  //
  matr = new AliSymMatrix(algTrack->GetNLocPar());
  rhs  = new TVectorD(algTrack->GetNLocPar());
  //
  for (int ip=algTrack->GetNPoints();ip--;) {
    AliAlgPoint* pnt = algTrack->GetPoint(ip);
    if (!pnt->ContainsMeasurement()) continue;
    for (int im=2;im--;) { // loop over 2 orthogonal measurements in each point
      double res = algTrack->GetResidual(im,ip);  // residual wrt current track
      double* deriv = algTrack->GetDerivative(im,ip); // its derivative wrt current params
      //
      for (int ipar=algTrack->GetNLocPar();ipar--;) { // loop over track parameters
	(*rhs)[ipar] -= deriv[ipar]*res/pnt->GetErrDiag(im);
	for (int jpar=ipar+1;jpar--;) {
	  (*matr)(ipar,jpar) += deriv[ipar]*deriv[jpar]/pnt->GetErrDiag(im);
	}
      }
    } // loop over 2 measurements per point
    //
  } // loop over points
  //
  // add constraints on material effect parameters
  for (int ip=algTrack->GetNPoints();ip--;) {
    AliAlgPoint* pnt = algTrack->GetPoint(ip);
    if (!pnt->ContainsMaterial()) continue;
    double sg2 = pnt->GetMSSigTheta2();
    if (IsZeroPos(sg2)) return kFALSE;
    int offs = pnt->GetMaxLocVarID();
    (*matr)(offs+AliAlgTrack::kMSTheta1,offs+AliAlgTrack::kMSTheta1) += 1./sg2;
    (*matr)(offs+AliAlgTrack::kMSTheta2,offs+AliAlgTrack::kMSTheta2) += 1./sg2;
    if (pnt->GetELossVaried()) {} // process eloss constraint
  }
  //
  return kTRUE;
}
