#include <TMath.h>
#include "AliAlgPoint.h"
#include "AliAlgAux.h"

using namespace AliAlgAux;

//_____________________________________
AliAlgPoint::AliAlgPoint() :
  fAlpha(0)
  ,fCosDiagErr(0)
  ,fSinDiagErr(0)
  ,fX2X0(0)
  ,fXTimesRho(0)
{
  // def c-tor
  for (int i=3;i--;) {
    fXYZTracking[i] = 0;
    fErrYZTracking[i] = 0;
  }
  //
}

//_____________________________________
void AliAlgPoint::GetResidualsDiag(const double* pos, double &resU, double &resV) const
{
  // calculate residuals in the frame where the errors are diagonal, given the position
  // of the track in the standard tracking frame
  double d0=pos[0]-fXYZTracking[1], d1=pos[1]-fXYZTracking[2];
  resU = fCosDiagErr*d0 - fSinDiagErr*d1;
  resV = fSinDiagErr*d0 + fCosDiagErr*d1;
  //
}

//_____________________________________
void AliAlgPoint::Init()
{
  // compute aux info
  const double kCorrToler = 1e-6;
  const double kDiagToler = 1e-14;  
  // 
  // compute parameters of tranformation to diagonal error matrix
  if (!IsZeroPos(fErrYZTracking[0]+fErrYZTracking[2])) { 
    // 
    SetContainsMeasurement();
    //
    // is there a correlation?
    if (SmallerAbs(fErrYZTracking[1]*fErrYZTracking[1]/(fErrYZTracking[0]*fErrYZTracking[2]),kCorrToler)) {
      fCosDiagErr = 1.;
      fSinDiagErr = 0.;
      fErrDiag[0] = fErrYZTracking[0];
      fErrDiag[1] = fErrYZTracking[2];
    }
    else {
      double dfd = 0.5*(fErrYZTracking[2] - fErrYZTracking[0]);
      double phi = 0;
      // special treatment if errors are equal
      if (TMath::Abs(dfd)<kDiagToler) phi = fErrYZTracking[1]>0 ? (TMath::Pi()*0.25) : (TMath::Pi()*0.75);
      else                            phi = 0.5*TMath::ATan2(fErrYZTracking[1],dfd);
      //
      fCosDiagErr = TMath::Cos(phi);
      fSinDiagErr = TMath::Sin(phi);
      //
      double det = dfd*dfd + fErrYZTracking[1]*fErrYZTracking[1];
      det = det>0 ? (TMath::Sqrt(det)) : 0;
      double smd = 0.5*(fErrYZTracking[0] + fErrYZTracking[2]);
      fErrDiag[0] = smd + det;
      fErrDiag[1] = smd - det;
    }
  }
  //
  if ( !(IsZeroPos(fX2X0) && IsZeroPos(fXTimesRho)) ) SetContainsMaterial();
  //
}

//_____________________________________
void AliAlgPoint::Print(Option_t* ) const
{
  // print
  printf("Alp:%+.3f Xtr: %+8.4f Ytr: %+8.4f Ztr: %+8.4f\n",
	 fAlpha,GetXTracking(),GetYTracking(),GetZTracking());
  if (ContainsMaterial()) {
    printf("X/X0: %.4f | X*rho: %.4f\n",fX2X0,fXTimesRho);
  }
}

//_____________________________________
void AliAlgPoint::Clear(Option_t* )
{
  // reset the point
  ResetBit(0xfffffff);
}
