#include <TMath.h>
#include "AliAlgPoint.h"
#include "AliAlgAux.h"

using namespace AliAlgAux;
using namespace TMath;

//_____________________________________
AliAlgPoint::AliAlgPoint()
  :fMaxLocVarID(0)
  ,fDetID(-1)
  ,fSID(-1)
  ,fAlphaSens(0)
  ,fXSens(0)
  ,fCosDiagErr(0)
  ,fSinDiagErr(0)
  ,fX2X0(0)
  ,fXTimesRho(0)
  ,fMSSigTheta2(0)
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
      if (Abs(dfd)<kDiagToler) phi = fErrYZTracking[1]>0 ? (Pi()*0.25) : (Pi()*0.75);
      else                            phi = 0.5*ATan2(fErrYZTracking[1],dfd);
      //
      fCosDiagErr = Cos(phi);
      fSinDiagErr = Sin(phi);
      //
      double det = dfd*dfd + fErrYZTracking[1]*fErrYZTracking[1];
      det = det>0 ? Sqrt(det) : 0;
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
  printf("Det%d SID:%4d Alp:%+.3f X:%+9.4f",GetDetID(),GetSID(),GetAlphaSens(),GetXSens());
  if (ContainsMeasurement()) {
    printf(" Meas: Xtr: %+8.4f Ytr: %+8.4f Ztr: %+8.4f",
	   GetXTracking(),GetYTracking(),GetZTracking());
  }
  if (ContainsMaterial()) {
    printf(" Mat: X/X0: %.4f | X*rho: %.4f\n",fX2X0,fXTimesRho);
  }
  //
  printf("\n");
}

//_____________________________________
void AliAlgPoint::Clear(Option_t* )
{
  // reset the point
  ResetBit(0xfffffff);
  fDetID = -1;
  fSID   = -1;
}
