#include <TMath.h>
#include <TString.h>
#include "AliAlgPoint.h"
#include "AliAlgAux.h"

using namespace AliAlgAux;
using namespace TMath;

//_____________________________________
AliAlgPoint::AliAlgPoint()
  :fMinLocVarID(0)
  ,fMaxLocVarID(0)
  ,fDetID(-1)
  ,fSID(-1)
  ,fAlphaSens(0)
  ,fXSens(0)
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
  memset(fMatCorrPar,0,5*sizeof(float));
  memset(fMatCorrCov,0,5*sizeof(float));
  memset(fMatDiag,0,5*5*sizeof(float));
  //
  memset(fTrParamWS,0,5*sizeof(double));
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
void AliAlgPoint::Print(Option_t* opt) const
{
  // print
  TString opts = opt;
  opts.ToLower();
  printf("%cDet%d SID:%4d Alp:%+.3f X:%+9.4f Meas:%s Mat: ",IsInvDir() ? '*':' ',
	 GetDetID(),GetSID(),GetAlphaSens(),GetXSens(),ContainsMeasurement() ? "ON":"OFF");
  if (!ContainsMaterial()) printf("OFF\n");
  else printf("x2X0: %.4f x*rho: %.4f | pars:[%3d:%3d)\n",GetX2X0(),GetXTimesRho(),GetMinLocVarID(),GetMaxLocVarID());
  //
  if (opts.Contains("meas") && ContainsMeasurement()) {
    printf("  MeasPnt: Xtr: %+9.4f Ytr: %+8.4f Ztr: %+9.4f | ErrYZ: %+e %+e %+e\n",
	   GetXTracking(),GetYTracking(),GetZTracking(),
	   fErrYZTracking[0],fErrYZTracking[1],fErrYZTracking[2]);
  }
  //
  if (opts.Contains("mat") && ContainsMaterial()) {
    printf("  MatCorr Par: %+.4e %+.4e %+.4e %+.4e %+.4e\n", 
	   fMatCorrPar[0], fMatCorrPar[1], fMatCorrPar[2], fMatCorrPar[3], fMatCorrPar[4]);
    printf("  MatCorr cov: %+.4e %+.4e %+.4e %+.4e %+.4e\n", 
	   fMatCorrCov[0], fMatCorrCov[1], fMatCorrCov[2], fMatCorrCov[3], fMatCorrCov[4]);
  }
  //
  if (opts.Contains("diag") && ContainsMaterial()) {
    printf("  Matrix for Mat.corr.errors diagonalization:\n");
    int npar = GetNMatPar();
    for (int i=0;i<npar;i++) {
      for (int j=0;j<npar;j++) printf("%+.4e ",fMatDiag[i][j]); 
      printf("\n");
    }
  }
  //
  if (opts.Contains("ws")) { // printf track state at this point stored during residuals calculation
    printf("Local Track: "); 
    for (int i=0;i<5;i++) printf("%+.3e ",fTrParamWS[i]); 
    printf("\n");
  }
  //
}

//_____________________________________
void AliAlgPoint::Clear(Option_t* )
{
  // reset the point
  ResetBit(0xfffffff);
  fMaxLocVarID = -1;
  fDetID = -1;
  fSID   = -1;
}

//__________________________________________________________________
Int_t AliAlgPoint::Compare(const TObject* b) const
{
  // sort points in direction opposite to track propagation, i.e.
  // 1) for tracks from collision: range in decreasing tracking X
  // 2) for cosmic tracks: upper leg (pnt->IsInvDir()==kTRUE) ranged in increasing X
  //                       lower leg - in decreasing X
  AliAlgPoint* pnt = (AliAlgPoint*)b;
  double x = GetXPoint();
  double xp = pnt->GetXPoint();
  if (!IsInvDir()) { // track propagates from low to large X via this point
    if (!pnt->IsInvDir()) { // via this one also
      return x>xp ? -1:1;   
    }
    else return -1; // range points of lower leg 1st
  }
  else { // this point is from upper cosmic leg: track propagates from large to low X
    if (pnt->IsInvDir()) { // this one also
      return x>xp ? 1:-1;
    }
    else return 1; // other point is from lower leg
  }
  //
}

//__________________________________________________________________
void AliAlgPoint::GetXYZGlo(Double_t r[3]) const
{
  // position in lab frame
  double cs=TMath::Cos(fAlphaSens);
  double sn=TMath::Sin(fAlphaSens);
  double x=GetXPoint(); 
  r[0] = x*cs - GetYTracking()*sn; 
  r[1] = x*sn + GetYTracking()*cs;
  r[2] = GetZTracking();
  //
}

//__________________________________________________________________
Double_t AliAlgPoint::GetPhiGlo() const
{
  // phi angle (-pi:pi) in global frame
  double xyz[3];
  GetXYZGlo(xyz);
  return ATan2(xyz[1],xyz[0]);
}

//__________________________________________________________________
Int_t AliAlgPoint::GetAliceSector() const
{
  // get global sector ID corresponding to this point phi
  return Phi2Sector(GetPhiGlo());  
}

//__________________________________________________________________
void AliAlgPoint::SetMatCovDiagonalizationMatrix(const TMatrixD& d)
{
  // save non-sym matrix for material corrections cov.matrix diagonalization
  // (actually, the eigenvectors are stored)
  int sz = d.GetNrows();
  for (int i=sz;i--;) for (int j=sz;j--;) fMatDiag[i][j] = d(i,j);
}

//__________________________________________________________________
void AliAlgPoint::SetMatCovDiag(const TVectorD& v)
{
  // save material correction diagonalized matrix 
  // (actually, the eigenvalues are stored w/o reordering them to correspond to the 
  // AliExternalTrackParam variables)
  for (int i=v.GetNrows();i--;) fMatCorrCov[i] = v(i);
}

//__________________________________________________________________
void AliAlgPoint::UnDiagMatCorr(const double* diag, double* nodiag) const
{
  // transform material corrections from the frame diagonalizing the errors to point frame
  // nodiag = fMatDiag * diag
  int np = GetNMatPar();
  for (int ip=np;ip--;) {
    double v = 0;
    for (int jp=np;jp--;) v += fMatDiag[ip][jp]*diag[jp];
    nodiag[ip] = v;
  }
  //
}

//__________________________________________________________________
void AliAlgPoint::DiagMatCorr(const double* nodiag, double* diag) const
{
  // transform material corrections from the AliExternalTrackParam frame to
  // the frame diagonalizing the errors
  // diag = fMatDiag^T * nodiag
  int np = GetNMatPar();
  for (int ip=np;ip--;) {
    double v = 0;
    for (int jp=np;jp--;) v += fMatDiag[jp][ip]*nodiag[jp];
    diag[ip] = v;
  }
  //
}
