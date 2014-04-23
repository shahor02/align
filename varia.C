#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliExternalTrackParam.h"
#include <TMath.h>
#endif


void ApplyMS(AliExternalTrackParam &trPar, double tms,double pms)
{
  //------------------------------------------------------------------------------
  // Modify track parameters in the tracking frame (dip angle lam, az. angle phi) 
  // by multiple scattering defined by polar and azumuthal scattering angles in 
  // the track collinear frame (tms and pms resp).
  // The updated direction vector in the tracking frame becomes
  //
  //  / Cos[lam]*Cos[phi] Cos[phi]*Sin[lam] -Sin[phi] \   / Cos[tms]         \
  //  | Cos[lam]*Sin[phi] Sin[lam]*Sin[phi]  Cos[phi] | x | Cos[pms]*Sin[tms]|
  //  \ Sin[lam]	       -Cos[lam]	0     /   \ Sin[pms]*Sin[tms]/
  //
  //------------------------------------------------------------------------------
  //
  double *par = (double*) trPar.GetParameter();
  //
  double snTms = TMath::Sin(tms), csTms = TMath::Cos(tms);
  double snPms = TMath::Sin(pms), csPms = TMath::Cos(pms);  
  double snPhi = par[2],  csPhi = TMath::Sqrt((1.-snPhi)*(1.+snPhi));
  double csLam = 1./TMath::Sqrt(1.+par[3]*par[3]), snLam = csLam*par[3];
  //
  double  r00 = csLam*csPhi, r01 = snLam*csPhi, &r02 = snPhi;
  double  r10 = csLam*snPhi, r11 = snLam*snPhi, &r12 = csPhi;
  double &r20 = snLam      ,&r21 = csLam;
  //
  double &v0 = csTms, v1 = snTms*csPms, v2 = snTms*snPms;
  //
  double px = r00*v0 + r01*v1 - r02*v2;
  double py = r10*v0 + r11*v1 + r12*v2;
  double pz = r20*v0 - r21*v1;
  //
  double pt = TMath::Sqrt(px*px + py*py);
  par[2] = py/pt;
  par[3] = pz/pt;
  par[4]*= csLam/pt;
  //
}
