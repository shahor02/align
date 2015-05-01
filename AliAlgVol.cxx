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

#include "AliAlgVol.h"
#include "AliAlgSens.h"
#include "AliAlignObjParams.h"
#include "AliGeomManager.h"
#include "AliAlgAux.h"
#include "AliLog.h"
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

ClassImp(AliAlgVol)

using namespace TMath;
using namespace AliAlgAux;

const char* AliAlgVol::fgkFrameName[AliAlgVol::kNVarFrames] = {"LOC","TRA"};
//
UInt_t      AliAlgVol::fgDefGeomFree = 
  BIT(AliAlgVol::kDOFTX)|BIT(AliAlgVol::kDOFTY)|BIT(AliAlgVol::kDOFTZ)|
  BIT(AliAlgVol::kDOFPH)|BIT(AliAlgVol::kDOFTH)|BIT(AliAlgVol::kDOFPS);
//
const char* AliAlgVol::fgkDOFName[AliAlgVol::kNDOFGeom]={"TX","TY","TZ","PSI","THT","PHI"};

//_________________________________________________________
AliAlgVol::AliAlgVol(const char* symname) :
  TNamed(symname,"")
  ,fVarFrame(kTRA)
  ,fX(0)
  ,fAlp(0)
  ,fNDOFs(0)
  ,fDOF(0)
  ,fNDOFGeomFree(0)
  ,fNDOFFree(0)
  ,fConstrChild(kDefChildConstr)
  //
  ,fParent(0)
  ,fChildren(0)
  //
  ,fNProcPoints(0)
  ,fFirstParGloID(-1)
  ,fParVals(0)
  ,fParErrs(0)
  //
  ,fMatL2G()
  ,fMatL2GOrig()
  ,fMatT2L()
{
  // def c-tor
  if (symname) { // real volumes have at least geometric degrees of freedom
    SetNDOFs(kNDOFGeom);
  }
  SetFreeDOFPattern(fgDefGeomFree);
}

//_________________________________________________________
AliAlgVol::~AliAlgVol()
{
  // d-tor
  delete fChildren;
}

//_________________________________________________________
void AliAlgVol::Delta2Matrix(TGeoHMatrix& deltaM, const Double_t *delta) const
{
  // prepare delta matrix for the volume from its
  // local delta vector (TGeoMatrix convension): dx,dy,dz,phi,theta,psi
  const double *tr=&delta[0],*rt=&delta[3]; // translation(cm) and rotation(degree) 
  AliAlignObjParams tempAlignObj;
  tempAlignObj.SetRotation(rt[0],rt[1],rt[2]);
  tempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  tempAlignObj.GetMatrix(deltaM);
  //  
}

//__________________________________________________________________
void AliAlgVol::GetModifiedMatrixL2G(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare for modified L2G matrix from its current matrix 
  // by applying local delta, i.e. glo = L2G*loc = L2G*delta*loc -> L2G = L2G*delta
  Delta2Matrix(matMod, delta);
  matMod.MultiplyLeft(&GetMatrixL2G());
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of LOCAL frame:
  // tra' = tau*tra = tau*T2L^-1*loc = T2L^-1*loc' = T2L^-1*delta*loc 
  // tau = T2L^-1*delta*T2L
  Delta2Matrix(matMod, delta);
  matMod.Multiply(&GetMatrixT2L());
  matMod.MultiplyLeft(&GetMatrixT2L().Inverse());
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of LOCAL frame of its PARENT; 
  // The relMat is matrix for transformation from child to parent frame: LOC = relMat*loc
  //
  // tra' = tau*tra = tau*T2L^-1*loc = T2L^-1*loc' = T2L^-1*relMat^-1*Delta*relMat*loc
  // tau = (relMat*T2L)^-1*Delta*(relMat*T2L)
  Delta2Matrix(matMod, delta);
  TGeoHMatrix tmp = relMat;
  tmp *= GetMatrixT2L();
  matMod.Multiply(&tmp);
  matMod.MultiplyLeft(&tmp.Inverse());
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of the same TRACKING frame:
  // tra' = tau*tra
  Delta2Matrix(matMod, delta);
}

//__________________________________________________________________
void AliAlgVol::GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const
{
  // prepare the variation matrix tau in volume TRACKING frame by applying 
  // local delta of modification of TRACKING frame of its PARENT; 
  // The relMat is matrix for transformation from child to parent frame: TRA = relMat*tra
  // (see DPosTraDParGeomTRA)
  //
  // tra' = tau*tra = tau*relMat^-1*TRA = relMat^-1*TAU*TRA = relMat^-1*TAU*relMat*tra
  // tau = relMat^-1*TAU*relMat
  Delta2Matrix(matMod, delta); // TAU
  matMod.Multiply(&relMat);
  matMod.MultiplyLeft(&relMat.Inverse());
}

//_________________________________________________________
Int_t AliAlgVol::CountParents() const
{
  // count parents in the chain
  int cnt = 0;
  const AliAlgVol* p = this;
  while( (p=p->GetParent()) ) cnt++;
  return cnt;
}

//____________________________________________
void AliAlgVol::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("Lev:%2d %s | %2d nodes | Effective X:%8.4f Alp:%+.4f | Used Points: %d\n",
	 CountParents(),GetSymName(),GetNChildren(),fX,fAlp,fNProcPoints);
  printf("     DOFs: Tot: %d Free: %d (offs: %5d) Geom: %d {",fNDOFs,fNDOFFree,fFirstParGloID,fNDOFGeomFree);
  for (int i=0;i<kNDOFGeom;i++) printf("%d",IsFreeDOF(i) ? 1:0); 
  printf("} in %s frame.",fgkFrameName[fVarFrame]);
  if (GetNChildren()) {
    printf(" Children constraints: {");
    for (int i=0;i<kNDOFGeom;i++) printf("%d",IsChildrenDOFConstrained(i) ? 1:0); 
    printf("}");
  }
  printf("\n");
  //
  if (opts.Contains("mat")) { // print matrices
    printf("L2G original: "); 
    GetMatrixL2GOrig().Print();
    printf("L2G misalign: "); 
    GetMatrixL2G().Print();
    printf("T2L (fake)  : "); 
    GetMatrixT2L().Print();
 }
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixL2G()
{
  // extract from geometry L2G matrix
  const char* path = GetSymName();
  if (gGeoManager->GetAlignableEntry(path)) {
    const TGeoHMatrix* l2g = AliGeomManager::GetMatrix(path);
    if (!l2g) AliFatalF("Failed to find L2G matrix for alignable %s",path);
    SetMatrixL2G(*l2g);
  }
  else { // extract from path
    if (!gGeoManager->CheckPath(path)) AliFatalF("Volume path %s not valid!",path);
    TGeoPhysicalNode* node = (TGeoPhysicalNode*)gGeoManager->GetListOfPhysicalNodes()->FindObject(path);
    TGeoHMatrix l2g;
    if (!node) {
      AliWarningF("Attention: volume %s was not misaligned, extracting original matrix",path);
      if (!AliGeomManager::GetOrigGlobalMatrix(path,l2g)) {
	AliFatalF("Failed to find ideal L2G matrix for %s",path);
      }      
    } 
    else {
      l2g = *node->GetMatrix();
    }
    SetMatrixL2G(l2g);
  }
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixL2GOrig()
{
  // extract from geometry ideal L2G matrix
  TGeoHMatrix mtmp;
  if (!AliGeomManager::GetOrigGlobalMatrix(GetSymName(),mtmp)) 
    AliFatalF("Failed to find ideal L2G matrix for %s",GetSymName());
  SetMatrixL2GOrig(mtmp);
  //
}

//____________________________________________
void AliAlgVol::PrepareMatrixT2L()
{
  // for non-sensors we define the fake tracking frame with the alpha angle being
  // the average angle of centers of its children
  //
  double tot[3]={0,0,0},loc[3]={0,0,0},glo[3];
  int nch = GetNChildren();
  for (int ich=nch;ich--;) {
    AliAlgVol* vol = GetChild(ich);
    vol->GetMatrixL2GOrig().LocalToMaster(loc,glo);
    for (int j=3;j--;) tot[j] += glo[j];
  }
  if (nch) for (int j=3;j--;) tot[j] /= nch;
  //
  fAlp = TMath::ATan2(tot[1],tot[0]);
  AliAlgAux::BringToPiPM(fAlp);
  //
  fX = TMath::Sqrt(tot[0]*tot[0]+tot[1]*tot[1]);
  //
  // 1st create Tracking to Global matrix
  fMatT2L.Clear();
  fMatT2L.RotateZ(fAlp*RadToDeg());
  fMatT2L.SetDx(fX);
  // then convert it to Tracking to Local  matrix
  fMatT2L.MultiplyLeft(&GetMatrixL2GOrig().Inverse());
  //
}

//____________________________________________
void AliAlgVol::SetMatrixT2L(const TGeoHMatrix& m)
{
  // set the T2L matrix and define tracking frame
  // Note that this method is used for the externally set matrices
  // (in case of sensors). For other volumes the tracking frame and matrix
  // is defined in the PrepareMatrixT2L method
  fMatT2L = m;
  SetTrackingFrame();
}

//__________________________________________________________________
void AliAlgVol::SetTrackingFrame()
{
  // Define tracking frame of the sensor
  // This method should be implemented for sensors, which receive the T2L
  // matrix from the geometry
  AliErrorF("Volume %s was supposed to implement its own method",GetName());
}

//__________________________________________________________________
void AliAlgVol::AssignDOFs(Int_t &cntDOFs, Float_t *pars, Float_t *errs)
{
  // Assigns offset of the DOFS of this volume in global array of DOFs, attaches arrays to volumes
  //
  fParVals = pars + cntDOFs;
  fParErrs = errs + cntDOFs;
  SetFirstParGloID(cntDOFs);
  cntDOFs += fNDOFs; // increment total DOFs count
  //
  int nch = GetNChildren();   // go over childs
  for (int ich=0;ich<nch;ich++) GetChild(ich)->AssignDOFs(cntDOFs,pars,errs);
  //
  return;
}

//__________________________________________________________________
void AliAlgVol::InitDOFs()
{
  // initialize degrees of freedom
  //
  // Do we need this strict condition?
  if (GetInitDOFsDone()) AliFatalF("Something is wrong, DOFs are already initialized for %s",GetName());
  for (int i=0;i<fNDOFs;i++) if (fParErrs[i]<0 && IsZeroAbs(fParVals[i])) FixDOF(i);
  CalcFree();
  //
  SetInitDOFsDone();
  //
}

//__________________________________________________________________
void AliAlgVol::CalcFree()
{
  // calculate free dofs
  fNDOFFree = fNDOFGeomFree = 0;
  for (int i=0;i<fNDOFs;i++) {
    if (!IsFreeDOF(i)) continue;
    fNDOFFree++;
    if (i<kNDOFGeom) fNDOFGeomFree++;
  }
  //
}

//__________________________________________________________________
void AliAlgVol::SetNDOFs(Int_t n)
{
  // book global degrees of freedom
  if (n<kNDOFGeom) n = kNDOFGeom;
  fNDOFs = n;
}

//__________________________________________________________________
void AliAlgVol::AddChild(AliAlgVol* ch) 
{
  // add child volume
  if (!fChildren) {
    fChildren = new TObjArray(); 
    fChildren->SetOwner(kFALSE);
  }
  fChildren->AddLast(ch);
}

//__________________________________________________________________
void AliAlgVol::SetParVals(Int_t npar,Double_t *vl,Double_t *er)
{
  // set parameters
  if (npar>fNDOFs) AliFatalF("Volume %s has %d dofs",GetName(),fNDOFs);
  for (int i=0;i<npar;i++) {
    fParVals[i] = vl[i];
    fParErrs[i] = er ? er[i] : 0;
  }
}

//__________________________________________________________________
Bool_t AliAlgVol::IsCondDOF(Int_t i) const
{
  // is DOF free and conditioned?
  return IsFreeDOF(i) && (!IsZeroAbs(GetParVal(i)) || !IsZeroAbs(GetParErr(i)));
}

//______________________________________________________
Int_t AliAlgVol::FinalizeStat()
{
  // finalize statistics on processed points
  if (IsSensor()) return fNProcPoints;
  fNProcPoints = 0;
  for (int ich=GetNChildren();ich--;) {
    AliAlgVol* child = GetChild(ich);
    fNProcPoints += child->FinalizeStat();
  }
  return fNProcPoints;
}

//______________________________________________________
void AliAlgVol::WritePedeParamFile(FILE* flOut, const Option_t *opt) const
{
  // contribute to params template file for PEDE
  enum {kOff,kOn};
  const char* comment[2] = {"  ","! "};
  const char* kKeyParam = "parameter";
  TString opts = opt;
  opts.ToLower();
  Bool_t showDef = opts.Contains("d"); // show free DOF even if not preconditioned
  Bool_t showFix = opts.Contains("f"); // show DOF even if fixed
  Bool_t showNam = opts.Contains("n"); // show volume name even if no nothing else is prinable
  //
  // is there something to print ?
  int nCond(0),nFix(0),nDef(0);
  for (int i=0;i<fNDOFs;i++) {
    if      (!IsFreeDOF(i)) nFix++;
    else if (IsCondDOF(i))  nCond++;
    else                    nDef++;
  }  
  //
  int cmt = nCond>0 ? kOff:kOn; // do we comment the "parameter" keyword for this volume
  if (!nFix) showFix = kFALSE;
  if (!nDef) showDef = kFALSE;
  //
  if (nCond || showDef || showFix || showNam) 
    fprintf(flOut,"%s%s %s\t\tDOF/Free: %d/%d (%s) %s\n",comment[cmt],kKeyParam,comment[kOn],
	    GetNDOFs(),GetNDOFFree(),fgkFrameName[fVarFrame],GetName());
  //
  if (nCond || showDef || showFix) {
    for (int i=0;i<fNDOFs;i++) {
      cmt = kOn;
      if      (IsCondDOF(i)) cmt = kOff;  // free-conditioned : MUST print
      else if (!IsFreeDOF(i)) {if (!showFix) continue;} // Fixed: print commented if asked
      else if (!showDef) continue;  // free-unconditioned: print commented if asked
      //
      fprintf(flOut,"%s %6d %+e %+e\t! %s p%d\n",comment[cmt],GetParGloID(i),GetParVal(i),GetParErr(i),
	      IsFreeDOF(i) ? "  ":"FX",i);
    }
    fprintf(flOut,"\n");
  }
  // children volume
  int nch = GetNChildren();
  if (nch && HasChildrenConstraint()) WriteChildrenConstraints(flOut);
  //
  for (int ich=0;ich<nch;ich++) GetChild(ich)->WritePedeParamFile(flOut,opt);
  //
}

//______________________________________________________
void AliAlgVol::WriteChildrenConstraints(FILE* flOut) const
{
  // write for PEDE eventual constraints on children movement in parent frame
  //
  enum {kOff,kOn};
  const char* comment[2] = {"  ","! "};
  const char* kKeyConstr = "constraint";
  //
  int nch = GetNChildren();
  float *cstrArr = new float[nch*kNDOFGeom*kNDOFGeom];
  // we need for each children the matrix for vector transformation from children frame 
  // (in which its DOFs are defined, LOC or TRA) to this parent variation frame
  // matRel = mPar^-1*mChild
  TGeoHMatrix mPar;
  if (IsFrameTRA()) GetMatrixT2G(mPar); // tracking to global
  else              mPar = GetMatrixL2G(); // local to global
  mPar = mPar.Inverse();
  //
  float *jac = cstrArr;
  for (int ich=0;ich<nch;ich++) {
    AliAlgVol* child = GetChild(ich);
    TGeoHMatrix matRel;
    if (child->IsFrameTRA()) child->GetMatrixT2G(matRel); // tracking to global
    else                     matRel = child->GetMatrixL2G(); // local to global
    matRel.MultiplyLeft(&mPar);
    ConstrCoefGeom(matRel,jac);
    jac += kNDOFGeom*kNDOFGeom; // matrix for next slot
  }
  //
  for (int ics=0;ics<kNDOFGeom;ics++) {
    if (!IsChildrenDOFConstrained(ics)) continue;
    fprintf(flOut,"\n%s%s\t%e\t%s %s on nodes of %s\n",comment[kOff],kKeyConstr,0.0,comment[kOn],fgkDOFName[ics],GetName());
    for (int ich=0;ich<nch;ich++) { // contribution from this children DOFs to constraint 
      AliAlgVol* child = GetChild(ich);
      jac = cstrArr + kNDOFGeom*kNDOFGeom*ich;
      for (int ip=0;ip<kNDOFGeom;ip++) {
	double jv = jac[ics*kNDOFGeom+ip];
	if (child->IsFreeDOF(ip)&&!IsZeroAbs(jv)) fprintf(flOut,"%6d %+.3e\t",child->GetParGloID(ip),jv);
      } // loop over DOF's of children contributing to this constraint
      fprintf(flOut,"%s from %s\n",comment[kOn],child->GetName());
    } // loop over children
  } // loop over constraints in parent volume
  //
  delete[] cstrArr;
}

//_________________________________________________________________
void AliAlgVol::ConstrCoefGeom(const TGeoHMatrix &matRD, float* jac/*[kNDOFGeom][kNDOFGeom]*/) const
{
  // If the transformation R brings the vector from "local" frame to "master" frame as V=R*v
  // then application of the small LOCAL correction tau to vector v is equivalent to
  // aplication of correction TAU in MASTER framce V' = R*tau*v = TAU*R*v
  // with TAU = R*tau*R^-1
  // Constraining the LOCAL modifications of child volumes to have 0 total movement in their parent
  // frame is equivalent to request that sum of all TAU matrices is unity matrix, or TAU-I = 0.
  //
  // This routine calculates derivatives of the TAU-I matrix over local corrections x,y,z, psi,tht,phi
  // defining matrix TAU. In small corrections approximation the constraint is equivalent to
  // Sum_over_child_volumes{ [dTAU/dParam]_ij * deltaParam } = 0
  // for all elements ij of derivative matrices. Since only 6 out of 16 matrix params are independent,
  // we request the constraint only for  [30](X), [31](Y), [32](Z), [12](psi), [02](tht), [01](phi)
  // Choice defined by convention of AliAlgObg::Angles2Matrix (need elements ~ linear in corrections)
  //
  TGeoHMatrix matRI = matRD.Inverse();
  const int ij[kNDOFGeom][2] = {{3,0},{3,1},{3,2},{1,2},{0,2},{0,1}};
  //
  const double *rd=matRD.GetRotationMatrix(), *ri=matRI.GetRotationMatrix();  
  const double /**td=matRD.GetTranslation(),*/*ti=matRI.GetTranslation();
  //
  // the angles are in degrees, while we use sinX->X approximation...
  const double cf[kNDOFGeom] = {1,1,1,DegToRad(),DegToRad(),DegToRad()};
  //
  double dDPar[kNDOFGeom][4][4] = {
    // dDX[4][4] 
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[0],rd[3],rd[6],0}},
    // dDY[4][4]
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[1],rd[4],rd[7],0}},
    // dDZ[4][4]
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[2],rd[5],rd[8],0}},
    // dDPSI[4][4]
    {{rd[2]*ri[3]-rd[1]*ri[6],rd[2]*ri[4]-rd[1]*ri[7],rd[2]*ri[5]-rd[1]*ri[8],0},
     {rd[5]*ri[3]-rd[4]*ri[6],rd[5]*ri[4]-rd[4]*ri[7],rd[5]*ri[5]-rd[4]*ri[8],0},
     {rd[8]*ri[3]-rd[7]*ri[6],rd[8]*ri[4]-rd[7]*ri[7],rd[8]*ri[5]-rd[7]*ri[8],0},
     {rd[2]*ti[1]-rd[1]*ti[2],rd[5]*ti[1]-rd[4]*ti[2],rd[8]*ti[1]-rd[7]*ti[2],0}},
    // dDTHT[4][4]
    {{rd[0]*ri[6]-rd[2]*ri[0], rd[0]*ri[7]-rd[2]*ri[1], rd[0]*ri[8]-rd[2]*ri[2],0},
     {rd[3]*ri[6]-rd[5]*ri[0], rd[3]*ri[7]-rd[5]*ri[1], rd[3]*ri[8]-rd[5]*ri[2],0},
     {rd[6]*ri[6]-rd[8]*ri[0], rd[6]*ri[7]-rd[8]*ri[1], rd[6]*ri[8]-rd[8]*ri[2],0},
     {rd[0]*ti[2]-rd[2]*ti[0], rd[3]*ti[2]-rd[5]*ti[0], rd[6]*ti[2]-rd[8]*ti[0],0}},
    // dDPHI[4][4]
    {{rd[1]*ri[0]-rd[0]*ri[3],rd[1]*ri[1]-rd[0]*ri[4],rd[1]*ri[2]-rd[0]*ri[5],0},
     {rd[4]*ri[0]-rd[3]*ri[3],rd[4]*ri[1]-rd[3]*ri[4],rd[4]*ri[2]-rd[3]*ri[5],0},
     {rd[7]*ri[0]-rd[6]*ri[3],rd[7]*ri[1]-rd[6]*ri[4],rd[7]*ri[2]-rd[6]*ri[5],0},
     {rd[1]*ti[0]-rd[0]*ti[1],rd[4]*ti[0]-rd[3]*ti[1],rd[7]*ti[0]-rd[6]*ti[1],0}},
  };
  //
  for (int cs=0;cs<kNDOFGeom;cs++) {
    int i=ij[cs][0],j=ij[cs][1];
    for (int ip=0;ip<kNDOFGeom;ip++) jac[cs*kNDOFGeom+ip] = dDPar[ip][i][j]*cf[ip]; // [cs][ip]
  }
}
