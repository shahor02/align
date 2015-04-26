#include "AliAlgMPRecord.h"
#include "AliAlgAux.h"
#include "AliAlgTrack.h"
#include "AliLog.h"
#include <TMath.h>


using namespace TMath;
using namespace AliAlgAux;

//_________________________________________________________
AliAlgMPRecord::AliAlgMPRecord()
  :fTrackID(0)
  ,fTimeStamp(0)
  ,fNResid(0)
  ,fNVarLoc(0)
  ,fNVarGlo(0)
  ,fNDLocTot(0)
  ,fNDGloTot(0)
  ,fNDLoc(0)
  ,fNDGlo(0)
  ,fResid(0)
  ,fResErr(0)
  ,fIDLoc(0)
  ,fIDGlo(0)
  ,fDLoc(0)
  ,fDGlo(0)
   //
  ,fNResidBook(0)
  ,fNDLocTotBook(0)
  ,fNDGloTotBook(0)
{
  // def c-tor
  
}

//_________________________________________________________
AliAlgMPRecord::~AliAlgMPRecord()
{
  // d-tor
  delete[] fNDLoc;
  delete[] fNDGlo;
  delete[] fResid;
  delete[] fResErr;
  delete[] fIDLoc;
  delete[] fIDGlo;
  delete[] fDLoc;
  delete[] fDGlo;
}

//_________________________________________________________
Bool_t AliAlgMPRecord::FillTrack(const AliAlgTrack* trc)
{
  // fill track info
  //
  if (!trc->GetDerivDone()) {
    AliError("Track derivatives are not yet evaluated");
    return kFALSE;
  }

  fNVarLoc = trc->GetNLocPar();    // number of local degrees of freedom in the track
  fNResid = 0;
  fNDLocTot = 0;
  fNDGloTot = 0;

  // 1) check sizes for buffers, expand if needed
  int np = trc->GetNPoints();
  int nres  = 0;
  int nlocd = 0;
  int nglod = 0;
  for (int ip=np;ip--;) {
    AliAlgPoint* pnt = trc->GetPoint(ip);
    if (pnt->ContainsMeasurement()) {
      nres    += 2;                     // every point has 2 residuals
      nlocd   += fNVarLoc+fNVarLoc;     // each residual has max fNVarLoc local derivatives
      nglod   += fNVarGlo+fNVarGlo;     // and max fNVarGlo global derivatives
    }
    if (pnt->ContainsMaterial()) {
      int nmatpar = pnt->GetNMatPar();
      nres    += nmatpar;               // each point with materials has nmatpar fake residuals 
      nlocd   += nmatpar;               // and nmatpar non-0 local derivatives (orthogonal)
    }
  }
  //
  Resize(nres,nlocd,nglod);
  int nParETP = trc->GetNLocExtPar();   // numnber of local parameters for reference track param
  //
  for (int ip=0;ip<np;ip++) {
    AliAlgPoint* pnt = trc->GetPoint(ip);
    if (pnt->ContainsMeasurement()) {
      for (int idim=0;idim<2;idim++) { // 2 dimensional orthogonal measurement
	fNDGlo[fNResid] = 0;
	fResid[fNResid]  = trc->GetResidual(idim,ip);
	fResErr[fNResid] = Sqrt(pnt->GetErrDiag(idim));
	double* deriv  = trc->GetDResDLoc(idim,ip);  // array of Dresidual/Dparams_loc
	int nnon0 = 0;
	for (int j=0;j<nParETP;j++) {       // derivatives over reference track parameters
	  if (SmallerAbs(deriv[j],kAlmostZeroD)) continue;
	  nnon0++;
	  fDLoc[fNDLocTot] = deriv[j];  // store non-0 derivative
	  fIDLoc[fNDLocTot] = j;          // and variable id
	  fNDLocTot++;
	}
	int lp0 = pnt->GetMinLocVarID();    // point may depend on material variables starting from this one
	int lp1 = pnt->GetMaxLocVarID();    // and up to this one (exclusive)
	for (int j=lp0;j<lp1;j++) {       // derivatives over material variables
	  if (SmallerAbs(deriv[j],kAlmostZeroD)) continue;
	  nnon0++;
	  fDLoc[fNDLocTot] = deriv[j];  // store non-0 derivative
	  fIDLoc[fNDLocTot] = j;          // and variable id
	  fNDLocTot++;
	}
	//
	fNDLoc[fNResid] = nnon0;          // local derivatives done, store thier number for this residual
	fNResid++;
      }
    }
    if (pnt->ContainsMaterial()) {     // material point can add 4 or 5 otrhogonal pseudo-measurements
      int nmatpar = pnt->GetNMatPar();  // residuals (correction expectation value)
      const float* expMatCorr = pnt->GetMatCorrExp(); // expected corrections (diagonalized)
      const float* expMatCov  = pnt->GetMatCorrCov(); // their diagonalized error matrix
      int offs  = pnt->GetMaxLocVarID() - nmatpar;    // start of material variables
      // here all derivatives are 1 = dx/dx
      for (int j=0;j<nmatpar;j++) {
	fNDGlo[fNResid]  = 0;
	fResid[fNResid]  = expMatCorr[j];
	fResErr[fNResid] = Sqrt(expMatCov[j]);
	fNDLoc[fNResid] = 1;                        // only 1 non-0 derivative
	fDLoc[fNDLocTot] = 1.0;
	fIDLoc[fNDLocTot]  = offs + j;              // variable id
	fNDLocTot++;
	fNResid++;
      }
    }
  }  
  //
  return kTRUE;
}

//________________________________________________
void AliAlgMPRecord::Resize(Int_t nresid, Int_t nloc, Int_t nglo)
{
  // resize container
  if (nresid>fNResidBook) {
    delete[] fNDLoc;
    delete[] fNDGlo;
    delete[] fResid;
    delete[] fResErr;
    fNDLoc = new Short_t[nresid];
    fNDGlo = new Int_t[nresid];
    fResid   = new Double32_t[nresid];
    fResErr  = new Double32_t[nresid];
    fNResidBook = nresid;
  }
  if (nloc>fNDLocTotBook) {
    delete[] fIDLoc;
    delete[] fDLoc;
    fIDLoc = new Short_t[nloc];
    fDLoc= new Double32_t[nloc];
    fNDLocTotBook = nloc;
  }
  if (nglo>fNDGloTotBook) {
    delete[] fIDGlo;
    delete[] fDGlo;
    fIDGlo = new Int_t[nglo];
    fDGlo= new Double32_t[nglo];
    fNDGloTotBook = nglo;
  }
  //
}

//____________________________________________
void AliAlgMPRecord::Print(const Option_t *) const
{
  // print info
  //
  printf("Track %d Event TimeStamp:%d Run:%d\n",fTrackID,fTimeStamp,GetRun());
  printf("NRes: %3d NLoc: %3d NGlo:%3d | Stored: Loc:%3d Glo:%5d\n",
	 fNResid,fNVarLoc,fNVarGlo,fNDLocTot,fNDGloTot);
  //
  int curLoc=0,curGlo=0;
  for (int ir=0;ir<fNResid;ir++) {
    int ndloc = fNDLoc[ir], ndglo = fNDGlo[ir];
    printf("Res:%3d %+e (%+e) | NDLoc:%3d NDGlo:%4d\n",ir,fResid[ir],fResErr[ir],ndloc,ndglo);
    //
    printf("Local Derivatives:\n");
    Bool_t eolOK = kTRUE;
    for (int id=0;id<ndloc;id++) {
      int jd = id+curLoc;
      printf("[%3d] %+.3e  ",fIDLoc[jd],fDLoc[jd]);
      if (((id+1)%5)==0) {printf("\n"); eolOK = kTRUE;}
      else eolOK = kFALSE;
    }
    if (!eolOK) printf("\n");
    curLoc += ndloc;
    //
    //
    printf("Global Derivatives:\n");
    eolOK = kTRUE;
    for (int id=0;id<ndglo;id++) {
      int jd = id+curGlo;
      printf("[%3d] %+.3e  ",fIDGlo[jd],fDGlo[jd]);
      if (((id+1)%5)==0) {printf("\n"); eolOK = kTRUE;}
      else eolOK = kFALSE;
    }
    if (!eolOK) printf("\n");
    curGlo += ndglo;
    //
  }
}
