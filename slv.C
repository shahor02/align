#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAlgMPRecord.h"
#include "AliAlgSteer.h"
#include "AliSymMatrix.h"
#include <TVectorD.h>
#include <TArrayI.h>
#include <TArrayD.h>
#include <TMath.h>
//
#endif

AliAlgSteer* algSteer = 0;


void slv(AliAlgMPRecord* rec)
{
  if (!algSteer) {printf("No AlgSteer\n"); return;}
  //
  int nres = rec->GetNResid();
  // 
  const float *recDGlo = rec->GetArrGlo();
  const float *recDLoc = rec->GetArrLoc();
  const short *recLabLoc = rec->GetArrLabLoc();
  const int   *recLabGlo = rec->GetArrLabGlo();
  int nvloc = rec->GetNVarLoc();
  int ngt = rec->GetNDGloTot();
  int nvglo = 0;
  //
  //
  AliSymMatrix* matpG = new AliSymMatrix(nvloc);
  TVectorD* vecpG = new TVectorD(nvloc);
  TVectorD& vecG = *vecpG;
  AliSymMatrix& matG = *matpG;
  //
  AliSymMatrix* matp = new AliSymMatrix(nvloc);
  TVectorD* vecp = new TVectorD(nvloc);
  TVectorD& vec = *vecp;
  AliSymMatrix& mat = *matp;
  //
  // residuals, accounting for global solution
  double *resid = new Double_t[nres];
  for (int irs=0;irs<nres;irs++) {
    double resOr = rec->GetResid(irs);
    resid[irs] = resOr;
    //
    int ndglo = rec->GetNDGlo(irs);
    int ndloc = rec->GetNDLoc(irs);
    for (int ig=0;ig<ndglo;ig++) {
      int lbI = recLabGlo[ig];
      int idP = algSteer->Label2ParID(lbI);
      if (idP<0) {printf("Did not find parameted for label %d\n",lbI);exit(1);}
      double parVal = algSteer->GetGloParVal()[idP];
      resid[irs] += parVal*recDGlo[ig];
    }
    //
    double  sg2inv = rec->GetResErr(irs);
    sg2inv = 1./sg2inv/sg2inv;
    //
    // Build matrix to solve local parameters
    for (int il=0;il<ndloc;il++) {
      int lbLI = recLabLoc[il]; // id of local variable
      vecG[lbLI] -= recDLoc[il]*resid[irs]*sg2inv;
      vec[lbLI]  -= recDLoc[il]*resOr*sg2inv;
      for (int jl=il+1;jl--;) {
	int lbLJ = recLabLoc[jl]; // id of local variable
	matG(lbLI,lbLJ) += recDLoc[il]*recDLoc[jl]*sg2inv;	  
	mat(lbLI,lbLJ)  += recDLoc[il]*recDLoc[jl]*sg2inv;	  
      }
    }
    //
    recLabGlo += ndglo; // prepare for next record
    recDGlo   += ndglo;
    recLabLoc += ndloc;
    recDLoc   += ndloc;
    //
  }
  //
  TVectorD vecSol(nvloc);
  TVectorD vecSolG(nvloc);

  if (!matp->SolveChol(vec,vecSol,kFALSE)) {
    printf("Failed to solve original track\n");
    delete matp;
    matp = 0;
  }
  if (!matpG->SolveChol(vecG,vecSolG,kFALSE)) {
    printf("Failed to solve track corrected for globals\n");
    delete matpG;
    matpG = 0;
  }
  // check
  recDGlo = rec->GetArrGlo();
  recDLoc = rec->GetArrLoc();
  recLabLoc = rec->GetArrLabLoc();
  recLabGlo = rec->GetArrLabGlo();
  //
  printf("Sol: L/LG:\n");
  int nExtP = (nvloc%4) ? 5:4;
  for (int i=0;i<nExtP;i++) printf("%+.3e/%+.3e ",vecSol[i],vecSolG[i]);    
  printf("\n");
  Bool_t nln = kTRUE;
  int cntL = 0;
  for (int i=nExtP;i<nvloc;i++) {
    nln = kTRUE;
    printf("%+.3e/%+.3e ",vecSol[i],vecSolG[i]);
    if ( ((++cntL)%4)==0 ) {printf("\n"); nln = kFALSE;}
  }
  if (!nln) printf("\n");
  printf("Check: \n");
  for (int irs=0;irs<nres;irs++) {
    double resOr = rec->GetResid(irs);
    double resL = resOr;
    double resLG = resid[irs];
    double  sg = rec->GetResErr(irs);
    //
    int ndglo = rec->GetNDGlo(irs);
    int ndloc = rec->GetNDLoc(irs);
    for (int il=0;il<ndloc;il++) {
      int lbLI = recLabLoc[il]; // id of local variable
      resL  += recDLoc[il]*vecSol[lbLI];
      resLG += recDLoc[il]*vecSolG[lbLI];
    }
    //
    printf("%3d (%9d) %6.4f | [%+.2e:%+7.2f] [%+.2e:%+7.2f] [%+.2e:%+7.2f]\n",
	   irs,ndglo ? recLabGlo[0]:-1,sg,resOr,resOr/sg,resL,resL/sg,resLG,resLG/sg);
    //
    recLabGlo += ndglo; // prepare for next record
    recDGlo   += ndglo;
    recLabLoc += ndloc;
    recDLoc   += ndloc;
  }  

}

