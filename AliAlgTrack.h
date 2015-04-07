#ifndef ALIALGTRACK_H
#define ALIALGTRACK_H

#include "AliExternalTrackParam.h"
#include "AliAlgPoint.h"
#include <TObjArray.h>
#include <TArrayD.h>

#define DEBUG 4

class AliAlgTrack: public AliExternalTrackParam
{
 public:

  enum {kCosmicBit=BIT(14),kFieldONBit=BIT(15),kResidDoneBit=BIT(16),kDerivDoneBit=BIT(17)};
  enum {kNKinParBOFF=4                       // N params for ExternalTrackParam part w/o field
	,kNKinParBON=5                       // N params for ExternalTrackParam part with field
	//
	,kNMSPar=4                           // N params per MS act       
	,kNELosPar=1                         // N params per e.loss act
	,kParY=0                             // Y parameter
	,kParZ=1                             // Z parameter
	,kParSnp=2                           // snp parameter
	,kParTgl=3                           // tgl parameter
	,kParq2Pt=4                          // q/pt parameter
	,kNMatDOFs                           // number of paremeters for material effects
  };
  AliAlgTrack();
  virtual ~AliAlgTrack();
  void         DefineDOFs();
  Double_t     GetMass()                         const {return fMass;}
  Double_t     GetMinX2X0Pt2Account()            const {return fMinX2X0Pt2Account;}
  Int_t        GetNPoints()                      const {return fPoints.GetEntriesFast();}
  AliAlgPoint* GetPoint(int i)                   const {return (AliAlgPoint*)fPoints[i];}
  void         AddPoint(AliAlgPoint* p)                {fPoints.AddLast(p);}
  void         SetMass(double m)                       {fMass = m;}
  void         SetMinX2X0Pt2Account(double v)          {fMinX2X0Pt2Account = v;}
  Int_t        GetNLocPar()                      const {return fNLocPar;}
  Int_t        GetNLocExtPar()                   const {return fNLocExtPar;}
  Int_t        GetInnerPointID()                 const {return fInnerPointID;}
  //
  virtual void Clear(Option_t *opt="");
  virtual void Print(Option_t *opt="")           const;
  //
  Bool_t PropagateToPoint(AliExternalTrackParam& tr, const AliAlgPoint* pnt, 
			  int minNSteps,double maxStep,Bool_t matCor, double* matPar=0);
  Bool_t PropagateParamToPoint(AliExternalTrackParam& tr, const AliAlgPoint* pnt); // param only
  Bool_t PropagateParamToPoint(AliExternalTrackParam* trSet, int nTr, const AliAlgPoint* pnt); // params only
  //
  Bool_t CalcResiduals(const double *params=0);
  Bool_t CalcResidDeriv(double *params=0);
  //
  Bool_t IsCosmic()                              const {return TestBit(kCosmicBit);}
  void   SetCosmic(Bool_t v=kTRUE)                     {SetBit(kCosmicBit);}
  Bool_t GetFieldON()                            const {return TestBit(kFieldONBit);}
  void   SetFieldON(Bool_t v=kTRUE)                    {SetBit(kFieldONBit,v);}
  Bool_t GetResidDone()                          const {return TestBit(kResidDoneBit);}
  void   SetResidDone(Bool_t v=kTRUE)                  {SetBit(kResidDoneBit,v);}
  Bool_t GetDerivDone()                          const {return TestBit(kDerivDoneBit);}
  void   SetDerivDone(Bool_t v=kTRUE)                  {SetBit(kDerivDoneBit,v);}
  //
  void   SortPoints();
  Bool_t IniFit();
  Bool_t ProcessMaterials();
  Bool_t CombineTracks(AliExternalTrackParam& trcL, const AliExternalTrackParam& trcU);
  //
  void     SetChi2(double c)                           {fChi2 = c;};
  Double_t GetChi2()                             const {return fChi2;}
  //
  // propagation methods
  void   CopyFrom(const AliExternalTrackParam* etp);
  Bool_t ApplyMatCorr(AliExternalTrackParam& trPar, const Double_t *corrPar, Bool_t eloss);
  Bool_t ApplyMatCorr(AliExternalTrackParam* trSet, int ntr, const Double_t *corrPar, Bool_t eloss);
  Bool_t ApplyELoss(AliExternalTrackParam& trPar, const AliAlgPoint* pnt);
  Bool_t ApplyELoss(AliExternalTrackParam* trSet, int ntr, const AliAlgPoint* pnt);
  //
  // misc methods, perhaps obsolete
  Bool_t ApplyMS(AliExternalTrackParam& trPar, double ms1,double ms2);
  Bool_t ApplyMS(AliExternalTrackParam* trSet, int ntr, double ms1,double ms2);
  Bool_t ApplyELoss(AliExternalTrackParam& trPar, double dE);
  Bool_t ApplyELoss(AliExternalTrackParam* trSet, int ntr, double dE);
  //
  Double_t  GetResidual(int dim, int pntID)       const {return  fResidA[dim][pntID];}
  Double_t *GetDerivative(int dim, int pntID)     const {return &fDerivA[dim][pntID*fNLocPar];}
  //
  void SetParams(AliExternalTrackParam& tr, double x, double alp, const double* par);
  void SetParams(AliExternalTrackParam* trSet, int ntr, double x, double alp, const double* par);
  void SetParam(AliExternalTrackParam& tr, int par, double val);
  void SetParam(AliExternalTrackParam* trSet, int ntr, int par, double val);
  void ModParam(AliExternalTrackParam& tr, int par, double delta);
  void ModParam(AliExternalTrackParam* trSet, int ntr, int par, double delta);
  //
  void RichardsonDeriv(const AliExternalTrackParam* trSet, const double *delta, 
		       const AliAlgPoint* pnt, double& derY, double& derZ);
  //
 protected: 
  //
  Bool_t CalcResidDeriv(double *params,Bool_t invert,int pFrom,int pTo);
  Bool_t CalcResiduals(const double *params,Bool_t invert,int pFrom,int pTo);
  Bool_t FitLeg(AliExternalTrackParam& trc, int pFrom,int pTo, Bool_t &inv);
  Bool_t ProcessMaterials(AliExternalTrackParam& trc, int pFrom,int pTo);
  //
  static Double_t RichardsonExtrap(double *val, int ord=1);
  static Double_t RichardsonExtrap(const double *val, int ord=1);
  //
 protected:

  //
  Int_t     fNLocPar;                    // number of local params
  Int_t     fNLocExtPar;                 // number of local params for the external track param
  Int_t     fInnerPointID;               // ID of inner point in sorted track. For 2-leg cosmics - innermost point of lower leg
  Bool_t    fNeedInv[2];                 // set if one of cosmic legs need inversion
  Double_t  fMinX2X0Pt2Account;          // minimum X2X0/pT accumulated between 2 points worth to account
  Double_t  fMass;                       // assumed mass
  Double_t  fChi2;                       // chi2 with current residuals
  TObjArray fPoints;                     // alignment points
  TArrayD   fResid[2];                   // residuals array
  TArrayD   fDeriv[2];                   // derivatives array  
  TArrayD   fLocPar;                     // local parameters array
  Double_t  *fResidA[2];                 //! fast access to residuals
  Double_t  *fDerivA[2];                 //! fast access to derivatives
  Double_t  *fLocParA;                   //! fast access to local params
  //
  ClassDef(AliAlgTrack,1)
};

//____________________________________________________________________________________________
inline void AliAlgTrack::SetParams(AliExternalTrackParam& tr, double x, double alp, const double* par)
{
  // set track params
  tr.SetParamOnly(x,alp,par);
  if (!GetFieldON()) ((double*)tr.GetParameter())[4] = 0.; // only 4 params are valid
  //
}

//____________________________________________________________________________________________
inline void AliAlgTrack::SetParams(AliExternalTrackParam* trSet, int ntr, double x, double alp, const double* par)
{
  // set parames for multiple tracks (VECTORIZE THIS)
  for (int itr=ntr;itr--;) {
    SetParams(trSet[itr],x,alp,par);
  }
}

//____________________________________________________________________________________________
inline void AliAlgTrack::SetParam(AliExternalTrackParam& tr, int par, double val)
{
  // set track parameter
  ((double*)tr.GetParameter())[par] = val;
}

//____________________________________________________________________________________________
inline void AliAlgTrack::SetParam(AliExternalTrackParam* trSet, int ntr, int par, double val)
{
  // set parames for multiple tracks (VECTORIZE THIS)
  for (int itr=ntr;itr--;) ((double*)trSet[itr].GetParameter())[par] = val;
}

//____________________________________________________________________________________________
inline void AliAlgTrack::ModParam(AliExternalTrackParam & tr, int par, double delta)
{
  // modify track parameter
  ((double*)tr.GetParameter())[par] += delta;
}

//____________________________________________________________________________________________
inline void AliAlgTrack::ModParam(AliExternalTrackParam* trSet, int ntr, int par, double delta)
{
  // modify track parameter (VECTORIZE THOS)
  for (int itr=ntr;itr--;) ModParam(trSet[itr],par,delta);
}

//______________________________________________
inline void AliAlgTrack::CopyFrom(const AliExternalTrackParam* etp)
{
  // assign kinematics
  Set(etp->GetX(),etp->GetAlpha(),etp->GetParameter(),etp->GetCovariance());
}


#endif
