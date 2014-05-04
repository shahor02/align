#ifndef ALIALGTRACK_H
#define ALIALGTRACK_H

#include <AliExternalTrackParam.h>
#include <TObjArray.h>
#include <TArrayD.h>
#include "AliAlgPoint.h"

class AliAlgTrack: public AliExternalTrackParam
{
 public:
  enum {kFieldONBit=BIT(14),kResidDoneBit=BIT(15),kDerivDoneBit=BIT(16)};
  enum {kNKinParBOFF=4                       // N params for ExternalTrackParam part w/o field
	,kNKinParBON=5                       // N params for ExternalTrackParam part with field
	,kNMSPar=2                           // N params per MS act       
	,kNELosPar=1                         // N params per e.loss act
  };
  enum {kMSTheta1,kMSTheta2,kELoss,kNMatDOFs};

  AliAlgTrack();
  virtual ~AliAlgTrack();
  void         DefineDOFs();
  Double_t     GetMass()                         const {return fMass;}
  Int_t        GetNPoints()                      const {return fPoints.GetEntriesFast();}
  AliAlgPoint* GetPoint(int i)                   const {return (AliAlgPoint*)fPoints[i];}
  void         AddPoint(AliAlgPoint* p)                {fPoints.AddLast(p);}
  void         SetMass(double m)                       {fMass = m;}
  Int_t        GetNLocPar()                      const {return fNLocPar;}
  Int_t        GetNLocExtPar()                   const {return fNLocExtPar;}
  //
  virtual void Clear(Option_t *opt="");
  virtual void Print(Option_t *opt="")           const;
  //
  Bool_t PropagateToPoint(AliExternalTrackParam& tr, const AliAlgPoint* pnt);
  Bool_t PropagateToPoint(AliExternalTrackParam* trSet, int nTr, const AliAlgPoint* pnt);
  //
  Bool_t CalcResiduals(const double *params);
  Bool_t CalcResidDeriv(const double *params);
  //
  Bool_t GetFieldON()                            const {return TestBit(kFieldONBit);}
  void   SetFieldON(Bool_t v=kTRUE)                    {SetBit(kFieldONBit,v);}
  Bool_t GetResidDone()                          const {return TestBit(kResidDoneBit);}
  void   SetResidDone(Bool_t v=kTRUE)                  {SetBit(kResidDoneBit,v);}
  Bool_t GetDerivDone()                          const {return TestBit(kDerivDoneBit);}
  void   SetDerivDone(Bool_t v=kTRUE)                  {SetBit(kDerivDoneBit,v);}
  // propagation methods
  //  Bool_t ApplyMS(AliExternalTrackParam& trPar, double tms,double pms);
  Bool_t ApplyMS(AliExternalTrackParam& trPar, double ms1,double ms2);
  Bool_t ApplyELoss(AliExternalTrackParam& trPar, double dE);
  Bool_t ApplyELoss(AliExternalTrackParam& trPar, const AliAlgPoint* pnt);
  //
  //  Bool_t ApplyMS(AliExternalTrackParam* trSet, int ntr, double tms,double pms);
  Bool_t ApplyMS(AliExternalTrackParam* trSet, int ntr, double ms1,double ms2);
  Bool_t ApplyELoss(AliExternalTrackParam* trSet, int ntr, double dE);
  Bool_t ApplyELoss(AliExternalTrackParam* trSet, int ntr, const AliAlgPoint* pnt);
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
  void   RichardsonDeriv(const AliExternalTrackParam* trSet, const double *delta, const AliAlgPoint* pnt, double& derY, double& derZ);
  //
 protected: 
  
  static Double_t RichardsonExtrap(double *val, int ord=1);
  static Double_t RichardsonExtrap(const double *val, int ord=1);

 protected:

  //
  Int_t     fNLocPar;                    // number of local params
  Int_t     fNLocExtPar;                 // number of local params for the external track param
  Double_t  fMass;                       // assumed mass
  TObjArray fPoints;                     // alignment points
  TArrayD   fResid[2];                   // residuals array
  TArrayD   fDeriv[2];                   // derivatives array  
  Double_t  *fResidA[2];                 // fast access to residuals
  Double_t  *fDerivA[2];                 // fast access to derivatives
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

#endif
