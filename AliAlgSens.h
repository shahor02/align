#ifndef ALIALGSENS_H
#define ALIALGSENS_H

#include "AliAlgVol.h"
#include <TMath.h>

class TObjArray;


/*--------------------------------------------------------
  Smallest alignment volume in detector branch, where the actual measurement is done.
  Apart from local2global matrix has also Tracking2Local matrix
  -------------------------------------------------------*/

class AliAlgSens : public AliAlgVol
{
 public:
  AliAlgSens(const char* name=0, Int_t vid=0, Int_t iid=0);
  virtual ~AliAlgSens();
  //
  virtual void AddChild(AliAlgVol*);
  //
  const TGeoHMatrix&  GetMatrixT2L()             const {return fMatT2L;}
  void  SetMatrixT2L(const TGeoHMatrix& m)             {fMatT2L = m;}
  //
  Int_t GetVolID()                              const  {return (Int_t)GetUniqueID();}
  void  SetVolID(Int_t v)                              {SetUniqueID(v);}
  //
  Int_t GetInternalID()                         const  {return fIntID;}
  void  SetInternalID(Int_t v)                         {fIntID = v;}
  //
  Double_t GetXTracking()                        const {return fX;}
  Double_t GetAlpTracking()                      const {return fAlp;}
  //
  virtual void DPosTraDParLoc(const double *tra, double* deriv) const;
  virtual void DPosTraDParLoc(const AliAlgVol* parent, const double *tra, double* deriv) const;
  //
  void GetModifiedMatrixL2G(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetModifiedMatrixT2L(TGeoHMatrix& matMod, const Double_t *delta) const;

  void GetDeltaMatrixTra(TGeoHMatrix& deltaM, const Double_t *delta) const;
  void DeltaTra2DeltaLoc(const TGeoHMatrix& deltaTra, TGeoHMatrix& deltaLoc) const;
  //
  virtual void   PrepareMatrixT2L();
  //
  virtual void   SetTrackingFrame();
  virtual Bool_t IsSensor()                       const {return kTRUE;}
  virtual void   Print(const Option_t *opt="")    const;
  //
  static Double_t SectAlpha(int i)                      {if (i>8)i-=18; return (i+0.5)*TMath::Pi()/9;}
  static Double_t SectDAlpha(int i)                     {return TMath::Pi()/9;}
  //
 protected:
  virtual Bool_t  IsSortable()                         const {return kTRUE;}
  virtual Int_t   Compare(const TObject* a)            const;

 protected:
  //
  Int_t    fIntID;                    // internal id within the detector
  Double_t fX;                        // tracking frame X offset
  Double_t fAlp;                      // tracking frame alpa
  TGeoHMatrix fMatT2L;                // tracking to local matrix
  //
  ClassDef(AliAlgSens,1)
};


#endif
