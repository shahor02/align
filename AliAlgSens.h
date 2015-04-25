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
  enum {kSkipBit=BIT(14)};

  AliAlgSens(const char* name=0, Int_t vid=0, Int_t iid=0);
  virtual ~AliAlgSens();
  //
  virtual void AddChild(AliAlgVol*);
  //
  Int_t GetVolID()                              const  {return (Int_t)GetUniqueID();}
  void  SetVolID(Int_t v)                              {SetUniqueID(v);}
  //
  Int_t GetInternalID()                         const  {return fIntID;}
  void  SetInternalID(Int_t v)                         {fIntID = v;}
  //
  // derivatives calculation
  virtual void DPosTraDParLOC(const double *tra, double* deriv) const;
  virtual void DPosTraDParTRA(const double *tra, double* deriv) const;
  virtual void DPosTraDParLOC(const AliAlgVol* parent, const double *tra, double* deriv) const;
  virtual void DPosTraDParTRA(const AliAlgVol* parent, const double *tra, double* deriv) const;
  //
  void GetModifiedMatrixT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetModifiedMatrixT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const;

  // void GetDeltaMatrixTra(TGeoHMatrix& deltaM, const Double_t *delta) const; ??
  //  void DeltaTra2DeltaLoc(const TGeoHMatrix& deltaTra, TGeoHMatrix& deltaLoc) const; ??
  //
  void            SetAddError(double y, double z)            {fAddError[0]=y;fAddError[1]=z;}
  const Double_t* GetAddError()                        const {return fAddError;} 
  //
  virtual void   PrepareMatrixT2L();
  //
  virtual void   SetTrackingFrame();
  virtual Bool_t IsSensor()                       const {return kTRUE;}
  virtual void   Print(const Option_t *opt="")    const;
  //
  void           SetSkip(Bool_t v=kTRUE)                {SetBit(kSkipBit,v);}
  Bool_t         GetSkip()                        const {return TestBit(kSkipBit);}
  //
 protected:
  //
  virtual Bool_t  IsSortable()                         const {return kTRUE;}
  virtual Int_t   Compare(const TObject* a)            const;
  //
 protected:
  //
  Int_t    fIntID;                    // internal id within the detector
  Double_t fAddError[2];              // additional error increment for measurement
  //
  ClassDef(AliAlgSens,1)
};


#endif
