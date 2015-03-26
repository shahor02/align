#ifndef ALIALGSENS_H
#define ALIALGSENS_H

#include "AliAlgVol.h"


class TObjArray;


/*--------------------------------------------------------
  Smallest alignment volume in detector branch, where the actual measurement is done.
  Apart from local2global matrix has also Tracking2Local matrix
  -------------------------------------------------------*/

class AliAlgSens : public AliAlgVol
{
 public:
  AliAlgSens(const char* name=0, Int_t vid=0);
  virtual ~AliAlgSens();
  //
  virtual void AddChild(AliAlgVol*);
  //
  const TGeoHMatrix&  GetMatrixT2L()             const {return fMatT2L;}
  void  SetMatrixT2L(const TGeoHMatrix& m)             {fMatT2L = m;}
  //
  Int_t GetVolID()                              const  {return (Int_t)GetUniqueID();}
  void  SetVolID(Int_t v    )                          {SetUniqueID(v);}
  //
  Double_t GetXTracking()                        const {return fXTracking;}
  Double_t GetAlpTracking()                      const {return fAlpTracking;}
  //
  virtual void DPosTraDParLoc(const double *tra, double* deriv) const;
  virtual void DPosTraDParLoc(const AliAlgVol* parent, const double *tra, double* deriv) const;
  //
  void GetModifiedMatrixG2L(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetModifiedMatrixT2L(TGeoHMatrix& matMod, const Double_t *delta) const;

  void GetDeltaMatrixTra(TGeoHMatrix& deltaM, const Double_t *delta) const;
  void DeltaTra2DeltaLoc(const TGeoHMatrix& deltaTra, TGeoHMatrix& deltaLoc) const;
  //
  virtual Bool_t IsSensor()                            const {return kTRUE;}
  //
 protected:
  virtual Bool_t  IsSortable()                         const {return kTRUE;}
  virtual Int_t   Compare(const TObject* a)            const;

 protected:
  //
  Double_t fXTracking;                // tracking frame X offset
  Double_t fAlpTracking;              // tracking frame alpa
  TGeoHMatrix fMatT2L;                // tracking to local matrix
  //
  ClassDef(AliAlgSens,1)
};


#endif
