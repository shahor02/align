#ifndef ALIALGVOL_H
#define ALIALGVOL_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TGeoMatrix.h>

class TObjArray;


/*--------------------------------------------------------
  Base class of alignable volume. Has at least geometric 
  degrees of freedom + user defined calibration DOFs.
  The name provided to constructor must be the SYMNAME which
  AliGeomManager can trace to geometry.
  -------------------------------------------------------*/

class AliAlgVol : public TNamed
{
 public:
  enum DOFGeom_t {kDOFTX,kDOFTY,kDOFTZ,kDOFPH,kDOFTH,kDOFPS,kNDOFGeom,kAllGeomDOF=0x3F};
  enum {kNDOFMax=32};
  enum Frame_t {kLOC,kTRA,kNVarFrames};  // variation frames defined
  enum {kInitDOFsDoneBit=BIT(14),kSkipBit=BIT(15)};
  //
  AliAlgVol(const char* symname=0);
  virtual ~AliAlgVol();
  //
  const char* GetSymName()                       const {return GetName();}
  //
  void       InitDOFs(Int_t &cntDOFs);
  //
  Frame_t    GetVarFrame()                       const {return fVarFrame;}
  void       SetVarFrame(Frame_t f)                    {fVarFrame = f;}
  Bool_t     IsFrameTRA()                        const {return fVarFrame == kTRA;}
  Bool_t     IsFrameLOC()                        const {return fVarFrame == kLOC;}
  //
  Bool_t     IsFreeDOFGeom(DOFGeom_t dof)        const {return (fDOF&(0x1<<dof))!=0;}
  void       SetFreeDOFGeom(DOFGeom_t dof)             {fDOF |= 0x1<<dof;}
  void       SetFreeDOFPattern(UInt_t pat)             {fDOF = pat;}
  UInt_t     GetFreeDOFPattern()                 const {return fDOF;}
  UInt_t     GetFreeDOFGeomPattern()             const {return fDOF&kAllGeomDOF;}

  //
  AliAlgVol* GetParent()                         const {return fParent;}
  void       SetParent(AliAlgVol* par)                 {fParent = par; if (par) par->AddChild(this);}
  Int_t      CountParents()                      const;
  //
  Int_t      GetNChildren()                      const {return fChildren ? fChildren->GetEntriesFast():0;}
  AliAlgVol* GetChild(int i)                     const {return fChildren ? (AliAlgVol*)fChildren->UncheckedAt(i):0;}
  virtual void AddChild(AliAlgVol* ch);
  //
  Double_t GetXTracking()                        const {return fX;}
  Double_t GetAlpTracking()                      const {return fAlp;}
  //
  Int_t      GetNProcessedPoints()               const {return fNProcPoints;}
  void       IncNProcessedPoints(Int_t step=1)         {fNProcPoints += step;}
  void       SetNProcessedPoints(Int_t v)              {fNProcPoints = v;}
  //
  Float_t*   GetParVals()                        const {return fParVals;}
  Double_t   GetParVal(int par)                  const {return fParVals[par];}
  Double_t   GetParErr(int par)                  const {return fParErrs[par];}
  Double_t   GetParConstraint(int par)           const {return fParCstr[par];}
  Int_t      GetParOffset(Int_t par)             const {return fParOffs[par];}
  //
  void       SetParVals(Double_t *vl,Int_t npar);          
  void       SetParVal(Int_t par,Double_t v=0)          {fParVals[par] = v;}
  void       SetParErr(Int_t par,Double_t e=0)          {fParErrs[par] = e;}
  void       SetParConstraint(Int_t par,Double_t s=1e6) {fParCstr[par] = s>0. ? s:0.0;}
  //
  void       SetNDOFTot(Int_t n=kNDOFGeom);
  Int_t      GetNDOFTot()                        const  {return fNDOFTot;}
  Int_t      GetNDOFFree()                       const  {return fNDOFFree;}
  Int_t      GetNDOFGeomFree()                   const  {return fNDOFGeomFree;}
  Int_t      GetParOffs(Int_t par)               const  {return fParOffs[par];}
  Int_t      GetParGloOffs(Int_t par)            const  {return fParOffs[par]<0?-1:fFirstParOffs+fParOffs[par];}
  void       SetFirstParOffs(Int_t id)                  {fFirstParOffs=id;}
  void       SetParOffs(Int_t par,Int_t offs)           {fParOffs[par]=offs;}
  //
  virtual void   PrepareMatrixT2L();
  virtual void   SetTrackingFrame();
  //
  const TGeoHMatrix&  GetMatrixL2G()             const {return fMatL2G;}
  const TGeoHMatrix&  GetMatrixL2GOrig()         const {return fMatL2GOrig;}
  void  SetMatrixL2G(const TGeoHMatrix& m)             {fMatL2G = m;}
  void  SetMatrixL2GOrig(const TGeoHMatrix& m)         {fMatL2GOrig = m;}
  virtual void   PrepareMatrixL2G();
  virtual void   PrepareMatrixL2GOrig();
  //
  void  GetMatrixT2G(TGeoHMatrix& m)             const;
  //
  const TGeoHMatrix&  GetMatrixT2L()             const {return fMatT2L;}
  void  SetMatrixT2L(const TGeoHMatrix& m);
  //
  void  Delta2Matrix(TGeoHMatrix& deltaM, const Double_t *delta)         const;
  void  GetModifiedMatrixL2G(TGeoHMatrix& matMod, const Double_t *delta) const;
  //
  // preparation of variation matrices
  void GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetDeltaT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const;
  void GetDeltaT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta, const TGeoHMatrix& relMat) const;
  //
  void    SetSkip(Bool_t v=kTRUE)                   {SetBit(kSkipBit,v);}
  Bool_t  GetSkip()                           const {return TestBit(kSkipBit);}
  //
  void    SetInitDOFsDone()                         {SetBit(kInitDOFsDoneBit);}
  Bool_t  GetInitDOFsDone()                   const {return TestBit(kInitDOFsDoneBit);}
  //
  virtual Bool_t IsSensor()                   const {return kFALSE;}
  //
  virtual void   Print(const Option_t *opt="")  const;
  //
  static void    SetDefGeomFree(UChar_t patt)    {fgDefGeomFree = patt;}
  static UChar_t GetDefGeomFree()                {return fgDefGeomFree;}
  //
 protected:
  //
  Frame_t    fVarFrame;               // Variation frame for this volume
  Double_t   fX;                      // tracking frame X offset
  Double_t   fAlp;                    // tracking frame alpa
  //
  Char_t     fNDOFTot;                // number of degrees of freedom, including fixed ones
  UInt_t     fDOF;                    // bitpattern degrees of freedom
  Char_t     fNDOFGeomFree;           // number of free geom degrees of freedom
  Char_t     fNDOFFree;               // number of all free degrees of freedom
  //
  AliAlgVol* fParent;             // parent volume
  TObjArray* fChildren;           // array of childrens
  //
  Int_t      fNProcPoints;        // n of processed points
  Int_t      fFirstParOffs;           // entry of the 1st free parameter in the global results array
  Char_t*    fParOffs;                //[fNDOFTot] offset for every parameters wrt the 1st free in global results array
  Float_t*   fParVals;            //[fNDOFTot]  values of the fitted params
  Float_t*   fParErrs;            //[fNDOFTot]  errors of the fitted params
  Float_t*   fParCstr;            //[fNDOFTot]  Gaussian type constraint on parameter, 0 means fixed param
  //
  TGeoHMatrix fMatL2G;            // local to global matrix, including current alignment
  TGeoHMatrix fMatL2GOrig;        // local to global matrix, ideal
  TGeoHMatrix fMatT2L;            // tracking to local matrix
  //
  static const char* fgkFrameName[kNVarFrames];
  static UInt_t      fgDefGeomFree;
  //
  ClassDef(AliAlgVol,1)
};

//-------------------------------------------------------------
inline void AliAlgVol::SetParVals(Double_t *vl,Int_t npar)
{
  // set parameters
  for (int i=npar;i--;) fParVals[i] = vl[i];
}

inline void AliAlgVol::GetMatrixT2G(TGeoHMatrix& m) const
{
  // compute tracking to global matrix, i.e. glo = T2G*tra = L2G*loc = L2G*T2L*tra
  m = GetMatrixL2G();
  m *= GetMatrixT2L();
}

#endif
