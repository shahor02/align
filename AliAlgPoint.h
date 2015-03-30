#ifndef ALIALGPOINT_H
#define ALIALGPOINT_H

#include <TObject.h>

class AliExternalTrackParam;


class AliAlgPoint : public TObject
{
 public:
  enum {kMaterialBit=BIT(14)     // point contains material
	,kMeasurementBit=BIT(15) // point contains measurement
	,kVaryELossBit=BIT(16)   // ELoss variation allowed
	,kUseBzOnly=BIT(17)      // use only Bz component (ITS)
	,kInvDir=BIT(18)         // propagation via this point is in decreasing X direction (upper cosmic leg)
  };
  
  AliAlgPoint();
  virtual   ~AliAlgPoint() {}
  //
  void       Init();
  //
  Double_t   GetAlphaSens()          const {return fAlphaSens;}
  Double_t   GetXSens()              const {return fXSens;}
  Double_t   GetXPoint()             const {return fXSens + GetXTracking();}
  Double_t   GetXTracking()          const {return fXYZTracking[0];}
  Double_t   GetYTracking()          const {return fXYZTracking[1];}
  Double_t   GetZTracking()          const {return fXYZTracking[2];}
  const Double_t* GetYZTracking()    const {return &fXYZTracking[1];}
  const Double_t* GetXYZTracking()   const {return fXYZTracking;}
  const Double_t* GetYZErrTracking() const {return fErrYZTracking;}
  //
  Int_t      GetDetID()              const {return fDetID;}
  Int_t      GetSID()                const {return fSID;}
  Int_t      GetMaxLocVarID()        const {return fMaxLocVarID;}
  Bool_t     ContainsMaterial()      const {return TestBit(kMaterialBit);}
  Bool_t     ContainsMeasurement()   const {return TestBit(kMeasurementBit);}
  Bool_t     GetELossVaried()        const {return TestBit(kVaryELossBit);}
  Bool_t     GetUseBzOnly()          const {return TestBit(kUseBzOnly);}
  Bool_t     IsInvDir()              const {return TestBit(kInvDir);}
  //
  Double_t   GetXTimesRho()          const {return fXTimesRho;}
  Double_t   GetX2X0()               const {return fX2X0;}
  void       SetXTimesRho(double v)        {fXTimesRho = v;}
  void       SetX2X0(double v)             {fX2X0 = v;}
  //
  void       SetDetID(Int_t id)                       {fDetID = (Char_t)id;}
  void       SetSID(Int_t id)                         {fSID = (Short_t)id;}
  //
  void       SetMaxLocVarID(Int_t id)                 {fMaxLocVarID = id;}
  void       SetELossVaried(Bool_t v=kTRUE)           {SetBit(kVaryELossBit,v);}
  void       SetContainsMaterial(Bool_t v=kTRUE)      {SetBit(kMaterialBit,v);}
  void       SetContainsMeasurement(Bool_t v=kTRUE)   {SetBit(kMeasurementBit,v);}
  void       SetUseBzOnly(Bool_t v=kTRUE)             {SetBit(kUseBzOnly,v);}
  void       SetInvDir(Bool_t v=kTRUE)                {SetBit(kInvDir,v);}
  //
  void       GetResidualsDiag(const double* pos, double &resU, double &resV) const;
  //
  void       SetAlphaSens(double a)                   {fAlphaSens = a;}
  void       SetXSens(double x)                       {fXSens = x;}
  void       SetXYZTracking(const double r[3])   {for (int i=3;i--;) fXYZTracking[i]=r[i];}
  void       SetXYZTracking(double x,double y,double z);
  void       SetYZErrTracking(double sy2, double syz, double sz2);
  void       SetYZErrTracking(const double *err) {for (int i=3;i--;) fErrYZTracking[i]=err[i];}
  Double_t   GetErrDiag(int i)             const {return fErrDiag[i];}
  //
  Double_t*  GetTrParamWS()                const {return (Double_t*)fTrParamWS;}
  void       SetTrParamWS(const double* param)   {for (int i=5;i--;) fTrParamWS[i] = param[i];}
  //
  Double_t   GetMSSigTheta2()              const {return fMSSigTheta2;}
  void       SetMSSigTheta2(double v)            {fMSSigTheta2 = v;}
  //
  virtual void Print(Option_t* option = "") const;
  virtual void Clear(Option_t* option = "");
  //
  protected:
  virtual Bool_t  IsSortable()                         const {return kTRUE;}
  virtual Int_t   Compare(const TObject* a)            const;
  //
 protected:
  //
  Int_t    fMaxLocVarID;                               // The residuals/derivatives depend on params<fMaxLocVarID.
                                                       // If the point contains materials, fMaxLocVarID also marks
                                                       // the parameters associated with this point
  Char_t   fDetID;                                     // DetectorID
  Short_t  fSID;                                       // sensor ID in the detector
  Double_t fAlphaSens;                                 // Alpha of tracking frame
  Double_t fXSens;                                     // X of tracking frame
  Double_t fXYZTracking[3];                            // X,Y,Z in tracking frame
  Double_t fCosDiagErr;                                // Cos of Phi of rotation in YZ plane which diagonalize errors
  Double_t fSinDiagErr;                                // Sin of Phi of rotation in YZ plane which diagonalize errors
  //
  Double_t fX2X0;                                      // X2X0 seen by the track (including inclination)
  Double_t fXTimesRho;                                 // Density*Length seen by the track (including inclination)
  Double_t fErrYZTracking[3];                          // errors in tracking frame
  Double_t fErrDiag[2];                                // diagonalized errors
  //
  Double_t fMSSigTheta2;                               //! sigma^2 of MS
  Double_t fTrParamWS[5];                              //! workspace for tracks params at this point
  //
  ClassDef(AliAlgPoint,1)
};

//____________________________________________________
inline void AliAlgPoint::SetXYZTracking(double x,double y,double z) 
{
  // assign tracking coordinates
  fXYZTracking[0]=x; fXYZTracking[1]=y; fXYZTracking[2]=z;
}

//____________________________________________________
inline void AliAlgPoint::SetYZErrTracking(double sy2,double syz,double sz2) 
{
  // assign tracking coordinates
  fErrYZTracking[0]=sy2; fErrYZTracking[1]=syz; fErrYZTracking[2]=sz2;
}

#endif
