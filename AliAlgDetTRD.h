#ifndef ALIALGDETTRD_H
#define ALIALGDETTRD_H

#include "AliAlgDet.h"

/*--------------------------------------------------------
  TRD detector wrapper
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgDetTRD : public AliAlgDet
{
 public:
  //
  enum {kCalibRCCorrDzDtgl,  // correction parameter for NonRC tracklets
	kNCalibParams};  // calibration parameters
  //
  AliAlgDetTRD(const char* title="");
  virtual ~AliAlgDetTRD();
  //
  virtual void DefineVolumes();  
  virtual void Print(const Option_t *opt="")              const;
  //
  Bool_t AcceptTrack(const AliESDtrack* trc,Int_t trtype) const;
  //
  virtual const char* GetCalibDOFName(int i)              const;
  //
  virtual void         WritePedeInfo(FILE* parOut,const Option_t *opt="") const;
  //
  Double_t GetNonRCCorrDzDtglWithCal()                    const {return GetNonRCCorrDzDtgl()+GetParVal(kCalibRCCorrDzDtgl);}
  Double_t GetNonRCCorrDzDtgl()                           const {return fNonRCCorrDzDtgl;}
  void     SetNonRCCorrDzDtgl(double v=1.055)                   {fNonRCCorrDzDtgl = v;}
  //
  const Double_t* GetExtraErrRC()                         const {return fExtraErrRC;} 
  void     SetExtraErrRC(double y=0.2, double z=1.0)            {fExtraErrRC[0]=y;fExtraErrRC[1]=z;}
  //  
 protected:
  //
  // -------- dummies --------
  AliAlgDetTRD(const AliAlgDetTRD&);
  AliAlgDetTRD& operator=(const AliAlgDetTRD&);
  //
 protected:
  //
  Double_t fNonRCCorrDzDtgl;     // correction in Z for non-crossing tracklets
  Double_t fExtraErrRC[2];       // extra errors for RC tracklets
  //
  static const char* fgkCalibDOFName[kNCalibParams];
  //
  ClassDef(AliAlgDetTRD,1);
};

#endif
