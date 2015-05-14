#ifndef ALIALGDETTRD_H
#define ALIALGDETTRD_H

#include "AliAlgDet.h"

/*--------------------------------------------------------
  TRD detector wrapper
  -------------------------------------------------------*/

class AliAlgDetTRD : public AliAlgDet
{
 public:
  AliAlgDetTRD(const char* title="");
  virtual ~AliAlgDetTRD();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc,Int_t trtype) const;
  //
 protected:
  //
  // -------- dummies --------
  AliAlgDetTRD(const AliAlgDetTRD&);
  AliAlgDetTRD& operator=(const AliAlgDetTRD&);
  //
 protected:

  ClassDef(AliAlgDetTRD,1);
};

#endif
