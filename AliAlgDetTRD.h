#ifndef ALIALGDETTRD_H
#define ALIALGDETTRD_H

#include "AliAlgDet.h"

class AliAlgDetTRD : public AliAlgDet
{
 public:
  AliAlgDetTRD(const char* title="");
  virtual ~AliAlgDetTRD();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc) const;
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
