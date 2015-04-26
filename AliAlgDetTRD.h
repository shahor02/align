#ifndef ALIALGDETTRD_H
#define ALIALGDETTRD_H

#include "AliAlgDet.h"

class AliAlgDetTRD : public AliAlgDet
{
 public:
  AliAlgDetTRD();
  AliAlgDetTRD(const char* name, const char* title="");
  virtual ~AliAlgDetTRD();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc) const;

 protected:

  ClassDef(AliAlgDetTRD,1);
};

#endif
