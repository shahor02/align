#ifndef ALIALGDETTOF_H
#define ALIALGDETTOF_H

#include "AliAlgDet.h"

class AliAlgDetTOF : public AliAlgDet
{
 public:
  AliAlgDetTOF();
  AliAlgDetTOF(const char* name, const char* title="");
  virtual ~AliAlgDetTOF();
  //
  virtual void DefineVolumes();  
  //
  Bool_t PresentInTrack(const AliESDtrack* trc) const;

 protected:

  ClassDef(AliAlgDetTOF,1);
};

#endif
