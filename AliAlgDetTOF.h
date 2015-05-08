#ifndef ALIALGDETTOF_H
#define ALIALGDETTOF_H

#include "AliAlgDet.h"

class AliAlgDetTOF : public AliAlgDet
{
 public:
  AliAlgDetTOF(const char* title="");
  virtual ~AliAlgDetTOF();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc) const;
  //
 protected:
  //
  // -------- dummies --------
  AliAlgDetTOF(const AliAlgDetTOF&);
  AliAlgDetTOF& operator=(const AliAlgDetTOF&);
  //
 protected:

  ClassDef(AliAlgDetTOF,1);
};

#endif
