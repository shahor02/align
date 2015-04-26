#ifndef ALIALGDETTPC_H
#define ALIALGDETTPC_H

#include "AliAlgDet.h"

class AliAlgDetTPC : public AliAlgDet
{
 public:
  AliAlgDetTPC();
  AliAlgDetTPC(const char* name, const char* title="");
  virtual ~AliAlgDetTPC();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc) const;

 protected:

  ClassDef(AliAlgDetTPC,1);
};

#endif
