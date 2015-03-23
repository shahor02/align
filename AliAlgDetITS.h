#ifndef ALIALGDETITS_H
#define ALIALGDETITS_H

#include "AliAlgDet.h"

class AliAlgDetITS : public AliAlgDet
{
 public:
  AliAlgDetITS();
  AliAlgDetITS(const char* name, const char* title="");
  virtual ~AliAlgDetITS();
  //
  virtual void DefineVolumes();  
  virtual void PrintHierarchy();
  //
 protected:

  ClassDef(AliAlgDetITS,1);
};

#endif
