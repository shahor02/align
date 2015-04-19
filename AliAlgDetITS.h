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
  //
  Bool_t PresentInTrack(const AliESDtrack* trc) const;

  void   SetAddErrorLr(int ilr, double sigY, double sigZ);
  void   SetSkipLr(int ilr);

 protected:

  ClassDef(AliAlgDetITS,1);
};

#endif
