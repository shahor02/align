#ifndef ALIALGSENSTRD_H
#define ALIALGSENSTRD_H

#include "AliAlgSens.h"


class TObjArray;


/*--------------------------------------------------------
  TRD sensor
  -------------------------------------------------------*/

class AliAlgSensTRD : public AliAlgSens
{
 public:
  AliAlgSensTRD(const char* name=0, Int_t vid=0, Int_t iid=0);
  virtual ~AliAlgSensTRD();
  //
  virtual void   SetTrackingFrame();
  //
 protected:
  //
  ClassDef(AliAlgSensTRD,1)
};


#endif
