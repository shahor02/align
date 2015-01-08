#ifndef ALIALGDET_H
#define ALIALGDET_H

#include <TNamed.h>
class AliESDtrack;
class AliAlgTrack;

class AliAlgDet : public TNamed
{
 public:
  AliAlgDet();
  AliAlgDet(const char* name, const char* title="");
  //
  Int_t   GetVolIDMin()              const {return fVolIDMin;}
  Int_t   GetVolIDMax()              const {return fVolIDMax;}

  void    SetVolIDMin(Int_t v)             {fVolIDMin = v;}
  void    SetVolIDMax(Int_t v)             {fVolIDMax = v;}
  //
  void    AcknowledgeNewRun(Int_t run);

  virtual Bool_t ProcessTrack(const AliESDtrack* esdTr, AliAlgTrack* fAlgTrack);

 protected:
  
  Int_t fVolIDMin;                   // min volID for this detector
  Int_t fVolIDMax;                   // max volID for this detector

  ClassDef(AliAlgDet,1)              // base class for detector global alignment
};


#endif
