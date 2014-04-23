#ifndef ALIALGDETECTOR_H
#define ALIALGDETECTOR_H

#include <TNamed.h>
class AliESDtrack;
class AliAlgTrack;

class AliAlgDetector : public TNamed
{
 public:
  AliAlgDetector();
  AliAlgDetector(const char* name, const char* title="");
  //
  Int_t   GetVolIDMin()              const {return fVolIDMin;}
  Int_t   GetVolIDMax()              const {return fVolIDMax;}

  void    SetVolIDMin(Int_t v)       const {return fVolIDMin = v;}
  void    SetVolIDMax(Int_t v)       const {return fVolIDMax = v;}
  //

  virtual Bool_t ProcessTrack(const AliESDtrack* esdTr, AliAlgTrack* fAlgTrack);

 protected:
  
  Int_t fVolIDMin;                   // min volID for this detector
  Int_t fVolIDMax;                   // max volID for this detector

  ClassDef(AliAlgDetector,1)
};


#endif
