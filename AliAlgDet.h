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

  Bool_t  VIDofDetector(Int_t id)    const {return id>=fVolIDMin && idM=<=fVolIDMax;}

  virtual Bool_t ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* fAlgTrack);
  virtual AliAlgPoint* TrackPoint2AlgPoint(int pntId, const TrackPointArray* trp);

 protected:
  
  Int_t fVolIDMin;                   // min volID for this detector
  Int_t fVolIDMax;                   // max volID for this detector
  //
  TClonesArray fSensorT2G;           // sensor tracking-to-global matrices
  //
  ClassDef(AliAlgDet,1)              // base class for detector global alignment
};


#endif
