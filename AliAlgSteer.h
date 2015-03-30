#ifndef ALIALGSTEER_H
#define ALIALGSTEER_H

#include "AliGeomManager.h"
#include "AliAlgTrack.h"
#include "AliESDtrack.h"

class AliAlgDet;

/*--------------------------------------------------------
  Steering class for the global alignment. Responsible for feeding the track data 
  to participating detectors and preparation of the millepede input.
  -------------------------------------------------------*/

class AliAlgSteer : public TObject
{
 public:
  enum {kNLrSkip=4};
  enum {kITS,kTPC,kTRD,kTOF,kHMPID,kNDetectors, kUndefined};
  //
  AliAlgSteer();
  virtual ~AliAlgSteer();
  void  Init();
  //
  void  AddDetector(UInt_t id, AliAlgDet* det=0);
  void  AddDetector(AliAlgDet* det);
  //
  void  AcknowledgeNewRun(Int_t run);
  void  SetRunNumber(Int_t run);
  Int_t GetRunNumber()                                    const {return fRunNumber;}
  //
  //----------------------------------------
  AliAlgTrack* GetAlgTrack()                              const {return (AliAlgTrack*)fAlgTrack;}
  Bool_t ProcessTrack(const AliESDtrack* esdTr);
  Bool_t AcceptTrack(const AliESDtrack* esdTr)            const;
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByDetID(Int_t i)                  const {return fDetPos[i]<0 ? 0:fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVolID(Int_t id)                 const;
  //----------------------------------------
  //
  virtual void Print(const Option_t *opt="")              const;
  //
 protected:


 protected:
  //
  Int_t         fNDet;                                    // number of deectors participating in the alignment
  Int_t         fRunNumber;                               // current run number
  AliAlgTrack*  fAlgTrack;                                // current alignment track 
  AliAlgDet*    fDetectors[kNDetectors];                  // detectors participating in the alignment
  Int_t         fDetPos[kNDetectors];                     // entry of detector in the fDetectors array
  //
  static const Int_t   fgkSkipLayers[kNLrSkip];           //  detector layers for which we don't need module matrices
  static const Char_t* fgkDetectorName[kNDetectors];      // names of detectors
  //
  ClassDef(AliAlgSteer,1)
};

#endif
