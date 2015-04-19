#ifndef ALIALGSTEER_H
#define ALIALGSTEER_H

#include "AliGeomManager.h"
#include "AliAlgTrack.h"
#include "AliAlgPoint.h"

#include <TMatrixDSym.h>
#include <TVectorD.h>

#include "AliSymMatrix.h"

class AliESDEvent;
class AliESDtrack;
class AliESDCosmicTrack;

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
  enum {kCosmLow,kCosmUp,kNCosmLegs};
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
  void  SetESDEvent(const AliESDEvent* ev)                      {fESDEvent = ev;}
  const AliESDEvent*  GetESDEvent()                       const {return fESDEvent;}
  //
  //----------------------------------------
  AliAlgPoint* GetRefPoint()                              const {return (AliAlgPoint*)&fRefPoint;}
  //
  AliAlgTrack* GetAlgTrack()                              const {return (AliAlgTrack*)fAlgTrack;}
  Bool_t ProcessTrack(const AliESDtrack* esdTr);
  Bool_t ProcessTrack(const AliESDCosmicTrack* esdCTr);
  Bool_t AcceptTrack(const AliESDtrack* esdTr)            const;
  Bool_t AcceptTrack(const AliESDtrack* esdPairCosm[kNCosmLegs]) const;
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByDetID(Int_t i)                  const {return fDetPos[i]<0 ? 0:fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVolID(Int_t id)                 const;
  void       ResetDetectors();
  //----------------------------------------
  //
  virtual void Print(const Option_t *opt="")              const;
  //
  static Char_t* GetDetNameByDetID(Int_t id)                    {return (Char_t*)fgkDetectorName[id];}
  //
  AliSymMatrix* BuildMatrix(TVectorD &vec);
  Bool_t        TestLocalSolution();
  //
 protected:
  AliAlgSteer(const AliAlgSteer&);
  //
 protected:
  //
  Int_t         fNDet;                                    // number of deectors participating in the alignment
  Int_t         fRunNumber;                               // current run number
  AliAlgTrack*  fAlgTrack;                                // current alignment track 
  AliAlgDet*    fDetectors[kNDetectors];                  // detectors participating in the alignment
  Int_t         fDetPos[kNDetectors];                     // entry of detector in the fDetectors array
  //
  AliAlgPoint   fRefPoint;                                //! reference point for track definition
  //
  const AliESDEvent* fESDEvent;                           //! externally set event
  //
  static const Int_t   fgkSkipLayers[kNLrSkip];           //  detector layers for which we don't need module matrices
  static const Char_t* fgkDetectorName[kNDetectors];      // names of detectors
  //
  ClassDef(AliAlgSteer,1)
};

#endif
