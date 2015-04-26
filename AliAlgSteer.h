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
class AliAlgMPRecord;

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
  void  InitDetectors();
  void  InitDOFs();
  
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
  Int_t GetMinDetAcc()                                    const {return fMinDetAcc;}
  void  SetMinDetAcc(Int_t n)                                   {fMinDetAcc=n;}
  //
  //----------------------------------------
  AliAlgPoint* GetRefPoint()                              const {return (AliAlgPoint*)&fRefPoint;}
  //
  AliAlgMPRecord* GetMPRecord()                           const {return (AliAlgMPRecord*)fMPRecord;}
  AliAlgTrack* GetAlgTrack()                              const {return (AliAlgTrack*)fAlgTrack;}
  Bool_t ProcessTrack(const AliESDtrack* esdTr);
  Bool_t ProcessTrack(const AliESDCosmicTrack* esdCTr);
  UInt_t  AcceptTrack(const AliESDtrack* esdTr)            const;
  UInt_t  AcceptTrack(const AliESDtrack* esdPairCosm[kNCosmLegs]) const;
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByDetID(Int_t i)                  const {return fDetPos[i]<0 ? 0:fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVolID(Int_t id)                 const;
  void       ResetDetectors();
  Int_t      GetNDOFs()                                   const {return fNDOFs;}
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
  Int_t         fNDOFs;                                   // number of degrees of freedom
  Int_t         fRunNumber;                               // current run number
  AliAlgMPRecord* fMPRecord;                              // MP record 
  AliAlgTrack*  fAlgTrack;                                // current alignment track 
  AliAlgDet*    fDetectors[kNDetectors];                  // detectors participating in the alignment
  Int_t         fDetPos[kNDetectors];                     // entry of detector in the fDetectors array
  //
  Int_t         fMinDetAcc;                               // min number of detector required in track
  Double_t*     fDOFPars;                                 //[fNDOFs] parameters for free DOFs
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
