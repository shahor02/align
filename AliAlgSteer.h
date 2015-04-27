#ifndef ALIALGSTEER_H
#define ALIALGSTEER_H

#include "AliGeomManager.h"
#include "AliAlgTrack.h"
#include "AliAlgPoint.h"
#include "AliSymMatrix.h"

#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TString.h>

class AliESDEvent;
class AliESDtrack;
class AliESDCosmicTrack;
class AliAlgDet;
class AliAlgMPRecord;
class TTree;
class TFile;


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
  enum {kInpStat,kAccStat,kNStatCl};
  enum {kRun,kEventColl,kEventCosm,kTrackColl,kTrackCosm, kMaxStat};
  //
  AliAlgSteer();
  virtual ~AliAlgSteer();
  void     InitDetectors();
  void     InitDOFs();  
  //
  void     AddDetector(UInt_t id, AliAlgDet* det=0);
  void     AddDetector(AliAlgDet* det);
  //
  void     AcknowledgeNewRun(Int_t run);
  void     SetRunNumber(Int_t run);
  Int_t    GetRunNumber()                                 const {return fRunNumber;}
  Bool_t   GetFieldOn()                                   const {return fFieldOn;}
  void     SetFieldOn(Bool_t v=kTRUE) {fFieldOn = v;}
  Bool_t   IsCosmicEvent()                                const {return fCosmicEvent;}
  void     SetCosmicEvent(Bool_t v=kTRUE) {fCosmicEvent = v;}
  Float_t  GetStat(int cls, int tp)                       const {return fStat[cls][tp];}
  //
  void     SetESDEvent(const AliESDEvent* ev)                   {fESDEvent = ev;}
  const    AliESDEvent* GetESDEvent()                     const {return fESDEvent;}
  void     SetESDtrack(const AliESDtrack* tr, int i=0)          {fESDTrack[i] = tr;}
  const    AliESDtrack* GetESDtrack(int i=0)              const {return fESDTrack[i];}
  //
  // Track selection
  Double_t GetPtMin()                                     const {return fPtMin;}
  void     SetPtMin(double pt=0.3)                              {fPtMin = pt;}  
  Double_t GetEtaMax()                                    const {return fEtaMax;}
  void     SetEtaMax(double eta=1.5)                            {fEtaMax = eta;}  
  Int_t    GetMinDetAcc()                                 const {return fMinDetAcc;}
  void     SetMinDetAcc(Int_t n)                                {fMinDetAcc=n;}
  //
  //----------------------------------------
  AliAlgPoint* GetRefPoint()                              const {return (AliAlgPoint*)&fRefPoint;}
  //
  AliAlgMPRecord* GetMPRecord()                           const {return (AliAlgMPRecord*)fMPRecord;}
  AliAlgTrack* GetAlgTrack()                              const {return (AliAlgTrack*)fAlgTrack;}
  Bool_t  ProcessEvent(const AliESDEvent* esdEv); 
  Bool_t  ProcessTrack(const AliESDtrack* esdTr);
  Bool_t  ProcessTrack(const AliESDCosmicTrack* esdCTr);
  UInt_t  AcceptTrack(const AliESDtrack* esdTr)            const;
  UInt_t  AcceptTrack(const AliESDtrack* esdPairCosm[kNCosmLegs]) const;
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByDetID(Int_t i)                  const {return fDetPos[i]<0 ? 0:fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVolID(Int_t id)                 const;
  void       ResetDetectors();
  Int_t      GetNDOFs()                                   const {return fNDOFs;}
  //----------------------------------------
  // output related
  void SetMPRecFileName(const char* name="mpRecord.root");
  const char* GetMPRecFileName()                          const {return fMPFileName.Data();}
  void   CloseMPOutput();
  void   InitMPOutput();
  Bool_t StoreProcessedTrack();
  void   PrintStatistics() const;
  //
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
  Bool_t        fFieldOn;                                 // field on flag
  Bool_t        fCosmicEvent;                             // cosmic event flag
  AliAlgTrack*  fAlgTrack;                                // current alignment track 
  AliAlgDet*    fDetectors[kNDetectors];                  // detectors participating in the alignment
  Int_t         fDetPos[kNDetectors];                     // entry of detector in the fDetectors array
  //  
  // Track selection
  Int_t         fMinDetAcc;                               // min number of detector required in track
  Double_t      fPtMin;                                   // min pT of tracks to consider
  Double_t      fEtaMax;                                  // eta cut on tracks
  //
  Double_t*     fDOFPars;                                 //[fNDOFs] parameters for free DOFs
  //
  AliAlgPoint   fRefPoint;                                //! reference point for track definition
  //
  const AliESDEvent* fESDEvent;                           //! externally set event
  const AliESDtrack* fESDTrack[kNCosmLegs];               //! externally set ESD tracks
  //
  // statistics
  Float_t fStat[kNStatCl][kMaxStat];                      // processing statistics
  static const Char_t* fgkStatClName[kNStatCl];           // stat classes names
  static const Char_t* fgkStatName[kMaxStat];             // stat type names  
  //
  AliAlgMPRecord* fMPRecord;                              // MP record 
  TTree*          fMPRecTree;                             //! tree to store MP record
  TFile*          fMPRecFile;                             //! file to store MP record tree
  TString         fMPFileName;                            //  file name for output
  //
  //
  static const Int_t   fgkSkipLayers[kNLrSkip];           //  detector layers for which we don't need module matrices
  static const Char_t* fgkDetectorName[kNDetectors];      // names of detectors
  //
  ClassDef(AliAlgSteer,1)
};

#endif
