#ifndef ALIALGSTEER_H
#define ALIALGSTEER_H

#include "AliGeomManager.h"
#include "AliAlgTrack.h"
#include "AliAlgPoint.h"
#include "AliSymMatrix.h"

#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TString.h>
#include <TArrayF.h>
#include <TArrayI.h>

class AliESDEvent;
class AliESDtrack;
class AliESDCosmicTrack;
class AliESDVertex;
class AliAlgDet;
class AliAlgMPRecord;
class TTree;
class TFile;
//
class Mille;


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
  enum MPOut_t {kMille,kMPRec,kMilleMPRec};
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
  Int_t    GetVtxMinCont()                                const {return fVtxMinCont;}
  void     SetVtxMinCont(int n)                                 {fVtxMinCont = n;}
  Int_t    GetVtxMaxCont()                                const {return fVtxMaxCont;}
  void     SetVtxMaxCont(int n)                                 {fVtxMaxCont = n;}
  Int_t    GetVtxMinContCS()                              const {return fVtxMinContCS;}
  void     SetVtxMinContCS(int n)                               {fVtxMinContCS = n;}
  //
  Bool_t   CheckDetectorPattern(UInt_t patt)              const;
  void     SetObligatoryDetector(Int_t detID, Bool_t v=kTRUE);
  void     SetEventSpeciiSelection(UInt_t sel)                  {fSelEventSpecii = sel;}
  UInt_t   GetEventSpeciiSelection()                      const {return fSelEventSpecii;}
  //
  void     SetVertex(const AliESDVertex* v)                     {fVertex = v;}
  const AliESDVertex* GetVertex()                         const {return fVertex;}
  //
  //----------------------------------------
  Float_t*   GetGloParVal()                               const {return (Float_t*)fGloParVal;}
  Float_t*   GetGloParErr()                               const {return (Float_t*)fGloParErr;}
  //
  AliAlgPoint* GetRefPoint()                              const {return (AliAlgPoint*)&fRefPoint;}
  //
  AliAlgMPRecord* GetMPRecord()                           const {return (AliAlgMPRecord*)fMPRecord;}
  AliAlgTrack* GetAlgTrack()                              const {return (AliAlgTrack*)fAlgTrack;}
  Bool_t     ProcessEvent(const AliESDEvent* esdEv); 
  Bool_t     ProcessTrack(const AliESDtrack* esdTr);
  Bool_t     ProcessTrack(const AliESDCosmicTrack* esdCTr);
  UInt_t     AcceptTrack(const AliESDtrack* esdTr, Bool_t strict=kTRUE)    const;
  UInt_t     AcceptTrackCosmic(const AliESDtrack* esdPairCosm[kNCosmLegs]) const;
  Bool_t     CheckSetVertex(const AliESDVertex* vtx);
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByDetID(Int_t i)                  const {return fDetPos[i]<0 ? 0:fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVolID(Int_t id)                 const;
  void       ResetDetectors();
  Int_t      GetNDOFs()                                   const {return fNDOFs;}
  //----------------------------------------
  // output related
  Bool_t FillMPRecData();
  Bool_t FillMilleData();
  Int_t  GetMPOutType()                                   const {return fMPOutType;}
  void   SetMPOutType(MPOut_t t)                                {fMPOutType = t;}
  void   SetMPDatFileName(const char* name="mpData");
  void   SetMPParFileName(const char* name="mpParam.txt");
  const  char* GetMPDatFileName()                          const {return fMPDatFileName.Data();}
  const  char* GetMPParFileName()                          const {return fMPParFileName.Data();}
  void   CloseMPRecOutput();
  void   CloseMilleOutput();
  void   InitMPRecOutput();
  void   InitMIlleOutput();
  Bool_t StoreProcessedTrack();
  void   PrintStatistics() const;
  //
  void   GenPedeParamFile(const Option_t *opt="") const;
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
  UInt_t        fSelEventSpecii;                          // consider only these event specii
  UInt_t        fObligatoryDetPattern;                    // pattern of obligatory detectors
  Bool_t        fCosmicSelStrict;                         // if true, each cosmic track leg selected like separate track
  Int_t         fMinDetAcc;                               // min number of detector required in track
  Double_t      fPtMin;                                   // min pT of tracks to consider
  Double_t      fEtaMax;                                  // eta cut on tracks
  Int_t         fVtxMinCont;                              // require min number of contributors in Vtx
  Int_t         fVtxMaxCont;                              // require max number of contributors in Vtx  
  Int_t         fVtxMinContCS;                            // min number of contributors to use as constraint
  //
  Float_t*      fGloParVal;                               //[fNDOFs] parameters for free DOFs
  Float_t*      fGloParErr;                               //[fNDOFs] parameters for free DOFs
  //
  AliAlgPoint   fRefPoint;                                //! reference point for track definition
  //
  const AliESDEvent* fESDEvent;                           //! externally set event
  const AliESDtrack* fESDTrack[kNCosmLegs];               //! externally set ESD tracks
  const AliESDVertex* fVertex;                            //! event vertex
  //
  // statistics
  Float_t fStat[kNStatCl][kMaxStat];                      // processing statistics
  static const Char_t* fgkStatClName[kNStatCl];           // stat classes names
  static const Char_t* fgkStatName[kMaxStat];             // stat type names  
  //
  MPOut_t         fMPOutType;                             // Format to store MP data
  Mille*          fMille;                                 // Mille interface
  AliAlgMPRecord* fMPRecord;                              // MP record 
  TTree*          fMPRecTree;                             //! tree to store MP record
  TFile*          fMPRecFile;                             //! file to store MP record tree
  TString         fMPDatFileName;                         //  file name for records binary data output
  TString         fMPParFileName;                         //  file name for MP steering params
  TArrayF         fMilleDBuffer;                          //! buffer for Mille Derivatives output
  TArrayI         fMilleIBuffer;                          //! buffer for Mille Indecis output
  //
  static const Int_t   fgkSkipLayers[kNLrSkip];           // detector layers for which we don't need module matrices
  static const Char_t* fgkDetectorName[kNDetectors];      // names of detectors
  static const Char_t* fgkMPDataExt[kMilleMPRec];         // extensions for MP2 binary data 
  //
  ClassDef(AliAlgSteer,1)
};

#endif
