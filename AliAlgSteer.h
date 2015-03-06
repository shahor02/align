#ifndef ALIALGSTEER_H
#define ALIALGSTEER_H

#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "AliAlgTrack.h"
#include "AliESDtrack.h"

class TGeoHMatrix;
class AliAlgDet;

/*--------------------------------------------------------
  Steering class for the global alignment. Responsible for feeding the track data 
  to participating detectors and preparation of the millepede input.
  -------------------------------------------------------*/

class AliAlgSteer : public TObject
{
 public:
  enum {kVIDmin,kVIDmax,kNVIDLim};
  enum {kMatStart,kNMat};
  enum {kNLrSkip=4};
  enum {kITS,kTPC,kTRD,kTOF,kHMPID,kNDetectors};
  enum {kGeomIdeal,    // ideal geometry matrices
	kGeomOld,      // matrices of aligned geometry used for the reco of data being used
	kGeomCurr,     // matrices of aligned geometry used for the alignment
	kNGeoms
  };  
  //
  AliAlgSteer();
  virtual ~AliAlgSteer();
  virtual void   LocalInit();
  //
  void  LoadMatrices(Int_t geomTyp);
  void  AddDetector(const char* name);
  
  void  AcknowledgeNewRun(Int_t run);
  void  SetRunNumber(Int_t run);
  Int_t GetRunNumber()                                    const {return fRunNumber;}

  Int_t GetMatrixID(Int_t volID)                          const;
  Int_t GetMatrixID(Int_t lr, Int_t mod)                  const;
  Int_t GetMatLrStart(Int_t lr)                           const {return fMatrixID[kMatStart][lr];}
  Int_t GetNMatLr(Int_t lr)                               const {return fMatrixID[kNMat][lr];}
  TGeoHMatrix* GetMatrixT2G(Int_t geomTyp, Int_t volID)   const;
  TGeoHMatrix* GetMatrixT2L(Int_t geomTyp, Int_t volID)   const;
  //
  //----------------------------------------
  Bool_t ProcessTrack(const AliESDtrack* esdTr);
  Bool_t AcceptTrack(const AliESDtrack* esdTr)            const;
  AliAlgDet* GetDetector(Int_t i)                         const {return fDetectors[i];}
  AliAlgDet* GetDetectorByType(Int_t i)                   const {return fDetPos[i]<0 ? 0 : fDetectors[fDetPos[i]];}
  AliAlgDet* GetDetectorByVID(Int_t id)                   const;
  Bool_t AddDetector(AliAlgDet* det);
  //----------------------------------------
  //
 protected:
  Bool_t        VIDofDetector(Int_t id)                   const;


 protected:
  //
  Int_t         fNDet;                                    // number of deectors participating in the alignment
  Int_t         fRunNumber;                               // current run number
  AliAlgTrack*  fAlgTrack;                                // current alignment track 
  AliAlgDet*    fDetectors[kNDetectors];                  // detectors participating in the alignment
  UShort_t      fDetVolMinMax[kNDetectors][kNVIDLim];     // min/max vol ID for each detector
  Int_t         fDetPos[kNDetectors];                     // entry of detector in the fDetectors array
  //
  Int_t fMatrixID[2][AliGeomManager::kLastLayer];         // start and N matrices for each layer in fMatrix
  TClonesArray* fMatrixT2G[kNGeoms];                      // tracking to global matrices for different geometries
  TClonesArray* fMatrixT2L[kNGeoms];                      // tracking to local matrices for different geometries
  //
  static const Int_t  fgkSkipLayers[kNLrSkip];           //  detector layers for which we don't need module matrices
  static const char*  fgkDetectorName[kNDetectors];      //  names of detectors
  //
  ClassDef(AliAlgSteer,1)
};

//__________________________________________________________________
inline Int_t AliAlgSteer::GetMatrixID(Int_t volID) const
{
  // get stored matrix ID from volID
  int mod=0, lr = AliGeomManager::VolUIDToLayer(volID,mod);
  return GetMatrixID(lr,mod);
}

//__________________________________________________________________
inline Int_t AliAlgSteer::GetMatrixID(Int_t lr, Int_t mod) const
{
  // get stored matrix ID from volID
  return fMatrixID[kNMat][lr]>0 ? fMatrixID[kMatStart][lr] + mod : -1;
}

//__________________________________________________________________
inline TGeoHMatrix* AliAlgSteer::GetMatrixT2G(Int_t geomTyp, Int_t volID)  const
{
  // get module tracking to global matrix
  int id = GetMatrixID(volID);
  return id<0 ? 0 : (TGeoHMatrix*)fMatrixT2G[geomTyp]->At(id);
}

//__________________________________________________________________
inline TGeoHMatrix* AliAlgSteer::GetMatrixT2L(Int_t geomTyp, Int_t volID)  const
{
  // get module tracking to local matrix
  int id = GetMatrixID(volID);
  return id<0 ? 0 : (TGeoHMatrix*)fMatrixT2L[geomTyp]->At(id);
}

//__________________________________________________________________
inline AliAlgDet* AliAlgSteer::GetDetectorByVID(Int_t id) const
{
  // get detector according to volume ID
  for (int i=fNDet;i--;) if (VIDofDetector(i)) return fDetectors[i];
  return 0;
}

//__________________________________________________________________
inline Bool_t AliAlgSteer::VIDofDetector(Int_t id) const
{
  // check if vid belongs to detector
  return id>=fDetVolMinMax[i][kVIDmin] && id<=fDetVolMinMax[i][kVIDmax]);
}

#endif
