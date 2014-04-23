#ifndef ALIGLOALG_H
#define ALIGLOALG_H

#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "AliAlgTrack.h"

class TGeoHMatrix;

class AliGloAlg : public TObject
{
 public:
  enum {kMatStart,kNMat};
  enum {kNLrSkip=4};
  enum {kITS,kTPC,kTRD,kTOF,kHMPID,kNDetectors};
  enum {kGeomIdeal,    // ideal geometry matrices
	kGeomOld,      // matrices of aligned geometry used for the reco of data being used
	kGeomCurr,     // matrices of aligned geometry used for the alignment
	kNGeoms
  };  
  //
  AliGloAlg();
  virtual ~AliGloAlg();
  virtual void   LocalInit();
  //
  void  LoadMatrices(Int_t geomTyp);
  void  AddDetector(const char* name);


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
  AliAlgDetector* GetDetector(Int_t i)                    const {return (AliAlgDetector*)fDetectors[i];}
  Bool_t AddDetector(AliAlgDetector* det);
  //----------------------------------------

 protected:
  //
  Int_t         fNDet;                            // number of deectors participating in the alignment
  AliAlgTrack*  fAlgTrack;                        // current alignment track 
  TObjArray     fDetectors;                       // detectors participating in the alignment



  Int_t fMatrixID[2][AliGeomManager::kLastLayer]; // start and N matrices for each layer in fMatrix
  TClonesArray* fMatrixT2G[kNGeoms];              // tracking to global matrices for different geometries
  TClonesArray* fMatrixT2L[kNGeoms];              // tracking to local matrices for different geometries
  //
  static const Int_t  fgkSkipLayers[kNLrSkip];      //! detector layers for which we don't need module matrices
  static const char*  fgkDetectorName[kNDetectors]; //! names of detectors
  //
  ClassDef(AliGloAlg,1)
};

//__________________________________________________________________
inline Int_t AliGloAlg::GetMatrixID(Int_t volID) const
{
  // get stored matrix ID from volID
  int mod=0, lr = AliGeomManager::VolUIDToLayer(volID,mod);
  return GetMatrixID(lr,mod);
}

//__________________________________________________________________
inline Int_t AliGloAlg::GetMatrixID(Int_t lr, Int_t mod) const
{
  // get stored matrix ID from volID
  return fMatrixID[kNMat][lr]>0 ? fMatrixID[kMatStart][lr] + mod : -1;
}

//__________________________________________________________________
inline TGeoHMatrix* AliGloAlg::GetMatrixT2G(Int_t geomTyp, Int_t volID)  const
{
  // get module tracking to global matrix
  int id = GetMatrixID(volID);
  return id<0 ? 0 : (TGeoHMatrix*)fMatrixT2G[geomTyp]->At(id);
}

//__________________________________________________________________
inline TGeoHMatrix* AliGloAlg::GetMatrixT2L(Int_t geomTyp, Int_t volID)  const
{
  // get module tracking to local matrix
  int id = GetMatrixID(volID);
  return id<0 ? 0 : (TGeoHMatrix*)fMatrixT2L[geomTyp]->At(id);
}

#endif
