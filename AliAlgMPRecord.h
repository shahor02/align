#ifndef ALIALGMPRECORD_H
#define ALIALGMPRECORD_H

#include <TObject.h>
class AliAlgTrack;

class AliAlgMPRecord : public TObject
{
 public:
  AliAlgMPRecord();
  virtual ~AliAlgMPRecord();
  //
  Int_t        GetRun()         const {return GetUniqueID();}
  void         SetRun(Int_t r)        {SetUniqueID(r);}
  UInt_t       GetTimeStamp()   const {return fTimeStamp;}
  void         SetTimeStamp(UInt_t t) {fTimeStamp = t;}
  Int_t        GetTrackID()     const {return fTrackID;}
  void         SetTrackID(UInt_t t)   {fTrackID = t;}
  //
  Int_t        GetNVarGlo()     const {return fNVarGlo;}
  void         SetNVarGlo(int n)      {fNVarGlo = n;}
  //
  Int_t        GetNResid()      const {return fNResid;}
  Int_t        GetNVarLoc()     const {return fNVarLoc;}
  //
  Int_t        GetNDLoc(int id) const {return fNDLoc[id];}
  Int_t        GetNDGlo(int id) const {return fNDGlo[id];}
  Double_t     GetResid(int id) const {return fResid[id];}
  Double_t     GetResErr(int id) const {return fResErr[id];}
  //
  //
  Bool_t       FillTrack(const AliAlgTrack* trc);
  //
  void         Resize(Int_t nresid, Int_t nloc, Int_t nglo);
  //
  virtual void Clear(const Option_t *opt="");
  virtual void Print(const Option_t *opt="") const;
 protected:
  //
  Int_t        fTrackID;         // track in the event
  UInt_t       fTimeStamp;       // event time stamp
  Short_t      fNResid;          // number of residuals for the track (=2 npoints)
  Short_t      fNVarLoc;         // number of local variables for the track
  Int_t        fNVarGlo;         // number of global variables defined
  Int_t        fNDLocTot;        // total number of non-zero local derivatives
  Int_t        fNDGloTot;        // total number of non-zero global derivatives
  //
  Short_t*     fNDLoc;           //[fNResid] number of non-0 local derivatives per residual
  Int_t*       fNDGlo;           //[fNResid] number of non-0 global derivatives per residual
  Double32_t*  fResid;           //[fNResid] residuals
  Double32_t*  fResErr;          //[fNResid] error associated to residual
  //
  Short_t*     fIDLoc;           //[fNDLocTot] ID of local variables for non-0 local derivatives
  Int_t*       fIDGlo;           //[fNDGloTot] ID of global variables for non-0 global derivatives  
  Double32_t*  fDLoc;            //[fNDLocTot] non-0 local derivatives
  Double32_t*  fDGlo;            //[fNDGloTot] non-0 global derivatives
  //
  // aux info
  Short_t      fNResidBook;      //! number of slots booked for residuals
  Short_t      fNDLocTotBook;    //! number of slots booked for local derivatives
  Int_t        fNDGloTotBook;    //! number of slots booked for global derivatives
  //
  ClassDef(AliAlgMPRecord,1);
};


#endif
