#ifndef ALIALGDET_H
#define ALIALGDET_H

#include <TNamed.h>
#include <TObjArray.h>
#include <stdio.h>
class AliESDtrack;
class AliAlgTrack;
class AliAlgPoint;
class AliAlgSens;
class AliAlgVol;
class AliAlgSteer;
class AliTrackPointArray;
class AliExternalTrackParam;

class AliAlgDet : public TNamed
{
 public:
  enum {kInitGeomDone=BIT(14),kInitDOFsDone=BIT(15),kDisableBit=BIT(16)};
  //
  AliAlgDet();
  AliAlgDet(const char* name, const char* title="");
  virtual ~AliAlgDet();
  Int_t   GetDetID()                            const {return GetUniqueID();}
  void    SetDetID(UInt_t tp);
  //
  virtual void  CacheReferenceOCDB();
  virtual void  AcknowledgeNewRun(Int_t run);
  virtual void  UpdateL2GRecoMatrices();
  //
  Int_t   VolID2SID(Int_t vid)                  const;
  Int_t   SID2VolID(Int_t sid)                  const {return sid<GetNSensors() ? fSID2VolID[sid] : -1;} //todo
  Int_t   GetNSensors()                         const {return fSensors.GetEntriesFast();}
  Int_t   GetNVolumes()                         const {return fVolumes.GetEntriesFast();}
  Int_t   GetVolIDMin()                         const {return fVolIDMin;}
  Int_t   GetVolIDMax()                         const {return fVolIDMax;}
  Bool_t  SensorOfDetector(Int_t vid)           const {return vid>=fVolIDMin && vid<=fVolIDMax;}
  void    SetAddError(double y, double z);
  const   Double_t* GetAddError()               const {return fAddError;} 
  //
  Int_t   GetNPoints()                          const {return fNPoints;}
  //
  void        SetAlgSteer(AliAlgSteer* s)             {fAlgSteer = s;}
  AliAlgSens* GetSensor(Int_t id)               const {return (AliAlgSens*)fSensors.UncheckedAt(id);}
  AliAlgSens* GetSensorByVolId(Int_t vid)       const {int sid=VolID2SID(vid); return sid<0 ? 0:GetSensor(sid);}
  AliAlgSens* GetSensor(const char* symname)    const {return (AliAlgSens*)fSensors.FindObject(symname);}
  AliAlgVol*  GetVolume(Int_t id)               const {return (AliAlgVol*)fVolumes.UncheckedAt(id);}
  AliAlgVol*  GetVolume(const char* symname)    const {return (AliAlgVol*)fVolumes.FindObject(symname);}
  //
  virtual Int_t InitGeom();
  virtual Int_t AssignDOFs();
  virtual void  InitDOFs();
  virtual void  Terminate();
  virtual void  AddVolume(AliAlgVol* vol);
  virtual void  DefineVolumes();
  virtual void  DefineMatrices();
  virtual void  Print(const Option_t *opt="")    const;
  virtual Int_t ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack,Bool_t inv=kFALSE);
  virtual AliAlgPoint* TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trp);
  virtual void  UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const;
  virtual void  SetUseErrorParam(Int_t v=0);
  Int_t         GetUseErrorParam()                   const {return fUseErrorParam;}
  //
  virtual Bool_t AcceptTrack(const AliESDtrack* trc) const = 0;
  //
  virtual AliAlgPoint* GetPointFromPool();
  virtual void ResetPool();
  //
  void      Disable()                                     {SetBit(kDisableBit);SetObligatory(kFALSE);}
  Bool_t    IsDisabled()                            const {return TestBit(kDisableBit);}
  //
  void      SetInitGeomDone()                             {SetBit(kInitGeomDone);}
  Bool_t    GetInitGeomDone()                       const {return TestBit(kInitGeomDone);}
  //
  void      SetInitDOFsDone()                             {SetBit(kInitDOFsDone);}
  Bool_t    GetInitDOFsDone()                       const {return TestBit(kInitDOFsDone);}
  //
  Int_t     GetNDOFs()                              const {return fNDOFs;}
  //
  void      SetTrackFlagSel(ULong_t f)                    {fTrackFlagSel = f;}
  ULong_t   GetTrackFlagSel()                       const {return fTrackFlagSel;}
  void      SetNPointsSel(Int_t n)                        {fNPointsSel = n;}
  Int_t     GetNPointsSel()                         const {return fNPointsSel;}
  Bool_t    IsObligatory()                          const {return fObligatory;}
  void      SetObligatory(Bool_t v=kTRUE);
  //
  void      WritePedeParamFile(FILE* flOut, const Option_t *opt="") const;
  void      WriteCalibrationResults()               const;
  void      WriteAlignmentResults()                 const;
  //
 protected:
  void     SortSensors();
  //
  // ------- dummies ---------
  AliAlgDet(const AliAlgDet&);
  AliAlgDet& operator=(const AliAlgDet&);
  //
 protected:
  //
  Int_t     fNDOFs;                      // number of DOFs free
  Int_t     fVolIDMin;                   // min volID for this detector (for sensors only)
  Int_t     fVolIDMax;                   // max volID for this detector (for sensors only)
  Int_t     fNSensors;                   // number of sensors (i.e. volID's)
  Int_t*    fSID2VolID;                    //[fNSensors] table of conversion from VolID to sid
  //
  // Track selection
  Bool_t    fObligatory;               // detector must be present in the track
  ULong_t   fTrackFlagSel;             // flag for track selection
  Int_t     fNPointsSel;               // min number of points to require                 
  //
  Int_t     fUseErrorParam;          // signal that points need to be updated using track info, 0 - no
  Double_t  fAddError[2];            // additional error increment for measurement
  TObjArray fSensors;                // all sensors of the detector
  TObjArray fVolumes;                // all volumes of the detector  
  //
  // this is transient info
  Int_t     fNPoints;                //! number of points from this detector
  Int_t     fPoolNPoints;            //! number of points in the pool
  Int_t     fPoolFreePointID;        //! id of the last free point in the pool
  TObjArray fPointsPool;             //! pool of aligment points
  //
  AliAlgSteer* fAlgSteer;            // pointer to alignment steering object
  //
  ClassDef(AliAlgDet,1);             // base class for detector global alignment
};


#endif
